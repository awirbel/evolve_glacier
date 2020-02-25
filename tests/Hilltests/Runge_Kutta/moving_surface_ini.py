"""
    This model creates initial files required to solve the free surface evolution using FEM
    and the FEniCS software framework https://fenicsproject.org/
    and PETSC's constrained snes solver https://www.mcs.anl.gov/petsc/

    Evolve
    Copyright (C) 2019  Anna Wirbel, Alexander Jarosch

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""


from dolfin import *
import pickle
import os
from cfg_datadir import data_dir, cmaxF, TS
import numpy as np
from utils import (ice_mass, get_timestep, save_pvd, set_velocity, mb_change,
                   write_surface, solve_SNES, solve_Netwon, ice_thickness,
                    set_mb, check_change)


#********************Read input data******************
# Load parameters
with open(os.path.join(data_dir, "data.pkl"), "r") as f:
    params = pickle.load(f)

# Define parameters from config file
fname = params['fname']
fname_old = params['fname_old']
t1 = params['t1']
dt_total_old = params['dt_total_old']
dt_ADV = params['dt_ADV']
mass_initial = params['mass_initial']
mass_previous = params['mass_previous']
bed_diff1 = params['bed_diff1']
tna = params['tna']
tnatext = params['tnatext']
dt_overall = params['dt_overall']

# FEniCS parameters
parameters['allow_extrapolation'] = True
# print "begin of surface move mesh: ", fname

# Load 3D geometry mesh from file
mesh = Mesh()
mainname = "%s" % (fname)
hdfin1 = HDF5File(mesh.mpi_comm(), os.path.join(data_dir, mainname+".hdf5"), "r")
hdfin1.read(mesh, "/mesh", False)

# Initialise boundaries on mesh
boundaries = FacetFunction("size_t", mesh)
hdfin1.read(boundaries, "/boundaries")

# Define mesh coordinates and FunctionSpaces on initial mesh
m_coors = mesh.coordinates().copy()
Q = FunctionSpace(mesh, "CG", 1)

#********************Separate surface and bed******************
# Define surface elevation as function on initial mesh
S = Function(Q)
surf_height = Expression('x[2]', degree=1)
S = interpolate(surf_height, Q)

# Create boundary mesh of initial mesh
bmesh = BoundaryMesh(mesh, "exterior")
D = bmesh.topology().dim()
bmesh.init(D-1,D)

#-----------Complex geometries---------
# Create new Mesh Function on boundary mesh to transfer boundary-info from full mesh on boundary mesh
# In glacier case to separate surface and bed geometry
bdim = bmesh.topology().dim()
boundary_boundaries = MeshFunction('size_t', bmesh, bdim)
boundary_boundaries.set_all(0)
for i, facet in enumerate(entities(bmesh, bdim)):
    parent_meshentity = bmesh.entity_map(bdim)[i]
    parent_boundarynumber = boundaries.array()[parent_meshentity]
    boundary_boundaries.array()[i] = parent_boundarynumber

# Create Submesh for surface and bed parts of the boundary mesh
submesh_surface = SubMesh(bmesh, boundary_boundaries, 101)
submesh_bed = SubMesh(bmesh, boundary_boundaries, 102)


#++++++++++++++SURFACE SUBMESH++++++++++++
# Define FunctionSpaces on submesh of the surface
V2 = VectorFunctionSpace(submesh_surface, 'CG', 2)
V2_1 = VectorFunctionSpace(submesh_surface, 'CG', 1)
Q2 = FunctionSpace(submesh_surface, 'CG', 1)

# Define surface elevation function on surface mesh
sd = Function(Q2)
surf_height = Expression('x[2]', degree=1)
sd = interpolate(surf_height, Q2)

# Define coordinates of surface mesh and set x[2] to 0.0 to get a flat mesh that is still 3D but x[2]=0 everywhere
s_coors = submesh_surface.coordinates()
s_old_coors = submesh_surface.coordinates().copy()
s_coors[:,2]= 0.0
sd_ini = sd.copy()


#++++++++++++++BED Submesh++++++++++++++++++++++
# Define FunctionSpace and bed elevation on bed mesh
Q_bed = FunctionSpace(submesh_bed, 'CG', 1)
s_bed = Function(Q_bed)
bed_height = Expression('x[2]', degree=1)
s_bed = interpolate(bed_height, Q_bed)

# Define bed mesh coordinates and set x[2] to 0.0 to get a flat mesh that is still 3D but x[2]=0 everywhere
bed_coors = submesh_bed.coordinates()
bed_old_coors = submesh_bed.coordinates().copy()
bed_coors[:,2]= 0.0
s_ini = s_bed.copy()


#+++++++++++++++++START OF SURFACE EVOLUTION+++++++++++++++++
# Read 2D computational mehs that covers entire domain
D = mesh.topology().dim()
mesh.init(D-1,D)
fname_surfmesh = "surface_2D"
Rmesh = Mesh(os.path.join(data_dir, "%s.xml" % (fname_surfmesh)))

# Define FunctionsSpaces on this mesh
Q3 = FunctionSpace(Rmesh, 'CG', 1)
V3 = VectorFunctionSpace(Rmesh, 'CG', 1)
# Extract dofs and corresponding coordinates of Function Spaces on computational mesh
gdim = Rmesh.geometry().dim()
dofmapQ3 = Q3.dofmap()
dofmapV3 = V3.dofmap()
# list of dofs in Q3, V3
dofsQ3 = dofmapQ3.dofs()
dofsV3 = dofmapV3.dofs()
# list of coordinates of the dofs in Q3, V3
dofsQ3_x = Q3.tabulate_dof_coordinates().reshape((-1, gdim))
dofV3_x = V3.tabulate_dof_coordinates().reshape((-1, gdim))

#*************Surface and bed elevation and MB on 2D computational mesh***********************
# Define surface and bed elevation as Function on 2D computational mesh
# therefore evaluate surface elevation on the flat submesh of surface on this mesh
# s_new is the surface elevation
# s_new_bed is the bed elevation + the bed_diff (difference of bed elevation and surface elevation where there is no ice)
# s_new_bed_ini is the actual bed elevation
s_new = Function(Q3)
s_new_bed = Function(Q3)
s_new_bed_ini = Function(Q3)


for dof, dofQ3_x in zip(dofsQ3, dofsQ3_x):
    s_new.vector()[dof] = sd(dofQ3_x[0], dofQ3_x[1], 0.0)
    s_new_bed.vector()[dof] = s_bed(dofQ3_x[0], dofQ3_x[1], 0.0) + bed_diff1    #this should actually be needed as we say the surface is always bed_diff1 higher than the bed!!
    s_new_bed_ini.vector()[dof] = s_bed(dofQ3_x[0], dofQ3_x[1], 0.0)


# Save functions to hdf5 files
# print "Save surface and bed from: ", fname_old
sfilename = "solution_%s" % (fname_old)
f1u = HDF5File(Rmesh.mpi_comm(),os.path.join(data_dir, sfilename+".hdf5"), 'w')
f1u.write(s_new, "/elevation")
sbfilename = "solutionbed_%s" % (fname_old)
f2u = HDF5File(Rmesh.mpi_comm(),os.path.join(data_dir, sbfilename+".hdf5"), 'w')
f2u.write(s_new_bed_ini, "/bed")
sbifilename = "solutionconst_%s" % (fname_old)
f3u = HDF5File(Rmesh.mpi_comm(),os.path.join(data_dir, sbifilename+".hdf5"), 'w')
f3u.write(s_new_bed, "/const")


with open(os.path.join(data_dir, 'data.pkl'), "w") as f:
    pickle.dump(params, f)


print 'surface move this is done now'
