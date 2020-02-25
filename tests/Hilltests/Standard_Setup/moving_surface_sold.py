"""
    This model solves the free surface evolution using FEM
    and the FEniCS software framework https://fenicsproject.org/
    and PETSC's constrained SNES solver https://www.mcs.anl.gov/petsc/

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
                    set_verticalvelocity,write_surface_v2, solve_SNES,
                     solve_Netwon, ice_thickness,
                    set_mb, check_change, subdomains_vels)


#********************Read input data******************
# Load parameters
with open(os.path.join(data_dir, "data.pkl"), "r") as f:
    params = pickle.load(f)

# Define parameters from configuration file
fname = params['fname']
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
print "begin of surface move mesh: ", fname

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
V = VectorFunctionSpace(mesh, 'CG', 2)
velocity = Function(V)

# Read 3D velocity from file
velfilename = "velocity_%s" % (fname)
hdfin81 = HDF5File(mesh.mpi_comm(), os.path.join(data_dir,velfilename+".hdf5"), 'r')
hdfin81.read(velocity, "/velocity")
hdfin81.close()
# In case user defined function has to be set for testing
#velocity.interpolate(Expression(("2.7e-7",  "1.e-8", "0.00"),  degree=1))

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

# Interpolate velocity from initial mesh onto surface mesh
# NOTE: polynomial degree=1 compared to degree=2 of initial mesh
vel2 = interpolate(velocity, V2_1)

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
# Read 2D computational mesh that covers entire domain
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

#*************Surface and bed elevation and MB on computational mesh***********************
# Define surface and bed elevation as Function on 2D computational mesh
# therefore evaluate surface elevation on the flat submesh of surface on this mesh
# s_new is the surface elevation
# s_new_bed is the bed elevation + the bed_diff (difference of bed elevation and surface elevation where there is no ice)
# s_new_bed_ini is the actual bed elevation
# h_ice1 is the ice thickness
s_new = Function(Q3)
s_new_bed = Function(Q3)
s_new_bed_ini = Function(Q3)
h_ice1 = Function(Q3)

for dof, dofQ3_x in zip(dofsQ3, dofsQ3_x):
    s_new.vector()[dof] = sd(dofQ3_x[0], dofQ3_x[1], 0.0)
    s_new_bed.vector()[dof] = s_bed(dofQ3_x[0], dofQ3_x[1], 0.0) + bed_diff1    #this should actually be needed as we say the surface is always bed_diff1 higher than the bed!!
    s_new_bed_ini.vector()[dof] = s_bed(dofQ3_x[0], dofQ3_x[1], 0.0)
    h_ice1.vector()[dof] = sd(dofQ3_x[0], dofQ3_x[1], 0.0) - (s_bed(dofQ3_x[0], dofQ3_x[1], 0.0) + bed_diff1)

# surf_ini is the initial surface elevation before moving but on the computational mesh
surf_ini = s_new.copy()
# Define ice thickness and ice mass at the start of the computations
hice_start = ice_mass(s_new, s_new_bed, Q3)
hice_t0 = assemble(hice_start * dx)
# Set mass balance rate
mb_func = set_mb(s_new, Q3, dofsQ3, dofsQ3_x)
# if constant value shall be set
#mb_func = Function(Q3)
#mb_func.assign(Constant(0.0))
# Set velocity on the computational mesh
v_new = set_velocity(vel2,s_new, s_new_bed_ini, Q3, V3, dofV3_x, dofsV3, bed_diff1)
v_vert = set_verticalvelocity(velocity, s_new, s_new_bed_ini,  Q3, dofsQ3, dofsQ3_x, bed_diff1)
# Compute exptected mass change due to mass balance rate field
[mb_start, mb0_start] = mb_change(mb_func, hice_start, dt_ADV, Q3)
mb_mass0 = hice_t0 - assemble(mb0_start * dx)

# Save fields for visualisation
# save_pvd("rmesh_hicestart", hice_start, pvd_number=tna)
# save_pvd("rmesh_ft0", mb_func)
# save_pvd("rmesh_snew", s_new, pvd_number=tna)
# save_pvd("rmesh_vnew", v_new, pvd_number=tna)
# save_pvd("rmesh_vert", v_vert, pvd_number=tna)
# save_pvd("rmesh_ice", h_ice1, pvd_number=tna)
# save_pvd("rmesh_bed", s_new_bed, pvd_number=tna)

# Compute surface gradient
DG = VectorFunctionSpace(Rmesh, 'DG', 0)
grad2 = project(grad(s_new), DG)
save_pvd("rmesh_grad2", grad2, pvd_number=tna)


# -------------------- START of SURFACE EVOLUTION -------------------------
# Perform computation on 2D computational mesh
# Required in SUPG stabilisation
h = CellSize(Rmesh)
n = FacetNormal(Rmesh)

# Set mass balance rate
f0 = Function(Q3)
f0.assign(mb_func)

# Initial condition: surface elevation at the beginning of the computations
u0 = Function(Q3)
u0 = s_new

#************************* Define subdomains **************************************
subdomains_Rnew3 = subdomains_vels(Rmesh, v_new, s_new, s_new_bed_ini, grad2, mb_start, dt_ADV, tna, bed_diff1)

# Define integration measures: this has to be defined as we have multiple meshes loaded in here
dx = Measure("dx", domain = Rmesh, subdomain_data = subdomains_Rnew3)
ds = Measure('ds', domain = Rmesh, subdomain_data = subdomains_Rnew3)

# ********************Set up FEM problem********************
# Test and trial functions
du, v = TrialFunction(Q3), TestFunction(Q3)
u = Function(Q3)

# Create domain boundaries
boundary_RDC = FacetFunction("size_t", Rmesh)
for f in facets(Rmesh):
    if f.exterior():
        boundary_RDC[f] = 1

# Set a Dirichlet boundary condition to the non-moving boundaries of the domain
# preC = u0
# bc = DirichletBC(Q3, preC, boundary_RDC, 1)
#bcs = [bc]
# Set natural (Neumann) boundary condition
bcs = []

# Set time step and time discretization
if TS == 'CN':
    print 'TIME DISCRETIZATION: Crank Nicholson'
    # Set time step - for Crank Nicholson
    # Define time step according to CFL
    [dt, velmax] = get_timestep(v_new, v_vert, mb_func, Rmesh, cmaxF)
    dt_form = Constant(dt)
    # Mid-point solution
    u_mid = 0.5 * (u0 + u)
    # Parameters for time stepping loop
    t_ref = 1
    T_ref = int(dt_ADV / dt) + 2

elif TS =='BE':
    print 'TIME DISCRETIZATION: Backward Euler'
    # Set time step - for Backward Euler
    dt_form = Constant(dt_ADV)
    dt = dt_ADV
    # Fully implict
    u_mid = u
    t_ref = 1
    T_ref = 1

# Residual for SUPG stabilisation
r = u - u0 + dt_form * (dot(v_new, grad(u_mid)) - v_vert - f0)

# Variational Form including different terms for different subdomains
F1 = v * (u - u0) * dx(1) + v * (u - u0) * dx(0) + dt_form * (v *(v_new[0] * grad(u_mid)[0] + v_new[1] * grad(u_mid)[1]) * dx(1) - v_vert*v*dx(1) - f0 * v * dx(1) - f0 * v * dx(0))

# Add SUPG stabilisation terms
vnorm = sqrt(dot(v_new, v_new))
F2 = F1 + (h / (2.0 * vnorm)) * dot(v_new, grad(v)) * r * dx(1)

# Add SOLD stabilisation terms following modified Burman approach in the SOLD 2 paper by John
gradunorm = sqrt(dot(grad(u_mid), grad(u_mid)))
rnorm = sqrt(dot(r, r))
tau = h / (2.0 *vnorm)
C = ((tau * vnorm * rnorm) / gradunorm) * ((vnorm * gradunorm) / (vnorm * gradunorm + rnorm))
F = F2 + C * dot(grad(u_mid), grad(v))*dx(1)

# Set constraint for surface elevation
# upper constraint
constraint_u = Expression(("Smax"), Smax=10000.0+DOLFIN_EPS, degree=1)
Smax = interpolate(constraint_u, Q3)
# lower constraint is the bed elevation + bed_diff1
Smin = s_new_bed

# # Save inequality constraints to File for visualisation
# save_pvd("Smax", Smax)
# save_pvd("Smin", Smin)

# create the derivative
dF = derivative(F, u, du)

# Define parameters for time stepping loop
dt_final = dt_ADV
flag_continue = 0
flag_stop = 0
dt_total_old = 0.0
dt_total = 0.0
# Time-stepping loop
while flag_stop < 0.5:

    print "Evolve time step: ", t1, "time step is: ", dt
    if dt >= dt_ADV:
        flag_continue = 1
        dt_new = dt_ADV
        dt_computed = dt_ADV

    # MB-Start------------Compute mass balance as a function of current surface elevation
    mb_func = set_mb(u0, Q3, dofsQ3, dofsQ3_x)
    # Set to Function
    f0.assign(mb_func)
    #save_pvd("rmesh_ft0", mb_func, pvd_number=t_ref)
    # MB-End---------------------------------------------------------------------------

    # ***********TIME STEP*****************
    if TS == 'CN':
        [dt_computed, velmax] = get_timestep(v_new, v_vert, mb_func, Rmesh, cmaxF)
        if (dt_total + dt_computed) > dt_final:
            dt_new = dt_final - dt_total
            dt_form.assign(dt_new)
            flag_continue = 1
            #print 'dt_new is set', dt_new
            dt_computed = dt_new
        else:
            #print 'dt_computed is set', dt_computed
            dt_form.assign(dt_computed)
    elif TS == 'BE':
        dt_computed = dt_ADV
        flag_stop = 1
    # ***********TIME STEP*****************

    # HICE-Start--------------Compute ice thickness second option--------------------------
    hice_new_computed = ice_thickness(Q3, dofsQ3_x, dofsQ3, u0, s_new_bed)
    # HICE-End------------------------------------------------------------------------------

    # Start++++++++++++++++++++++++++++Perform computations++++++++++++++++++++++
    u = solve_Netwon(F, u, bcs, dF, solv=True)
    u = solve_SNES(F, u, du, bcs, Smin, Smax, solv=True)
    u.rename("usnes", "unses")

    # Save solution to file
    save_pvd("Solution_%s_%04d_new" % (tnatext, tna), u, pvd_number=t_ref)

    # MBC-Start----------------Compute mass change-----------------------------
    # Compute actual mass change as diff of hice_previous and hice_new
    hice_new = ice_mass(u, s_new_bed, Q3)
    hice_previous = ice_mass(u0, s_new_bed, Q3)
    hice_new_mass = assemble(hice_new * dx)
    hice_previous_mass = assemble(hice_previous * dx)
    print "This is the ice mass at time", t_ref, hice_new_mass, "at previous time", hice_previous_mass
    # Actual ice mass change within computational time step
    mb_change_solution = assemble(hice_previous * dx) - assemble(hice_new * dx)
    # MBC-End----------------------------------------------------------------------------

    # Compute ice thickness change
    hice_change = check_change(s_new, s_new_bed_ini, u, Q3, dofsQ3_x, dofsQ3)
    # save_pvd("rmesh_hicechange", hice_change, pvd_number=tna)

    # Compute total time step
    dt_total = dt_computed + dt_total

    # Copy solution from previous interval to provide initial condition
    u0.assign(u)

    t_ref = t_ref + 1
    print "Total time step is: ", dt_total
    print int(dt_total), int(dt_final)

    # Check if dt_ADV exceeded within next computation, if TRUE change dt to fill up dt_ADV
    if int(dt_total) == int(dt_final):
        print "stop iterations computations done"
        flag_stop = 1.0


#----------------- Write result to stl to create new 3D geometry START---------------------
write_surface_v2(u, s_new_bed_ini)
#----------------- Write result to stl create new 3D geometry START---------------------


# Compute ice thickness
hice_solution = ice_thickness(Q3, dofsQ3_x, dofsQ3, u, s_new_bed)
save_pvd("rmesh_hicesolution", hice_solution, pvd_number=tna)


# Update parameters and save to file
params['t1'] = t1 + 1
params['tna'] = tna + 1
params['dt_total_old'] = dt_total
params['fname'] = 'glacier_%04d' % (tna)
params['fnameold'] = 'glaciernew%d' % (tna)
params['dt_overall'] = dt_overall

with open(os.path.join(data_dir, 'data.pkl'), "w") as f:
    pickle.dump(params, f)

print 'surface move this is done now'
