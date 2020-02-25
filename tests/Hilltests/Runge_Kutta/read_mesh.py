"""
    Here, mesh boundary are initialised and the mesh is stored to an hdf5 file

    Evolve
    Copyright (C) 2019  Anna Wirbel

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
from cfg_datadir import data_dir
import numpy as np
from utils import save_pvd

# Load parameters
with open(os.path.join(data_dir, "data.pkl"), "r") as f:
    params = pickle.load(f)

# Define parameters from file
fname = params['fname']
nx1 = params['nx1']
nx2 = params['nx2']
ny1 = params['ny1']
ny2 = params['ny2']
t1 = params['t1']

# Print start info
print "begin read and save mesh: ", fname

# Load mesh from file
mesh = Mesh(os.path.join(data_dir, "%s.xml" % (fname)))
# if boundaries and subdomains are defined as physical groups in gmsh
#subdomains = MeshFunction("size_t", mesh, os.path.join(data_dir, "%s_physical_region.xml" % (fname)))
#boundaries = MeshFunction("size_t", mesh, os.path.join(data_dir, "%s_facet_region.xml" % (fname)))

# Define boundaries on 3D mesh using facet normals
tdim = mesh.topology().dim()
mesh.init(tdim-1, tdim)
boundaries = FacetFunction("size_t", mesh)
for f in facets(mesh):
    # all exterior facets
    if f.exterior() and f.normal().z() < 0.0:
        boundaries[f] = 102
    if f.exterior() and f.normal().z() > 0.0:
        boundaries[f] = 101

class BoundaryNorth(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (near(x[1], ny2))

class BoundarySouth(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (near(x[1], ny1))

class BoundaryWest(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (near(x[0], nx1))

class BoundaryEast(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (near(x[0], nx2))

north = BoundaryNorth()
south = BoundarySouth()
west = BoundaryWest()
east = BoundaryEast()

# Mark boundaries
north.mark(boundaries, 106)
south.mark(boundaries, 104)
west.mark(boundaries, 105)
east.mark(boundaries, 103)


# Save mesh to hdf5
m3filename = fname
fu10 = HDF5File(mesh.mpi_comm(), os.path.join(data_dir, m3filename+".hdf5"), 'w')
fu10.write(mesh, "/mesh")
#fu10.write(subdomains, "/subdomains")
fu10.write(boundaries, "/boundaries")


# Save for visualisation
save_pvd("3D_mesh_%s" % (fname), mesh)
#save_pvd("3D_boundaries_%s" % (fname), boundaries)

print "read mesh done - t1: ", t1, mesh.hmin()
