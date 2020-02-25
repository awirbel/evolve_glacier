"""
    This is the utilities file with all functions required from the model to solve the free surface evolution using FEM
    and the FEniCS software framework https://fenicsproject.org/

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

from __future__ import division, print_function
import numpy as np
from dolfin import *
import os
from cfg_datadir import (data_dir, ab_params, acc_params, mb_height,
                         yins, siny, cmaxF)


def readmeshvelocity(fname_vel, nx1, nx2, ny1, ny2, order):

    # Load mesh from file
    mesh = Mesh(os.path.join(data_dir, "%s.xml" % (fname_vel)))

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
    m3filename = fname_vel
    fu10 = HDF5File(mesh.mpi_comm(), os.path.join(data_dir, m3filename+".hdf5"), 'w')
    fu10.write(mesh, "/mesh")
    #fu10.write(subdomains, "/subdomains")
    fu10.write(boundaries, "/boundaries")
