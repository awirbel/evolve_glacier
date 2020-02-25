"""
    Create refined mesh geometry

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
from cfg_datadir import data_dir, cmaxF
import numpy as np
from utils import (write_surface)


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

# Load 2D computational mesh from previous time step
fname_surfmesh = "surface_2D"
mesh = Mesh(os.path.join(data_dir, "%s_%04d.xml" % (fname_surfmesh, t1)))

# Read solution from previous time step to define new surface elvevation
Q = FunctionSpace(mesh, "CG", 1)
# print "Read now in MAKE MESH surface and bed from: ", fname
s = Function(Q)
b = Function(Q)
sfilename = "solution_%s" % (fname)
sin1 = HDF5File(mesh.mpi_comm(), os.path.join(data_dir,sfilename+".hdf5"), 'r')
sin1.read(s, "/elevation")
sbfilename = "solutionbed_%s" % (fname)
sinb = HDF5File(mesh.mpi_comm(), os.path.join(data_dir,sbfilename+".hdf5"), 'r')
sinb.read(b, "/bed")

# Make new 3D mesh
write_surface(s, b)

# Update parameters and save to file
params['t1'] = t1 + 1
params['tna'] = tna + 1
params['fname'] = 'glacier_%04d' % (tna)

with open(os.path.join(data_dir, 'data.pkl'), "w") as f:
    pickle.dump(params, f)


print "Update :) "
