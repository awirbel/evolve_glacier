"""
    icetools code Newton iteration version
    original version can be found: https://sourceforge.net/projects/icetools/

    Copyright 2015-2017 Alexander H. Jarosch, modified by Anna Wirbel

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>

"""

from dolfin import *
import pickle
import os
from cfg_datadir import data_dir
from utils import save_pvd

# Load  parameters
with open(os.path.join(data_dir, "data.pkl"), "r") as f:
    params = pickle.load(f)

# Define parameters from file
fname = params['fname']
nx1 = params['nx1']
nx2 = params['nx2']
ny1 = params['ny1']
ny2 = params['ny2']

parameters["form_compiler"]["quadrature_degree"] = 2
print "begin of ice dynamics, glacier mesh name: ", fname

# Load mesh from file and boundaries
mesh = Mesh()
mname = "%s" % (fname)
hdfin3 = HDF5File(mesh.mpi_comm(), os.path.join(data_dir, mname+".hdf5"), "r")
hdfin3.read(mesh, "/mesh", False)
#subdomains = CellFunction("size_t", mesh)
#hdfin3.read(subdomains, "/subdomains")
boundaries = FacetFunction("size_t", mesh)
hdfin3.read(boundaries, "/boundaries")

# define some general constants
g = -9.81               # gravitational constant
rho = 917.0             # fluid density
Aglen = 2.4e-24         # Glen flow parameter for temperate ice (Cuffey & Paterson,2010 p. 73)
nglen = 3.0             # Glen's n
glen_fact = 0.5 * Aglen**(-1.0/nglen)

# Define the body force f, i.e. gravity, driving the flow
f_x0 = 0.0
f_x1 = 0.0
f_x2 = g*rho
f = Constant((f_x0, f_x1, f_x2))

# Define function spaces
P2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
TH = P2 * P1
W = FunctionSpace(mesh, TH)

# Define solution function
w = Function(W)
(u, p) = split(w)
(v, q) = TestFunctions(W)

# Apply a no-slip boundary condition for velocity
noslip = Constant((0.0,0.0, 0.0))
bc0 = DirichletBC(W.sub(0), noslip, boundaries, 102)

# Apply a no flow through domain boundary
bc1 = DirichletBC(W.sub(0).sub(0), 0.0, boundaries, 103)  # for x component
bc2 = DirichletBC(W.sub(0).sub(0), 0.0, boundaries, 105)  # for x component
bc3 = DirichletBC(W.sub(0).sub(1), 0.0, boundaries, 104)  # for y component
bc4 = DirichletBC(W.sub(0).sub(1), 0.0, boundaries, 106)  # for y component

# Collect boundary conditions
bcs = [bc0, bc1, bc2, bc3, bc4]


##### New Code ######
nu = 8e13 # linear viscosity for for first linear guess

epsilon = sym(grad(u))
F = (2*nu*inner(epsilon, grad(v)) - div(u)*q - div(v)*p)*dx - inner(v, f)*dx

dF = derivative(F, w)
pde = NonlinearVariationalProblem(F, w, bcs, dF)
solver = NonlinearVariationalSolver(pde)
solver.parameters["symmetric"] = True
solver.parameters["newton_solver"]["maximum_iterations"] = 5
solver.parameters["newton_solver"]["linear_solver"] = 'superlu_dist'
solver.solve()

# define the nonlinear stokes equation directly:
def visc(u):

    eps_dot = sqrt(0.5*inner(sym(grad(u)),sym(grad(u)))) # second invariant of strain
    nu_out = glen_fact * eps_dot**((1.0 - nglen)/nglen)
#     return nu_out
    return Min(nu_out, 2e15)    # would introduce a viscosity limit

nu = visc(u)

epsilon = sym(grad(u))
F = (2*nu*inner(epsilon, grad(v)) - div(u)*q - div(v)*p)*dx - inner(v, f)*dx

dF = derivative(F, w)
pde = NonlinearVariationalProblem(F, w, bcs, dF)
solver = NonlinearVariationalSolver(pde)
solver.parameters["symmetric"] = True
solver.parameters["newton_solver"]["maximum_iterations"] = 80
solver.parameters["newton_solver"]["error_on_nonconvergence"]  = False
solver.parameters["newton_solver"]["relaxation_parameter"] = 0.6
solver.parameters["newton_solver"]["relative_tolerance"] = 1E-4
solver.parameters["newton_solver"]["linear_solver"] = 'superlu_dist'
solver.solve()

(u, p) = w.split(deepcopy=True)


# Save the final solution
ufilename = "velocity_%s" % (fname)
fu = HDF5File(mesh.mpi_comm(),os.path.join(data_dir, ufilename+".hdf5"), 'w')
fu.write(u, "/velocity")

# Save the final solution for visualisation
save_pvd("3D_velocity_%s" % (fname), u)

print "end of ice dynamics"
