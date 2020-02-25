"""
    This model solves the free surface evolution with a 2nd order Runge Kutta scheme using FEM
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
                   write_surface_v2, solve_SNES, solve_Netwon, ice_thickness,
                    set_mb, check_change, set_verticalvelocity, subdomains_vels)
from utils_RK import (readmeshvelocity)

#********************Read input data******************
# Load parameters
with open(os.path.join(data_dir, "data.pkl"), "r") as f:
    params = pickle.load(f)

# Define parameters from configuration file
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
nx1 = params['nx1']
nx2 = params['nx2']
ny1 = params['ny1']
ny2 = params['ny2']

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


# print "Read surface and bed from: ", fname_old
# Read functions from previous time step to initialise the problem
s_new = Function(Q3)
s_new_bed = Function(Q3)
s_new_bed_ini = Function(Q3)
sfilename = "solution_%s" % (fname_old)
sin1 = HDF5File(Rmesh.mpi_comm(), os.path.join(data_dir,sfilename+".hdf5"), 'r')
sin1.read(s_new, "/elevation")
sbfilename = "solutionbed_%s" % (fname_old)
sinb = HDF5File(Rmesh.mpi_comm(), os.path.join(data_dir,sbfilename+".hdf5"), 'r')
sinb.read(s_new_bed_ini, "/bed")
sb2filename = "solutionconst_%s" % (fname_old)
sinb2 = HDF5File(Rmesh.mpi_comm(), os.path.join(data_dir,sb2filename+".hdf5"), 'r')
sinb2.read(s_new_bed, "/const")
# surf_ini is the initial surface elevation before moving
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
v_new = set_velocity(velocity, s_new, s_new_bed_ini, Q3, V3, dofV3_x, dofsV3, bed_diff1)
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

# Define integration measures: this has to be defined as we have multiple meshes loaded in this function
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
print 'TIME DISCRETIZATION: Trapezoidal Rule, Runge Kutta 2nd order'
# Set time step with CFL condition
# Define time step according to CFL
[dt, velmax] = get_timestep(v_new, v_vert, mb_func, Rmesh, cmaxF)
dt_form = Constant(dt)
# Parameters for time stepping loop
t_ref = 1
T_ref = int(dt_ADV / dt) + 2

# Set constraint for surface elevation
# upper constraint
constraint_u = Expression(("Smax"), Smax=10000.0+DOLFIN_EPS, degree=1)
Smax = interpolate(constraint_u, Q3)
# lower constraint is the bed elevation + bed_diff1
Smin = s_new_bed

# Save inequality constraints to File for visualisation
# save_pvd("Smax", Smax)
# save_pvd("Smin", Smin)

# Define parameters for time stepping loop
dt_final = dt_ADV
flag_continue = 0
flag_stop = 0
dt_total_old = 0.0
dt_total = 0.0
# Time-stepping loop
while flag_stop < 0.5:

    print "Evolve time step: ", t1, " in computation: ", t_ref, "of", T_ref
    if dt >= dt_ADV:
        flag_continue = 1
        dt_new = dt_ADV
        dt_computed = dt_ADV

    # MB-Start------------Compute mass balance as a function of current surface elevation
    mb_func = set_mb(u0, Q3, dofsQ3, dofsQ3_x)
    # Set to Function
    f0.assign(mb_func)
    save_pvd("rmesh_ft0", mb_func, pvd_number=t_ref)
    # MB-End---------------------------------------------------------------------------

    [dt_computed, velmax] = get_timestep(v_new, v_vert, mb_func, Rmesh, cmaxF)
    if (dt_total + dt_computed) > dt_final:
        dt_new = dt_final - dt_total
        dt_form.assign(dt_new)
        flag_continue = 1
        # print 'dt_new is set', dt_new
        dt_computed = dt_new
    else:
        # print 'dt_computed is set', dt_computed
        dt_form.assign(dt_computed)
    # ***********TIME STEP*****************

    # HICE-Start--------------Compute ice thickness second option--------------------------
    hice_new_computed = ice_thickness(Q3, dofsQ3_x, dofsQ3, u0, s_new_bed)
    # HICE-End------------------------------------------------------------------------------

    # Start++++++++++++++++++++++++++++Perform computations++++++++++++++++++++++
    # ---------Variational form for explicit Euler---------
    # Test and trial functions
    du, v = TrialFunction(Q3), TestFunction(Q3)
    u = Function(Q3)
    # Residual
    r = u - u0 + dt_form * (dot(v_new, grad(u0)) - v_vert - f0)

    # equation solved on subdomains for glacier and outside
    F1 = v * (u - u0) * dx(1) + v * (u - u0) * dx(0) + dt_form * (
                v * (v_new[0] * grad(u0)[0] + v_new[1] * grad(u0)[1]) * dx(1) - v_vert * v * dx(1) - f0 * v * dx(
            1) - f0 * v * dx(0))

    # Add SUPG stabilisation terms
    vnorm = sqrt(dot(v_new, v_new))
    F11 = F1 + (h / (2.0 * vnorm)) * dot(v_new, grad(v)) * r * dx(1)

    # Add SOLD stabilisation terms following modified Burman approach in John and Knobloch (2008)
    gradunorm = sqrt(dot(grad(u0), grad(u0)))
    rnorm = sqrt(dot(r, r))
    vnorm = sqrt(dot(v_new, v_new))
    tau = h / (2.0 *vnorm)
    C = ((tau * vnorm * rnorm) / gradunorm) * ((vnorm * gradunorm) / (vnorm * gradunorm + rnorm))
    F = F11 + C * dot(grad(u0), grad(v))*dx(1)

    # create the derivative
    dF = derivative(F, u, du)

    # solve for u with Newton solver but not constrained to receive initial guess
    u = solve_Netwon(F, u, bcs, dF, solv=True)

    # solve for u with SNES solver for first step in Runge Kutta 2nd order
    u = solve_SNES(F, u, du, bcs, Smin, Smax, solv=True)
    #save_pvd("SolutionT0_%s_%04d_new" % (tnatext, tna), u, pvd_number=t_ref)

    # Derive new velocity field with this solution_
    #----------------- Write result to stl to create new 3D geometry
    #print "RK step 1: write surface"
    write_surface_v2(u, s_new_bed_ini)
    # Create new 3D mesh
    # print "RK step 1: make new mesh"
    os.system('dolfin-convert %s/volmesh_block.msh %s/%s_vel.xml' % (data_dir, data_dir, fname))
    # print "RK step 1: mesh created"
    # Create new name
    fname_vel = "%s_vel" % (fname)
    order = 1
    # Save 3D mesh for intermediate time step
    readmeshvelocity(fname_vel, nx1, nx2, ny1, ny2, order)
    # print "RK step 1: read 3D velocity mesh , fname_vel", fname_vel, fname
    # print "RK step 1: compute velocity_RK"
    # Compute 3D velocity field for intermediate time step
    os.system('mpirun -np 4 python ice_dynamics3D_RK.py')
    mesh_vel = Mesh()
    mainname_vel = "%s" % (fname_vel)
    hdfin1_vel = HDF5File(mesh_vel.mpi_comm(), os.path.join(data_dir, mainname_vel+".hdf5"), "r")
    hdfin1_vel.read(mesh_vel, "/mesh", False)
    hdfin1_vel.close()
    # Define mesh coordinates and FunctionSpaces on initial mesh
    V = VectorFunctionSpace(mesh_vel, 'CG', 2)
    velocity_RK = Function(V)
    # Read 3D velocity from file
    velfilename_vel = "velocity_%s" % (fname_vel)
    hdfin_vel = HDF5File(mesh_vel.mpi_comm(), os.path.join(data_dir,velfilename_vel+".hdf5"), 'r')
    hdfin_vel.read(velocity_RK, "/velocity")
    hdfin_vel.close()
    # print "RK step 1: velocities computed, start solving"
    # Set mass balance rate
    mb_func2 = set_mb(u, Q3, dofsQ3, dofsQ3_x)
    f02 = Function(Q3)
    f02.assign(mb_func2)
    # Set velocity
    v_new2 = set_velocity(velocity_RK, u, s_new_bed_ini, Q3, V3, dofV3_x, dofsV3, bed_diff1)
    v_vert2 = set_verticalvelocity(velocity_RK, u, s_new_bed_ini,  Q3, dofsQ3, dofsQ3_x, bed_diff1)

    # ---------Variational form for 2nd order---------
    # Test and trial functions
    du2, v2 = TrialFunction(Q3), TestFunction(Q3)
    u2 = Function(Q3)
    # Residual
    r2 = (u2 - u0) + (0.5 * dt_form) * ((dot(v_new2, grad(u)) - v_vert2 - f02) + (dot(v_new, grad(u0)) - v_vert - f0))

    # equation solved on subdomains for glacier and outside
    F12 = v2 * (u2 - u0) * dx(1) + v2 * (u2 - u0) * dx(0) + (0.5 * dt_form) * ((v2 * (v_new2[0] * grad(u)[0] + v_new2[1] * grad(u)[1]) * dx(1) - v_vert2 * v2 * dx(1) - f02 * v2 * dx(1) - f02 * v2 * dx(0)
          ) + (v2 * (v_new[0] * grad(u0)[0] + v_new[1] * grad(u0)[1]) * dx(1) - v_vert * v2 * dx(1) - f0 * v2 * dx(1) - f0 * v2 * dx(0)))

    # Add SUPG stabilisation terms
    vnorm2 = sqrt(dot(v_new2, v_new2))
    F22 = F12 + 0.5 * (((h / (2.0 * vnorm2)) * dot(v_new2, grad(v2)) * r2 * dx(1)) + ((h / (2.0 * vnorm)) * dot(v_new, grad(v2)) * r * dx(1)))

    # Add SOLD stabilisation terms following modified Burman approach in the SOLD 2 paper by John
    gradunorm2 = sqrt(dot(grad(u), grad(u)))
    rnorm2 = sqrt(dot(r2, r2))
    vnorm2 = sqrt(dot(v_new2, v_new2))
    tau2 = h / (2.0 *vnorm2)
    C2 = ((tau2 * vnorm2 * rnorm2) / gradunorm2) * ((vnorm2 * gradunorm2) / (vnorm2 * gradunorm2 + rnorm2))
    F2 = F22 + 0.5 * ((C2 * dot(grad(u), grad(v2))*dx(1) + (C * dot(grad(u0), grad(v2))*dx(1))))

    u2 = solve_SNES(F2, u2, du2, bcs, Smin, Smax, solv=True)
    u2.rename("u", "u")
    # Save solution to file
    save_pvd("Solution_%s_%04d_new" % (tnatext, tna), u2, pvd_number=t_ref)

    # MBC-Start----------------Compute mass change due to MB-----------------------------
    # Compute actual mass change as diff of hice_previous and hice_new
    hice_new = ice_mass(u2, s_new_bed, Q3)
    hice_previous = ice_mass(u0, s_new_bed, Q3)
    hice_new_mass = assemble(hice_new * dx)
    hice_previous_mass = assemble(hice_previous * dx)
    print "This is the ice mass at time", t_ref, hice_new_mass, "at previous time", hice_previous_mass
    # Actual ice mass change within computational time step
    mb_change_solution = assemble(hice_previous * dx) - assemble(hice_new * dx)
    # MBC-End----------------------------------------------------------------------------

    # Compute ice thickness change
    hice_change = check_change(s_new, s_new_bed_ini, u2, Q3, dofsQ3_x, dofsQ3)
    # save_pvd("rmesh_hicechange", hice_change, pvd_number=tna)

    # Compute total time step
    dt_total = dt_computed + dt_total

    # Copy solution from previous interval to provide initial condition
    u0.assign(u2)
    v_new = v_new2
    v_vert = v_vert2

    t_ref = t_ref + 1
    print "Total time step is: ", dt_total
    print int(dt_total), int(dt_final)

    # Check if dt_ADV exceeded within next computation, if TRUE change dt to fill up dt_ADV
    if int(dt_total) == int(dt_final):
        print "stop iterations computations done"
        flag_stop = 1.0


#----------------- Write result to stl to create new 3D geometry START---------------------
write_surface_v2(u2, s_new_bed_ini)
#----------------- Write result to stl create new 3D geometry START---------------------


# Compute ice thickness
hice_solution = ice_thickness(Q3, dofsQ3_x, dofsQ3, u, s_new_bed)
save_pvd("rmesh_hicesolution", hice_solution, pvd_number=tna)

# Save the final functions for the next time step
# print "Save surface and bed from: ", fname
sfilename = "solution_%s" % (fname)
f1u = HDF5File(Rmesh.mpi_comm(),os.path.join(data_dir, sfilename+".hdf5"), 'w')
f1u.write(u2, "/elevation")
sbfilename = "solutionbed_%s" % (fname)
f2u = HDF5File(Rmesh.mpi_comm(),os.path.join(data_dir, sbfilename+".hdf5"), 'w')
f2u.write(s_new_bed_ini, "/bed")
sbifilename = "solutionconst_%s" % (fname)
f3u = HDF5File(Rmesh.mpi_comm(),os.path.join(data_dir, sbifilename+".hdf5"), 'w')
f3u.write(s_new_bed, "/const")


# Update parameters and save to file
params['t1'] = t1 + 1
params['tna'] = tna + 1
params['dt_total_old'] = dt_total
params['fname'] = 'glacier_%04d' % (tna)
params['fname_old'] = fname
params['dt_overall'] = dt_overall

with open(os.path.join(data_dir, 'data.pkl'), "w") as f:
    pickle.dump(params, f)

print 'surface move this is done now'
