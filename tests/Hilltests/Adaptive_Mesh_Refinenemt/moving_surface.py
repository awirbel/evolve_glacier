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
                    set_verticalvelocity,write_surface, solve_SNES,
                     solve_Netwon, ice_thickness, refine_mesh,
                    set_mb, check_change, subdomains_vels)


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



#+++++++++++++++++START OF SURFACE EVOLUTION+++++++++++++++++
# Read old 2D computational mesh that covers entire domain
fname_surfmesh = "surface_2D"
Rmesh_old = Mesh(os.path.join(data_dir, "%s_%04d.xml" % (fname_surfmesh, t1-1)))
# Define FunctionsSpaces on this mesh
Q3_old = FunctionSpace(Rmesh_old, 'CG', 1)
# print "Read surface and bed from: ", fname_old
s_new_old = Function(Q3_old)
s_new_bed_old = Function(Q3_old)
s_new_bed_ini_old = Function(Q3_old)
sfilename = "solution_%s" % (fname_old)
sin1 = HDF5File(Rmesh_old.mpi_comm(), os.path.join(data_dir,sfilename+".hdf5"), 'r')
sin1.read(s_new_old, "/elevation")
sbfilename = "solutionbed_%s" % (fname_old)
sinb = HDF5File(Rmesh_old.mpi_comm(), os.path.join(data_dir,sbfilename+".hdf5"), 'r')
sinb.read(s_new_bed_ini_old, "/bed")
sb2filename = "solutionconst_%s" % (fname_old)
sinb2 = HDF5File(Rmesh_old.mpi_comm(), os.path.join(data_dir,sb2filename+".hdf5"), 'r')
sinb2.read(s_new_bed_old, "/const")

# Read 2D computational mesh that covers entire domain
D = mesh.topology().dim()
mesh.init(D-1,D)
fname_surfmesh = "surface_2D"
Rmesh = Mesh(os.path.join(data_dir, "%s_%04d.xml" % (fname_surfmesh, t1)))
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
s_new = interpolate(s_new_old, Q3)
s_new_bed = interpolate(s_new_bed_old, Q3)
s_new_bed_ini = interpolate(s_new_bed_ini_old, Q3)
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
v_new = set_velocity(velocity, s_new, s_new_bed_ini, Q3, V3, dofV3_x, dofsV3, bed_diff1)
v_vert = set_verticalvelocity(velocity, s_new, s_new_bed_ini,  Q3, dofsQ3, dofsQ3_x, bed_diff1)
# Compute exptected mass change due to mass balance rate field
[mb_start, mb0_start] = mb_change(mb_func, hice_start, dt_ADV, Q3)
mb_mass0 = hice_t0 - assemble(mb0_start * dx)

# Save fields for visualisation
# save_pvd("rmesh_hicestart", hice_start, pvd_number=tna)
save_pvd("rmesh_ft0", mb_func)
# save_pvd("rmesh_snew", s_new, pvd_number=tna)
save_pvd("rmesh_vnew", v_new, pvd_number=tna)
save_pvd("rmesh_vert", v_vert, pvd_number=tna)
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
F = F1 + (h / (2.0 * vnorm)) * dot(v_new, grad(v)) * r * dx(1)

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

grads = project(grad(u), DG)
# save_pvd("rmesh_grads", grads, pvd_number=tna)

# Compute ice thickness
hice_solution = ice_thickness(Q3, dofsQ3_x, dofsQ3, u, s_new_bed)
save_pvd("rmesh_hicesolution", hice_solution, pvd_number=tna)

sp = np.zeros((len(dofsQ3), 2))
sp_func = Function(Q3)
count = 0
for dof, dofQ3_x in zip(dofsQ3, dofsQ3_x):
    if hice_solution(dofQ3_x[0], dofQ3_x[1]) > 1.0:
        if np.sqrt(grads(dofQ3_x[0], dofQ3_x[1])[0]**2 + grads(dofQ3_x[0], dofQ3_x[1])[1]**2) > 1.25:
            if (nx1 + 100.0 < dofQ3_x[0] < nx2 - 100.0) and (ny1 + 100.0 < dofQ3_x[1] < ny2 - 100.0):
                sp[count, 0] = dofQ3_x[0]
                sp[count, 1] = dofQ3_x[1]
                sp_func.vector()[dof] = 1.0
                count = count + 1
                if (u(dofQ3_x[0], dofQ3_x[1]) - s_new_bed_ini(dofQ3_x[0], dofQ3_x[1])) < 1.0:
                    print "smaller "

sp = sp[0:count, :]
save_pvd("rmesh_spfunc", sp_func, pvd_number=tna)
print "sp", len(sp), count
refine_mesh(sp, t1+1)


# Save the final solution
print "Save surface and bed from: ", fname
sfilename = "solution_%s" % (fname)
f1u = HDF5File(Rmesh.mpi_comm(),os.path.join(data_dir, sfilename+".hdf5"), 'w')
f1u.write(u, "/elevation")
sbfilename = "solutionbed_%s" % (fname)
f2u = HDF5File(Rmesh.mpi_comm(),os.path.join(data_dir, sbfilename+".hdf5"), 'w')
f2u.write(s_new_bed_ini, "/bed")
sbifilename = "solutionconst_%s" % (fname)
f3u = HDF5File(Rmesh.mpi_comm(),os.path.join(data_dir, sbifilename+".hdf5"), 'w')
f3u.write(s_new_bed, "/const")
# Update parameters and save to file

params['dt_total_old'] = dt_total
params['fname_old'] = fname
params['dt_overall'] = dt_overall

with open(os.path.join(data_dir, 'data.pkl'), "w") as f:
    pickle.dump(params, f)

print 'surface move this is done now'
