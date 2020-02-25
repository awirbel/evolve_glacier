"""
    This is the utilities file with all functions required from the model to solve the free surface evolution using FEM
    and the FEniCS software framework https://fenicsproject.org/

    Evolve
    Copyright (C) 2019  Anna Wirbel, Alexander Jarosch

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.



"""

from __future__ import division, print_function
import numpy as np
from dolfin import *
import os
from cfg_datadir import (data_dir, ab_params, acc_params, mb_height,
                         yins, siny, cmaxF)


def ice_mass(s_elev, b_elev, Qice):

    """Compute ice thickness function as surface - bed constraint!"""

    hice = Function(Qice)
    s_vec = s_elev.vector()[:]
    b_vec = b_elev.vector()[:]
    hice_vec = s_vec - b_vec
    hice.vector()[:] = hice_vec

    return hice


def ice_thickness(Qice, dofsQice_x, dofsQice, s_elev, b_elev):

    """ Compute ice thickness second option :) """

    hice_new_computed = Function(Qice)
    for dof, dofQice_x in zip(dofsQice, dofsQice_x):
        hice_new_computed.vector()[dof] = s_elev(dofQice_x[0], dofQice_x[1]) - b_elev(dofQice_x[0], dofQice_x[1])

    return hice_new_computed


def mb_change(mb_timestep, h_actual, dt, Qmb):

    """Compute expected mass change as a function of:
     mb rate, vertical velocity, ice thickness (due to constraint) and dt"""

    mb_mass = Function(Qmb)
    mb0_mass = Function(Qmb)
    hice_vec = h_actual.vector()[:]
    mb_tvec = mb_timestep.vector()[:] * dt
    mb_vec = hice_vec + mb_tvec
    mb0_vec = mb_vec.copy()
    mb0_vec[mb_vec < 0.0] = 0.0
    mb_mass.vector()[:] = mb_vec
    mb0_mass.vector()[:] = mb0_vec

    return mb_mass, mb0_mass


def check_change(s_elev, b_elev, s_actual, Qice, dofsQice_x, dofsQice):

    """ Compute change in ice thickness in percentage of previous ice thickness """

    hice_change = Function(Qice)
    for dof, dofQice_x in zip(dofsQice, dofsQice_x):
        hice_change.vector()[dof] = (100. / (s_elev(dofQice_x[0], dofQice_x[1]) - b_elev(dofQice_x[0], dofQice_x[1]))) *  (s_actual(dofQice_x[0], dofQice_x[1]) - b_elev(dofQice_x[0], dofQice_x[1]))

    return hice_change


def set_mb(s_elev, Qmb, dofsQmb, dofsQmb_x):

    """ Set mass balance on computational 2D mesh """

    # mb functions as a polynom
    mb_func_acc = np.poly1d(acc_params)
    mb_func_ab = np.poly1d(ab_params)

    mb1_func = Function(Qmb)
    for dof, dofQmb_x in zip(dofsQmb, dofsQmb_x):
        if (s_elev(dofQmb_x[0], dofQmb_x[1]) >= mb_height):
            mb1_func.vector()[dof] = mb_func_acc(s_elev(dofQmb_x[0], dofQmb_x[1])) * siny
        else:
            mb1_func.vector()[dof] = mb_func_ab(s_elev(dofQmb_x[0], dofQmb_x[1])) * siny

    return mb1_func


def get_timestep(v_hor, v_ver, mb, meshC, cmaxF):

    """Compute time step with CFL condition"""

    # ***********TIME STEP*****************
    # Velocity components and mass balance rate as numpy arrays
    vx_max11 = np.asarray(v_hor.vector()[:])
    vx_max1 = vx_max11[::2]
    vx_max = np.amax(np.absolute(vx_max1))
    vy_max11 = np.asarray(v_hor.vector()[:])
    vy_max1 = vy_max11[1::2]
    vy_max = np.amax(np.absolute(vy_max1))
    vz_max = np.amax(np.absolute(np.asarray(v_ver.vector()[:])))
    mb_max = np.amax(np.absolute(np.asarray(mb.vector()[:])))

    # Compute max of vels and mb rate
    velmax = max(vx_max, vy_max, vz_max, mb_max)
    # Computation time step following CFL condition
    dt_computed = ((cmaxF * meshC.hmin()) / (velmax))
    # ***********TIME STEP*****************

    return dt_computed, velmax


def save_pvd(pvdfilename, dataset, pvd_number=0):

    """ Save function to pvd file """

    data_dir1 = data_dir
    if pvd_number > 0:
        file_pvd = File(os.path.join(data_dir1, ("%s%04d.pvd" % (pvdfilename, pvd_number))))
        file_pvd << dataset
    else:
        file_pvd = File(os.path.join(data_dir1, ("%s.pvd" % (pvdfilename))))
        file_pvd << dataset



def set_velocity(velocity, s_elev, b_elev, Qv, Vv, dofsVv_x, dofsVv, bed_diff1):

    """Set velocity onto 2D computational mesh, within artificial ice layer set velocity to zero"""

    # Velocity on 2D computational mesh, has to be separated into vertical and horizontal components
    # Define velocity on 2D computational mesh x and y component
    v_new1 = Function(Vv)
    count = 0

    # Set 3D velocity onto dofs of 2D mesh, which is surface of 3D mesh
    for dof, dof_x in zip(dofsVv, dofsVv_x):
        if count % 2 == 0:
            if s_elev(dof_x[0], dof_x[1]) - b_elev(dof_x[0], dof_x[1]) > (bed_diff1 + 1.e-6):
                v_new1.vector()[dof] = velocity(dof_x[0], dof_x[1], s_elev(dof_x[0], dof_x[1]))[0]
            else:
                v_new1.vector()[dof] = 1.e-50

        elif count % 2 != 0:
            if s_elev(dof_x[0], dof_x[1]) - b_elev(dof_x[0], dof_x[1]) > (bed_diff1 + 1.e-6):
                v_new1.vector()[dof] = velocity(dof_x[0], dof_x[1], s_elev(dof_x[0], dof_x[1]))[1]
            else:
                v_new1.vector()[dof] = 1.e-50
        count = count + 1

    return v_new1


def set_verticalvelocity(velocity, s_elev, b_elev, Qmb, dofsQmb, dofsQmb_x, bed_diff1):

    """Set vertical velocity onto 2D computational mesh, within artificial ice layer set velocity to zero"""

    # Define velocity on 2D computational mesh z component
    v_vert1 = Function(Qmb)
    for dof, dofQmb_x in zip(dofsQmb, dofsQmb_x):
        if s_elev(dofQmb_x[0], dofQmb_x[1]) - b_elev(dofQmb_x[0], dofQmb_x[1]) > (bed_diff1 + 1.e-6):
            v_vert1.vector()[dof] = velocity(dofQmb_x[0], dofQmb_x[1], s_elev(dofQmb_x[0], dofQmb_x[1]))[2]

    return v_vert1


def subdomains_vels(Rmesh, velocity_h, s_elev, b_elev, surfgrad, mb_start, dt_ADV, tna, bed_diff1):

    """ Define subdomains that separate ice-covered and ice-free regions and include a velocity dependent buffer zone """

    #************************* Define subdomains **************************************
    # Mark all cells for given ice thickness and velocity threshold
    subdomains_Rnew1 = CellFunction("size_t", Rmesh)
    subdomains_Rnew1.set_all(0)
    cell_dofs1 = []
    grid_s = Rmesh.hmin() * 0.25
    flag_rcells = 0
    vels = np.zeros((Rmesh.num_cells()))
    countvel = 0
    for cell in cells(Rmesh):
        p = cell.midpoint()
        if ((s_elev(p) - b_elev(p)) > bed_diff1 + 0.5):
            if abs((np.sqrt(velocity_h(p)[0]**2 + velocity_h(p)[1]**2)) * dt_ADV) > grid_s:
                subdomains_Rnew1[cell] = 1
                vels[countvel] = abs((np.sqrt(velocity_h(p)[0]**2 + velocity_h(p)[1]**2)) * dt_ADV)
                flag_rcells = 1
                countvel = countvel + 1
    #save_pvd("rmesh_markersvel", subdomains_Rnew1, pvd_number=tna)

    # Compute all cell neighbours
    tdim = Rmesh.topology().dim()
    Rmesh.init(tdim - 1, tdim)
    cell_neighbors = {cell.index(): sum((filter(lambda ci: ci != cell.index(),
                                                facet.entities(tdim))
                                        for facet in facets(cell)), [])
                    for cell in cells(Rmesh)}


    # Compute how big the velocity dependent buffer zone should be
    if flag_rcells == 1:
        n_its = int(np.ceil(np.amax(vels) /grid_s))
        # print("glacier in advancing subdomain regime")
    else:
        n_its = 1
        # print("glacier in retreating subdomain regime")

    # Mark all neighbours of cells, repeat according to ice velocities
    marked_dofs = []
    for m in range(n_its+1):
        print("iteration: ", m, "of ", n_its)
        rcells11 = np.where(subdomains_Rnew1.array()==1)
        rcells1 = np.asarray(rcells11[0])

        subdomains_Rnew1 =  CellFunction("size_t", Rmesh)
        subdomains_Rnew1.set_all(0)
        for k in rcells1:
            for i in range(len(cell_neighbors[int(k)])):
                subdomains_Rnew1[int(cell_neighbors[int(k)][i])] = 1
    # save_pvd("rmesh_2markers", subdomains_Rnew1, pvd_number=tna)

    # Mark all cells where ice thickness is more than bed_diff + threshold
    # and include those from the velocity dependent buffer
    subdomains_Rnew3 = CellFunction("size_t", Rmesh)
    subdomains_Rnew3.set_all(0)
    for cell in cells(Rmesh):
        p = cell.midpoint()
        if s_elev(p) - b_elev(p) >= bed_diff1+ 0.1:
            subdomains_Rnew3[cell] = 1
        elif subdomains_Rnew1[cell] == 1:
            subdomains_Rnew3[cell] = 1
    # Save final subdomains for visualisation
    save_pvd("rmesh_markers", subdomains_Rnew3, pvd_number=tna)

    return subdomains_Rnew3


def solve_Netwon(F, u, bcs, dF, solv=False):

    """ Set up non-linear problem to be solved with Newton solver
        without any constraints first to provide initial guess,
        if solv=True also solve """

    # Define the solver parameters
    newton_solver_parameters = {"nonlinear_solver": "newton",
                            "newton_solver"     : { "linear_solver"   : "superlu_dist",
                                                "maximum_iterations": 50,
                                                "report": True,
                                                "error_on_nonconvergence": False,
                                                "relative_tolerance": 6.E-9,
                                                "absolute_tolerance": 6.e-9,
                                                }}

    # Set up the non-linear problem without constraints
    problem = NonlinearVariationalProblem(F, u, bcs, dF)
    solver = NonlinearVariationalSolver(problem)
    solver.parameters.update(newton_solver_parameters)
    if solv == True:
        (iter, converged) = solver.solve()

    return u


def solve_SNES(F, u, du, bcs, Smin, Smax, solv=False):


    """ Set up non-linear problem to be solved with PETSc's SNES solver
        and constraints, if solv=True also solve"""

    # Define the solver parameters
    snes_solver_parameters = {"nonlinear_solver": "snes",
                              "snes_solver": {"linear_solver": "superlu_dist",
                                              "method": "vinewtonrsls",
                                              "maximum_iterations": 50,
                                              "report": True,
                                              "error_on_nonconvergence": False,
                                              "relative_tolerance": 1.E-8,
                                              "absolute_tolerance": 1.e-8,
                                              }}

    # create the derivative
    dF = derivative(F, u, du)
    # Set up the non-linear problem
    problem = NonlinearVariationalProblem(F, u, bcs, dF)
    # Set the non-linear solver and parameters
    solver = NonlinearVariationalSolver(problem)
    solver.parameters.update(snes_solver_parameters)
    # info(solver.parameters, True)
    # Set also constraints
    problem.set_bounds(Smin, Smax)
    if solv == True:
        (iter, converged) = solver.solve()

    return u

def refine_mesh(sp, tna):


    """ Create attractor fields to refine mesh for 2D and 3D meshes """

    # Parameters to define attractor field, depend on number of points in initial geometry
    len_p = len(sp)
    m1  = 30

    # print("2D mesh write", "surface_2D_%04d.geo" % (tna))
    with open(os.path.join(data_dir, "surface_2D_%04d.geo" % (tna)), "a") as myfile:

        myfile.write("lc2=20;\n")
        for m in range(m1, m1 + len_p):
            myfile.write("Point(%d) = {%.3f, %.3f, %.3f, lc2};\n" % (m, sp[m - m1,0], sp[m - m1,1], 0.0))

        myfile.write("Field[1] = Attractor; \n")
        myfile.write("Field[1].NodesList = {%d" % (m1))
        for m in range(m1 + 1, m1 + len_p):
            myfile.write(",%d" % (m))
        myfile.write("};\n")
        myfile.write("Field[2] = Threshold;\nField[2].IField = 1;\nField[2].LcMin = %.2f;\n" % (5))
        myfile.write("Field[2].LcMax = lc2;\nField[2].DistMin = %.2f;\nField[2].DistMax = %.2f;\nField[3] = Min;\n" % (0.1, 400))
        myfile.write("Field[3].FieldsList ={2};\nBackground Field = 3;\n")


    with open(os.path.join(data_dir, "surface_shell_domain.geo"), "a") as myfile:

        myfile.write("Extrude {0, 0, layext} {\n")
        myfile.write("Surface{8}; Layers{layno};\n")
        myfile.write("}\n")

        myfile.write("lc2=20;\n")
        for m in range(m1, m1 + len_p):
            myfile.write("Point(%d) = {%.3f, %.3f, %.3f, lc2};\n" % (m, sp[m - m1,0], sp[m - m1,1], 0.0))

        myfile.write("Field[1] = Attractor; \n")
        myfile.write("Field[1].NodesList = {%d" % (m1))
        for m in range(m1 + 1, m1 + len_p):
            myfile.write(",%d" % (m))
        myfile.write("};\n")
        myfile.write("Field[2] = Threshold;\nField[2].IField = 1;\nField[2].LcMin = %.2f;\n" % (5))
        myfile.write("Field[2].LcMax = lc2;\nField[2].DistMin = %.2f;\nField[2].DistMax = %.2f;\nField[3] = Min;\n" % (0.1, 400))
        myfile.write("Field[3].FieldsList ={2};\nBackground Field = 3;\n")

        myfile.write("Physical Surface(101) = {8};\n")
        myfile.write("Physical Surface(102) = {30};\n")
        myfile.write("Physical Surface(103) = {25};\n")
        myfile.write("Physical Surface(104) = {29};\n")
        myfile.write("Physical Surface(105) = {17};\n")
        myfile.write("Physical Surface(106) = {21};\n")


def write_surface(s_elev, b_elev):

    """  Write new surface elevation into vol_domain.stl  """

    # Open the input and output files
    file_in = 'surface_shell_domain.stl'
    file_out = 'vol_domain.stl'
    fp = open(os.path.join(data_dir, file_in), 'r')
    fp_out = open(os.path.join(data_dir, file_out), 'w')

    for line in fp:
        linesplit = line.split()
        if linesplit[0] == 'vertex':
            if linesplit[3] == '0.000':
                xcoord = float(linesplit[1])
                ycoord = float(linesplit[2])

                # get the z coord
                zcoord = s_elev(xcoord, ycoord)
                fp_out.write("    vertex %.03f %.03f %.03f\n" % (xcoord, ycoord,
                                                           zcoord))

            elif linesplit[3] == '-1000.000':
                xcoord = float(linesplit[1])
                ycoord = float(linesplit[2])

                # get the z coord
                zcoord = b_elev(xcoord, ycoord)
                fp_out.write("    vertex %.03f %.03f %.03f\n" % (xcoord, ycoord,
                                                                 zcoord))
        else:
            fp_out.write(line)

    fp_out.close()
    fp.close()


def write_surface_v2(s_elev, b_elev):

    """ create new meshfile directly, based on block mesh """

    # Open the input and output files
    file_in = 'volmesh_initblock.msh'
    file_out = 'volmesh_block.msh'
    fp = open(os.path.join(data_dir, file_in), 'r')
    fp_out = open(os.path.join(data_dir, file_out), 'w')

    for line in fp:
        linesplit = line.split()
        if len(linesplit) == 4:
            xcoord = float(linesplit[1])
            ycoord = float(linesplit[2])
            z_mult = float(linesplit[3]) / 1000.


            s_point = s_elev(xcoord, ycoord)
            b_point = b_elev(xcoord, ycoord)
            h_point = s_point - b_point
            z_coord = b_point + (h_point * z_mult)
            fp_out.write(linesplit[0]+" %.03f %.03f %.03f\n" %
                         (xcoord, ycoord, z_coord))

        else:
            fp_out.write(line)


    fp_out.close()
    fp.close()
