"""
    configuration file

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

import pickle
import os
from cfg_datadir import data_dir
import numpy as np

"""

    Save parameters that have to be passed between functions in a pickle as dictionary

    :param
    fname               name prefix of gmsh mesh files
    mass_initial        mass of initial concentration at t=0
    mass_previous       mass of previous time step
    velmax              max velocity in refined area
    dt_total_old        total time step in previous computation step
    dt_ADV              surface evolution time step
    nx1, nx2, ny1, ny2  x/y corner coordinates of computational domain
    t1                  number of actual surface evolution time step
    bed_diff1           thickness of artificial ice layer
    tnatext             name for solution files
    tna                 name suffix
    dt_overall          overall time step

"""


anim_name_def = "anim01.pvd"


run_params = dict(
    fname = "glacier_0000",
    mass_initial = 0.0,
    mass_previous = 0.0,
    velmax = 0.,
    dt_total_old = 0.0,
    dt_ADV = 1.0 * 3600. *24. *365.,
    nx1=599040.00,
    ny1=6741040.0,
    nx2=603920.0,
    ny2=6746960.0,
    t1 = 0,
    bed_diff1 = 1.0,
    t_all = 0,
    tnatext = "A",
    tna = 01,
    dt_overall = 0.0
    )


with open(os.path.join(data_dir, "data.pkl"), "wb" ) as f:
    pickle.dump(run_params, f)
