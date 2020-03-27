"""
    Data directory and constant parameters

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

import numpy as np

# general parameters
data_dir = '/workingdirectory/data/'
# Courant number
cmaxF = 0.1
# year in seconds, seconds in year
yins = (3600. * 365. * 24.)
siny = 1. / (3600. * 365. * 24.)
# mass balance parameters
# here accumulation and ablation, mb_height defines where to use which gradient
y = [0.0, 5]
x = [0.0, 100]
coefficients = np.polyfit(x, y, 1)
acc_params = coefficients
ab_params = coefficients
mb_height = 10.0
# TIME DISCRETIZATION
TS = 'CN'
