#
# Copyright 2024-2026 Fabian L. Seidler
#
# This file is part of volbear.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
"""
Collection of natural (and non-natural) constants used throught this project.
"""

import numpy as np

M_EARTH_KG = 5.974e24
R_EARTH_M = 6.371e6

M_JUP_KG = 1.898e27
R_JUP_M = 6.9911e7

M_SUN_KG = 1.989e30
R_SUN_M = 6.95508e8

VMF_BSE: float = 0.00028485
""" Palme & O'Neill 2013 """

SOLAR_VOLATILE_MRATIOS: np.ndarray = np.array(
    [0.3925180671673516, 4.315e-3, 1.205e-3, 0.552e-3]
)
""" Lodders 2021, He/H, C/H, N/H and S/H mass ratio """

BSE_VOLATILE_MRATIOS: np.ndarray = np.array([0.0, 0.83, 0.0166667, 1.6666667])
""" BSE, Palme O'Neil 2013, mass ratio """

DIFF_BSE_SOL: np.ndarray = BSE_VOLATILE_MRATIOS - SOLAR_VOLATILE_MRATIOS

BSE_WT_COMP: dict = {
    "SiO2": 45.4,
    "MgO": 36.77,
    "FeO": 8.10,
    "Al2O3": 4.49,
    "CaO": 3.65,
}
""" BSE, Palme O'Neil 2013 """

MASSFRAC_OF_OXY_IN_BSE: float = 0.4433
""" Mass fraction of oxygen in BSE, Palme & O'Neill 2013, Table 4 """
