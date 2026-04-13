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
Basic interior structure model to compute radius as function of mass and core mass fraction.
"""

from importlib import resources
import numpy as np
from scipy.interpolate import interp1d

# load precomputed mass-radius curves
_zeng_mr_earthlike_file = resources.files("volbear.data") / "mr_Earth_zeng16.txt"
_data = np.genfromtxt(_zeng_mr_earthlike_file)
_fit = interp1d(_data[:, 0], _data[:, 1])


def interior_structure(planet_mass_in_mearth: float, cmf: float) -> float:
    """
    Evaluates the interior structure model for the condensed part of the planet,
    i.e. a silicate mantle and an iron core with the massfraction cmf.

    Parameters
    ----------
        planet_mass_in_mearth : float
            Mass of the planet, in Earth masses.
        cmf : float
            Core mass fraction. Not actually used, but is a placeholder in case a future model
            has to use it.

    Returns
    -------
        planet_radius_in_rearth : float
            Radius of the planet, in Earth radii.
    """
    planet_radius_in_rearth: float = _fit(planet_mass_in_mearth)[()]

    return planet_radius_in_rearth
