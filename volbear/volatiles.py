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
Utility functions to get volatile abundances (ratios or masses) from the parametrisation
of mixing_line_param + VMF.
"""

from typing import List
import numpy as np
import pandas as pd

from volbear.utility import kwargs_are_float_or_equilong_arrays
from volbear.nat_const import (
    M_EARTH_KG,
    VMF_BSE,
    SOLAR_VOLATILE_MRATIOS,
    DIFF_BSE_SOL,
)


def get_elem_ratios(mixing_line_param: float | np.ndarray) -> pd.Series | pd.DataFrame:
    """
    Get mass ratios of He/H, C/H, N/H, and S/H from `mixing_line_param`.

    Parameters
    ----------
        mixing_line_param : float or np.ndarray
            `mixing_line_param` (ξ) from Seidler et al. 2025.
    Returns
    -------
        ratios : pd.Series or pd.DataFrame
            pandas Series (for scalar input) or pandas DataFrame (for array input)
            containing ratios He/H, C/H, N/H, and S/H.
    """
    # Compute the ratios
    if np.isscalar(mixing_line_param):
        ratios = SOLAR_VOLATILE_MRATIOS + mixing_line_param * DIFF_BSE_SOL
        # Modify the He/H ratio for mixing_line_param=1 to avoid issues with zeros
        ratios[0] = max(ratios[0], 1e-12)
        return pd.Series(
            data=ratios,
            index=["He/H", "C/H", "N/H", "S/H"],
        )
    else:
        # For array inputs, perform broadcasting
        mixing_line_param = np.asarray(
            mixing_line_param
        )  # Ensure mixing_line_param is an ndarray
        ratios = (
            SOLAR_VOLATILE_MRATIOS + mixing_line_param[:, np.newaxis] * DIFF_BSE_SOL
        )
        # Modify the He/H ratio for mixing_line_param=1 to avoid issues with zeros
        ratios[:, 0] = np.maximum(ratios[:, 0], 1e-12)
        return pd.DataFrame(
            data=ratios,
            columns=["He/H", "C/H", "N/H", "S/H"],
        )


def get_elem_massfracs(mixing_line_param: float | np.ndarray) -> pd.Series | np.ndarray:
    """
    Get mass fractions of elements H, He, C, N & S from `mixing_line_param`.

    Parameters
    ----------
        mixing_line_param : float or np.ndarray
            `mixing_line_param` from Seidler et al. 2025.
    Returns
    -------
        elem_massfracs : pd.Series or np.ndarray
            pandas Series or numpy array containing mass fractions of H, He, C, N, S.
    """
    elem_ratios: pd.Series | pd.DataFrame = get_elem_ratios(
        mixing_line_param=mixing_line_param
    )
    elem_names: List[str] = ["H", "He", "C", "N", "S"]

    # If `mixing_line_param` is a scalar, return a pandas Series
    if np.isscalar(mixing_line_param):
        s = pd.Series()
        for elem in elem_names:
            if elem == "H":
                s[elem] = 1.0
            else:
                s[elem] = elem_ratios[elem + "/H"]
        return s / s.sum()
    # else, return DataFrame
    else:
        df = pd.DataFrame()
        for elem in elem_names:
            if elem == "H":
                df[elem] = np.ones(len(elem_ratios))
            else:
                df[elem] = elem_ratios[elem + "/H"]
        return df.div(df.sum(axis=1), axis=0)


@kwargs_are_float_or_equilong_arrays()
def get_elem_masses(
    planet_mass_in_mearth: float | np.ndarray,
    log_vmf_factor: float,
    mixing_line_param: float | np.ndarray,
) -> pd.Series | pd.DataFrame:
    """
    Get masses (kg) of H, He, C, N, and S from `mixing_line_param`.

    Parameters
    ----------
        mixing_line_param : float or np.ndarray
            `mixing_line_param` from Seidler et al. 2025.
    Returns
    -------
        ratios : pd.Series or pd.DataFrame
            pandas Series or numpy array containing masses (in kg) of H, He, C, N, S.
    """

    # Volatile constraints, in total mass
    vmf: float | np.ndarray = 10**log_vmf_factor * VMF_BSE
    total_mass_of_volatiles: float | np.ndarray = (
        vmf * planet_mass_in_mearth * M_EARTH_KG
    )

    elem_mass_fracs: pd.Series | pd.DataFrame = get_elem_massfracs(mixing_line_param)

    elem_masses_df: pd.DataFrame = pd.DataFrame()
    for elem in ["H", "He", "C", "N", "S"]:
        _entry: float | pd.Series = total_mass_of_volatiles * elem_mass_fracs[elem]
        if isinstance(_entry, float):
            _entry = [_entry]
        elem_masses_df[elem] = _entry

    # all elements should have the same length; use it to determine return type
    if len(elem_masses_df) == 1:
        return elem_masses_df.iloc[0]
    else:
        return elem_masses_df
