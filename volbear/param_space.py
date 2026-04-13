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
A `ParameterSpace` object stores the values that describe a atmosphere-interior model, either a
single point (all entries floats), or combinations of floats and arrays.
"""

from typing import Union, Tuple, Optional, List, Self, Dict
from dataclasses import dataclass, asdict
import numpy as np
import pandas as pd
from molmass import Formula

# atmodeller
from atmodeller import (
    Planet,
)
from atmodeller.containers import FixedFugacityConstraint
from atmodeller.thermodata import IronWustiteBuffer

# muspell
from muspell import Melt

# volbear
from volbear.interior_structure_model import interior_structure
from volbear.volatiles import (
    get_elem_masses,
)
from volbear.nat_const import M_EARTH_KG, R_EARTH_M, BSE_WT_COMP
from volbear.utility import to_float, kwargs_are_float_or_equilong_arrays, ensure_tuple


@dataclass(frozen=True)
class ParameterSpace:
    """
    A `ParameterSpace` object stores the values that describe a atmosphere-interior model, either
    a single point (all entries floats), or combinations of floats and arrays.
    """

    h_mass_kg: Union[float, np.ndarray]
    he_mass_kg: Union[float, np.ndarray]
    c_mass_kg: Union[float, np.ndarray]
    n_mass_kg: Union[float, np.ndarray]
    s_mass_kg: Union[float, np.ndarray]
    delta_iw: Union[float, np.ndarray]
    temperature: Union[float, np.ndarray]
    planet_mass_mearth: Union[float, np.ndarray]
    planet_radius_rearth: Union[float, np.ndarray]
    core_mass_frac: Union[float, np.ndarray]
    mantle_melt_frac: Union[float, np.ndarray]
    melt_comp_wt: Dict[str, float]
    """
    Tuple containing the mineral gas species that `atmodeller` should use as fugacity constraints.

    NOTE: For each element, include only ONE species! Otherwise, `atmodeller` cannot build the 
    fugacity constraints as they are no longer unique.
    """

    def is_singular_point(self) -> bool:
        """
        Checks if this parameter space is a singular point, i.e., every field is float. Used by
        the `outgassing` routine that couples to `phaethon`, which does only work for a singular
        point and not with arrays.
        """
        for _, value in self.__dict__.items():
            if isinstance(value, (list, np.ndarray)):
                if len(value) > 1:
                    return False
        return True

    def get_df(self) -> pd.DataFrame:
        """
        Cast sample into a pandas DataFrame.
        """

        # ---- Single simulation ---- #
        if self.is_singular_point():
            as_dict = self.asdict()
            for key in "melt_comp_wt":
                as_dict.pop(key)
            return pd.DataFrame(pd.Series(as_dict)).T

        # ---- Batch simulation ---- #
        frame: pd.DataFrame = pd.DataFrame()
        for key, data in self.__dict__.items():
            if key not in ("melt_comp_wt"):
                if isinstance(data, float):
                    # skip for now, because we want the pandas frame to be build based on the
                    # arrays
                    frame[key] = np.nan  # np.array([data])
                else:
                    frame[key] = np.array(data)

        # second iteration, fill in whatever was a float before
        for key, data in self.__dict__.items():
            if key not in ("melt_comp_wt"):
                if isinstance(data, float):
                    # skip for now, because we want the pandas frame to be build based on the
                    # arrays
                    frame[key] = np.full(len(frame), data)
        return frame

    def get_param_names(self) -> List[str]:
        """
        List of names of parameters of this dataclass, e.g. 'he_mass_kg', 'delta_iw', etc.

        Returns
        -------
            param_names : List[str]
                List of param names.
        """
        return list(self.__dataclass_fields__.keys())

    def asdict(self) -> dict:
        """
        Returns itself as dictionary.
        """
        return asdict(self)

    def get_atmodeller_constraints(
        self,
        elems_with_fugconstr: Tuple[str],
        mineral_vapour_pressures: Optional[pd.DataFrame] = None,
        progbar: bool = False,
    ) -> Tuple[Planet, dict, dict]:
        """
        Construct atmodeller constraints from sample.

        Parameters
        ----------
            elems_with_fugconstr: Tuple[str]
                A list/tuple of elements or molecular species that have a fugacity constraint,
                excluding oxygen (O2).
            mineral_vapour_pressures : Optional[pd.DataFrame]
                A pandas DataFrame containing the vapour pressures of mineral gases (SiO, Mg, etc).
                Pressures must be in bar. Optional. Default is `None`, in which case it will be
                computed on the fly.
            progbar: bool
                Should a progressbar be shown when computing mineral vapour pressures? Default
                is `False`.
        Returns
        -------
            planet : atmodeller.Planet
                An atmodeller planet object.
            mass_constraints : dict
                Dictionary with the constraints on mass of different volatile species.
            fugacity_constraints : dict
                Dictionary with the constraints on fugacity. This typically is for O2 and the
                mineral gas vapour pressures. O2 is added automatically.
        """

        # load or calculate mineral vapour pressures
        _mineral_vapour_pressures: pd.DataFrame = pd.DataFrame()
        if mineral_vapour_pressures is not None:
            _mineral_vapour_pressures = mineral_vapour_pressures
        else:
            melt: Melt = Melt(buffer="IW")
            melt.set_chemistry(wt=self.melt_comp_wt)
            logp_mineral_gases, _, _, _, _ = melt.run_series(
                temperature=self.temperature, dlogfO2=self.delta_iw, progbar=progbar
            )
            _mineral_vapour_pressures = 10**logp_mineral_gases

        # init planet object
        planet = Planet(
            temperature=self.temperature,
            planet_mass=self.planet_mass_mearth * M_EARTH_KG,
            surface_radius=self.planet_radius_rearth * R_EARTH_M,
            mantle_melt_fraction=self.mantle_melt_frac,
            core_mass_fraction=self.core_mass_frac,
        )

        # mass constraints
        mass_constraints = {
            "H": self.h_mass_kg,
            "He": self.he_mass_kg,
            "C": self.c_mass_kg,
            "N": self.n_mass_kg,
            "S": self.s_mass_kg,
        }

        # fugacity constraints
        fugacity_constraints: list = {"O2_g": IronWustiteBuffer(self.delta_iw)}
        for species_name in ensure_tuple(elems_with_fugconstr):
            if species_name == "O2":
                continue
            atm_name: str = (
                Formula(species_name).formula + "_g"
            )  # we assume that its a gas
            if _mineral_vapour_pressures.shape[0] > 1:
                fugacity = to_float(_mineral_vapour_pressures[species_name].to_numpy())
            else:
                fugacity = to_float(_mineral_vapour_pressures[species_name].iloc[0])
            fugacity_constraints[atm_name] = FixedFugacityConstraint(fugacity)

        return planet, mass_constraints, fugacity_constraints


@dataclass(frozen=True)
class LinearMixingParameterSpace(ParameterSpace):
    """
    A param space whichs volatile species are mixed on a linear mixing line between a fully
    VIBSE-like or SOLAR-like mixing line, see Seidler et al. 2026 (TODO: insert DOI).
    """

    mixing_line_param: Union[float, np.ndarray]
    log_vmf_factor: Union[float, np.ndarray]

    # pylint: disable=too-many-arguments
    @classmethod
    @kwargs_are_float_or_equilong_arrays(exceptions=("cls", "melt_comp_wt"))
    def new(
        cls,
        *,
        planet_mass_mearth: Union[float, int, np.ndarray],
        cmf: Union[float, int, np.ndarray],
        mantle_meltfrac: Union[float, int, np.ndarray],
        temperature: Union[float, int, np.ndarray],
        log_vmf_factor: Union[float, int, np.ndarray],
        mixing_line_param: Union[float, int, np.ndarray],
        delta_iw: Union[float, int, np.ndarray],
        melt_comp_wt: Optional[Union[float, int, np.ndarray]] = None,
    ) -> Self:
        """
        Generate a uniformly sampled parameter space to explore for atmosphere-interior
        interaction at the volbear.

        Parameters
        ----------
            planet_mass_mearth : float | int | np.ndarray
                Planet mass in Earth masses.
                Can be a scalar or array.
            cmf : float | int | np.ndarray
                Core mass fraction. Values should be in the range [0, 1].
                Can be a scalar or array.
            mantle_meltfrac : float | int | np.ndarray
                Mantle melt fraction. Values should be in the range [0, 1].
                Can be a scalar or array.
            temperature : float | int | np.ndarray
                Temperature at the magma-ocean-atmosphere interface, in kelvin.
                Can be a scalar or array.
            log_vmf_factor : float | int | np.ndarray
                Logarithmic vapor mass-fraction scaling factor used to compute volatile
                abundances in the atmosphere (base-10 log). Can be a scalar or array.
            mixing_line_param : float | int | np.ndarray
                Parameter controlling the bulk composition along a compositional mixing
                line (unitless). Interpreted as the mixing fraction between endmember
                compositions VIBSE and SOLAR; typically in the range [0, 1].
                Can be a scalar or array.
            delta_iw : float | int | np.ndarray
                Oxygen fugacity offset relative to the iron-wüstite (IW) buffer, in log10
                units (dex). Can be positive or negative; may be a scalar or array.
            melt_comp_wt : float | int | np.ndarray, optional
                Mantle melt composition given as weight fractions for the tracked oxide or
                element species. Defaults to BSE_WT_COMP (bulk silicate Earth weight
                composition). Accepts a single composition vector or an array of vectors
                matching other input dimensions.

        Returns
        -------
        Sample
            A Sample object (or equivalent data structure) containing `n_samples` sets of
            parameters sampled uniformly from the given limits for each planetary property.
        """

        if melt_comp_wt is None:
            melt_comp_wt = BSE_WT_COMP

        # Compute the volatile masses
        volatile_masses: pd.Series = get_elem_masses(
            planet_mass_in_mearth=to_float(planet_mass_mearth),
            log_vmf_factor=log_vmf_factor,
            mixing_line_param=to_float(mixing_line_param),
        )

        # Finish
        return LinearMixingParameterSpace(
            planet_mass_mearth=to_float(planet_mass_mearth),
            planet_radius_rearth=to_float(
                interior_structure(planet_mass_in_mearth=planet_mass_mearth, cmf=cmf)
            ),
            core_mass_frac=to_float(cmf),
            mantle_melt_frac=to_float(mantle_meltfrac),
            temperature=to_float(temperature),
            delta_iw=to_float(delta_iw),
            h_mass_kg=to_float(volatile_masses["H"]),
            he_mass_kg=to_float(volatile_masses["He"]),
            c_mass_kg=to_float(volatile_masses["C"]),
            n_mass_kg=to_float(volatile_masses["N"]),
            s_mass_kg=to_float(volatile_masses["S"]),
            melt_comp_wt=melt_comp_wt,
            mixing_line_param=to_float(mixing_line_param),
            log_vmf_factor=log_vmf_factor,
        )
