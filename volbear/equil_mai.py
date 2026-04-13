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
Equilibrator for the magma-ocean-atmosphere interface (MAI).
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Union, Callable, Tuple, Optional, List, Dict, Self
import numpy as np
import pandas as pd

# atmodeller
from atmodeller import EquilibriumModel, Planet
from atmodeller.output import Output

# phaethon ecosystem
from muspell import Melt
from phaethon.gas_mixture import IdealGasMixture

# volbear
from volbear.param_space import ParameterSpace
from volbear.chemical_network import ChemicalNetwork


def _get_gas_csv(
    interior_atmosphere: EquilibriumModel,
    param_space: ParameterSpace,
    func: Callable,
) -> None:
    """
    Write a simulation output in a pandas data frame. The output
    is constructed from the `func` argument output.

    Parameters
    ----------
        func : Callable
            Function that takes `gas_name` as parameter (gas names are found in `species_names`,
            e.g. CO2, CO, ...), performs any calculation and returns the respective parameter.
    Returns
    -------
        df : pd.DataFrame
            pandas DataFrame with the result.
    """
    df = param_space.get_df()
    for gas_species in interior_atmosphere.species.get_species_names():
        if not gas_species.endswith("_cr"):
            _name: str = gas_species[:-2]
            df[_name] = func(gas_species)

    return df


class MagmaOceanAtmosphereEquilibrator:
    """
    Equilibrate an atmosphere with an underlying magma ocean.
    """

    interior_atmosphere: EquilibriumModel
    """
    Engine; solves for the melt-atmosphere equilibrium, using the species defined in `species`.
    """
    _chemical_network: ChemicalNetwork
    """ Contains `atmodeller` species and identity of mineral species fugacity constraints."""

    _atm_planet: Planet
    """ Atmodeller planet object """

    _mass_constraints: Tuple
    """ Constraints on the total mass of volatiles, i.e. H, He, C, N & S"""

    _fugacity_constraints: Tuple
    """ Constraints on fugacity of species, i.e. O2 and the mineral gases (SiO, Mg, ...)"""

    output: Output
    """ Full atmodeller output. """

    __param_space: ParameterSpace
    """ Parameter space to be explored (oxygen fugacity, volatile compostion & amount, etc.). """

    def __init__(self, chemical_network: ChemicalNetwork) -> None:
        """
        Inits the runner. Only depends on species, as all chemistry is defined there.

        Parameters
        ---------
            chemical_network : ChemicalNetwork
                A chemical network object which encapsulates solubility laws.
        """
        self._chemical_network = chemical_network

        # filled later
        self.interior_atmosphere = None
        self.__param_space = None

    @staticmethod
    def _get_mineral_vapour_pressures(param_space: ParameterSpace) -> pd.DataFrame:
        """
        Computes the mineral vapor pressures based on the chemical composition, temperature and
        redox state of the melt.

        Parameters
        ---------
            sample : ParameterSpace
                A sample object that contains the parameters required for the calculation,
                including the melt composition (`melt_comp_wt`), temperature (`temperature`),
                and oxygen fugacity difference (`delta_iw`).

        Returns
        -------
            mineral_vapour_pressures : pd.DataFrame:
                A DataFrame containing the mineral vapor pressures, which are the antilogs of the
                calculated log vapor pressures (in Pa). Each row corresponds to a mineral species.
        """
        melt: Melt = Melt(buffer="IW")
        melt.set_chemistry(wt=param_space.melt_comp_wt)
        logp_mineral_gases, _, _, _, _ = melt.run_series(
            temperature=param_space.temperature,
            dlogfO2=param_space.delta_iw,
            progbar=False,
        )

        return 10**logp_mineral_gases

    def equil_melt_and_atmosphere(
        self,
        param_space: ParameterSpace,
        mineral_vapour: Optional[Union[os.PathLike, pd.DataFrame]] = None,
    ) -> Output:
        """
        Equilibrate the melt-atmosphere chemistry in all points defined by the paramter space.

        Parameters
        ----------
            sample : ParameterSpace
                The ParameterSpace on which to operate on.
            mineral_vapour : Optional[Union[os.PathLike, pd.DataFrame]]
                An object containing mineral vapour pressures, in bar. If `None` provided, compute
                them based on the information in `sample`.
        """

        # store sample for later
        self.__param_space = param_space

        # Check if mineral vapour is None, and compute the vapour pressures if so
        if mineral_vapour is None:
            mineral_vapour_pressures: pd.DataFrame = self._get_mineral_vapour_pressures(
                param_space=param_space
            )

        # If mineral vapour is already a DataFrame, use it directly
        elif isinstance(mineral_vapour, pd.DataFrame):
            mineral_vapour_pressures: pd.DataFrame = mineral_vapour

        # If mineral vapour is a string (path), but the file doesn't exist, compute and save
        # pressures to file. Generate parent directory if it does not exist yet.
        elif not os.path.isfile(mineral_vapour) and isinstance(mineral_vapour, str):
            if not mineral_vapour.endswith(".csv"):
                # Raise an error if the file doesn't have a .csv extension
                raise FileNotFoundError(
                    "Mineral vapour file should have a `.csv` ending."
                )

            # Compute mineral vapour pressures
            mineral_vapour_pressures: pd.DataFrame = self._get_mineral_vapour_pressures(
                param_space=param_space
            )

            # Ensure the directory for the file exists
            _filepath = Path(mineral_vapour)
            _filepath.parent.mkdir(parents=True, exist_ok=True)

            # Save the computed pressures to the CSV file
            mineral_vapour_pressures.to_csv(mineral_vapour)

        # If the file exists, read the mineral vapour pressures from the CSV file
        elif os.path.isfile(mineral_vapour):
            mineral_vapour_pressures: pd.DataFrame = pd.read_csv(
                mineral_vapour, index_col=0
            )

        # Raise an error if none of the above conditions are met
        else:
            err_msg: str = (
                "`mineral_vapour` is not a valid DataFrame,"
                + " and does not point to a valid file"
            )
            raise TypeError(err_msg)

        # get atmodeller constraints
        self._atm_planet, self._mass_constraints, self._fugacity_constraints = (
            param_space.get_atmodeller_constraints(
                mineral_vapour_pressures=mineral_vapour_pressures,
                elems_with_fugconstr=self._chemical_network.elems_with_fugconstr,
            )
        )

        # solve
        self.interior_atmosphere = EquilibriumModel(self._chemical_network.species)
        self.interior_atmosphere.solve(
            state=self._atm_planet,
            mass_constraints=self._mass_constraints,
            fugacity_constraints=self._fugacity_constraints,
        )
        self.output = self.interior_atmosphere.output

        return self.output

    def write_result(self, outpath: os.PathLike) -> None:
        """
        Dump all interesting results into a folder.

        Parameters
        ----------
            outpath : os.PathLike
                The folder where to store the files.
        """

        os.makedirs(outpath, exist_ok=True)

        df = self.__param_space.get_df()
        df.to_csv(outpath + "sample.csv")

        self.output.to_excel(outpath + "atmodeller_output")

    def get_ideal_gas(self) -> IdealGasMixture:
        """
        Turns the gas composition returned by `atmodeller` into an ideal gas object processable by
        `phaethon`.

        This function has hardly any use outside of the phaethon main pipeline.
        """

        # assert that self.__sample is a singular point in a parameter space, i.e., check dim.
        if not self.__param_space.is_singular_point():
            raise ValueError("The parameterspace is not a singular point.")

        # retrieve pressures from atmodeller
        p_bar: dict = {}
        sol_dict: dict = self.interior_atmosphere.output.asdict()
        for gas_name in self.interior_atmosphere.species_network.gas_species_names:
            gas_pressure: float = float(sol_dict[gas_name]["pressure"][0])
            p_bar = p_bar | {gas_name.removesuffix("_g"): gas_pressure}

        # return ideal gas (HELIOS & FastChem assume ideal gas)
        return IdealGasMixture.new_from_pressure(p_bar=p_bar)


# pylint: disable = too-many-instance-attributes
@dataclass(frozen=True)
class SimulationResult:
    """
    Loads and stores the output of a full simulation, and presents tools for quick and easy access
    of the derived atmospheric and planetary parameters.
    """

    param_space: pd.DataFrame
    gas_species: List[str]
    cond_species: List[str]
    elems: List[str]
    volume_mixing_ratios: pd.DataFrame
    dissolved_number_gas: pd.DataFrame
    atmo_number_gas: pd.DataFrame
    dissolved_number_elems: pd.DataFrame
    atmo_number_elems: pd.DataFrame
    log_pressures: pd.DataFrame
    dissolved_ppmw_gas: pd.DataFrame
    atmo_massfrac: pd.Series
    full_dict: Dict[str, pd.DataFrame]

    @classmethod
    def load(cls, path: str | os.PathLike) -> Self:
        """
        Initializes the frames; reads the raw atmodeller input and extracts the relevant info for
        postprocessing.
        """

        path_to_csv = Path(path) / "sample.csv"
        param_space = pd.read_csv(path_to_csv, index_col=0)

        full_dict: Dict[str, pd.DataFrame] = pd.read_excel(
            path + "atmodeller_output.xlsx", sheet_name=None
        )

        # names: gases, condensates & elements
        gas_species: List[str] = []
        cond_species: List[str] = []
        elems: List[str] = []
        for key in full_dict.keys():
            if key.endswith("_g"):
                gas_species.append(key.split("_g")[0])
            elif key.endswith("_cd"):
                cond_species.append(key.split("_cr")[0])
            elif key.startswith("element_"):
                elems.append(key.split("element_")[1])

        # atmospheric mass fraction
        atmo_massfrac = full_dict["gas"]["mass"] / (
            full_dict["gas"]["mass"] + full_dict["state"]["planet_mass"]
        )

        # helper functions for below
        def from_gases(source_var) -> pd.DataFrame:
            return _load_gases_from_atmodeller_output(
                full_dict, gases=gas_species, source_var=source_var
            )

        def from_elems(source_var) -> pd.DataFrame:
            return _load_elements_from_atmodeller_output(
                full_dict, elems=elems, source_var=source_var
            )

        return cls(
            param_space=param_space,
            gas_species=gas_species,
            cond_species=cond_species,
            elems=elems,
            volume_mixing_ratios=from_gases("volume_mixing_ratio"),
            dissolved_number_gas=from_gases("dissolved_number"),
            atmo_number_gas=from_gases("gas_number"),
            dissolved_number_elems=from_elems("dissolved_number"),
            atmo_number_elems=from_elems("gas_number"),
            log_pressures=np.log10(from_gases("pressure")),
            dissolved_ppmw_gas=from_gases("dissolved_ppmw"),
            atmo_massfrac=atmo_massfrac,
            full_dict=full_dict,
        )


def _load_elements_from_atmodeller_output(
    atmodeller_dict, elems: List[str], source_var: str
) -> pd.DataFrame:
    """
    Loads a variable associated to an element from an atmodeller output dict.
    """
    return pd.DataFrame(
        {elem: atmodeller_dict[f"element_{elem}"][source_var] for elem in elems}
    )


def _load_gases_from_atmodeller_output(
    atmodeller_dict, gases: List[str], source_var: str
) -> pd.DataFrame:
    """
    Loads a variable associated to a gas species from an atmodeller output dict.
    """
    return pd.DataFrame({gas: atmodeller_dict[f"{gas}_g"][source_var] for gas in gases})
