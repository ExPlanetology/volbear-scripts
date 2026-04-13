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
Chemical network representing the atmosphere and soluble melt species. Used by `atmodeller`.
"""

from os import PathLike
from importlib import resources
from copy import copy
from typing import Tuple, Self
from dataclasses import dataclass

import yaml
from atmodeller import ChemicalSpecies, SpeciesNetwork
from atmodeller.solubility import get_solubility_models, NoSolubility

SOLUBILITY_MODELS = get_solubility_models()
DEFAULT_CHEMICAL_NETWORK_FILE = resources.files("volbear.data") / "chem_network.yaml"


@dataclass(frozen=True)
class ChemicalNetwork:
    """
    Stores the chemical network of atmospheric gas species, exsolved mineral vapour, graphite
    condensates and dissolved volatiles. Basically, a convient wrapper around 
    `atmodeller.SpeciesNetwork`.
    """

    name: str
    """ Name/identifier of the network. """

    species: SpeciesNetwork
    """ Network of interacting chemical species. An atmodeller object. """

    elems_with_fugconstr: Tuple[str]
    """ Elements other than oxygen that have fugacity constraints, i.e. the mineral gases Si, Mg,
    and Fe. Not necessarily an 'element' in the strict sense, because Si-gas is replaced by SiO for
    stability reasons. """

    @classmethod
    def new(cls, config_file: str | PathLike = DEFAULT_CHEMICAL_NETWORK_FILE) -> Self:
        """
        Builds a new chemical network from a configuration file.

        Parameters
        ----------
            config_file : str | PathLike
                Configuration file containing the desired gas species, condensates and solubility
                laws for volatiles.
        """

        with open(config_file, "r", encoding="utf8") as cfg_file:
            config = yaml.safe_load(cfg_file)

        # gas species; they might have a solubility law attached.
        chemical_species = []
        for gas_formula, gas_props in config["gases"].items():
            if gas_props.get("solubility") is None:
                solubility = NoSolubility()
            else:
                solubility = SOLUBILITY_MODELS[gas_props.get("solubility")]

            gas = ChemicalSpecies.create_gas(
                formula=gas_formula,
                solubility=solubility,
            )
            chemical_species.append(gas)

        # condensates
        for cond_formula in config["condensates"]:
            condensate = ChemicalSpecies.create_condensed(formula=cond_formula)
            chemical_species.append(condensate)

        # assemble full reaction network
        species_network = SpeciesNetwork(chemical_species)

        # mineral gas fugacity constraints
        unique_elements = list(species_network.unique_elements)
        mineral_gases = copy(unique_elements)
        for required_element in ["H", "He", "C", "N", "O", "S"]:
            assert (
                required_element in unique_elements
            ), f"At least one gas/condensate species of element '{required_element}' is required"
            mineral_gases.remove(required_element)

        # SiO seems to be the more stable constraint, based on earlier versions of atmodeller
        # (becase more abundant? TODO: verify!)
        mineral_gases = [
            "SiO" if rocky_gas == "Si" else rocky_gas for rocky_gas in mineral_gases
        ]

        return cls(
            name=config["name"],
            species=species_network,
            elems_with_fugconstr=mineral_gases,
        )
