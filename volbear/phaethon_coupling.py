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

# pylint: disable=E1101

"""
Unified set-up for the radiative transfer problem (phaethon) in this study.
"""

import os
from importlib import resources
from dataclasses import dataclass, field, asdict, replace
from typing import Dict, Set, Self, Optional, List

# atmodeller

from phaethon import (
    Planet,
    Star,
    PlanetarySystem,
    FastChemCoupler,
    CondensationMode,
    OutgassingProtocol,
    IdealGasMixture,
)
from phaethon.pipeline import PhaethonPipeline
from phaethon.iterator import MeltTemperatureIterator
from phaethon.postradtrans.petitradtrans_coupling import PetitRadtransCoupler

from volbear.param_space import ParameterSpace, LinearMixingParameterSpace
from volbear.chemical_network import ChemicalNetwork
from volbear.interior_structure_model import interior_structure
from volbear.equil_mai import MagmaOceanAtmosphereEquilibrator

DEFAULT_FASTCHEM_GASES = resources.files("volbear.data") / "FastChem_logK.dat"
DEFAULT_FASTCHEM_CONDENSATES = (
    resources.files("volbear.data") / "FastChem_logK_condensates.dat"
)

DEFAULT_OPAC_SPECIES_HELIOS: Set[str] = {
    "H",
    "H2",
    "He",
    "H2O",
    "OH",
    "CO",
    "CO2",
    "SO",
    "SO2",
    "SH",
    "H2S",
    "SiO",
    "MgO",
    "Mg",
    "Fe",
    "N2",
    "NO",
    "CH4",
    "O2",
}
""" Default opacity species used for the radiative transfer (in phaethon) """

DEFAULT_RAYLEIGH_SCATTERERS_HELIOS: Set[str] = {
    "H",
    "H2",
    "He",
    "H2O",
    "CO",
    "CO2",
    "N2",
    "O2",
    "e-",
}
""" Default rayleigh species used for the radiative transfer (in phaethon) """

DEFAULT_OPAC_SPECIES_PRT: List[str] = [
    "Fe",
    "Mg",
    "Si",
    "SiO",
    "MgO",
    "CO2",
    "CO",
    "H2O",
    "SO2",
    "CH4",
    "HCN",
    "H2S",
    "SO",
    "SH",
]
""" Default opacity species used by petitRADTRANS. The spectrally active molecules and atoms are
the same as for the HELIOS opacities, but some auxillary species requried by HELIOS for the MMW
(e.g., H, He) are not required here. """

DEFAULT_RAYLEIGH_SCATTERERS_PRT = ["CO2", "H2", "H", "He", "H2O", "N2", "O2"]
""" Rayleigh scatterers for petitRADTRANS """

DEFAULT_HELIOS_CONFIG_FILE = resources.files("volbear.data") / "helios_config.dat"


class AtmospellError(Exception):
    """
    Error in the execution of `atmodeller`.
    """

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)


class AtmospellEngine(OutgassingProtocol):
    """
    `OutgassingProtocol` engine that couples `atmodeller` and `muspell`, preparing them for
    `phaethon`. De facto a wrapper around `MagmaOceanAtmosphereEquilibrator`.
    """

    param_space: ParameterSpace
    """ Stores the parameter space of magma-ocean atmosphere interface interactions. """

    runner: MagmaOceanAtmosphereEquilibrator
    """ Runs the simulation. """

    vapour: IdealGasMixture
    """ 
    Stores the output as ideal gas, based on the computed pressures.
    
    NOTE: Output pressures are  ALWAYS converted to an ideal gas mixture since `FastChem` and 
    `HELIOS` make this assumption, even if `species` use a non-ideal EOS! 
    """

    def __init__(
        self,
        param_space: ParameterSpace,
        chemical_network: ChemicalNetwork,
    ) -> None:
        self.param_space = param_space
        self.equilibrator = MagmaOceanAtmosphereEquilibrator(
            chemical_network=chemical_network
        )

    def get_info(self) -> dict:
        """
        Return information about this simulation that will be dumped in a json-file.
        """
        return asdict(self.param_space)

    def equilibriate(self, temperature: float) -> IdealGasMixture:
        """
        Calculates the vapour composition & pressure above the planets surface.
        """

        # reset temperature of frozen dataclass
        self.param_space = replace(self.param_space, temperature=temperature)

        # run simulation
        try:
            self.equilibrator.equil_melt_and_atmosphere(param_space=self.param_space)
        except Exception as e:
            raise AtmospellError(message=e) from e

        # return ideal gas (HELIOS & FastChem assume ideal gas)
        self.vapour = self.equilibrator.get_ideal_gas()

        return self.vapour


@dataclass
class FullyCoupledPlanet:
    """
    Pipeline to streamline the full atmosphere-interior exchange with subsequent radiative
    transfer to obtain atmospheric pressure-temperature profiles and spectra.
    """

    planetary_system: PlanetarySystem
    chemical_network: ChemicalNetwork
    parameter_space: ParameterSpace
    cond_mode: CondensationMode

    opac_path: str
    opac_species: Set[str]
    scatterers: Set[str]

    nvcc_kws: Optional[Dict[str, str]] = field(default_factory=lambda: None)

    @classmethod
    def from_mixingline_and_irradtemp(
        cls,
        *,
        chemical_network: ChemicalNetwork,
        mass: float,
        cmf: float,
        melt_frac: float,
        log_vmf_factor: float,
        mixing_line_param: float,
        delta_iw: float,
        cond_mode: CondensationMode,
        intrinsic_temp: float,
        dilution_factor: float,
        irrad_temp: float,
        star: Star,
        opac_species_helios: Set[str],
        scatterers_helios: Set[str],
        planet_name: str = "unnamed planet",
        melt_albedo: float = 0.0,
        opac_path: Optional[str] = None,
        nvcc_kws: Optional[Dict[str, str]] = None,
    ) -> Self:

        """
        Build a full planet simulation; mildly hardcoded for the cases shown in Seidler et al. 2026

        Parameters
        ----------
            chemical_network: ChemicalNetwork
                Chemical network of atmosphere species (dissolvable and non-dissolvable gases, 
                condensates). 
            mass: float
                Mass of the planet, in Earth-masses.
            cmf: float
                Core mass fraction.
            melt_frac: float
                Mantle melt fraction.
            log_vmf_factor: float
                Logarithm of volatile mass fraction factor (Ψ)
            mixing_line_param: float
                Mixing line paramter (ξ)
            delta_iw: float
                Oxygen fugacity shift relative to the iron-wüstite buffer (IW); in dex.
            cond_mode: CondensationMode
                FastChem condensation mode. Must NEVER be rainout!
            intrinsic_temp: float
                Intrinsic temperature, in Kelvin.
            dilution_factor: float
                Dilution factor accounting for atmospheric heat redistribtuion.
            irrad_temp: float
                Irradiation temperature of the planet (as defined by Seidler et al. 2026)
            star: Star
                A `phaethon`-Star object.
            opac_species_helios: Set[str]
                A set of opacity species that HELIOS should use (both molecular absorbers and
                continuum absorbers, i.e., CIA)
            scatterers_helios: Set[str]
                A set of Rayleigh species that HELIOS should use.
            melt_albedo:
                Albedo of the silicate melt.
            opac_path:
                Path to the opacities; default is None, and it will be overwritten from the 
                environment varialbe `OPAC_PATH` if not specified.
            nvcc_kws: Optional[Dict[str, str]] = None
                Optional arguments for the CUDA compiler.
        """

        if nvcc_kws is None:
            nvcc_kws = {}

        # evaluate interior structure model
        radius: float = interior_structure(mass, cmf)

        # phaethon planet
        planet = Planet(
            name=planet_name,
            mass=mass,
            radius=radius,
            bond_albedo=melt_albedo,
            dilution_factor=dilution_factor,
            intrinsic_temperature=intrinsic_temp,
        )

        # planetary system
        planetary_system = PlanetarySystem.build_from_irrad_temp(
            irrad_temp=irrad_temp, planet=planet, star=star
        )

        # parameter space
        parameter_space = LinearMixingParameterSpace.new(
            planet_mass_mearth=mass,
            cmf=cmf,
            mantle_meltfrac=melt_frac,
            temperature=2500.0,
            log_vmf_factor=log_vmf_factor,
            mixing_line_param=mixing_line_param,
            delta_iw=delta_iw,
        )

        return cls(
            planetary_system=planetary_system,
            chemical_network=chemical_network,
            parameter_space=parameter_space,
            cond_mode=cond_mode,
            opac_path=os.environ.get("OPAC_PATH") if opac_path is None else opac_path,
            opac_species=opac_species_helios,
            scatterers=scatterers_helios,
            nvcc_kws=nvcc_kws,
        )

    def run(
        self,
        outdir: str,
        *,
        helios_config_file: str | os.PathLike = DEFAULT_HELIOS_CONFIG_FILE,
        skip_exists: bool = False,
        t_melt_init: Optional[float] = None,
        logfile_name: str = "phaethon.log",
        prt_radtrans: PetitRadtransCoupler,
    ):
        """
        Runs the full simulation for a planet: chemical equilibrium at MAI and radiative transfer
        in the atmosphere.

        Parameters
        ----------
            helios_param_file : str
                The path to the parameterfile of HELIOS
            nlayer : int
                Layers in the atmosphere. Default is 38.
            skip_exists : bool
                If target already exists, skip computation
        """

        _outdir = outdir
        if not _outdir.endswith("/"):
            _outdir += "/"

        # skip target if it already exists?
        if skip_exists:
            # TODO: improve this - maybe deposit a "failure" flag in the output if it failed?
            # Also, check if HELIOS has properly deposited its stuff.
            if os.path.isdir(_outdir + "HELIOS_iterative/"):
                # TODO: Inform logger
                return

        # vapour engine object for phaethon
        outgassing: AtmospellEngine = AtmospellEngine(
            param_space=self.parameter_space,
            chemical_network=self.chemical_network,
        )

        # setting up FastChem
        fastchem_coupler = FastChemCoupler(
            path_to_eqconst=DEFAULT_FASTCHEM_GASES,
            path_to_condconst=DEFAULT_FASTCHEM_CONDENSATES,
            verbosity_level=0,  # silent
            ref_elem="O",  # oxygen is always present and is a better reference elem than H
            cond_mode=self.cond_mode,
        )

        # setting-up phaethon
        pipeline = PhaethonPipeline(
            planetary_system=self.planetary_system,
            outgassing=outgassing,
            fastchem_coupler=fastchem_coupler,
            outdir=_outdir,
            opac_species=self.opac_species,
            scatterers=self.scatterers,
            opacity_path=self.opac_path,
            iterator=MeltTemperatureIterator(
                delta_temp_abstol=35.0,
                max_iter=15,
                tmelt_limits=(50, 5000),
            ),
            postradtrans=prt_radtrans,
        )

        # run phaethon
        pipeline.run(
            config_file=helios_config_file,
            nvcc_kws=self.nvcc_kws,
            t_melt_init=t_melt_init,
            logfile_name=logfile_name,
        )
