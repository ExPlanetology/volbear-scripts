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
The main script; generates fully equilibrated atmospheres for magma-ocean planets.
"""

import os
import numpy as np
from astropy import units
from phaethon import Star, CondensationMode
from phaethon.postradtrans.petitradtrans_coupling import PetitRadtransCoupler

from volbear.chemical_network import ChemicalNetwork, DEFAULT_CHEMICAL_NETWORK_FILE
from volbear.phaethon_coupling import (
    FullyCoupledPlanet,
    DEFAULT_OPAC_SPECIES_HELIOS,
    DEFAULT_RAYLEIGH_SCATTERERS_HELIOS,
    DEFAULT_OPAC_SPECIES_PRT,
    DEFAULT_RAYLEIGH_SCATTERERS_PRT,
    DEFAULT_HELIOS_CONFIG_FILE,
)

SUN = Star(
    name="Sun",
    mass=1.0 * units.M_sun,
    radius=1.0 * units.R_sun,
    t_eff=5770.0 * units.K,
    distance=10.0 * units.pc,
    metallicity=0.0 * units.dex,
)
SUN.get_spectrum_from_file(
    outdir="output/grid/",
    source_file="volbear/data/sun_gueymard_2003_modified.txt",
    opac_file_for_lambdagrid=os.environ.get("OPAC_PATH") + "SiO_opac_ip_kdistr.h5",
    skiprows=9,
    plot_and_tweak=False,
    w_conversion_factor=1e-7,
    flux_conversion_factor=1e10,
)

prt_radtrans = PetitRadtransCoupler(
    line_species=list(DEFAULT_OPAC_SPECIES_PRT),
    wlen_bords_micron=(0.2, 60),
    rayleigh_species=list(DEFAULT_RAYLEIGH_SCATTERERS_PRT),
)

# =================================================================================================
# GRID
# =================================================================================================

for irrad_temp in [2500]:
    for Ψ in [0]:
        for ξ in [1.0]:
            for delta_iw in [3]:
                planet = FullyCoupledPlanet.from_mixingline_and_irradtemp(
                    chemical_network=ChemicalNetwork.new(DEFAULT_CHEMICAL_NETWORK_FILE),
                    mass=8.0,
                    cmf=0.325,
                    melt_frac=1.0,
                    log_vmf_factor=Ψ,
                    mixing_line_param=ξ,
                    delta_iw=delta_iw,
                    cond_mode=CondensationMode.EQ_COND,
                    intrinsic_temp=0.0,
                    dilution_factor=2 / 3,
                    irrad_temp=irrad_temp,
                    star=SUN,
                    opac_species_helios={"H2O"}, #DEFAULT_OPAC_SPECIES_HELIOS,
                    scatterers_helios=DEFAULT_RAYLEIGH_SCATTERERS_HELIOS,
                    planet_name="HRE",
                )

                planet.run(
                    helios_config_file=DEFAULT_HELIOS_CONFIG_FILE,
                    outdir=f"output/grid/tirr_{irrad_temp}/Ψ{Ψ}/ξ{ξ}/IW{delta_iw}/",
                    t_melt_init=3000.0,
                    skip_exists=False,
                    prt_radtrans=prt_radtrans,
                )

# =================================================================================================
# MASS SERIE
# =================================================================================================

MASS_ARR = np.logspace(np.log10(0.8), 1, 12)

for irrad_temp in [1250, 2500]:
    for ξ in [0.0, 1.0]:
        for delta_iw in [-3, 3]:
            for mass in MASS_ARR:
                planet = FullyCoupledPlanet.from_mixingline_and_irradtemp(
                    chemical_network=ChemicalNetwork.new(DEFAULT_CHEMICAL_NETWORK_FILE),
                    mass=8.0,
                    cmf=0.325,
                    melt_frac=1.0,
                    log_vmf_factor=0,
                    mixing_line_param=ξ,
                    delta_iw=delta_iw,
                    cond_mode=CondensationMode.EQ_COND,
                    intrinsic_temp=0.0,
                    dilution_factor=2 / 3,
                    irrad_temp=irrad_temp,
                    star=SUN,
                    opac_species_helios=DEFAULT_OPAC_SPECIES_HELIOS,
                    scatterers_helios=DEFAULT_RAYLEIGH_SCATTERERS_HELIOS,
                    planet_name="HRE",
                )

                planet.run(
                    helios_config_file=DEFAULT_HELIOS_CONFIG_FILE,
                    outdir=f"output/mass_series/tirr_{irrad_temp}/Ψ{0}/ξ{ξ}/IW{delta_iw}/mass_{mass}/",
                    t_melt_init=3000.0,
                    skip_exists=False,
                    prt_radtrans=prt_radtrans,
                )

# =================================================================================================
# EFFECT OF INTRINSIC TEMPERATURE
# =================================================================================================

IRRAD_TEMP = 2500

for intrinsic_temperature in [100, 200, 300]:
    for ξ in [0.0, 1.0]:
        for delta_iw in [-3, 3]:
            planet = FullyCoupledPlanet.from_mixingline_and_irradtemp(
                chemical_network=ChemicalNetwork.new(DEFAULT_CHEMICAL_NETWORK_FILE),
                mass=8.0,
                cmf=0.325,
                melt_frac=1.0,
                log_vmf_factor=0,
                mixing_line_param=ξ,
                delta_iw=delta_iw,
                cond_mode=CondensationMode.EQ_COND,
                intrinsic_temp=intrinsic_temperature,
                dilution_factor=2 / 3,
                irrad_temp=IRRAD_TEMP,
                star=SUN,
                opac_species_helios=DEFAULT_OPAC_SPECIES_HELIOS,
                scatterers_helios=DEFAULT_RAYLEIGH_SCATTERERS_HELIOS,
                planet_name="HRE",
            )

            planet.run(
                helios_config_file=DEFAULT_HELIOS_CONFIG_FILE,
                outdir=f"output/tint_series/tirr_{IRRAD_TEMP}/tint_{intrinsic_temperature}/Ψ{0}/ξ{ξ}/IW{delta_iw}/",
                t_melt_init=3000.0,
                skip_exists=False,
                prt_radtrans=prt_radtrans,
            )
