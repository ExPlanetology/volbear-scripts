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
Determine vapour pressure evolution along a `series`, i.e. a single parameter (or a combination of
related paramters) evolve within a certain range while all others are kept constant. Shown in Fig.
2, 3, G.1 and G.2
"""

import logging
import numpy as np

# atmodeller
from atmodeller import debug_logger

# volbear
from volbear.param_space import LinearMixingParameterSpace, ParameterSpace
from volbear.equil_mai import MagmaOceanAtmosphereEquilibrator
from volbear.chemical_network import ChemicalNetwork, DEFAULT_CHEMICAL_NETWORK_FILE
from volbear.nat_const import BSE_WT_COMP

# =================================================================================================
# SETTING UP 
# =================================================================================================

logger = debug_logger()
logger.setLevel(logging.INFO)

chemical_network = ChemicalNetwork.new(config_file=DEFAULT_CHEMICAL_NETWORK_FILE)
mai_equilibrator = MagmaOceanAtmosphereEquilibrator(chemical_network=chemical_network)


BASE_OUTPATH: str = "output/MAI_only/"

# =================================================================================================
# FUGACITY
# =================================================================================================

delta_iw: np.ndarray = np.linspace(-6.5, 6.5, 100)
temperature: float = 3000.
OUTPATH: str = BASE_OUTPATH + "fugacity_series/"

for ξ in [0.0, 1.0]:
    for Ψ in [-1, 0, 1]:
        param_space: ParameterSpace = LinearMixingParameterSpace.new(
            planet_mass_mearth=8.,
            cmf=0.325,
            mantle_meltfrac=1.,
            temperature=temperature,
            log_vmf_factor=Ψ,
            mixing_line_param=ξ,
            delta_iw=delta_iw,
            melt_comp_wt=BSE_WT_COMP,
        )

        mai_equilibrator.equil_melt_and_atmosphere(param_space=param_space)
        mai_equilibrator.write_result(
            outpath=OUTPATH
            + f"T{temperature}_ξ{ξ}_Ψ{Ψ}/"
        )

# =================================================================================================
# BATCH (for appendix)
# =================================================================================================

delta_iw_arr = np.linspace(-6.0, 6.0, 100)
OUTPATH: str = BASE_OUTPATH + "batch/"

# different figure for different temperatures
for temperature in [3000]:
    for ξ in [0.0, 0.33, 0.66, 1.0]:
        for Ψ in [-1, 0, 1]:
            param_space: ParameterSpace = LinearMixingParameterSpace.new(
                planet_mass_mearth=8.,
                cmf=0.325,
                mantle_meltfrac=1.,
                temperature=temperature,
                log_vmf_factor=Ψ,
                mixing_line_param=ξ,
                delta_iw=delta_iw_arr,
                melt_comp_wt=BSE_WT_COMP,
            )

            mai_equilibrator.equil_melt_and_atmosphere(param_space=param_space)
            mai_equilibrator.write_result(
                outpath=OUTPATH
                + f"T{temperature}_ξ{ξ}_Ψ{Ψ}/"
            )