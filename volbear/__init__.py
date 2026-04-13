#
# Copyright 2024-2025 Fabian L. Seidler
#
# This file is part of moatif.
#
# moatif is free software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# moatif is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with moatif. If not,
# see <https://www.gnu.org/licenses/>.
#
"""
Do **sp**ectra of h**o**t **r**ocky **e**xoplanets constrain its **t**hermod**y**namics at the
**m**agma-**o**cean **at**mosphere **i**nter**f**ace?

Scripts for running models and generating plots for our project about hot rocky Super-Earths and
the expected atmosphere-interior exchange. Which spectral signatures can be related to oxygen
fugacity and co? What atmospheric structures do we expect?
"""
import dotenv
import warnings

if not dotenv.load_dotenv():
    warnings.warn("No .env file was found!")


# from moatif.rt_utils import (
#     SporetyPlanet,
#     EnergyBudgetFromIrradiationTemperature,
#     EnergyBudgetFromOrbit,
#     AtmospellError
# )
