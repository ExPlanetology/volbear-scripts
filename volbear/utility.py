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
Module for all things considered utility.
"""

from typing import Union, Optional, Tuple, List
import numpy as np
import pandas as pd

# from jaxlib.xla_extension import ArrayImpl


def kwargs_are_float_or_equilong_arrays(exceptions: Optional[Union[str, Tuple]] = None):
    """
    Decorator to test if a function has compatible args, i.e. arguments are either all floats or
    floats mixed with multiple arrays of equal length.
    """

    if exceptions is None:
        exceptions = ()

    def decorator_func(function):
        def wrapper(*args, **kwargs):
            # Separate floats and arrays
            floats: List[str] = []
            arrays: List[str] = {}
            for arg, value in kwargs.items():
                if not arg in exceptions:
                    if isinstance(value, (float, int)):
                        floats.append(arg)
                    elif isinstance(value, (list, np.ndarray)):
                        arrays = arrays | {arg: value}
                    else:
                        raise TypeError(
                            f"Argument {arg} must be either float, int, or array."
                        )

            # Check if all arrays have the same length
            if arrays:
                lengths = [len(arr) for arg, arr in arrays.items()]
                if not all(length == lengths[0] for length in lengths):
                    err_msg: str = (
                        "All array arguments must have the same length. Lenghts: "
                    )
                    for arg, arr in arrays.items():
                        err_msg += f"{arg}: {len(arr)}, "
                    raise ValueError(err_msg[:-2])

            # Call the original function
            return function(*args, **kwargs)

        return wrapper

    return decorator_func


def to_float(
    x: Union[float, int, str, np.ndarray, pd.Series],
) -> Union[float, np.ndarray]:
    """
    Unified method to convert scalars and arrays to float.
    """
    if isinstance(x, np.ndarray):
        return x.astype(float)
    if isinstance(x, pd.Series):
        return x.to_numpy(dtype=float)
    if isinstance(x, list):
        return np.array(x, dtype=float)
    # if isinstance(x, ArrayImpl):
    #     if len(x) == 1:
    #         return float(x[0])
    #     else:
    #         return np.array(x)
    return float(x)


def ensure_tuple(unkown_type: Union[str, Tuple[str]]) -> Tuple[str]:
    """Ensures that a Tuple type is generated, even if a string (e.g. ("MyString")) is passed."""
    if isinstance(unkown_type, str):
        return (unkown_type,)
    return unkown_type
