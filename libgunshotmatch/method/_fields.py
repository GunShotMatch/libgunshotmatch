#!/usr/bin/env python3
#
#  _fields.py
"""
Internal attrs field helpers.
"""
#
#  Copyright © 2023 Dominic Davis-Foster <dominic@davis-foster.co.uk>
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in all
#  copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
#  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
#  DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
#  OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
#  OR OTHER DEALINGS IN THE SOFTWARE.
#

# stdlib
from typing import Iterable, Optional, Set, Tuple, TypeVar, Union, overload

# 3rd party
from pyms.Utils.Time import time_str_secs

__all__ = ["convert_crop_mass_range", "convert_rt_range", "convert_sg_window", "default_base_peak_filter"]

_FT = TypeVar("_FT")


def convert_crop_mass_range(value: Optional[Iterable]) -> Optional[Tuple[int, int]]:
	"""
	Convert a value to a tuple of integers for the ``crop_mass_range`` option.

	:param value:
	"""

	if not value:
		return None

	value_as_tuple = tuple([int(x) for x in value])
	if len(value_as_tuple) != 2:
		err_msg = f"crop_mass_range must be a 2-element tuple. Got {len(value_as_tuple)} values: {value_as_tuple}"
		raise ValueError(err_msg)

	return value_as_tuple


@overload
def convert_sg_window(value: str) -> str: ...


@overload
def convert_sg_window(value: int) -> int: ...


def convert_sg_window(value: Union[str, int]) -> Union[str, int]:
	"""
	Convert a Savitzky-Golay window size parameter value.

	:param value:
	"""

	if isinstance(value, int):
		return value
	elif isinstance(value, str):
		# Check the string can be correctly parsed
		time_str_secs(value)
		return value
	else:
		raise TypeError("'savitzky_golay_window' must be either an integer or a string")


def convert_rt_range(value: Optional[Iterable]) -> Optional[Tuple[float, float]]:
	"""
	Convert a value to a tuple of floats for the ``rt_range`` option.

	:param value:
	"""

	if not value:
		return None

	value_as_tuple = tuple([float(x) for x in value])
	if len(value_as_tuple) != 2:
		raise ValueError(f"rt_range must be a 2-element tuple. Got {len(value_as_tuple)} values: {value_as_tuple}")

	return value_as_tuple


def default_base_peak_filter() -> Set[int]:
	"""
	Returns the default value for the ``base_peak_filter`` option.
	"""

	return {73, 147}
