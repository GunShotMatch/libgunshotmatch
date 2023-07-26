#!/usr/bin/env python3
#
#  _fields.py
"""
Internal attrs field helpers.
"""
#
#  Copyright © 2020-2023 Dominic Davis-Foster <dominic@davis-foster.co.uk>
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
from typing import Dict, Mapping, Optional, Union

# 3rd party
import pandas  # type: ignore[import]
from pyms_nist_search import ReferenceData


def _attrs_convert_cas(cas: str) -> str:
	if cas == "0-00-0":
		cas = "---"

	return cas


_reference_data_error_msg = ''.join([
		"'reference_data' must be a `pyms_nist_search.ReferenceData` object,",
		"a dictionary representing a `ReferenceData` object,",
		"or `None`",
		])


def _attrs_convert_reference_data(reference_data: Union[Dict, ReferenceData, None]) -> Optional[ReferenceData]:

	if reference_data is None:
		return None

	elif isinstance(reference_data, ReferenceData):
		return reference_data

	elif isinstance(reference_data, dict):
		expected_keys = {"name", "cas", "formula", "contributor", "nist_no", "id", "mw", "mass_spec", "synonyms"}
		if set(reference_data.keys()) != expected_keys:
			# print(set(reference_data.keys()))
			raise TypeError(_reference_data_error_msg)
		else:
			return ReferenceData(**reference_data)

	else:
		raise TypeError(_reference_data_error_msg)


def _attrs_convert_ms_comparison(
		ms_comparison: Union[Mapping[str, float], pandas.Series, None],
		) -> pandas.Series:

	if ms_comparison is None:
		return pandas.Series()
	elif isinstance(ms_comparison, pandas.Series):
		return ms_comparison
	elif isinstance(ms_comparison, Mapping):
		return pandas.Series(ms_comparison)
	else:
		raise TypeError("'ms_comparison' must be a mapping or a pandas.Series")
