#!/usr/bin/env python3
#
#  _spectra.py
"""
Internal spectrum comparison utilities.
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
from typing import List, Optional, Sequence, Tuple

# 3rd party
from pyms.Spectrum import MassSpectrum  # type: ignore[import]

# this package
from libgunshotmatch.utils import ms_comparison

PermsListType = List[Tuple[int, int]]


def _row_ms_comparison(ms_data: Sequence[MassSpectrum], perms: PermsListType) -> List[Optional[float]]:
	"""
	Performs pairwise mass spectral comparisons of the mass spectra in ``ms_data``.

	:param ms_data:
	:param perms: List of tuples giving the indices of the specta in ``ms_data`` to compare.

	:returns: List of similarity scores for each permutation.
	"""

	similarity_list = []

	for perm in perms:

		top_ms = ms_data[perm[0]]
		bottom_ms = ms_data[perm[1]]

		similarity_list.append(ms_comparison(top_ms, bottom_ms))

	return similarity_list


def _map_func(
		peak_number: int,
		ms_data: Sequence[MassSpectrum],
		perms: PermsListType,
		) -> Tuple[int, List[Optional[float]]]:
	"""
	Multiprocessing helper function.

	Wraps :func:`~.row_ms_comparison` and also returns the peak number.
	"""

	return peak_number, _row_ms_comparison(ms_data, perms)
