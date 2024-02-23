#!/usr/bin/env python3
#
#  unknown.py
"""
Comparison between project(s) and unknowns.
"""
#
#  Copyright © 2024 Dominic Davis-Foster <dominic@davis-foster.co.uk>
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

# 3rd party
import pandas  # type: ignore[import-untyped]
from pyms.DPA.Alignment import Alignment

# this package
from libgunshotmatch.project import Project
from libgunshotmatch.utils import create_alignment

# this package
from ._utils import _get_padded_peak_list, _PaddedPeakList

__all__ = ("filter_alignment_to_consolidate", "get_padded_peak_list")


def filter_alignment_to_consolidate(unknown: Project) -> Alignment:
	"""
	Filter peaks from an unknown's ``alignment`` to only those which survived the ``consolidate`` process.

	:param project:
	"""

	assert unknown.consolidated_peaks is not None

	cp_rts = [cp.rt for cp in unknown.consolidated_peaks]

	aligned_peaks = [peak for peak in unknown.alignment.peakpos[0] if peak.rt in cp_rts]

	# Sanity check
	for c_expr_peaks, expr_peaks in zip([aligned_peaks], unknown.alignment.peakpos):
		for peak in c_expr_peaks:
			assert peak in expr_peaks

	# Create new Alignment object
	return create_alignment([aligned_peaks], [unknown.name])


def get_padded_peak_list(unknown: Project, alignment_rts: pandas.DataFrame) -> _PaddedPeakList:
	"""
	Returns a list of consolidated peaks for the unknown sample, based on the between-project alignment.

	:param unknown:
	:param alignment_rts: Pandas DataFrame giving retention times for the peak alignment.
		The output of :meth:`~pyms.DPA.Alignment.Alignment.get_peak_alignment`.
	"""

	rts = list(alignment_rts[unknown.name])

	return _get_padded_peak_list(unknown, rts)
