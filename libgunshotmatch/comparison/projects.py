#!/usr/bin/env python3
#
#  project.py
"""
Comparison between projects.
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
import numpy
import pandas  # type: ignore[import-untyped]
from pyms.DPA.Alignment import Alignment

# this package
from libgunshotmatch.project import Project
from libgunshotmatch.utils import create_alignment

# this package
from ._utils import _get_padded_peak_list, _PaddedPeakList

__all__ = ("filter_alignment_to_consolidate", "get_padded_peak_list")


def filter_alignment_to_consolidate(project: Project) -> Alignment:
	"""
	Filter peaks from a project's ``alignment`` to only those which survived the ``consolidate`` process.

	:param project:
	"""

	assert project.consolidated_peaks is not None

	# Sort expr_code and peakpos into order from datafile_data
	desired_order = list(project.datafile_data)[::-1]
	sort_map = [desired_order.index(code) for code in project.alignment.expr_code]
	expr_code = [project.alignment.expr_code[idx] for idx in sort_map]
	peakpos = [project.alignment.peakpos[idx] for idx in sort_map]

	consolidated_peak_retention_times = [cp.rt_list for cp in project.consolidated_peaks]

	aligned_peaks_surviving_consolidate = []
	for aligned_peaks in zip(*peakpos):
		aprt = [p.rt for p in reversed(aligned_peaks)]
		if aprt in consolidated_peak_retention_times:
			aligned_peaks_surviving_consolidate.append(aligned_peaks)

	alignment_surviving_consolidate = list(zip(*aligned_peaks_surviving_consolidate))

	# Sanity check
	for c_expr_peaks, expr_peaks in zip(alignment_surviving_consolidate, peakpos):
		for peak in c_expr_peaks:
			assert peak in expr_peaks

	# Create new Alignment object
	return create_alignment(alignment_surviving_consolidate, expr_code)


def get_padded_peak_list(project: Project, alignment_rts: pandas.DataFrame) -> _PaddedPeakList:
	"""
	Returns a list of consolidated peaks for the project, based on the between-project alignment.

	:param project:
	:param alignment_rts: Pandas DataFrame giving retention times for the peak alignment.
		The output of :meth:`~pyms.DPA.Alignment.Alignment.get_peak_alignment`.

	"""

	rts = [numpy.mean(row[1:]) for row in alignment_rts[list(project.datafile_data)].itertuples()]

	return _get_padded_peak_list(project, rts)
