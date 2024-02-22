#!/usr/bin/env python3
#
#  _utils.py
"""
Utility functions.
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

# stdlib
from typing import List, MutableSequence, Optional

# 3rd party
import numpy

# this package
from libgunshotmatch import project
from libgunshotmatch.consolidate import ConsolidatedPeak

_PaddedPeakList = MutableSequence[Optional[ConsolidatedPeak]]


def _get_padded_peak_list(project: project.Project, rts: List[float]) -> _PaddedPeakList:

	assert project.consolidated_peaks is not None

	padded_cp_list: _PaddedPeakList = []

	for cp in project.consolidated_peaks:
		if rts:
			top_rt = rts.pop(0)
			while numpy.isnan(top_rt):
				padded_cp_list.append(None)
				if rts:
					top_rt = rts.pop(0)
				else:
					break
		padded_cp_list.append(cp)

	return padded_cp_list
