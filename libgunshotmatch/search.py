#!/usr/bin/env python3
#
#  search.py
"""
Library search functions.
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
from typing import Iterable, List

# 3rd party
import pandas  # type: ignore[import]
import pyms_nist_search
from pyms.Peak.Class import Peak  # type: ignore[import]

# this package
from libgunshotmatch.peak import QualifiedPeak
from libgunshotmatch.utils import round_rt

__all__ = ["identify_peaks"]


def identify_peaks(
		engine: pyms_nist_search.Engine,
		peaks_to_identify: Iterable[float],
		peak_list: List[Peak],
		n_hits: int = 10,
		verbose: bool = False,
		) -> List[QualifiedPeak]:
	"""
	Identify the peaks in ``peak_list`` where their retention times are in ``peaks_to_identify``.

	:param engine:
	:param peaks_to_identify:
	:param peak_list:
	:param n_hits: The number of hits to return for each peak.
	:param verbose: Enable debug logging
	"""

	# TODO: Shared engine between multiple calls to identify_peaks
	# (perhaps wrap this function in a class)

	# Convert float retention times to Decimal
	# rt_list = [rounders(rt, "0.0000000000") for rt in target_times]
	target_times = pandas.Series(peaks_to_identify).apply(round_rt)

	# Remove NaN values
	rt_list = [rt for rt in target_times if not rt.is_nan()]

	# Sort smallest to largest
	rt_list.sort()

	# # Obtain area for each peak
	# peak_area_list = get_area_list(self.peak_list)
	peaks = []

	# Filter to those peaks present in all samples, by UID
	for peak in peak_list:

		rounded_rt = round_rt(peak.rt / 60)

		if rounded_rt in rt_list:
			qualified_peak = QualifiedPeak.from_peak(peak)
			qualified_peak.peak_number = target_times[target_times == rounded_rt].index[0]

			ms = qualified_peak.mass_spectrum

			if verbose:
				print(f"Identifying peak at rt {rounded_rt} minutes...")

			hit_list = engine.full_spectrum_search(ms, n_hits)

			# Add search results to peak
			for hit in hit_list:
				qualified_peak.hits.append(hit)

			peaks.append(qualified_peak)

	return peaks
