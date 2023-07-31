#!/usr/bin/env python3
#
#  project.py
"""
Represents a collection of repeat analyses.
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
import os
from typing import Any, Dict, List, Mapping, Optional, Type

# 3rd party
import attr
import numpy
import pandas  # type: ignore[import]
import pyms_nist_search
from domdf_python_tools.typing import PathLike
from pyms.DPA.Alignment import Alignment  # type: ignore[import]
from pyms.Peak.Class import Peak  # type: ignore[import]

# this package
from libgunshotmatch import gzip_util
from libgunshotmatch.consolidate import (
		ConsolidatedPeak,
		ConsolidatedPeakFilter,
		match_counter,
		pairwise_ms_comparisons
		)
from libgunshotmatch.datafile import Repeat
from libgunshotmatch.peak import PeakList, QualifiedPeak, peak_from_dict

__all__ = ["Project"]


@attr.define
class Project:
	"""
	A project represents the aligned peaks from multiple datafiles.
	"""

	#: The name of the project.
	name: str

	#: Peak alignment for the repeats in this project.
	alignment: Alignment

	# datafile_data: Dict[str, DatafileDataElement]

	#: Mapping of repeat names to :class:`~.Repeat` objects.
	datafile_data: Dict[str, Repeat]

	#: List of peaks after :meth:`~.consolidate` is performed. :py:obj:`None` initially.
	consolidated_peaks: Optional[List[ConsolidatedPeak]] = attr.field(default=None)

	def to_dict(self) -> Dict[str, Any]:
		"""
		Returns a dictionary representation of this :class:`~.Project`.

		All keys are native, JSON-serializable, Python objects.
		"""

		alignment_as_dict = {
				"peaks": [PeakList(x).to_list() for x in self.alignment.peakpos],
				"expr_code": self.alignment.expr_code,
				"similarity": self.alignment.similarity,
				}

		if self.consolidated_peaks is None:
			consolidated_peaks_as_list = None
		else:
			consolidated_peaks_as_list = [cp.to_dict() for cp in self.consolidated_peaks]

		datafile_data_as_dict = {k: v.to_dict() for k, v in self.datafile_data.items()}

		return {
				"name": self.name,
				"alignment": alignment_as_dict,
				"datafile_data": datafile_data_as_dict,  # "datafile_data": list(self.datafile_data.keys()),
				"consolidated_peaks": consolidated_peaks_as_list,
				}

	def export(self, output_dir: PathLike) -> str:
		"""
		Export as a ``gsmp`` file.

		:returns: The output filename.
		"""

		export_filename = os.path.join(output_dir, f"{self.name}.gsmp")
		gzip_util.write_gzip_json(export_filename, self.to_dict(), indent=0)
		return export_filename

	@classmethod
	def from_file(cls: Type["Project"], filename: PathLike) -> "Project":
		"""
		Parse a ``gsmp`` file.

		:param filename: The input filename.
		"""

		as_dict: Dict[str, Any] = gzip_util.read_gzip_json(filename)  # type: ignore[assignment]
		return cls.from_dict(as_dict)

	@classmethod
	def from_dict(cls: Type["Project"], d: Mapping[str, Any]) -> "Project":
		"""
		Construct a :class:`~.Project` from a dictionary.

		:param d:
		"""

		alignment_as_dict = d["alignment"]
		alignment = Alignment(None)
		alignment_peaks: List[List[Peak]] = []
		for row in alignment_as_dict["peaks"]:
			alignment_peaks.append([])
			for peak in row:
				# print(peak)
				peak_obj = peak_from_dict(peak)

				alignment_peaks[-1].append(peak_obj)

		alignment.peakpos = alignment_peaks
		alignment.peakalgt = numpy.transpose(alignment.peakpos).tolist()
		alignment.expr_code = alignment_as_dict["expr_code"]
		alignment.similarity = alignment_as_dict["similarity"]

		consolidated_peaks_as_list = d["consolidated_peaks"]
		if consolidated_peaks_as_list is None:
			consolidated_peaks = None
		else:
			consolidated_peaks = [ConsolidatedPeak.from_dict(cp) for cp in consolidated_peaks_as_list]

		datafile_data = {k: Repeat.from_dict(v) for k, v in d["datafile_data"].items()}

		return cls(
				name=d["name"],
				alignment=alignment,
				consolidated_peaks=consolidated_peaks,
				datafile_data=datafile_data,
				)

	def consolidate(
			self,
			engine: pyms_nist_search.Engine,
			peak_filter: Optional[ConsolidatedPeakFilter] = None,
			) -> pandas.DataFrame:
		"""
		Consolidate the compound identification from the experiments into a single dataset.

		:param engine:
		:param peak_filter: Filter for the consolidated peaks.

		:returns: :class:`pandas.DataFrame` giving the results of pairwise mass spectral comparisons
			between the repeats for each aligned peak.
		"""

		ms_comparison_df = pairwise_ms_comparisons(self.alignment)

		peak_numbers: List[int] = []
		peak: Optional[QualifiedPeak]

		qualified_peak_array = []

		# for experiment in project.alignment.expr_code:
		for experiment in self.datafile_data:
			qualified_peaks = self.datafile_data[experiment].qualified_peaks
			assert qualified_peaks is not None
			for peak in qualified_peaks:
				assert peak.peak_number is not None
				peak_numbers.append(peak.peak_number)
			qualified_peak_array.append(qualified_peaks)

		# Convert peak_numbers to a set and sort smallest to largest
		peak_numbers = sorted(set(peak_numbers))

		consolidated_peaks = match_counter(
				engine=engine,
				peak_numbers=peak_numbers,
				qualified_peaks=qualified_peak_array,
				ms_comp_data=ms_comparison_df,
				)

		if peak_filter is None:
			self.consolidated_peaks = consolidated_peaks
		else:
			self.consolidated_peaks = peak_filter.filter(consolidated_peaks)

		return ms_comparison_df

		# chart_data = make_chart_data(self)
