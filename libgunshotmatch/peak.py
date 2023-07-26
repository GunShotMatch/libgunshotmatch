#!/usr/bin/env python3
#
#  peak.py
"""
Classes representing peaks, and functions for peak filtering.
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
from typing import Any, Dict, List, Mapping, Optional, Sequence, Type

# 3rd party
import numpy
import pandas  # type: ignore[import]
import sdjson
from domdf_python_tools.paths import PathPlus
from domdf_python_tools.typing import PathLike
from pyms.BillerBiemann import num_ions_threshold  # type: ignore[import]
from pyms.DPA.Alignment import exprl2alignment  # type: ignore[import]
from pyms.DPA.PairwiseAlignment import Alignment, PairwiseAlignment, align_with_tree  # type: ignore[import]
from pyms.Experiment import Experiment  # type: ignore[import]
from pyms.IonChromatogram import IonChromatogram  # type: ignore[import]
from pyms.Noise.Analysis import window_analyzer  # type: ignore[import]
from pyms.Peak.Class import AbstractPeak, ICPeak, Peak  # type: ignore[import]
from pyms.Peak.List.Function import sele_peaks_by_rt  # type: ignore[import]
from pyms.Spectrum import MassSpectrum  # type: ignore[import]
from pyms_nist_search import SearchResult

__all__ = [
		"PeakList",
		"QualifiedPeak",
		"QualifiedPeakList",
		"align_peaks",
		"filter_aligned_peaks",
		"filter_peaks",
		"peak_from_dict",
		"write_alignment"
		]


class QualifiedPeak(Peak):
	"""
	A Peak that has been identified using NIST MS Search and contains a list of possible identities.

	:param rt: Retention time.
	:param ms: The mass spectrum at the apex of the peak.
	:param minutes: Retention time units flag. If :py:obj:`True`, retention time
		is in minutes; if :py:obj:`False` retention time is in seconds.
	:param outlier: Whether the peak is an outlier.
	:param hits: List of possible identities for the peak.
	:param peak_number: Optional numerical identifier for the :class:`~pyms.Peak.Class.Peak`, such as in an :class:`~.pyms.DPA.Alignment.Alignment`.
	"""

	#: List of possible identities for the peak.
	hits: List[SearchResult]
	#: Optional numerical identifier for the peak, such as in an :class:`~.pyms.DPA.Alignment.Alignment`.
	peak_number: Optional[int]

	def __init__(
			self,
			rt: float = 0.0,
			ms: Optional[MassSpectrum] = None,
			minutes: bool = False,
			outlier: bool = False,
			hits: Optional[List[SearchResult]] = None,
			peak_number: Optional[int] = None
			):

		Peak.__init__(self, rt, ms, minutes, outlier)

		if hits is not None:
			if not isinstance(hits, list) or not isinstance(hits[0], SearchResult):
				raise TypeError("'hits' must be a list of SearchResult objects")

		if hits is None:
			self.hits = []
		else:
			self.hits = hits

		self.peak_number = peak_number

	def __new__(  # noqa: D102
		cls,
		rt: float = 0.0,
		ms: Optional[MassSpectrum] = None,
		minutes: bool = False,
		outlier: bool = False,
		hits: Optional[List[SearchResult]] = None,
		peak_number: Optional[int] = None
		):
		# Overrides __new__ method which warns if the Peak is for an IC (not applicable here)

		obj = object.__new__(cls)
		obj.__init__(rt, ms, minutes, outlier, hits, peak_number)  # type: ignore[misc]
		return obj

	def __eq__(self, other) -> bool:  # noqa: MAN001
		"""
		Return whether this QualifiedPeak object is equal to another object.

		:param other: The other object to test equality with.
		"""

		if isinstance(other, self.__class__):
			return (
					self.UID == other.UID and self.bounds == other.bounds and self.rt == other.rt
					and self.mass_spectrum == other.mass_spectrum and self.area == other.area
					and self.hits == other.hits and self.peak_number == other.peak_number
					)

		return NotImplemented

	@classmethod
	def from_peak(cls: Type["QualifiedPeak"], peak: Peak) -> "QualifiedPeak":
		"""
		Construct :class:`~.QualifiedPeak` from a :class:`~pyms.Peak.Class.Peak`.

		The resulting :class:`~.QualifiedPeak` will not have :attr:`~.hits` or :attr:`~.peak_number` set,
		but those attributes can be set after calling this method.

		:param peak:
		"""

		if not isinstance(peak, AbstractPeak):
			raise TypeError("'peak' must be a Peak object")

		if isinstance(peak, ICPeak):
			new_peak = cls(peak.rt, peak.ic_mass, False, peak.is_outlier)
		else:
			new_peak = cls(peak.rt, peak.mass_spectrum, False, peak.is_outlier)

		new_peak.area = peak.area
		new_peak.bounds = peak.bounds
		new_peak._UID = peak.UID

		return new_peak

	def __repr__(self) -> str:
		return f"<Qualified Peak: {self.rt}>"

	def __str__(self) -> str:
		return self.__repr__()

	def to_dict(self) -> Dict[str, Any]:
		"""
		Returns a dictionary representation of this peak.

		All keys are native, JSON-serializable, Python objects.
		"""

		try:
			ion_areas = self.ion_areas
		except ValueError:
			ion_areas = 0

		return {
				"UID": self.UID,
				"area": self.area,
				"bounds": self.bounds,
				"ion_areas": ion_areas,
				"is_outlier": self.is_outlier,
				"mass_spectrum": {
						"intensity_list": self.mass_spectrum.intensity_list,
						"mass_list": self.mass_spectrum.mass_list,
						},
				"rt": self.rt,
				"hits": [hit.to_dict() for hit in self.hits],
				"peak_number": self.peak_number,
				}

	@classmethod
	def from_dict(cls: Type["QualifiedPeak"], d: Mapping[str, Any]) -> "QualifiedPeak":
		"""
		Construct a :class:`~.QualifiedPeak` from a dictionary.

		:param d:
		"""

		hits = [SearchResult.from_dict(hit) for hit in d["hits"]]
		peak_obj = QualifiedPeak(
				d["rt"],
				MassSpectrum(**d["mass_spectrum"]),
				outlier=d["is_outlier"],
				hits=hits,
				peak_number=d["peak_number"]
				)
		peak_obj.bounds = d["bounds"]
		peak_obj._UID = d["UID"]
		peak_obj.area = d["area"]
		if d["ion_areas"]:
			peak_obj.ion_areas = d["ion_areas"]

		return peak_obj


# @prettify_docstrings
class PeakList(List[Peak]):
	"""
	Represents a list of peaks.
	"""

	#: String identifier for the datafile the peaks were detected in.
	datafile_name: Optional[str] = None

	def __repr__(self) -> str:
		if self.datafile_name:
			return f"{self.__class__.__name__}(datafile={self.datafile_name}; <{len(self)} peaks>)"
		else:
			return f"{self.__class__.__name__}(<{len(self)} peaks>)"

	def __str__(self) -> str:
		return self.__repr__()

	def to_list(self) -> List[Dict[str, Any]]:
		"""
		Return a list of pure-Python dictionaries representing the peaks and their mass spectra.
		"""

		peaks_as_pure_list: List[Dict[str, Any]] = []
		for peak in self:
			try:
				ion_areas = peak.ion_areas
			except ValueError:
				ion_areas = None

			peaks_as_pure_list.append({
					"UID": peak.UID,
					"area": peak.area,
					"bounds": peak.bounds,
					"ion_areas": ion_areas,
					"is_outlier": peak.is_outlier,
					"mass_spectrum": {
							"intensity_list": peak.mass_spectrum.intensity_list,
							"mass_list": peak.mass_spectrum.mass_list,
							},
					"rt": peak.rt,
					})

		return peaks_as_pure_list


# @prettify_docstrings
class QualifiedPeakList(List[QualifiedPeak]):
	"""
	Represents a list of qualified peaks.
	"""

	#: String identifier for the datafile the peaks were detected in.
	datafile_name: Optional[str] = None

	def __repr__(self) -> str:
		if self.datafile_name:
			return f"{self.__class__.__name__}(datafile={self.datafile_name}; <{len(self)} peaks>)"
		else:
			return f"{self.__class__.__name__}(<{len(self)} peaks>)"

	def __str__(self) -> str:
		return self.__repr__()

	def to_list(self) -> List[Dict[str, Any]]:
		"""
		Return a list of pure-Python dictionaries representing the peaks and their mass spectra.
		"""

		return [qp.to_dict() for qp in self]


def filter_peaks(
		peak_list: List[Peak],
		tic: IonChromatogram,
		noise_filter: bool = True,
		noise_threshold: int = 2,
		base_peak_filter: Sequence[int] = (73, 147),
		rt_range: Optional[Sequence[float]] = None,
		) -> PeakList:
	"""
	Filter a list of peaks to remove noise and peaks due to e.g. column bleed.

	:param peak_list:
	:param tic: The TIC of the GC-MS data from which these peaks were identified.
	:param noise_filter: Whether to perform automatic noise filtering of the peak list.
	:param noise_threshold: The minimum number of ions that must have intensities above the noise floor, otherwise the peak is excluded.
	:param base_peak_filter: Peaks whose base peak is at one of the listed masses (*m/z*) are excluded.
	:param rt_range: Optional retention time range (in minutes) to filter the peak list to.
	"""

	if noise_filter:
		# Filtering peak lists with automatic noise filtering
		noise_level = window_analyzer(tic)
		# should we also do rel_threshold() here?
		# https://pymassspec.readthedocs.io/en/master/pyms/BillerBiemann.html#pyms.BillerBiemann.rel_threshold
		peak_list = num_ions_threshold(peak_list, noise_threshold, noise_level)

	if rt_range is not None:
		# Crop time range
		peak_list = PeakList(sele_peaks_by_rt(peak_list, [f"{rt_range[0]}m", f"{rt_range[1]}m"]))

	final_peak_list = PeakList()

	for peak in peak_list:
		# Get mass and intensity lists for the mass spectrum at the apex of the peak
		apex_mass_list = peak.mass_spectrum.mass_list
		apex_mass_spec = peak.mass_spectrum.mass_spec

		# Determine the intensity of the base peak in the mass spectrum
		base_peak_intensity = max(apex_mass_spec)

		# Determine the index of the base peak in the mass spectrum
		base_peak_index = [
				index for index, intensity in enumerate(apex_mass_spec) if intensity == base_peak_intensity
				][0]

		# Finally, determine the mass of the base peak
		base_peak_mass = apex_mass_list[base_peak_index]

		# skip the peak if the base peak is at e.g. m/z 73, i.e. septum bleed
		if base_peak_mass in base_peak_filter:  # method.base_peak_filter:
			continue

		final_peak_list.append(peak)

	return final_peak_list


@sdjson.register_encoder(numpy.int64)
def _encode_numpy_int64(obj: numpy.int64) -> int:
	# Necessary as to_dict() sometimes includes int64 which isn't natively serializable.
	return int(obj)


def align_peaks(
		peaks: List[PeakList],
		rt_modulation: float = 2.5,
		gap_penalty: float = 0.3,
		min_peaks: int = 1,
		) -> Alignment:
	"""
	Perform peak alignment.

	:param peaks: List of list of identified peaks. Each :class:`~.PeakList` must have its :attr:`~.PeakList.datafile_name` attribute set.
	:param rt_modulation: Retention time tolerance parameter for pairwise alignments.
	:param gap_penalty: Gap parameter for pairwise alignments.
	:param min_peaks: Minimum number of peaks required for the alignment position to survive filtering.
	"""

	print("\nAligning\n")

	expr_list = []
	for peak_list in peaks:
		if peak_list.datafile_name is None:
			raise ValueError("Cannot align peaks with PeakList.datafile_name unset")
		else:
			expr_list.append(Experiment(peak_list.datafile_name, peak_list))

	F1: List[Alignment] = exprl2alignment(expr_list)
	# F1: List[Alignment] = exprl2alignment([Experiment(d["datafile"].name, d["peak_list"]) for d in expr_list])

	T1 = PairwiseAlignment(
			F1,
			rt_modulation,
			gap_penalty,
			)

	A1: Alignment = align_with_tree(T1, min_peaks=min_peaks)

	return A1


def write_alignment(
		alignment: Alignment,
		project_name: str,
		output_dir: PathLike,
		require_all_datafiles: bool = False,
		) -> None:
	"""
	Write the alignment data (retention times, peak areas, mass spectra) to disk.

	The output files are as follows:


	* :file:`{{project_name}}_alignment_rt.csv`, containing the aligned retention times.
	* :file:`{{project_name}}_alignment_area.csv`, containing the peak areas for the corresponding aligned retention times.
	* :file:`{{project_name}}_alignment_rt.json`, containing the aligned retention times.
	* :file:`{{project_name}}_alignment_area.json`, containing the peak areas for the corresponding aligned retention times.
	* :file:`{{project_name}}_alignment_ms.json`, containing the mass spectra for the corresponding aligned retention times.

	:param alignment:
	:param project_name: The name of the project. Prefixed to all filenames.
	:param output_dir: Directory to store the output files in.
	:param require_all_datafiles: Whether the peak must be present in all experiments to be included in the data frame.
	"""

	output_dir_p = PathPlus(output_dir)

	alignment.write_csv(
			output_dir_p / f"{project_name}_alignment_rt.csv",
			output_dir_p / f"{project_name}_alignment_area.csv",
			)

	rt_alignment = alignment.get_peak_alignment(require_all_expr=require_all_datafiles)
	rt_alignment.to_json(output_dir_p / f"{project_name}_alignment_rt.json")

	area_alignment = alignment.get_area_alignment(require_all_expr=require_all_datafiles)
	area_alignment.to_json(output_dir_p / f"{project_name}_alignment_area.json")

	ms_alignment = alignment.get_ms_alignment(require_all_expr=require_all_datafiles)
	# ms_alignment.to_json(output_dir_p / 'alignment_ms.json')
	with (output_dir_p / f"{project_name}_alignment_ms.json").open('w') as fp:
		sdjson.dump(ms_alignment.to_dict(), fp)


def filter_aligned_peaks(
		alignment: Alignment,
		top_n_peaks: int = 80,
		min_peak_area: float = 0,
		) -> pandas.DataFrame:
	"""
	Filter aligned peaks by minimum average peak area, and to the top ``n`` largest peaks.

	:param alignment:
	:param top_n_peaks: Filter to the largest ``n`` peaks.
	:param min_peak_area: Exclude aligned peaks with an average peak area below this threshold.

	:returns: :class:`pandas.DataFrame` giving the retention times of the aligned peaks.
	"""

	# Get peak area and retention times from Alignment
	area_alignment: pandas.DataFrame = alignment.get_area_alignment(require_all_expr=False)
	rt_alignment: pandas.DataFrame = alignment.get_peak_alignment(require_all_expr=False)

	# Calculate average peak area for each of the aligned peaks
	area_alignment["mean"] = area_alignment[alignment.expr_code].mean(axis=1)
	area_alignment["stdev"] = area_alignment[alignment.expr_code].std(axis=1)

	# print(area_alignment[["mean", "stdev"]])

	# Sort mean_peak_areas from largest (top) to smallest
	area_alignment = area_alignment.sort_values(by="mean")

	########

	print(f"Filtering to the largest {top_n_peaks} peaks with an average peak area above {min_peak_area}")

	# Get indices of largest n peaks based on `ident_top_peaks`
	top_peaks_indices = []
	# print("tail of area_alignment=", area_alignment.tail(top_n_peaks))

	# Limit to the largest `ident_top_peaks` peaks
	for peak_no, areas in area_alignment.tail(top_n_peaks).iterrows():
		# Ignore peak if average peak area is less then min_peak_area
		if areas["mean"] >= min_peak_area:
			top_peaks_indices.append(peak_no)

	# Remove peaks from rt_alignment if they are not in top_peaks_indices,
	# i.e. they are one of the largest n largest peaks
	rt_alignment = rt_alignment.filter(top_peaks_indices, axis=0)

	return rt_alignment


def peak_from_dict(d: Dict[str, Any]) -> Peak:
	"""
	Construct a :class:`~pyms.Peak.Class.Peak` from a dictionary.

	:param d:
	"""

	peak_obj = Peak(d["rt"], MassSpectrum(**d["mass_spectrum"]), outlier=d["is_outlier"])
	peak_obj.bounds = d["bounds"]
	peak_obj._UID = d["UID"]
	peak_obj.area = d["area"]
	if d["ion_areas"]:
		peak_obj.ion_areas = d["ion_areas"]

	return peak_obj


def _to_peak_list(a_list: List[Peak]) -> PeakList:
	"""
	Internal utility to coerce a list of peaks to an actual :class:`~.PeakList`.
	"""

	if isinstance(a_list, PeakList):
		return a_list
	else:
		return PeakList(a_list)
