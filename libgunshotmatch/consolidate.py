#!/usr/bin/env python3
#
#  consolidate.py
"""
Functions for combining peak identifications across aligned peaks into a single set of results.
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
from typing import Any, Dict, Iterator, List, Mapping, Optional, Tuple, Type, Union

# 3rd party
import attr
import numpy
import pandas  # type: ignore[import]
from pyms.Spectrum import MassSpectrum  # type: ignore[import]
from pyms_nist_search import ReferenceData

__all__ = ["ConsolidatedPeak", "ConsolidatedSearchResult"]


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


@attr.define(frozen=True)
class ConsolidatedSearchResult:
	"""
	Represents a candidate compound for a peak.

	This is determined from a set of :class:`SearchResults <pyms_nist_search.search_result.SearchResult>` for a set of aligned peaks.
	"""

	# TODO:  Not currently copying hit_prob from SearchResult

	#: The name of the candidate compound.
	name: str

	#: The CAS number of the compound.
	cas: str = attr.field(converter=_attrs_convert_cas)

	mf_list: List[int] = attr.field(default=attr.Factory(list))
	"""
	List of Match Factors comparing the mass spectrum of the peak with the reference spectrum in each aligned peak.

	Will contain NaN where the compound was not in the hit list for a peak.
	"""

	rmf_list: List[int] = attr.field(default=attr.Factory(list))
	"""
	List of Reverse Match Factors comparing the reference spectrum with the spectrum for each aligned peak.

	Will contain NaN where the compound was not in the hit list for a peak.
	"""

	hit_numbers: List[int] = attr.field(default=attr.Factory(list))
	"""
	List of "hit" numbers from NIST MS Search.

	Lower is better. Will contain NaN where the compound was not in the hit list for a peak.
	"""

	#: The reference mass spectrum for the compound from the NIST library.
	reference_data: Optional[ReferenceData] = attr.field(converter=_attrs_convert_reference_data, default=None)

	@property
	def match_factor(self) -> float:
		"""
		The average match factor.

		Missing values (where the compound was not in the hit list for a peak) are excluded from the calculation.
		"""

		return numpy.nanmean(self.mf_list)  # type: ignore[return-value]  # mypy thinks type is "floating[Any]"

	@property
	def match_factor_stdev(self) -> float:
		"""
		The standard deviation of the match factors.

		Missing values (where the compound was not in the hit list for a peak) are excluded from the calculation.
		"""

		return numpy.nanstd(self.mf_list)  # type: ignore[return-value]  # mypy thinks type is "floating[Any]"

	@property
	def reverse_match_factor(self) -> float:
		"""
		The average reverse match factor.

		Missing values (where the compound was not in the hit list for a peak) are excluded from the calculation.
		"""

		return numpy.nanmean(self.rmf_list)  # type: ignore[return-value]  # mypy thinks type is "floating[Any]"

	@property
	def reverse_match_factor_stdev(self) -> float:
		"""
		The standard deviation of the reverse match factors.

		Missing values (where the compound was not in the hit list for a peak) are excluded from the calculation.
		"""

		return numpy.nanstd(self.rmf_list)  # type: ignore[return-value]  # mypy thinks type is "floating[Any]"

	@property
	def average_hit_number(self) -> float:
		"""
		The average hit number.

		Missing values (where the compound was not in the hit list for a peak) are excluded from the calculation.
		"""

		return numpy.nanmean(self.hit_numbers)  # type: ignore[return-value]  # mypy thinks type is "floating[Any]"

	@property
	def hit_number_stdev(self) -> float:
		"""
		The standard deviation of the hit numbers.

		Missing values (where the compound was not in the hit list for a peak) are excluded from the calculation.
		"""

		return numpy.nanstd(self.hit_numbers)  # type: ignore[return-value]  # mypy thinks type is "floating[Any]"

	def __len__(self) -> int:
		"""
		The number of aligned peaks the compound appeared in the hit list for.
		"""

		return numpy.count_nonzero(~numpy.isnan(self.hit_numbers))

	def __repr__(self) -> str:
		return f"<Consolidated Search Result: {self.name} \tmf={self.match_factor}\tn={len(self)}>"

	def __str__(self) -> str:
		return self.__repr__()

	def to_dict(self) -> Dict[str, Any]:
		"""
		Returns a dictionary representation of this :class:`~.ConsolidatedSearchResult`.

		All keys are native, JSON-serializable, Python objects.
		"""

		if self.reference_data is None:
			reference_data_as_dict = None
		else:
			reference_data_as_dict = self.reference_data.to_dict()

		return {
				"name": self.name,
				"cas": self.cas,
				"mf_list": self.mf_list,
				"rmf_list": self.rmf_list,
				"hit_numbers": self.hit_numbers,
				"reference_data": reference_data_as_dict,
				}

	@classmethod
	def from_dict(cls: Type["ConsolidatedSearchResult"], d: Dict[str, Any]) -> "ConsolidatedSearchResult":
		"""
		Construct a :class:`~.ConsolidatedSearchResult` from a dictionary.

		:param d:
		"""

		if d["reference_data"] is None:
			reference_data = None
		else:
			reference_data = ReferenceData.from_dict(d["reference_data"])

		return cls(
				name=d["name"],
				cas=d["cas"],
				mf_list=d["mf_list"],
				rmf_list=d["rmf_list"],
				hit_numbers=d["hit_numbers"],
				reference_data=reference_data,
				)


def _attrs_convert_hits(hits: Optional[List[ConsolidatedSearchResult]]) -> List[ConsolidatedSearchResult]:
	if hits is None:
		hits = []

	hits = list(hits)
	for hit in hits:
		if not isinstance(hit, ConsolidatedSearchResult):
			raise TypeError("'hits' must be a list of ConsolidatedSearchResult objects")

	return hits


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


@attr.define(init=False)
class ConsolidatedPeak:
	"""
	A Peak that has been produced by consolidating the properties and search results of several qualified peaks.

	:param rt_list: List of retention times of the aligned peaks.
	:param area_list: List of peak areas for the aligned peaks.
	:param ms_list: List of mass spectra for the aligned peaks.
	:param minutes: Retention time units flag.
		If :py:obj:`True`, retention time is in minutes;
		if :py:obj:`False` retention time is in seconds
	:param hits: Optional list of possible identities for this peak.
	:param meta: Optional dictionary for storing e.g. peak number or whether the peak should be hidden.
	:param ms_comparison: Mapping or Pandas :class:`~pandas.Series` giving pairwise mass spectral comparison scores.
	"""

	#: List of retention times of the aligned peaks.
	rt_list: List[float] = attr.field(converter=list)

	#: List of peak areas for the aligned peaks.
	area_list: List[float] = attr.field(converter=list)

	#: List of mass spectra for the aligned peaks.
	ms_list: List[MassSpectrum] = attr.field(converter=list)

	#: Optional list of possible identities for this peak.
	hits: List[ConsolidatedSearchResult] = attr.field(converter=_attrs_convert_hits)

	#: Pairwise mass spectral comparison scores.
	ms_comparison: pandas.Series = attr.field(converter=_attrs_convert_ms_comparison)

	#: Optional dictionary for storing e.g. peak number or whether the peak should be hidden.
	meta: Dict[str, Any]

	def __init__(
			self,
			rt_list: List[float],
			area_list: List[float],
			ms_list: List[MassSpectrum],
			*,
			minutes: bool = False,
			hits: Optional[List[ConsolidatedSearchResult]] = None,
			ms_comparison: Union[Mapping[str, float], pandas.Series, None] = None,
			meta: Optional[Dict[str, Any]] = None
			):

		# TODO: Type check rt_list and ms_list
		# 		if not isinstance(rt, (int, float)):
		# 			raise TypeError("'rt' must be a number")
		#
		# 			if ms and not isinstance(ms, MassSpectrum) and not isinstance(ms, (int, float)):
		# 				raise TypeError("'ms' must be a number or a MassSpectrum object")
		#

		if minutes:
			rt_list = [rt * 60.0 for rt in rt_list]

		self.rt_list = rt_list
		self.area_list = area_list
		self.ms_list = ms_list
		self.meta = meta or {}
		self.ms_comparison = ms_comparison

		if hits is None:
			self.hits = []
		else:
			self.hits = hits

	@property
	def rt(self) -> float:
		"""
		The average retention time across the aligned peaks.
		"""

		return numpy.nanmean(self.rt_list)  # type: ignore[return-value]  # mypy thinks type is "floating[Any]"

	@property
	def rt_stdev(self) -> float:
		"""
		The standard deviation of the retention time across the aligned peaks.
		"""

		return numpy.nanstd(self.rt_list)  # type: ignore[return-value]  # mypy thinks type is "floating[Any]"

	@property
	def area(self) -> float:
		"""
		The average peak area across the aligned peaks.
		"""

		return numpy.nanmean(self.area_list)  # type: ignore[return-value]  # mypy thinks type is "floating[Any]"

	@property
	def area_stdev(self) -> float:
		"""
		The standard deviation of the peak area across the aligned peaks.
		"""

		return numpy.nanstd(self.area_list)  # type: ignore[return-value]  # mypy thinks type is "floating[Any]"

	@property
	def average_ms_comparison(self) -> float:
		"""
		The average of the pairwise mass spectral comparison scores.
		"""

		# return self._average_ms_comparison

		if self.ms_comparison.empty:
			return 0
		else:
			return numpy.nanmean(self.ms_comparison)

	@property
	def ms_comparison_stdev(self) -> float:
		"""
		The standard deviation of the pairwise mass spectral comparison scores.
		"""

		# return self._ms_comparison_stdev

		if self.ms_comparison.empty:
			return 0
		else:
			return numpy.nanstd(self.ms_comparison)

	# def _calculate_spectra(self):
	# 	"""
	# 	Calculate Combined and Averaged spectra
	# 	"""

	# 	mass_lists = []
	# 	intensity_lists = []

	# 	for spec in self._ms_list:

	# 		if spec:
	# 			# print(spec.mass_list)
	# 			mass_lists.append(spec.mass_list)
	# 			intensity_lists.append(spec.intensity_list)
	# 		else:
	# 			# print()
	# 			pass

	# 	if all_equal(mass_lists):
	# 		mass_list = mass_lists[0]
	# 		# print(intensity_lists)
	# 		combined_intensity_list = list(sum(map(numpy.array, intensity_lists)))
	# 		self._combined_mass_spectrum = MassSpectrum(
	# 				mass_list=mass_list, intensity_list=combined_intensity_list
	# 				)

	# 		# averaged_intensity_list = [intensity / len(mass_lists) for intensity in combined_intensity_list]

	# 		averaged_intensity_list = []
	# 		avg_intensity_array = numpy.array(intensity_lists)
	# 		for column in avg_intensity_array.T:
	# 			if sum(column) == 0 or numpy.count_nonzero(column) == 0:
	# 				averaged_intensity_list.append(0)
	# 			else:
	# 				averaged_intensity_list.append(sum(column) / numpy.count_nonzero(column))

	# 		self._averaged_mass_spectrum = MassSpectrum(
	# 				mass_list=mass_list, intensity_list=averaged_intensity_list
	# 				)

	# 	else:
	# 		warnings.warn("Mass Ranges Differ. Unable to process")
	# 		self._combined_mass_spectrum = None
	# 		self._averaged_mass_spectrum = None

	# @property
	# def combined_mass_spectrum(self):
	# 	return self._combined_mass_spectrum

	# @property
	# def averaged_mass_spectrum(self):
	# 	return self._averaged_mass_spectrum

	def __repr__(self) -> str:
		return f"<Consolidated Peak: {self.rt}>"

	def __str__(self) -> str:
		return self.__repr__()

	# def __eq__(self, other):
	# 	"""
	# 	Return whether this ConsolidatedPeak object is equal to another object

	# 	:param other: The other object to test equality with
	# 	:type other: object

	# 	:rtype: bool
	# 	"""

	# 	if isinstance(other, self.__class__):
	# 		if self.rt_list == other.rt_list and self.area_list == other.area_list:
	# 			return self._ms_list == other._ms_list
	# 		#: TODO: compare hits, meta and ms_comparison
	# 		return False
	# 	else:
	# 		return NotImplemented

	def to_dict(self) -> Dict[str, Any]:
		"""
		Returns a dictionary representation of this :class:`~.ConsolidatedPeak`.

		All keys are native, JSON-serializable, Python objects.
		"""

		return {
				"rt_list": self.rt_list,
				"area_list": self.area_list,
				"meta": self.meta,
				"hits": [hit.to_dict() for hit in self.hits],
				"ms_list": [dict(ms) if ms else None for ms in self.ms_list],
				"ms_comparison": self.ms_comparison.to_dict(),
				}

	def __iter__(self) -> Iterator[Tuple[str, Any]]:
		yield from self.to_dict().items()

	def __len__(self) -> int:
		return numpy.count_nonzero(~numpy.isnan(self.rt_list))

	@classmethod
	def from_dict(cls: Type["ConsolidatedPeak"], d: Dict[str, Any]) -> "ConsolidatedPeak":
		"""
		Construct a :class:`~.ConsolidatedPeak` from a dictionary.

		:param d:
		"""

		hits = []
		for hit in d["hits"]:
			hits.append(ConsolidatedSearchResult.from_dict(hit))

		return cls(
				rt_list=d["rt_list"],
				area_list=d["area_list"],
				ms_list=d["ms_list"],
				meta=d["meta"],
				ms_comparison=d["ms_comparison"],
				hits=hits,
				)
