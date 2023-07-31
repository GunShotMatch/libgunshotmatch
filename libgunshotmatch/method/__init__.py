#!/usr/bin/env python3
#
#  __init__.py
"""
Methods for GunShotMatch analysis.
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

try:
	# 3rd party
	import tomllib  # type: ignore[import]
except ImportError:
	import tomli as tomllib  # type: ignore[no-redef]

# stdlib
import json
from typing import Any, Dict, List, Mapping, Optional, Set, Tuple, Type, TypeVar

# 3rd party
import attr
import tomli_w

# this package
from libgunshotmatch.method._fields import (
		Boolean,
		Integer,
		Number,
		String,
		convert_crop_mass_range,
		convert_rt_range,
		default_base_peak_filter
		)
from libgunshotmatch.utils import _fix_init_annotations, _to_list

__all__ = [
		"MethodBase",
		"Method",
		"IntensityMatrixMethod",
		"PeakDetectionMethod",
		"PeakFilterMethod",
		"AlignmentMethod",
		"ConsolidateMethod",
		]

_MB = TypeVar("_MB", bound="MethodBase")


class MethodBase:
	"""
	Base class for methods.
	"""

	@classmethod
	def from_dict(cls: Type["MethodBase"], method: Mapping[str, Any]) -> "MethodBase":
		"""
		Construct a :class:`~.MethodBase` from a dictionary.

		:param method:
		"""

		return cls(**method)

	def to_dict(self) -> Dict[str, Any]:
		"""
		Convert a :class:`~.MethodBase` to a dictionary.
		"""

		return attr.asdict(self, recurse=True)  # type: ignore[arg-type]  # This is an ABC but not an attrs class.

	@classmethod
	def _coerce(cls: Type[_MB], method: Any) -> _MB:
		if isinstance(method, cls):
			return method
		elif isinstance(method, Mapping):
			return cls(**method)
		else:

			class_name = cls.__name__
			if class_name.startswith("aeiouAEIOU"):
				raise TypeError(f"Cannot convert {type(method).__name__} to an {class_name}")
			# TODO: edge cases
			else:
				raise TypeError(f"Cannot convert {type(method).__name__} to a {class_name}")


@_fix_init_annotations
@attr.define
class IntensityMatrixMethod(MethodBase):
	"""
	Method used for constructing an intensity matrix from a datafile.
	"""

	#: The range of masses to which the GC-MS data should be limited to.
	crop_mass_range: Optional[Tuple[int, int]] = attr.field(default=(50, 500), converter=convert_crop_mass_range)

	#: Whether to perform Savitzky-Golay smoothing.
	savitzky_golay: bool = Boolean.field(default=True)

	#: Whether to perform Tophat baseline correction.
	tophat: bool = Boolean.field(default=True)

	#: The structure size for Tophat baseline correction.
	tophat_structure_size: str = String.field(default="1.5m")


@_fix_init_annotations
@attr.define
class PeakDetectionMethod(MethodBase):
	"""
	Method used for Biller-Biemann peak detection.
	"""

	#: Number of scans over which to consider a maxima to be a peak.
	points: int = Integer.field(default=10)

	#: Number of scans to combine peaks from to compensate for spectra skewing.
	scans: int = Integer.field(default=1)


@_fix_init_annotations
@attr.define
class PeakFilterMethod(MethodBase):
	"""
	Method used for peak filtering.
	"""

	#: Whether to perform automatic noise filtering of the peak list.
	noise_filter: bool = Boolean.field(default=True)

	#: The minimum number of ions that must have intensities above the noise floor, otherwise the peak is excluded.
	noise_threshold: int = Integer.field(default=2)

	# TODO: non-integer binned data
	#: Peaks whose base peak is at one of the listed masses (m/z) are excluded.
	base_peak_filter: Set[int] = attr.field(
			default=attr.Factory(default_base_peak_filter),
			converter=set,
			validator=attr.validators.instance_of(set)
			)

	#: Optional retention time range (in minutes) to filter the peak list to.
	rt_range: Optional[Tuple[float, float]] = attr.field(default=None, converter=convert_rt_range)


@_fix_init_annotations
@attr.define
class AlignmentMethod(MethodBase):
	"""
	Method used for peak alignment.
	"""

	#: Retention time tolerance parameter for pairwise alignments.
	rt_modulation: bool = Number.field(default=2.5)

	#: Gap parameter for pairwise alignments.
	gap_penalty: float = Number.field(default=0.3)

	#: Minimum number of peaks required for the alignment position to survive filtering.
	min_peaks: int = Integer.field(default=1)

	#: Number of peaks (starting with the largest) to include in the output.
	top_n_peaks: int = Integer.field(default=80)

	#: Minimum area of peaks to include in the output.
	min_peak_area: float = Number.field(default=0.0)


@_fix_init_annotations
@attr.define
class ConsolidateMethod(MethodBase):
	"""
	Method used for consolidation (finding most likely identity for aligned peaks).

	:param min_appearances: Number of times the hit must appear across individual the aligned peaks.
		Consolidated peaks where the most common hit appears fewer times than this will be excluded.
		If set to ``-1`` this number of repeats in the project are used.

	.. versionchanged:: 0.2.0  Added the ``min_appearances`` argument.
	"""

	name_filter: List[str] = attr.field(converter=_to_list, default=attr.Factory(list))
	"""
	List of glob-style matches for compound names.

	Consolidated peaks matching any of these will be excluded.
	"""

	min_match_factor: int = Integer.field(default=600)
	"""
	Minimum average match factor.

	Consolidated peaks with an average match factor below this will be excluded.
	"""

	min_appearances: int = Integer.field(default=-1)
	"""
	Number of times the hit must appear across individual the aligned peaks.

	Consolidated peaks where the most common hit appears fewer times than this will be excluded.

	If set to ``-1`` this number of repeats in the project are used.

	.. versionadded:: 0.2.0
	"""


# target_range = 4.0,37.0


def _submethod_field(submethod_type: Type[_MB]) -> _MB:
	# Actually returns attr.Attribute, but mypy doesn't like it
	return attr.field(default=attr.Factory(submethod_type), converter=submethod_type._coerce)


@_fix_init_annotations
@attr.define
class Method(MethodBase):
	"""
	Overall GunShotMatch method.
	"""

	#: Method used for constructing an intensity matrix from a datafile.
	intensity_matrix: IntensityMatrixMethod = _submethod_field(IntensityMatrixMethod)

	#: Method used for Biller-Biemann peak detection.
	peak_detection: PeakDetectionMethod = _submethod_field(PeakDetectionMethod)

	#: Method used for peak filtering.
	peak_filter: PeakFilterMethod = _submethod_field(PeakFilterMethod)

	#: Method used for peak alignment.
	alignment: AlignmentMethod = _submethod_field(AlignmentMethod)

	#: Method used for consolidation (finding most likely identity for aligned peaks).
	consolidate: ConsolidateMethod = _submethod_field(ConsolidateMethod)

	@classmethod
	def from_toml(cls: Type["Method"], toml_string: str) -> "Method":
		"""
		Parse a :class:`~.Method` from a TOML string.

		:param toml_string:
		"""

		parsed_toml = tomllib.loads(toml_string)
		return cls(**parsed_toml["method"])

	@classmethod
	def from_json(cls: Type["Method"], json_string: str) -> "Method":
		"""
		Parse a :class:`~.Method` from a JSON string.

		:param json_string:
		"""

		parsed_json = json.loads(json_string)
		return cls(**parsed_json["method"])

	def to_toml(self) -> str:
		"""
		Convert a :class:`~.Method` to a TOML string.
		"""

		return tomli_w.dumps({"method": self.to_dict()})
