#!/usr/bin/env python3
#
#  utils.py
"""
Utility functions.
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
from decimal import Decimal
from typing import Any, Iterable, List, Optional, Sequence, Type, TypeVar, Union

# 3rd party
import numpy
from chemistry_tools.spectrum_similarity import SpectrumSimilarity
from mathematical.utils import rounders
from pyms.Spectrum import MassSpectrum  # type: ignore[import]
from scipy.stats import truncnorm  # type: ignore[import]

__all__ = ("round_rt", "get_truncated_normal", "ms_comparison")


def round_rt(rt: Union[str, float, Decimal]) -> Decimal:
	"""
	Truncate precision of retention time to 10 decimal places.

	:param rt:
	"""

	# Limit to 10 decimal places as that's what Pandas writes JSON data as;
	# no need for greater precision.

	return rounders(rt, "0.0000000000")


def get_truncated_normal(
		mean: float,
		sd: float,
		low: float = 0,
		upp: float = 10,
		count: int = 10,
		random_state: Optional[int] = None,
		) -> Sequence[float]:
	"""
	Returns ``count`` values from a truncated normal distrubition.

	:param mean: The midpoint of the normal distribution.
	:param sd: The spread of the normal distribution (the standard deviation).
	:param low: The lower bound.
	:param upp: The upper bound.
	:param count:
	:param random_state: Optional seed for the random number generator.
	"""

	# From https://stackoverflow.com/a/74448424
	# By toco_tico https://stackoverflow.com/users/1060349/toto-tico
	# CC BY-SA 4.0

	dist = truncnorm((low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)
	return dist.rvs(count, random_state=random_state)


def ms_comparison(top_ms: MassSpectrum, bottom_ms: MassSpectrum) -> Optional[float]:
	"""
	Performs a Mass Spectrum similarity calculation two mass spectra.

	:param top_ms:
	:param bottom_ms:

	If either of ``top_ms`` or ``bottom_ms`` is :py:obj:`None` then :py:obj:`None` is returned,
	otherwise a comparison score is returned.
	"""

	if top_ms is None or bottom_ms is None:
		return None

	top_spec = numpy.column_stack((top_ms.mass_list, top_ms.mass_spec))
	bottom_spec = numpy.column_stack((bottom_ms.mass_list, bottom_ms.mass_spec))

	sim = SpectrumSimilarity(
			top_spec,
			bottom_spec,
			b=1,
			xlim=(45, 500),  # TODO: configurable or taken from spectra
			)

	match, rmatch = sim.score()
	return match * 1000


_O = TypeVar("_O", bound=object)


def _fix_init_annotations(method: Type[_O]) -> Type[_O]:
	init_annotations = method.__init__.__annotations__
	cls_annotations = method.__annotations__

	for k, v in cls_annotations.items():
		if k in init_annotations:
			if init_annotations[k] is Any:
				init_annotations[k] = v
		else:
			init_annotations[k] = v

	return method


def _to_list(l: Iterable[str]) -> List[str]:
	"""
	Attrs type hint helper for converting to a list.

	Otherwise the errors are:

	libgunshotmatch/consolidate/__init__.py:701: error: Argument "name_filter" to "ConsolidatedPeakFilter" has incompatible type "List[str]"; expected "Iterable[_T]"  [arg-type]
	libgunshotmatch/project.py:202: error: List item 0 has incompatible type "str"; expected "_T"  [list-item]
	"""

	return list(l)
