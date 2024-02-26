#!/usr/bin/env python3
#
#  __init__.py
"""
Comparison between projects and unknowns.

The two submodules, :mod:`~.comparison.projects` and :mod:`~.comparison.unknowns`,
provide identical APIs with different internals to handle reference projects (containing two or more repeats)
and unknown samples (from a single datafile).

.. versionadded:: 0.8.0
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
from typing import List, Sequence, Tuple, Union

# 3rd party
from pyms.DPA.Alignment import Alignment
from pyms.DPA.PairwiseAlignment import PairwiseAlignment, align_with_tree

# this package
from libgunshotmatch.project import Project

# this package
from . import projects, unknowns
from ._utils import _PaddedPeakList

__all__ = ("align_projects", "get_padded_peak_lists")

# Aliases to prevent clashes with argument names
_projects_mod = projects
_unknowns_mod = unknowns


def align_projects(
		projects: Union[Sequence[Project], Project] = (),
		unknowns: Union[Sequence[Project], Project] = (),
		D: float = 2.5,
		gap: float = 0.3,
		) -> Alignment:
	"""
	Align multiple projects and/or unknowns.

	:param projects:
	:param unknowns:
	:param D: Retention time tolerance for pairwise alignments (in seconds).
	:param gap: Gap penalty for pairwise alignments.

	:rtype:

	.. versionchanged:: 0.9.0

		* Added ``D`` and ``gap`` arguments.
		* ``projects`` and ``unknowns`` can now be a single :class:`~.Project`.
	"""

	if isinstance(projects, Project):
		projects = [projects]

	if isinstance(unknowns, Project):
		unknowns = [unknowns]

	project_alignments = map(_projects_mod.filter_alignment_to_consolidate, projects)
	unknown_alignments = map(_unknowns_mod.filter_alignment_to_consolidate, unknowns)

	pwa = PairwiseAlignment([*project_alignments, *unknown_alignments], D=float(D), gap=float(gap))
	return align_with_tree(pwa)


def get_padded_peak_lists(
		alignment: Alignment,
		projects: Union[Sequence[Project], Project] = (),
		unknowns: Union[Sequence[Project], Project] = (),
		) -> Tuple[List[_PaddedPeakList], List[_PaddedPeakList]]:
	"""
	Pad the consolidated peak lists in each project/unknown, from the given between-project alignment.

	:param alignment:
	:param projects:
	:param unknowns:

	.. versionchanged:: 0.9.0  ``projects`` and ``unknowns`` can now be a single :class:`~.Project`.
	"""

	if isinstance(projects, Project):
		projects = [projects]

	if isinstance(unknowns, Project):
		unknowns = [unknowns]

	data = alignment.get_peak_alignment(require_all_expr=False, minutes=False)

	projects_padded_cp = [_projects_mod.get_padded_peak_list(p, data) for p in projects]
	unknowns_padded_cp = [_unknowns_mod.get_padded_peak_list(p, data) for p in unknowns]
	return projects_padded_cp, unknowns_padded_cp
