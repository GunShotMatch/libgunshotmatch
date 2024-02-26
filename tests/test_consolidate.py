# stdlib
from operator import attrgetter
from typing import Optional

# 3rd party
import numpy
import pyms_nist_search
import pytest
from coincidence.regressions import AdvancedDataRegressionFixture, _representer_for
from domdf_python_tools.paths import PathPlus
from domdf_python_tools.typing import PathLike
from pytest_regressions.data_regression import RegressionYamlDumper
from pytest_regressions.dataframe_regression import DataFrameRegressionFixture

# this package
from libgunshotmatch.consolidate import ConsolidatedPeakFilter, pairwise_ms_comparisons
from libgunshotmatch.project import Project

# Test consolidate process from gsmp file


@_representer_for(
		numpy.int64,
		numpy.int32,
		numpy.float64,
		)
def _represent_numpy(dumper: RegressionYamlDumper, data: int):  # noqa: MAN002
	return dumper.represent_data(int(data))


class MockEngine(pyms_nist_search.Engine):
	"""
	Engine that returns :py:obj:`None` for the reference data.
	"""

	def __init__(
			self,
			lib_path: PathLike,
			lib_type: int = ...,
			work_dir: Optional[PathLike] = None,
			debug: bool = False,
			):
		pass

	def get_reference_data(self, spec_loc: int) -> pyms_nist_search.ReferenceData:
		return None  # type: ignore[return-value]


def test_consolidate(
		advanced_data_regression: AdvancedDataRegressionFixture,
		dataframe_regression: DataFrameRegressionFixture,
		):
	project_file = PathPlus(__file__).parent / "ELEY .22LR.gsmp"
	project = Project.from_file(project_file)
	# advanced_data_regression.check(sdjson.dumps(project.to_dict(), indent=2))

	cp_filter = ConsolidatedPeakFilter(
			name_filter=["*silane*", "*silyl*", "*siloxy*"],
			min_match_factor=600,
			min_appearances=5,
			verbose=True,
			)

	ms_comparison_df = project.consolidate(MockEngine(''), cp_filter).astype(int)

	assert project.consolidated_peaks is not None

	consolidated_peaks = []
	for cp in project.consolidated_peaks:
		cp_as_dict = cp.to_dict()
		cp_as_dict["ms_comparison"] = cp.ms_comparison.astype(int).to_dict()
		consolidated_peaks.append(cp_as_dict)
	advanced_data_regression.check(consolidated_peaks)
	dataframe_regression.check(ms_comparison_df, basename="test_consolidate_ms_comparison_df")


def test_consolidated_peak_attributes():

	project = Project.from_file("tests/Eley Hymax.gsmp")
	assert project.consolidated_peaks is not None

	peaks = sorted(project.consolidated_peaks, key=attrgetter("area"), reverse=True)
	largest_peak = peaks[0]
	assert largest_peak.rt == 1276.2899780273438
	assert largest_peak.rt_stdev == 0.8302141496232701
	assert largest_peak.area == 7313500780.940475
	assert largest_peak.area_stdev == 509223188.0045349
	assert largest_peak.average_ms_comparison == 999.775888623915
	assert largest_peak.ms_comparison_stdev == 0.16859417266292104

	top_hit = largest_peak.hits[0]
	assert top_hit.match_factor == 927.25
	assert top_hit.match_factor_stdev == 6.139014578904337
	assert top_hit.reverse_match_factor == 929.0
	assert top_hit.reverse_match_factor_stdev == 6.519202405202649
	assert top_hit.average_hit_number == 1.0
	assert len(top_hit) == 4

	seventh_hit = largest_peak.hits[6]
	assert seventh_hit.match_factor == 643.75
	assert seventh_hit.match_factor_stdev == 3.491060010942235
	assert seventh_hit.reverse_match_factor == 652.5
	assert seventh_hit.reverse_match_factor_stdev == 3.640054944640259
	assert seventh_hit.average_hit_number == 7.75
	assert seventh_hit.hit_number_stdev == 0.82915619758885
	assert len(seventh_hit) == 4

	bottom_hit = largest_peak.hits[-1]
	assert bottom_hit.match_factor == 626.0
	assert bottom_hit.match_factor_stdev == 0.0
	assert bottom_hit.reverse_match_factor == 647.0
	assert bottom_hit.reverse_match_factor_stdev == 0.0
	assert bottom_hit.average_hit_number == 10.0
	assert bottom_hit.hit_number_stdev == 0.0
	assert len(bottom_hit) == 1


def test_pairwise_ms_comparison(dataframe_regression: DataFrameRegressionFixture):
	project_file = PathPlus(__file__).parent / "ELEY .22LR.gsmp"
	project = Project.from_file(project_file)

	ms_comparison_df = pairwise_ms_comparisons(project.alignment, parallel=False).astype(int)
	dataframe_regression.check(ms_comparison_df)

	ms_comparison_df = pairwise_ms_comparisons(project.alignment).astype(int)
	dataframe_regression.check(ms_comparison_df)


# TODO: creation etc.
