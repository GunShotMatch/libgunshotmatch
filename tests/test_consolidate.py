# stdlib
from typing import Optional

# 3rd party
import numpy
import pyms_nist_search
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
		numpy.float64,
		)
def _represent_mappings(dumper: RegressionYamlDumper, data: int):  # noqa: MAN002
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
		advanced_data_regression: AdvancedDataRegressionFixture, dataframe_regression: DataFrameRegressionFixture
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


def test_pairwise_ms_comparison(dataframe_regression: DataFrameRegressionFixture):
	project_file = PathPlus(__file__).parent / "ELEY .22LR.gsmp"
	project = Project.from_file(project_file)

	ms_comparison_df = pairwise_ms_comparisons(project.alignment, parallel=False).astype(int)
	dataframe_regression.check(ms_comparison_df)

	ms_comparison_df = pairwise_ms_comparisons(project.alignment).astype(int)
	dataframe_regression.check(ms_comparison_df)


# TODO: creation etc.
