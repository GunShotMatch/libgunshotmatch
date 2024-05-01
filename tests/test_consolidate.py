# stdlib
from operator import attrgetter
from typing import Optional

# 3rd party
import pyms_nist_search
import pytest
from coincidence.regressions import AdvancedDataRegressionFixture, AdvancedFileRegressionFixture
from domdf_python_tools.paths import PathPlus
from domdf_python_tools.typing import PathLike
from pytest_regressions.dataframe_regression import DataFrameRegressionFixture

# this package
from libgunshotmatch.consolidate import (
		ConsolidatedPeakFilter,
		InvertedFilter,
		combine_spectra,
		pairwise_ms_comparisons
		)
from libgunshotmatch.consolidate._fields import _attrs_convert_reference_data
from libgunshotmatch.project import Project

# Test consolidate process from gsmp file


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
		advanced_file_regression: AdvancedFileRegressionFixture,
		dataframe_regression: DataFrameRegressionFixture,
		capsys,
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

	advanced_file_regression.check(capsys.readouterr().out)

	assert project.consolidated_peaks is not None

	consolidated_peaks = []
	for cp in project.consolidated_peaks:
		cp_as_dict = cp.to_dict()
		cp_as_dict["ms_comparison"] = cp.ms_comparison.astype(int).to_dict()
		consolidated_peaks.append(cp_as_dict)
	advanced_data_regression.check(consolidated_peaks)
	dataframe_regression.check(ms_comparison_df, basename="test_consolidate_ms_comparison_df")


def test_inverted_consolidate(
		advanced_data_regression: AdvancedDataRegressionFixture,
		advanced_file_regression: AdvancedFileRegressionFixture,
		dataframe_regression: DataFrameRegressionFixture,
		capsys,
		):
	project_file = PathPlus(__file__).parent / "ELEY .22LR.gsmp"
	project = Project.from_file(project_file)
	# advanced_data_regression.check(sdjson.dumps(project.to_dict(), indent=2))

	cp_filter = InvertedFilter(
			name_filter=["*silane*", "*silyl*", "*siloxy*"],
			min_match_factor=600,
			min_appearances=5,
			verbose=True,
			)

	ms_comparison_df = project.consolidate(MockEngine(''), cp_filter).astype(int)

	advanced_file_regression.check(capsys.readouterr().out)

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


def test__attrs_convert_reference_data():
	project = Project.from_file("tests/Eley Hymax.gsmp")
	assert project.consolidated_peaks is not None

	peaks = sorted(project.consolidated_peaks, key=attrgetter("area"), reverse=True)
	largest_peak = peaks[0]
	reference_data = largest_peak.hits[0].reference_data
	assert reference_data is not None
	via_converter = _attrs_convert_reference_data(reference_data.to_dict())
	assert via_converter is not None

	assert reference_data.name == via_converter.name
	assert reference_data.cas == via_converter.cas
	assert reference_data.formula == via_converter.formula
	assert reference_data.contributor == via_converter.contributor
	assert reference_data.nist_no == via_converter.nist_no
	assert reference_data.id == via_converter.id
	assert reference_data.mw == via_converter.mw
	assert reference_data.exact_mass == via_converter.exact_mass
	assert reference_data.synonyms == via_converter.synonyms
	assert reference_data.mass_spec == via_converter.mass_spec

	assert _attrs_convert_reference_data(reference_data) is reference_data

	with pytest.raises(TypeError, match="'reference_data' must be a `pyms_nist_search.ReferenceData` object,"):
		# Wrong keys
		_attrs_convert_reference_data({})

	with pytest.raises(TypeError, match="'reference_data' must be a `pyms_nist_search.ReferenceData` object,"):
		# Wrong type
		_attrs_convert_reference_data([])  # type: ignore[arg-type]


def test_combine_spectra(advanced_data_regression: AdvancedDataRegressionFixture):
	project = Project.from_file("tests/Eley Super Game.gsmp")
	assert project.consolidated_peaks is not None

	spectra = []
	for peak in project.consolidated_peaks:
		spectra.append(combine_spectra(peak))
	advanced_data_regression.check(spectra)
