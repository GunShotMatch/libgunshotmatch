# stdlib
from typing import Any, Dict, Type

# 3rd party
import pytest
from coincidence.regressions import AdvancedDataRegressionFixture

# this package
from libgunshotmatch.method import (
		AlignmentMethod,
		ConsolidateMethod,
		IntensityMatrixMethod,
		Method,
		MethodBase,
		PeakDetectionMethod,
		PeakFilterMethod,
		SavitzkyGolayMethod
		)


@pytest.mark.parametrize(
		"cls",
		[
				pytest.param(AlignmentMethod, id="AlignmentMethod"),
				pytest.param(PeakFilterMethod, id="PeakFilterMethod"),
				pytest.param(ConsolidateMethod, id="ConsolidateMethod"),
				pytest.param(PeakDetectionMethod, id="PeakDetectionMethod"),
				pytest.param(SavitzkyGolayMethod, id="SavitzkyGolayMethod"),
				pytest.param(IntensityMatrixMethod, id="IntensityMatrixMethod"),
				pytest.param(Method, id="Method"),
				]
		)
def test_default_methods(advanced_data_regression: AdvancedDataRegressionFixture, cls: Type[MethodBase]):
	method = cls()
	advanced_data_regression.check(method.to_dict())


@pytest.mark.parametrize(
		"cls, kwargs",
		[
				pytest.param(SavitzkyGolayMethod, {"enable": False}, id="savgol_false"),
				pytest.param(SavitzkyGolayMethod, {"window": 8}, id="savgol_window"),
				pytest.param(SavitzkyGolayMethod, {"window": "1m"}, id="savgol_window_str_min"),
				pytest.param(SavitzkyGolayMethod, {"window": "20s"}, id="savgol_window_str_sec"),
				pytest.param(SavitzkyGolayMethod, {"degree": 5}, id="savgol_degree"),
				pytest.param(SavitzkyGolayMethod, {"window": 8, "degree": 5}, id="savgol_window_degree"),
				pytest.param(IntensityMatrixMethod, {"crop_mass_range": [100, 200]}, id="im_mass_range"),
				pytest.param(
						IntensityMatrixMethod, {"tophat_structure_size": "2m"}, id="im_tophat_structure_size"
						),
				pytest.param(
						IntensityMatrixMethod, {"savitzky_golay": {"enable": False}}, id="im_savgol_dict_false"
						),
				pytest.param(
						IntensityMatrixMethod, {"savitzky_golay": {"window": 10}}, id="im_savgol_dict_window"
						),
				pytest.param(IntensityMatrixMethod, {"savitzky_golay": False}, id="im_savgol_false"),
				pytest.param(PeakDetectionMethod, {"points": 8}, id="peak_detection_points"),
				pytest.param(PeakDetectionMethod, {"scans": 3}, id="peak_detection_scans"),
				pytest.param(PeakFilterMethod, {"noise_filter": False}, id="peak_filter_noise_filter"),
				pytest.param(PeakFilterMethod, {"noise_threshold": 5}, id="peak_filter_noise_threshold"),
				pytest.param(
						PeakFilterMethod, {"base_peak_filter": [1, 2, 3, 4, 5]},
						id="peak_filter_base_peak_filter_list"
						),
				pytest.param(
						PeakFilterMethod, {"base_peak_filter": {1, 2, 3, 4, 5}},
						id="peak_filter_base_peak_filter_set"
						),
				pytest.param(PeakFilterMethod, {"rt_range": [5, 10]}, id="peak_filter_rt_range_list"),
				pytest.param(PeakFilterMethod, {"rt_range": (15, 20)}, id="peak_filter_rt_range_tuple"),
				pytest.param(AlignmentMethod, {"rt_modulation": 10}, id="alignment_rt_modulation"),
				pytest.param(AlignmentMethod, {"gap_penalty": 0.2}, id="alignment_gap_penalty"),
				pytest.param(AlignmentMethod, {"min_peaks": 5}, id="alignment_min_peaks"),
				pytest.param(AlignmentMethod, {"top_n_peaks": 50}, id="alignment_top_n_peaks"),
				pytest.param(AlignmentMethod, {"min_peak_area": 1.5e6}, id="alignment_min_peak_area"),
				# ConsolidateMethod
				# Method
				]
		)
def test_partial_arguments(
		advanced_data_regression: AdvancedDataRegressionFixture, cls: Type[MethodBase], kwargs: Dict[str, Any]
		):
	method = cls(**kwargs)
	advanced_data_regression.check(method.to_dict())

	method = cls.from_dict(kwargs)
	advanced_data_regression.check(method.to_dict())


# TODO: test errors
# TODO: test from dict and to/from toml
