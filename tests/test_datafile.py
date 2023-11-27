# stdlib
import datetime
import os

# 3rd party
import pytest
import sdjson
from coincidence.regressions import AdvancedDataRegressionFixture, AdvancedFileRegressionFixture
from domdf_python_tools.paths import PathPlus
from pyms.BillerBiemann import BillerBiemann  # type: ignore[import]
from pyms.GCMS.Class import GCMS_data  # type: ignore[import]
from pyms.Peak.Function import peak_sum_area  # type: ignore[import]

# this package
from libgunshotmatch.datafile import Datafile, Repeat
from libgunshotmatch.method import Method
from libgunshotmatch.peak import PeakList, filter_peaks


@pytest.mark.parametrize(
		"filename",
		[
				"ELEY_1_SUBTRACT.JDX",
				"ELEY_2_SUBTRACT.JDX",
				"ELEY_3_SUBTRACT.JDX",
				"ELEY_4_SUBTRACT.JDX",
				"ELEY_5_SUBTRACT.JDX",
				]
		)
def test_datafile_from_jdx(
		filename: str,
		advanced_file_regression: AdvancedFileRegressionFixture,
		monkeypatch,
		tmp_pathplus: PathPlus
		):
	monkeypatch.setenv("USERNAME", "test-user")

	method = Method()

	path = PathPlus(__file__).parent / filename
	datafile = Datafile.new(path.stem, path)

	gcms_data: GCMS_data = datafile.load_gcms_data()

	datafile.prepare_intensity_matrix(
			gcms_data,
			savitzky_golay=method.intensity_matrix.savitzky_golay,
			tophat=method.intensity_matrix.tophat,
			tophat_structure_size=method.intensity_matrix.tophat_structure_size,
			crop_mass_range=method.intensity_matrix.crop_mass_range,
			)

	assert datafile.intensity_matrix is not None

	as_dict = {
			"name": datafile.name,
			"original_filename": os.path.split(datafile.original_filename)[1],
			"original_filetype": int(datafile.original_filetype),
			"description": datafile.description,
			"version": datafile.version,
			}

	advanced_file_regression.check(sdjson.dumps(as_dict), extension=".json")

	advanced_file_regression.check(
			sdjson.dumps(datafile.intensity_matrix.time_list, indent=2),
			extension=".im-times.json",
			)
	advanced_file_regression.check(
			sdjson.dumps(datafile.intensity_matrix.mass_list, indent=2),
			extension=".im-masses.json",
			)
	# advanced_data_regression.check([[str(round_rt(intensity)) for intensity in row] for row in datafile.intensity_matrix.intensity_array])

	datafile.export(tmp_pathplus)
	assert (tmp_pathplus / (datafile.name + ".gsmd")).is_file()


@pytest.mark.parametrize(
		"name", [
				"ELEY_1_SUBTRACT",
				"ELEY_2_SUBTRACT",
				"ELEY_3_SUBTRACT",
				"ELEY_4_SUBTRACT",
				"ELEY_5_SUBTRACT",
				]
		)
def test_load_datafile(name: str, advanced_file_regression: AdvancedFileRegressionFixture):
	path = PathPlus(__file__).parent / f"{name}.gsmd"
	datafile = Datafile.from_file(path)
	advanced_file_regression.check(sdjson.dumps(datafile.to_dict()))


def prepare_peak_list(datafile: Datafile, gcms_data: GCMS_data, method: Method) -> PeakList:
	"""
	Construct and filter the peak list.

	:param datafile:
	:param gcms_data:
	:param method:
	"""

	im = datafile.intensity_matrix

	peak_list: PeakList = PeakList(
			BillerBiemann(
					im,
					points=method.peak_detection.points,
					scans=method.peak_detection.scans,
					)
			)
	print(" Peak list before filtering:", peak_list)

	filtered_peak_list: PeakList = filter_peaks(
			peak_list,
			gcms_data.tic,
			noise_filter=method.peak_filter.noise_filter,
			noise_threshold=method.peak_filter.noise_filter,
			base_peak_filter=method.peak_filter.base_peak_filter,
			rt_range=method.peak_filter.rt_range,
			)
	print(" Peak list after filtering:", filtered_peak_list)

	for peak in filtered_peak_list:
		peak.area = peak_sum_area(im, peak)

	filtered_peak_list.datafile_name = datafile.name
	return filtered_peak_list


@pytest.mark.flaky(reruns=3)
def test_peaks(advanced_data_regression: AdvancedDataRegressionFixture, monkeypatch):
	monkeypatch.setenv("USERNAME", "test-user")

	method = Method()
	path = PathPlus(__file__).parent / f"ELEY_4_SUBTRACT.JDX"
	datafile = Datafile.new(path.stem, path)
	gcms_data: GCMS_data = datafile.load_gcms_data()
	datafile.prepare_intensity_matrix(
			gcms_data,
			savitzky_golay=method.intensity_matrix.savitzky_golay,
			tophat=method.intensity_matrix.tophat,
			tophat_structure_size=method.intensity_matrix.tophat_structure_size,
			crop_mass_range=method.intensity_matrix.crop_mass_range,
			)
	peak_list = prepare_peak_list(datafile, gcms_data, method)
	assert peak_list.datafile_name == "ELEY_4_SUBTRACT"
	# TODO: advanced_data_regression.check(peak_list.to_list())
	# repeat = Repeat(datafile, peak_list)


@pytest.mark.parametrize(
		"name", [
				"ELEY_1_SUBTRACT",
				"ELEY_2_SUBTRACT",
				"ELEY_3_SUBTRACT",
				"ELEY_4_SUBTRACT",
				"ELEY_5_SUBTRACT",
				]
		)
def test_create_repeat(name: str, tmp_pathplus: PathPlus):
	path = PathPlus(__file__).parent / f"{name}.gsmd"
	jdx_path = PathPlus(__file__).parent / f"{name}.JDX"
	datafile = Datafile.from_file(path)
	gcms_data: GCMS_data = datafile.load_gcms_data(jdx_path)
	method = Method()
	peak_list = prepare_peak_list(datafile, gcms_data, method)
	repeat = Repeat(
			datafile,
			peak_list,
			user="foo",
			device="bar",
			date_created=datetime.datetime.fromtimestamp(0),
			date_modified=datetime.datetime.fromtimestamp(0),
			)
	repeat.export(tmp_pathplus)
	assert (tmp_pathplus / (datafile.name + ".gsmr")).is_file()
	Repeat.from_file(tmp_pathplus / (datafile.name + ".gsmr"))
