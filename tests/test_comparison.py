# 3rd party
from coincidence.regressions import AdvancedDataRegressionFixture, AdvancedFileRegressionFixture
from pytest_regressions.dataframe_regression import DataFrameRegressionFixture

# this package
from libgunshotmatch import comparison
from libgunshotmatch.project import Project


def test_align_two_projects_and_unknown(
		capsys,
		advanced_file_regression: AdvancedFileRegressionFixture,
		dataframe_regression: DataFrameRegressionFixture,
		):
	project1 = Project.from_file("tests/Eley Super Game.gsmp")
	project2 = Project.from_file("tests/Eley Hymax.gsmp")
	unknown = Project.from_file("tests/Eley Super Game 2.gsmp")

	print("Loaded")

	A1 = comparison.align_projects([project1, project2], [unknown])
	(padded_p1_cp, padded_p2_cp), (padded_unkn_cp, ) = comparison.get_padded_peak_lists(
			A1,
			[project1, project2],
			[unknown],
			)

	advanced_file_regression.check(capsys.readouterr().out)
	dataframe_regression.check(A1.get_peak_alignment())


def test_align_project_and_unknown(
		capsys,
		advanced_file_regression: AdvancedFileRegressionFixture,
		advanced_data_regression: AdvancedDataRegressionFixture,
		dataframe_regression: DataFrameRegressionFixture,
		):
	project1 = Project.from_file("tests/Eley Super Game.gsmp")
	unknown = Project.from_file("tests/Eley Super Game 2.gsmp")

	print("Loaded")

	A1 = comparison.align_projects(project1, unknown)
	(padded_p1_cp, ), (padded_unkn_cp, ) = comparison.get_padded_peak_lists(
			A1,
			project1,
			unknown,
			)

	advanced_file_regression.check(capsys.readouterr().out)
	dataframe_regression.check(A1.get_peak_alignment())

	paired_retention_times = []
	for p_cp, u_cp in zip(padded_p1_cp, padded_unkn_cp):
		p_rt = p_cp.rt if p_cp else None
		u_rt = u_cp.rt if u_cp else None
		paired_retention_times.append((p_rt, u_rt))
	advanced_data_regression.check(paired_retention_times)
