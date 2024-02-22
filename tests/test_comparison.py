# this package
from libgunshotmatch import comparison
from libgunshotmatch.project import Project


def test_align_two_projects_and_unknown():
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
	# TODO: check output and capsys.readouterr
