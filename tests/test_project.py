# 3rd party
import sdjson
from coincidence.regressions import AdvancedFileRegressionFixture
from domdf_python_tools.paths import PathPlus

# this package
from libgunshotmatch.project import Project


def test_load(advanced_file_regression: AdvancedFileRegressionFixture):
	project_file = PathPlus(__file__).parent / "ELEY .22LR.gsmp"
	project = Project.from_file(project_file)
	advanced_file_regression.check(sdjson.dumps(project.to_dict()))


# TODO: creation etc.
