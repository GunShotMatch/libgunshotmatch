# 3rd party
import numpy
from coincidence.regressions import _representer_for
from pytest_regressions.data_regression import RegressionYamlDumper
from yaml.representer import RepresenterError

pytest_plugins = ("coincidence", )


def represent_undefined(self, data):  # noqa: MAN001,MAN002
	raise RepresenterError("cannot represent an object", data, type(data))


RegressionYamlDumper.represent_undefined = represent_undefined  # type: ignore[method-assign]


@_representer_for(
		numpy.int64,
		numpy.int32,
		numpy.float64,
		)
def _represent_numpy(dumper: RegressionYamlDumper, data: int):  # noqa: MAN002
	return dumper.represent_data(int(data))


# @_representer_for(
# 		numpy.float64,
# 		)
# def _represent_mappings(dumper: RegressionYamlDumper, data: float):  # noqa: MAN002
# 	if data == 0:
# 		return dumper.represent_data(0)
# 	return dumper.represent_data(str(round_rt(data)))
