# 3rd party
from pytest_regressions.data_regression import RegressionYamlDumper
from yaml.representer import RepresenterError

pytest_plugins = ("coincidence", )


def represent_undefined(self, data):  # noqa: MAN001,MAN002
	raise RepresenterError("cannot represent an object", data, type(data))


RegressionYamlDumper.represent_undefined = represent_undefined  # type: ignore[method-assign]
