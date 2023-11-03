# 3rd party
from domdf_python_tools.paths import PathPlus

# this package
from libgunshotmatch.gzip_util import read_gzip_json, write_gzip_json


def test_roundtrip(tmp_pathplus: PathPlus):
	data = {"hello": "world", "number": 1234, "list": [1, 2, 3, 'a', 'b', 'c']}

	write_gzip_json(tmp_pathplus / "test.json.gz", data)
	assert read_gzip_json(tmp_pathplus / "test.json.gz") == data

	write_gzip_json(tmp_pathplus / "test.gsmp", data)
	assert read_gzip_json(tmp_pathplus / "test.gsmp") == data
