#!/usr/bin/env python3
#
#  gzip_util.py
"""
Read and write gzipped JSON.
"""
#
#  Copyright © 2023 Dominic Davis-Foster <dominic@davis-foster.co.uk>
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in all
#  copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
#  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
#  DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
#  OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
#  OR OTHER DEALINGS IN THE SOFTWARE.
#

# stdlib
import gzip
from typing import Any, Dict, List, OrderedDict, Tuple, Union

# 3rd party
import sdjson
from domdf_python_tools.paths import PathPlus
from domdf_python_tools.typing import PathLike

__all__ = ["read_gzip_json", "write_gzip_json"]

JSONOutput = Union[str, int, float, bool, None, Dict[str, Any], List[Any]]
JSONInput = Union[JSONOutput, Tuple[Any, ...], OrderedDict]


def read_gzip_json(path: PathLike) -> JSONOutput:
	"""
	Load JSON from a gzipped file.

	:param path: The filename to read from.

	:returns: The loaded JSON content.
	"""

	with gzip.open(PathPlus(path), 'r') as f:
		return sdjson.load(f)


def write_gzip_json(path: PathLike, data: JSONInput, indent: int = 2) -> None:
	"""
	Write JSON to a gzip file.

	:param path: The filename to write to.
	:param data: The JSON-serializable data to output.
	:param indent: Number of spaces used to indent JSON.
	"""

	json_data = sdjson.dumps(data, indent=indent)

	with gzip.open(PathPlus(path), 'w') as f:
		f.write(json_data.encode("utf-8"))
