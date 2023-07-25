#!/usr/bin/env python3
#
#  _fields.py
"""
Internal attrs field helpers.
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
from typing import Any, Generic, Iterable, Optional, Set, Tuple, Type, TypeVar

# 3rd party
import attr

__all__ = [
		"Boolean",
		"FieldType",
		"Integer",
		"Number",
		"String",
		"convert_crop_mass_range",
		"convert_rt_range",
		"default_base_peak_filter"
		]

_FT = TypeVar("_FT")


class FieldType(Generic[_FT]):
	"""
	Customisable method field type.
	"""

	#: The Python type of the object.
	field_type: Type[_FT]

	#: String name of the field.
	field_name: str

	def __new__(cls: Type["FieldType"]) -> Type["FieldType"]:  # type: ignore[misc]
		"""
		Make the ``FieldType`` class unitializable.

		Mypy doesn't like it, hence the ``type: ignore``.
		"""

		return cls

	@classmethod
	def validator(cls: Type["FieldType"], inst: object, attr: attr.Attribute, value: Any) -> None:
		"""
		Check if a value conforms to this field's expected datatype.

		This method is called by attrs.
		"""

		if not isinstance(value, cls.field_type):
			err_msg = f"'{attr.name}' must be {cls.field_name} (got {value!r} that is a {value.__class__!r})."
			raise TypeError(err_msg)

	@classmethod
	def field(cls: Type["FieldType"], default: _FT) -> _FT:
		"""
		Construct an attrs field.
		"""
		# Actually returns attr.Attribute, but mypy doesn't like it

		return attr.field(
				default=default,
				converter=cls.field_type,
				validator=cls.validator,
				)


class Boolean(FieldType):
	"""
	Boolean method field type.
	"""

	field_type = bool
	field_name = "boolean"


class String(FieldType):
	"""
	String method field type.
	"""

	field_type = str
	field_name = "a string"


class Integer(FieldType):
	"""
	Integer method field type.
	"""

	field_type = int
	field_name = "an integer"


class Number(FieldType):
	"""
	Numerical method field type.
	"""

	field_type = float
	field_name = "a number"

	@classmethod
	def validator(cls: Type["FieldType"], inst: object, attr: attr.Attribute, value: Any) -> None:
		"""
		Check if a value is a number.

		This method is called by attrs.
		"""

		if not isinstance(value, (int, float)):
			err_msg = f"'{attr.name}' must be {cls.field_name} (got {value!r} that is a {value.__class__!r})."
			raise TypeError(err_msg)


def convert_crop_mass_range(value: Optional[Iterable]) -> Optional[Tuple[int, int]]:
	"""
	Convert a value to a tuple of integers for the ``crop_mass_range`` option.
	"""

	if not value:
		return None

	value_as_tuple = tuple([int(x) for x in value])
	if len(value_as_tuple) != 2:
		err_msg = f"crop_mass_range must be a 2-element tuple. Got {len(value_as_tuple)} values: {value_as_tuple}"
		raise ValueError(err_msg)

	return value_as_tuple  # type: ignore[return-value]


def convert_rt_range(value: Optional[Iterable]) -> Optional[Tuple[float, float]]:
	"""
	Convert a value to a tuple of floats for the ``rt_range`` option.
	"""

	if not value:
		return None

	value_as_tuple = tuple([float(x) for x in value])
	if len(value_as_tuple) != 2:
		raise ValueError(f"rt_range must be a 2-element tuple. Got {len(value_as_tuple)} values: {value_as_tuple}")

	return value_as_tuple  # type: ignore[return-value]


def default_base_peak_filter() -> Set[int]:
	"""
	Returns the default value for the ``base_peak_filter`` option.
	"""

	return {73, 147}
