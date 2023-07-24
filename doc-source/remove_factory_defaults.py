# stdlib
import inspect
from typing import Any, Callable, List, Optional, Tuple

# 3rd party
import attr
import sphinx_toolbox.more_autodoc.typehints
from sphinx.application import Sphinx
from sphinx.util.inspect import signature as Signature
from sphinx.util.inspect import stringify_signature
from sphinx_toolbox.utils import is_namedtuple, singleton


def preprocess_class_defaults(
		obj: Callable
		) -> Tuple[Optional[Callable], Optional[inspect.Signature], List[inspect.Parameter]]:
	"""
	Pre-processes the default values for the arguments of a class.

	.. versionadded:: 0.8.0

	:param obj: The class.

	:return: The class signature and a list of arguments/parameters.
	"""

	init: Optional[Callable[..., Any]] = getattr(obj, "__init__", getattr(obj, "__new__", None))

	if is_namedtuple(obj):
		init = getattr(obj, "__new__")

	try:
		signature = Signature(inspect.unwrap(init))  # type: ignore[arg-type]
	except ValueError:  # pragma: no cover
		return init, None, []

	parameters = []
	preprocessor_list = sphinx_toolbox.more_autodoc.typehints.default_preprocessors

	for argname, param in signature.parameters.items():
		default = param.default

		if default is not inspect.Parameter.empty:
			for check, preprocessor in preprocessor_list:
				if check(default):
					default = preprocessor(default)
					break

			else:
				if hasattr(obj, "__attrs_attrs__") and default is attr.NOTHING:
					if argname in {"date_created", "date_modified"}:
						default = singleton("datetime.datetime.now()")
					elif argname == "user":
						default = singleton("getpass.getuser()")
					elif argname == "device":
						default = singleton("socket.gethostname()")

		parameters.append(param.replace(annotation=inspect.Parameter.empty, default=default))

	return init, signature, parameters


def process_signature(  # noqa: MAN001
		app: Sphinx,
		what: str,
		name: str,
		obj,
		options,
		signature,
		return_annotation: Any,
		) -> Optional[Tuple[str, None]]:

	if name not in {"libgunshotmatch.datafile.Datafile", "libgunshotmatch.datafile.Repeat"}:
		return sphinx_toolbox.more_autodoc.typehints.process_signature(
				app, what, name, obj, options, signature, return_annotation
				)

	if not callable(obj):
		return None

	original_obj = obj

	obj, signature, parameters = preprocess_class_defaults(obj)

	obj = inspect.unwrap(obj)

	if not getattr(obj, "__annotations__", None):
		return None

	# The generated dataclass __init__() and class are weird and need extra checks
	# This helper function operates on the generated class and methods
	# of a dataclass, not an instantiated dataclass object. As such,
	# it cannot be replaced by a call to `dataclasses.is_dataclass()`.
	def _is_dataclass(name: str, what: str, qualname: str) -> bool:
		if what == "method" and name.endswith(".__init__"):
			# generated __init__()
			return True
		if what == "class" and qualname.endswith(".__init__"):
			# generated class
			return True
		return False

	# The generated dataclass __init__() is weird and needs the second condition
	if (
			hasattr(obj, "__qualname__") and "<locals>" in obj.__qualname__
			and not _is_dataclass(name, what, obj.__qualname__)
			):
		sphinx_toolbox.more_autodoc.typehints.sat_logger.warning(
				"Cannot treat a function defined as a local function: '%s'  (use @functools.wraps)", name
				)
		return None

	if parameters:
		if inspect.isclass(original_obj) or (what == "method" and name.endswith(".__init__")):
			del parameters[0]
		elif what == "method":

			try:
				outer = inspect.getmodule(obj)
				if outer is not None:
					for clsname in obj.__qualname__.split('.')[:-1]:
						outer = getattr(outer, clsname)
			except AttributeError:
				outer = None

			method_name = obj.__name__
			if method_name.startswith("__") and not method_name.endswith("__"):
				# If the method starts with double underscore (dunder)
				# Python applies mangling so we need to prepend the class name.
				# This doesn't happen if it always ends with double underscore.
				class_name = obj.__qualname__.split('.')[-2]
				method_name = f"_{class_name}{method_name}"

			if outer is not None:
				method_object = outer.__dict__[method_name] if outer else obj
				if not isinstance(method_object, (classmethod, staticmethod)):
					del parameters[0]

			else:
				if not inspect.ismethod(obj) and parameters[0].name in {"self", "cls", "_cls"}:
					del parameters[0]

	signature = signature.replace(parameters=parameters, return_annotation=inspect.Signature.empty)

	return stringify_signature(signature), None  # .replace('\\', '\\\\')


def setup(app: Sphinx):
	app.connect("autodoc-process-signature", process_signature, priority=0)
