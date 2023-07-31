#!/usr/bin/env python3
#
#  datafile.py
"""
Represents a parsed datafile.
"""
#
#  Copyright © 2020-2023 Dominic Davis-Foster <dominic@davis-foster.co.uk>
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
import getpass
import os
import socket
from datetime import datetime
from statistics import mean, median
from typing import Any, Dict, List, Mapping, NamedTuple, Optional, Tuple, Type

# 3rd party
import attr
import pyms.Noise.SavitzkyGolay  # type: ignore[import]
import pyms.TopHat  # type: ignore[import]
from domdf_python_tools.paths import PathPlus
from domdf_python_tools.typing import PathLike
from enum_tools import IntEnum
from pyms.GCMS.Class import GCMS_data  # type: ignore[import]
from pyms.IntensityMatrix import IntensityMatrix, build_intensity_matrix_i  # type: ignore[import]

# this package
from libgunshotmatch import gzip_util
from libgunshotmatch.peak import PeakList, QualifiedPeak, _to_peak_list, peak_from_dict

__all__ = ["Datafile", "FileType", "Repeat", "get_info_from_gcms_data", "GCMSDataInfo"]


class FileType(IntEnum):
	"""
	Represents the input datafile types supported by PyMassSpec.
	"""

	#: JCAMP-DX (https://iupac.org/wp-content/uploads/2021/08/JCAMP-DX_MS_1994.pdf)
	JDX = 0

	#: mzML (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3073315/)
	MZML = 1

	#: Analytical Data Interchange Format for Mass Spectrometry (https://www.astm.org/e1947-98r22.html)
	ANDI = 2

	@classmethod
	def _attrs_convert(cls: Type["FileType"], value: int) -> "FileType":
		if value in iter(cls):  # type: ignore[attr-defined, arg-type]
			return cls(value)

		raise ValueError(f"Unrecognised file format {value}")


@attr.define
class Datafile:
	"""
	Represents a single datafile in a project.

	:default user: taken from the currently logged-in user
	:default device: taken from the current device's hostname
	:default date_created: is the current date and time
	:default date_modified: is the current date and time
	"""

	#: The name of the :class:`~.Datafile`.
	name: str
	#: The filename of the file the :class:`~.Datafile` was created from.
	original_filename: str
	#: The filetype of the file the :class:`~.Datafile` was created from.
	original_filetype: FileType = attr.field(
			converter=FileType._attrs_convert,  # type: ignore[misc]  # doesn't like the converter
			)
	#: A description of the :class:`~.Datafile`.
	description: str = attr.field(default='')
	#: PyMassSpec :class:`~pyms.IntensityMatrix.IntensityMatrix` object.
	intensity_matrix: Optional[IntensityMatrix] = attr.field(default=None)

	#: The user who created the :class:`~.Datafile`.
	user: str = attr.field(factory=getpass.getuser)
	#: The device that created the :class:`~.Datafile`.
	device: str = attr.field(factory=socket.gethostname)
	#: The date and time the :class:`~.Datafile` was created.
	date_created: datetime = attr.field(factory=datetime.now)
	#: The date and time the :class:`~.Datafile` was last modified.
	date_modified: datetime = attr.field(factory=datetime.now)
	#: File format version
	version: int = attr.field(default=1)

	@classmethod
	def new(cls: Type["Datafile"], name: str, filename: PathLike) -> "Datafile":
		"""
		Construct a new :class:`~.Datafile` from a file.

		:param name: The name of the Datafile
		:param filename:
		"""

		date_modified = date_created = datetime.now()

		filename_p = PathPlus(filename)

		if filename_p.suffix.lower() == ".jdx":
			original_filetype = FileType.JDX
		elif filename_p.suffix.lower() == ".cdf":
			original_filetype = FileType.ANDI
		elif filename_p.suffix.lower() == ".mzml":
			original_filetype = FileType.MZML
		else:
			raise ValueError(f"Unrecognised file format for file {filename}")

		return cls(
				name=name,
				date_created=date_created,
				date_modified=date_modified,
				original_filename=os.fspath(filename),
				original_filetype=original_filetype,
				)

	def load_gcms_data(self) -> GCMS_data:
		"""
		Load GC-MS data from the datafile.
		"""

		if self.original_filetype == FileType.JDX:
			# 3rd party
			from pyms.GCMS.IO.JCAMP import JCAMP_reader  # type: ignore[import]
			gcms_data = JCAMP_reader(self.original_filename)

		elif self.original_filetype == FileType.MZML:
			# 3rd party
			from pyms.GCMS.IO.MZML import mzML_reader  # type: ignore[import]
			gcms_data = mzML_reader(self.original_filename)

		elif self.original_filetype == FileType.ANDI:
			# 3rd party
			from pyms.GCMS.IO.ANDI import ANDI_reader  # type: ignore[import]
			gcms_data = ANDI_reader(self.original_filename)

		else:
			# Shouldn't get here due to filetype validation at attrs' end
			raise ValueError(f"Unrecognised file format {self.original_filetype._name_}")

		return gcms_data

	def prepare_intensity_matrix(
			self,
			gcms_data: GCMS_data,
			savitzky_golay: bool = True,
			tophat: bool = True,
			tophat_structure_size: str = "1.5m",  # Ignored if tophat=False
			crop_mass_range: Optional[Tuple[int, int]] = None,
			) -> IntensityMatrix:
		"""
		Build an :class:`~pyms.IntensityMatrix.IntensityMatrix` for the datafile.

		:param gcms_data:
		:param savitzky_golay: Whether to perform Savitzky-Golay smoothing.
		:param tophat: Whether to perform Tophat baseline correction.
		:param tophat_structure_size: The structure size for Tophat baseline correction.
		:param crop_mass_range: The range of masses to which the GC-MS data should be limited to.

		"""
		intensity_matrix = build_intensity_matrix_i(gcms_data)

		# Show the m/z of the maximum and minimum bins
		print(f" Minimum m/z bin: {intensity_matrix.min_mass}")
		print(f" Maximum m/z bin: {intensity_matrix.max_mass}")

		# Crop masses
		if crop_mass_range is not None:
			# None means don't crop

			# min_mass, max_mass = 50, 400
			min_mass, max_mass = crop_mass_range

			# Catch case where numbers are flipped
			if min_mass >= max_mass:
				raise ValueError(
						'\n'.join([
								"Cannot crop mass range when `max mass` is less than `min mass`.\n'",
								"Did you put the numbers the wrong way around? The expected order is (<min>, <max>)",
								])
						)

			if min_mass < intensity_matrix.min_mass:
				min_mass = intensity_matrix.min_mass
			if max_mass > intensity_matrix.max_mass:
				max_mass = intensity_matrix.max_mass

			intensity_matrix.crop_mass(min_mass, max_mass)

		# Perform Data filtering
		n_mz = intensity_matrix.size[1]

		# Iterate over each IC in the intensity matrix
		for index in range(n_mz):
			# print("\rWorking on IC#", index+1, '  ',end='')
			ic = intensity_matrix.get_ic_at_index(index)

			if savitzky_golay:
				# Perform Savitzky-Golay smoothing.
				# Note that Turbomass does not use smoothing for qualitative method.
				ic = pyms.Noise.SavitzkyGolay.savitzky_golay(ic)

			if tophat:
				# Perform Tophat baseline correction
				# Top-hat baseline Correction seems to bring down noise,
				#  		retaining shapes, but keeps points on actual peaks
				# ic = tophat(ic, struct=method.tophat_struct)
				ic = pyms.TopHat.tophat(ic, struct=tophat_structure_size)

			# Set the IC in the intensity matrix to the filtered one
			intensity_matrix.set_ic_at_index(index, ic)

		self.intensity_matrix = intensity_matrix
		return intensity_matrix

	def to_dict(self) -> Dict[str, Any]:
		"""
		Returns a dictionary representation of this :class:`~.Datafile`.

		All keys are native, JSON-serializable, Python objects.
		"""

		if self.intensity_matrix is None:
			im_as_dict = None
		else:
			im_as_dict = {
					"times": self.intensity_matrix.time_list,
					"masses": self.intensity_matrix.mass_list,
					"intensities": self.intensity_matrix.intensity_array.tolist(),
					}

		return {
				"name": self.name,
				"original_filename": self.original_filename,
				"original_filetype": int(self.original_filetype),
				"description": self.description,
				"intensity_matrix": im_as_dict,
				"user": self.user,
				"device": self.device,
				"date_created": self.date_created.isoformat(),
				"date_modified": self.date_modified.isoformat(),
				"version": self.version,
				}

	@classmethod
	def from_dict(cls: Type["Datafile"], d: Mapping[str, Any]) -> "Datafile":
		"""
		Construct a :class:`~.Datafile` from a dictionary.

		:param d:
		"""

		im_as_dict = d["intensity_matrix"]
		if im_as_dict is None:
			intensity_matrix = None
		else:
			intensity_matrix = IntensityMatrix(
					time_list=im_as_dict["times"],
					mass_list=im_as_dict["masses"],
					intensity_array=im_as_dict["intensities"],
					)

		optional_keys = {}
		if "user" in d:
			optional_keys["user"] = d["user"]
		if "device" in d:
			optional_keys["device"] = d["device"]
		if "date_created" in d:
			optional_keys["date_created"] = datetime.fromisoformat(d["date_created"])
		if "date_modified" in d:
			optional_keys["date_modified"] = datetime.fromisoformat(d["date_modified"])
		if "version" in d:
			optional_keys["version"] = d["version"]

		return cls(
				name=d["name"],
				original_filename=d["original_filename"],
				original_filetype=FileType(d["original_filetype"]),
				description=d["description"],
				intensity_matrix=intensity_matrix,
				**optional_keys,
				)

	def export(self, output_dir: PathLike) -> str:
		"""
		Export as a ``.gsmd`` file.

		:returns: The output filename.
		"""

		export_filename = os.path.join(output_dir, f"{self.name}.gsmd")
		gzip_util.write_gzip_json(export_filename, self.to_dict(), indent=0)
		return export_filename

	@classmethod
	def from_file(cls: Type["Datafile"], filename: PathLike) -> "Datafile":
		"""
		Parse a ``gsmd`` file.

		:param filename: The input filename.
		"""

		as_dict: Dict[str, Any] = gzip_util.read_gzip_json(filename)  # type: ignore[assignment]
		return Datafile.from_dict(as_dict)


class GCMSDataInfo(NamedTuple):
	"""
	Represents information about a :class:`pyms.GCMS.Class.GCMS_data` object returned by :func:`get_info_from_gcms_data`.
	"""

	#: The minimum and maximum retention times.
	rt_range: Tuple[float, float]

	#: The average time step between scans.
	time_step: float

	#: The standard deviation of the time steps between scans.
	time_step_stdev: float

	#: The total number of scans.
	n_scans: int

	#: The minimum and maximum mass (*m/z*) values.
	mz_range: Tuple[float, float]

	#: The mean average number of masses per scan.
	num_mz_mean: float

	#: The median number of masses per scan.
	num_mz_median: float


def get_info_from_gcms_data(gcms_data: GCMS_data) -> GCMSDataInfo:
	"""
	Returns a information about the data in a :class:`pyms.GCMS.Class.GCMS_data` object.
	"""

	# TODO: within pyms make read only properties for these private attributes

	# calculate median number of m/z values measured per scan
	scan_list = gcms_data.scan_list

	n_list = [len(scan) for scan in scan_list]
	n_mz_mean = mean(n_list)
	n_mz_median = median(n_list)

	info = GCMSDataInfo(
			rt_range=(gcms_data._min_rt / 60, gcms_data._max_rt / 60),
			time_step=gcms_data._time_step,
			time_step_stdev=gcms_data._time_step_std,
			n_scans=len(scan_list),
			mz_range=(gcms_data._min_mass, gcms_data._max_mass),
			num_mz_mean=n_mz_mean,
			num_mz_median=n_mz_median,
			)

	return info


@attr.define
class Repeat:
	"""
	Represents a repeat sample in a project, constructed from a datafile.

	:default user: taken from the currently logged-in user.
	:default device: taken from the current device's hostname.
	:default date_created: is the current date and time.
	:default date_modified: is the current date and time.
	"""

	#: The :class:`~.Datafile` for this repeat.
	datafile: Datafile

	peaks: PeakList = attr.field(converter=_to_peak_list)

	#: Peaks containing identities from library search. This is usually populated after peak alignment.
	qualified_peaks: Optional[List[QualifiedPeak]] = attr.field(default=None)

	#: The user who created the :class:`~.Repeat`.
	user: str = attr.field(factory=getpass.getuser)
	#: The device that created the :class:`~.Repeat`.
	device: str = attr.field(factory=socket.gethostname)
	#: The date and time the :class:`~.Repeat` was created.
	date_created: datetime = attr.field(factory=datetime.now)
	#: The date and time the :class:`~.Repeat` was last modified.
	date_modified: datetime = attr.field(factory=datetime.now)
	#: File format version
	version: int = attr.field(default=1)

	def to_dict(self) -> Dict[str, Any]:
		"""
		Returns a dictionary representation of this :class:`~.Repeat`.

		All keys are native, JSON-serializable, Python objects.
		"""

		if self.qualified_peaks is None:
			qualified_peaks_as_list = None
		else:
			qualified_peaks_as_list = [qp.to_dict() for qp in self.qualified_peaks]

		return {
				"datafile": self.datafile.to_dict(),
				"peaks": self.peaks.to_list(),
				"qualified_peaks": qualified_peaks_as_list,
				"user": self.user,
				"device": self.device,
				"date_created": self.date_created.isoformat(),
				"date_modified": self.date_modified.isoformat(),
				"version": self.version,
				}

	@classmethod
	def from_dict(cls: Type["Repeat"], d: Mapping[str, Any]) -> "Repeat":
		"""
		Construct a :class:`~.Repeat` from a dictionary.

		:param d:
		"""

		datafile = Datafile.from_dict(d["datafile"])

		peaks = PeakList(peak_from_dict(peak) for peak in d["peaks"])

		qualified_peaks_as_list = d["qualified_peaks"]
		if qualified_peaks_as_list is None:
			qualified_peaks = None
		else:
			qualified_peaks = [QualifiedPeak.from_dict(peak) for peak in d["qualified_peaks"]]

		optional_keys = {}
		if "user" in d:
			optional_keys["user"] = d["user"]
		if "device" in d:
			optional_keys["device"] = d["device"]
		if "date_created" in d:
			optional_keys["date_created"] = datetime.fromisoformat(d["date_created"])
		if "date_modified" in d:
			optional_keys["date_modified"] = datetime.fromisoformat(d["date_modified"])
		if "version" in d:
			optional_keys["version"] = d["version"]

		return cls(
				datafile=datafile,
				peaks=peaks,
				qualified_peaks=qualified_peaks,
				**optional_keys,
				)
