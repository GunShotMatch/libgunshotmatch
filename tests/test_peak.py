# 3rd party
from pyms.Peak import Peak
from pyms.Spectrum import MassSpectrum

# this package
from libgunshotmatch.peak import QualifiedPeak


def test_qualified_from_peak():
	peak = Peak(
			rt=123.45,
			ms=MassSpectrum([1, 2, 3, 4, 5], [500, 400, 300, 200, 100]),
			)

	peak.area = 12345678.9
	peak.bounds = (7, 8, 9)

	qp = QualifiedPeak.from_peak(peak)
	assert qp.rt == 123.45
	assert qp.mass_spectrum == MassSpectrum([1, 2, 3, 4, 5], [500, 400, 300, 200, 100])
	assert qp.area == 12345678.9
	assert qp.bounds == (7, 8, 9)
	assert qp.UID == peak.UID
