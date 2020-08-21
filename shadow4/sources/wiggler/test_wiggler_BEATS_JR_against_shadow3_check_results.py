import numpy
from srxraylib.plot.gol import plot, set_qt

set_qt()

a = numpy.loadtxt("test_wiggler_BEATS_JR_against_shadow3.dat")
print(a.shape)
plot(a[:, 0], 1e6 * a[:, 1],
     a[:, 0], 1e6 * a[:, 2],
     a[:, 0], 1e6 * a[:, 3],legend=["shadow3","shadow4","diff"],
     xtitle="Photon energy [eV]",ytitle="FWHM [urad]")

for i in range(a.shape[0]):
     print("Checking Energy: %f  diff: %f  urad; shadow3 * 20 percent: %f urad" % (a[i, 0],1e6 * numpy.abs(a[i, 3]), 1e6 * 0.2 * numpy.abs(a[i, 1])))
     assert (numpy.abs(a[i, 3]) < 0.2 * numpy.abs(a[i, 1]))