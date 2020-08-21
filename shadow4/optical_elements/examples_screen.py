
import numpy
from shadow4.beam.beam import Beam
from shadow4.sources.source_geometrical.gaussian import SourceGaussian
from shadow4.optical_elements.screen import Screen


from shadow4.compatibility.beam3 import Beam3
from Shadow.ShadowTools import plotxy

from syned.beamline.optical_elements.ideal_elements.lens import IdealLens as SyIdelLens
from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.beamline_element import BeamlineElement

from syned.beamline.optical_elements.ideal_elements.screen import Screen as SyScreen
from syned.beamline.optical_elements.absorbers.absorber import Absorber as SyAbsorber
from syned.beamline.optical_elements.absorbers.beam_stopper import BeamStopper as SyBeamStopper
from syned.beamline.optical_elements.absorbers.filter import Filter as SyFilter
from syned.beamline.optical_elements.absorbers.slit import Slit as SySlit



def get_sigmas_radiation(photon_energy,undulator_length):
    import scipy.constants as codata
    lambdan = 1e-10 * codata.h*codata.c/codata.e*1e10 / photon_energy # in m
    print("wavelength in m",lambdan)
    return 1e6*2.740/4/numpy.pi*numpy.sqrt(lambdan*undulator_length),1e6*0.69*numpy.sqrt(lambdan/undulator_length)


def test_with_collimated_beam(do_plot=True):

    #
    # collimated source
    #
    src = SourceGaussian.initialize_collimated_source(number_of_rays=10000,sigmaX=1e-6,sigmaZ=1e-6)

    interface = 'old' # 'new'

    if interface == 'new':
        beam = src.get_beam()
    elif interface == "old":
        beam = Beam()
        beam.genSource(src)
    print(beam.info())
    SX, SZ = (1e6*beam.get_standard_deviation(1),1e6*beam.get_standard_deviation(3))

    if do_plot:
        beam3 = Beam3.initialize_from_shadow4_beam(beam)
        plotxy(beam3,1,3,nbins=100,title="SOURCE")

    #
    # screen definition
    #

    syscreen1 = SyScreen(name="Undefined")

    coordinates_syned = ElementCoordinates(p=100.0, q=0.0)

    beamline_element_syned = BeamlineElement(optical_element=syscreen1, coordinates=coordinates_syned)

    screen1 = Screen(beamline_element_syned=beamline_element_syned)

    print(screen1.info())

    #
    # trace
    #
    interface = 'new'
    if interface == 'new':
        beam2 = screen1.trace_beam(beam)
    elif interface == 'old':
        beam2 = beam.traceOE(screen1, 1, overwrite=False)

    #
    if do_plot:
        beam3 = Beam3.initialize_from_shadow4_beam(beam2)
        plotxy(beam3,1,3,nbins=100,title="FOCAL PLANE")




if __name__ == "__main__":
    test_with_collimated_beam(do_plot=False)
