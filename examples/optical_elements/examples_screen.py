
import numpy

from shadow4.sources.source_geometrical.source_gaussian import SourceGaussian
from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen, S4ScreenElement

from shadow4.tools.graphics import plotxy

from syned.beamline.element_coordinates import ElementCoordinates

from syned.beamline.shape import Rectangle, Ellipse


def get_sigmas_radiation(photon_energy,undulator_length):
    import scipy.constants as codata
    lambdan = 1e-10 * codata.h*codata.c/codata.e*1e10 / photon_energy # in m
    print("wavelength in m",lambdan)
    return 1e6*2.740/4/numpy.pi*numpy.sqrt(lambdan*undulator_length),1e6*0.69*numpy.sqrt(lambdan/undulator_length)


def example_screen(do_plot=True):

    #
    # collimated source
    #
    src = SourceGaussian.initialize_collimated_source(nrays=10000,sigmaX=1e-6,sigmaZ=1e-6)

    beam = src.get_beam()

    print(beam.info())

    #
    # screen definition
    #

    screen1 = S4ScreenElement(optical_element=S4Screen(), coordinates=ElementCoordinates(p=100.0, q=0.0), input_beam=beam)

    print(screen1.info())

    beam2, tmp = screen1.trace_beam()

    #
    if do_plot:
        plotxy(beam2,1,3,nbins=100,title="SCREEN")


def example_slit(do_plot=True):

    src = SourceGaussian.initialize_collimated_source(nrays=10000,sigmaX=1e-6,sigmaZ=1e-6)
    beam = src.get_beam()


    #
    # slit definition
    #
    boundary_shape = Rectangle(x_left=-0.5e-6, x_right=0.5e-6, y_bottom=-0.5e-6, y_top=0.5e-6)
    coordinates = ElementCoordinates(p=100.0, q=0.0)
    optical_element = S4Screen(name="slit1", boundary_shape=boundary_shape,
                                i_abs=False, i_stop=False, thick=0.0, file_abs="")

    slit1 = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    print(slit1.info())

    #
    # trace
    #

    beam2, tmp = slit1.trace_beam()


    #
    if do_plot:
        plotxy(beam2,1,3,nbins=100,title="SLIT", nolost=True)



def example_beam_stopper(do_plot=True):

    src = SourceGaussian.initialize_collimated_source(nrays=10000,sigmaX=1e-6,sigmaZ=1e-6)
    beam = src.get_beam()


    #
    # slit definition
    #
    boundary_shape = Ellipse(-1e-6, 1e-6, -0.5e-6, 0.5e-6)
    coordinates = ElementCoordinates(p=100.0, q=0.0)
    optical_element = S4Screen(name="slit1", boundary_shape=boundary_shape,
                                i_abs=False, i_stop=True, thick=0.0, file_abs="")

    screen1 = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    print(screen1.info())

    #
    # trace
    #

    beam2, tmp = screen1.trace_beam()


    #
    if do_plot:
        plotxy(beam2,1,3,nbins=100,title="BEAM STOPPER", nolost=True)


def example_filter(do_plot=True):

    src = SourceGaussian.initialize_collimated_source(nrays=10000,sigmaX=1e-6,sigmaZ=1e-6)
    beam = src.get_beam()

    #
    # slit definition
    #


    from shadow4.physical_models.prerefl.prerefl import PreRefl
    PreRefl.prerefl(interactive=False, SYMBOL="Be", DENSITY=1.848, FILE="Be.dat", E_MIN=100.0, E_MAX=20000.0, E_STEP=100.0)



    optical_element = S4Screen(name="filter1", boundary_shape=None,
                                i_abs=True, i_stop=False, thick=10e-6, file_abs="Be.dat")
    coordinates = ElementCoordinates(p=100.0, q=0.0)
    filter1 = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)


    print(filter1.info())

    #
    # trace
    #

    beam2, tmp = filter1.trace_beam()

    #
    if do_plot:
        plotxy(beam2,1,3,nbins=100,title="FILTER", nolost=True)

    print("Intensity: ", beam2.intensity())

def example_holed_filter(do_plot=True):

    src = SourceGaussian.initialize_collimated_source(nrays=10000,sigmaX=1e-6,sigmaZ=1e-6)
    beam = src.get_beam()


    #
    # slit definition
    #
    boundary_shape = Rectangle(x_left=-0.5e-6, x_right=0.5e-6, y_bottom=-0.5e-6, y_top=0.5e-6)

    optical_element = S4Screen(name="filter1", boundary_shape=boundary_shape,
                                i_abs=True, i_stop=True, thick=10e-6, file_abs="Be.dat")

    coordinates = ElementCoordinates(p=100.0, q=0.0)
    filter1 = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    print(filter1.info())

    #
    # trace
    #

    beam2, tmp = filter1.trace_beam()


    #
    if do_plot:
        plotxy(beam2,1,3,nbins=100,title="HOLED FILTER", nolost=True)


if __name__ == "__main__":
    from srxraylib.plot.gol import set_qt
    set_qt()

    do_plot = True

    example_screen(do_plot=do_plot)
    example_slit(do_plot=do_plot)
    example_beam_stopper(do_plot=do_plot)
    example_filter(do_plot=do_plot)
    example_holed_filter(do_plot=do_plot)
