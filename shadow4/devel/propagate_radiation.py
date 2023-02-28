"""
TODO: remove this file
"""

from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D
from srxraylib.plot.gol import plot

#
# Import section
#
import numpy
from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters
from syned.beamline.beamline_element import BeamlineElement
from syned.beamline.element_coordinates import ElementCoordinates
from wofryimpl.propagator.propagators1D.fresnel_zoom import FresnelZoom1D


input_wavefront = GenericWavefront1D()
input_wavefront = input_wavefront.load_h5_file("/users/srio/OASYS1.1/minishadow/minishadow/undulator/tmp.h5","wfr")

# input_wavefront.rescale_amplitude( 1.0/numpy.sqrt(input_wavefront.get_intensity().max()))

# input_wavefront.set_spherical_wave(radius=100,complex_amplitude=numpy.abs(input_wavefront.get_intensity()) )
# plot(input_wavefront.get_abscissas()*1e6,input_wavefront.get_intensity())

from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D

optical_element = WOScreen1D()

#
# propagating
#
#
propagation_elements = PropagationElements()
beamline_element = BeamlineElement(optical_element=optical_element,
                coordinates=ElementCoordinates(p=-100.000000,q=0.000000,
                angle_radial=numpy.radians(0.000000),
                angle_azimuthal=numpy.radians(0.000000)))
propagation_elements.add_beamline_element(beamline_element)
propagation_parameters = PropagationParameters(wavefront=input_wavefront.duplicate(),propagation_elements = propagation_elements)
propagation_parameters.set_additional_parameters('magnification_x', 0.01 ) #0.010000)

#
propagator = PropagationManager.Instance()
try:
    propagator.add_propagator(FresnelZoom1D())
except:
    pass
output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,handler_name='FRESNEL_ZOOM_1D')


plot(output_wavefront.get_abscissas()*1e6,output_wavefront.get_intensity())
