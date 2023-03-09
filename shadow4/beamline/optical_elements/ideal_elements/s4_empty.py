import numpy
from syned.beamline.optical_elements.ideal_elements.screen import Screen
from syned.beamline.element_coordinates import ElementCoordinates

from shadow4.beamline.s4_optical_element import S4OpticalElementDecorator
from shadow4.beamline.s4_beamline_element import S4BeamlineElement
from shadow4.beam.s4_beam import S4Beam

class S4Empty(Screen, S4OpticalElementDecorator):
    def __init__(self, name="Undefined"):
        super().__init__(name=name)

class S4EmptyElement(S4BeamlineElement):
    def __init__(self,
                 optical_element : S4Empty = None,
                 coordinates : ElementCoordinates = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4Empty(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         input_beam=input_beam)

    def trace_beam(self, **params):
        p, q, angle_radial, angle_radial_out, angle_azimuthal = self.get_coordinates().get_positions()

        theta_grazing1 = numpy.pi / 2 - angle_radial
        theta_grazing2 = numpy.pi / 2 - angle_radial_out
        alpha1 = angle_azimuthal

        #
        input_beam = self.get_input_beam().duplicate()

        #
        # put beam in mirror reference system
        #
        input_beam.rotate(alpha1, axis=2)
        input_beam.rotate(theta_grazing1, axis=1)
        input_beam.translation([0.0, -p * numpy.cos(theta_grazing1), p * numpy.sin(theta_grazing1)])

        #
        # oe does nothing
        #
        pass

        #
        # from oe reference system to image plane
        #
        output_beam = input_beam.duplicate()
        output_beam.change_to_image_reference_system(theta_grazing2, q)


        return output_beam, input_beam

