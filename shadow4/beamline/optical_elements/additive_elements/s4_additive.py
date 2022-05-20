import numpy

from shadow4.syned.element_coordinates import ElementCoordinates
from syned.beamline.optical_element import OpticalElement
from shadow4.beamline.s4_optical_element import S4OpticalElement
from shadow4.beamline.s4_beamline_element import S4BeamlineElement

from shadow4.beamline.optical_elements.mirrors.s4_surface_data_mirror import S4SurfaceDataMirror, S4SurfaceDataMirrorElement


class S4Additive(OpticalElement, S4OpticalElement):
    def __init__(self, name="Undefined",optical_elements_list=[]):
        super().__init__(name=name)

        self._optical_elements_list = optical_elements_list


    def get_optical_elements_list(self):
        return self._optical_elements_list


class S4AdditiveElement(S4BeamlineElement):

    def __init__(self, optical_element=None, coordinates=None):
        super().__init__(optical_element if optical_element is not None else S4Additive(),
                         coordinates if coordinates is not None else ElementCoordinates())


    def trace_beam(self,beam1):

        oe_list = self.get_optical_element().get_optical_elements_list()

        if len(oe_list) == 1:
            print(">>> only 1 oe in additive")

        if not isinstance(oe_list[0], S4SurfaceDataMirror):
            raise Exception("First element in additive list must be S4SurfaceDataMirror")

        p, q, angle_radial, angle_radial_out, angle_azimuthal = self.get_coordinates().get_positions()

        theta_grazing1 = numpy.pi / 2 - angle_radial
        theta_grazing2 = numpy.pi / 2 - angle_radial_out
        alpha1 = angle_azimuthal

        #
        beam = beam1.duplicate()

        #
        # put beam in mirror reference system
        #
        beam.rotate(alpha1, axis=2)
        beam.rotate(theta_grazing1, axis=1)
        beam.translation([0.0, -p * numpy.cos(theta_grazing1), p * numpy.sin(theta_grazing1)])

        #
        # recurrent call
        #

        # TODO: this is NOT CORRECT....  TO BE REDONE....
        # we should first calculate the added surface and then apply surface_data to get the intercept
        for oe in oe_list:
            mirr, normal = oe.apply_geometrical_model(beam)


        #
        # from oe reference system to image plane
        #
        beam_out = beam.duplicate()
        beam_out.change_to_image_reference_system(theta_grazing2, q)


        return beam_out, beam




