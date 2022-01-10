import numpy

from shadow4.syned.shape import Rectangle

from shadow4.syned.element_coordinates import ElementCoordinates
from shadow4.syned.refractors.interface import Interface

from shadow4.physical_models.prerefl.prerefl import PreRefl

from shadow4.beamline.s4_beamline_element import S4BeamlineElement


class S4Interface(Interface):

    def __init__(self,
                 name="Undefined",
                 boundary_shape=None,
                 surface_shape=None,
                 material_object=None,
                 material_image=None, ):

        Interface.__init__(self,
                        name=name,
                        surface_shape=surface_shape,
                        boundary_shape=boundary_shape,
                        material_object=material_object,
                        material_image=material_image,
                        )

        # self._f_refl = f_refl
        # self._file_refl = file_refl
        # self._refraction_index = refraction_index



class S4InterfaceElement(S4BeamlineElement):
    
    def __init__(self, optical_element=None, coordinates=None):
        super().__init__(optical_element if optical_element is not None else S4Interface(),
                         coordinates if coordinates is not None else ElementCoordinates())
    
    def trace_beam(self, beam_in, flag_lost_value=-1):

        p = self.get_coordinates().p()
        q = self.get_coordinates().q()
        theta_grazing1 = numpy.pi / 2 - self.get_coordinates().angle_radial()
        alpha1 = self.get_coordinates().angle_azimuthal()

        #
        beam = beam_in.duplicate()

        #
        # put beam in mirror reference system
        #
        beam.rotate(alpha1, axis=2)
        beam.rotate(theta_grazing1, axis=1)
        beam.translation([0.0, -p * numpy.cos(theta_grazing1), p * numpy.sin(theta_grazing1)])

        #
        # reflect beam in the mirror surface
        #
        soe = self.get_optical_element() #._optical_element_syned
        # print(">>> CCC", soe.get_surface_shape().get_conic_coefficients())

        # v_in = beam.get_columns([4,5,6])
        if not isinstance(soe, Interface): # undefined
            raise Exception("Undefined refractive interface")
        else:
            print(">>>>>> self.apply_local_refraction(beam)")
            beam_mirr, normal = self.apply_local_refraction(beam)
            print(">>>>>> BACK FROM self.apply_local_refraction(beam)")
        #
        # apply mirror boundaries
        #
        beam_mirr.apply_boundaries_syned(soe.get_boundary_shape(), flag_lost_value=flag_lost_value)

        #
        # TODO" apply lens absorption
        #

        #
        # TODO: write angle.xx for comparison
        #


        #
        # from element reference system to image plane
        #

        beam_out = beam_mirr.duplicate()
        beam_out.change_to_image_reference_system(theta_grazing1, q)

        return beam_out, beam_mirr

        return beam, beam

    def apply_local_refraction(self, beam):
        raise NotImplementedError()


if __name__ == "__main__":
    pass
