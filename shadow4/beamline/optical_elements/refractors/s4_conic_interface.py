from shadow4.syned.shape import Conic
from shadow4.beamline.optical_elements.refractors.s4_interface import S4InterfaceElement, S4Interface
from shadow4.syned.element_coordinates import ElementCoordinates
from shadow4.optical_surfaces.s4_conic import S4Conic

from shadow4.beamline.s4_optical_element import S4ConicOpticalElement

class S4ConicInterface(S4Interface, S4ConicOpticalElement):
    def __init__(self,
                 name="Conic Refractive Interface",
                 boundary_shape=None,
                 material_object=None,
                 material_image=None,
                 f_r_ind=0,
                 r_ind_obj=1.0,
                 r_ind_ima=1.0,
                 r_attenuation_obj=0.0,
                 r_attenuation_ima=0.0,
                 file_r_ind_obj="",
                 file_r_ind_ima="",
                 conic_coefficients=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                 ):

        S4ConicOpticalElement.__init__(self, conic_coefficients)
        S4Interface.__init__(self,
                             name,
                             boundary_shape,
                             self._conic_surface_shape,
                             material_object=material_object,
                             material_image=material_image,
                             f_r_ind=f_r_ind,
                             r_ind_obj=r_ind_obj,
                             r_ind_ima=r_ind_ima,
                             r_attenuation_obj=r_attenuation_obj,
                             r_attenuation_ima=r_attenuation_ima,
                             file_r_ind_obj=file_r_ind_obj,
                             file_r_ind_ima=file_r_ind_ima,
                          )

class S4ConicInterfaceElement(S4InterfaceElement):
    def __init__(self, optical_element=None, coordinates=None):
        super().__init__(optical_element if optical_element is not None else S4ConicInterface(),
                         coordinates if coordinates is not None else ElementCoordinates())
        if not isinstance(self.get_optical_element().get_surface_shape(), Conic):
            raise ValueError("Wrong Optical Element: only Conic shape is accepted")

    def apply_local_refraction(self, beam):
        surface_shape = self.get_optical_element().get_surface_shape()

        ccc = S4Conic.initialize_from_coefficients(surface_shape.get_conic_coefficients())

        oe = self.get_optical_element()
        refraction_index_object, refraction_index_image = oe.get_refraction_indices()


        beam_mirr, normal = ccc.apply_refraction_on_beam(beam, refraction_index_object, refraction_index_image)

        return beam_mirr, normal