from syned.beamline.shape import Conic
from syned.beamline.element_coordinates import ElementCoordinates

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.optical_elements.refractors.s4_interface import S4InterfaceElement, S4Interface
from shadow4.optical_surfaces.s4_conic import S4Conic
from shadow4.beamline.s4_optical_element import S4ConicOpticalElementDecorator

class S4ConicInterface(S4Interface, S4ConicOpticalElementDecorator):
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

        S4ConicOpticalElementDecorator.__init__(self, conic_coefficients)
        S4Interface.__init__(self,
                             name,
                             boundary_shape,
                             self._curved_surface_shape,
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
    def __init__(self,
                 optical_element : S4ConicInterface = None,
                 coordinates : ElementCoordinates = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4ConicInterface(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         input_beam=input_beam)
        if not isinstance(self.get_optical_element().get_surface_shape(), Conic):
            raise ValueError("Wrong Optical Element: only Conic shape is accepted")

    def apply_local_refraction(self, beam):
        surface_shape = self.get_optical_element().get_surface_shape()

        ccc = S4Conic.initialize_from_coefficients(surface_shape.get_conic_coefficients())

        oe = self.get_optical_element()
        refraction_index_object, refraction_index_image = oe.get_refraction_indices()

        footprint, normal = ccc.apply_refraction_on_beam(beam, refraction_index_object, refraction_index_image)

        return footprint, normal
