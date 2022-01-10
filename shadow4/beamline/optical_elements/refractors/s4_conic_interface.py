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
                 conic_coefficients=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                 ):

        S4ConicOpticalElement.__init__(self, conic_coefficients)
        S4Interface.__init__(self, name, boundary_shape, self._conic_surface_shape,
                          material_object=material_object, material_image=material_image)

class S4ConicInterfaceElement(S4InterfaceElement):
    def __init__(self, optical_element=None, coordinates=None):
        super().__init__(optical_element if optical_element is not None else S4ConicInterface(),
                         coordinates if coordinates is not None else ElementCoordinates())
        if not isinstance(self.get_optical_element().get_surface_shape(), Conic):
            raise ValueError("Wrong Optical Element: only Conic shape is accepted")

    def apply_local_refraction(self, beam):
        surface_shape = self.get_optical_element().get_surface_shape()

        print(">>>>> Conic refractive interface", surface_shape)

        ccc = S4Conic.initialize_from_coefficients(surface_shape.get_conic_coefficients())

        mat1 = self.get_optical_element().get_material_object()
        mat2 = self.get_optical_element().get_material_image()
        print(">>>>>>> materials: ", mat1, mat2, type(mat1))

        if isinstance(mat1,float) or isinstance(mat1, int):
            refraction_index_object = mat1
        else:
            raise Exception("Not implemented")

        if isinstance(mat2,float) or isinstance(mat1, int):
            refraction_index_image = mat2
        else:
            raise Exception("Not implemented")

        beam_mirr, normal = ccc.apply_refraction_on_beam(beam, refraction_index_object, refraction_index_image)

        return beam_mirr, normal
