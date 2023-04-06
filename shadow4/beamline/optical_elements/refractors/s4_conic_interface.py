from syned.beamline.shape import Conic
from syned.beamline.element_coordinates import ElementCoordinates

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.optical_elements.refractors.s4_interface import S4InterfaceElement, S4Interface
from shadow4.optical_surfaces.s4_conic import S4Conic
from shadow4.beamline.s4_optical_element_decorators import S4ConicOpticalElementDecorator

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

        self.__inputs = {
             "name": name,
             "boundary_shape": boundary_shape,
             "material_object": material_object,
             "material_image": material_image,
             "f_r_ind": f_r_ind,
             "r_ind_obj": r_ind_obj,
             "r_ind_ima": r_ind_ima,
             "r_attenuation_obj": r_attenuation_obj,
             "r_attenuation_ima": r_attenuation_ima,
             "file_r_ind_obj": file_r_ind_obj,
             "file_r_ind_ima": file_r_ind_ima,
             "conic_coefficients": repr(conic_coefficients),
        }
    def to_python_code(self, **kwargs):
        txt = self.to_python_code_boundary_shape()
        txt_pre = """

from shadow4.beamline.optical_elements.refractors.s4_conic_interface import S4ConicInterface
optical_element = S4ConicInterface(name='{name:s}', boundary_shape=boundary_shape,
    material_object='{material_object:s}', material_image='{material_image:s}',
    f_r_ind={f_r_ind:g}, r_ind_obj={r_ind_obj:g}, r_ind_ima={r_ind_ima:g},
    r_attenuation_obj={r_attenuation_obj:g}, r_attenuation_ima={r_attenuation_ima:g},
    file_r_ind_obj='{file_r_ind_obj:s}', file_r_ind_ima='{file_r_ind_ima:s}',
    conic_coefficients={conic_coefficients:s})
"""
        txt += txt_pre.format(**self.__inputs)
        return txt


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

    def to_python_code(self, **kwargs):
        txt = "\n\n# optical element number XX"
        txt += self.get_optical_element().to_python_code()
        coordinates = self.get_coordinates()
        txt += "\nfrom syned.beamline.element_coordinates import ElementCoordinates"
        txt += "\ncoordinates = ElementCoordinates(p=%g, q=%g, angle_radial=%g, angle_azimuthal=%g, angle_radial_out=%g)" % \
               (coordinates.p(), coordinates.q(), coordinates.angle_radial(), coordinates.angle_azimuthal(), coordinates.angle_radial_out())
        txt += "\nfrom shadow4.beamline.optical_elements.refractors.s4_conic_interface import S4ConicInterfaceElement"
        txt += "\nbeamline_element = S4ConicInterfaceElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)"
        txt += "\n\nbeam, mirr = beamline_element.trace_beam()"
        return txt

if __name__ == "__main__":
    import numpy

    #
    # collimated source
    #
    from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical
    src = SourceGeometrical()
    src.set_energy_distribution_singleline(value=5000, unit='A')
    src.set_spatial_type_rectangle(width=1e-3, height=1e-3)
    src.set_angular_distribution_uniform(0, 0, 0, 0)

    beam = src.get_beam()

    print(beam.info())
    SX, SZ = (1e6 * beam.get_standard_deviation(1), 1e6 * beam.get_standard_deviation(3))

    #
    # interface definition
    #

    optical_element = S4ConicInterface(
        name="Conic Refractive Interface",
        boundary_shape=None,
        material_object="vacuum",
        material_image="glass",
        f_r_ind=0,
        r_ind_obj=1.0,
        r_ind_ima=1.5,
        conic_coefficients=[1.0, 1.0, 1.0, 0.0, -0.0, -0.0, 0.0, 0.0, 3350.0e-3, 0.0],
    )

    coordinates = ElementCoordinates(p=0.0, q=5000.0e-3,
                                       angle_radial=0.0, angle_azimuthal=0.0, angle_radial_out=numpy.pi)

    interface_element = S4ConicInterfaceElement(optical_element=optical_element,
                                         coordinates=coordinates,
                                         input_beam=beam)

    print(interface_element.info())
    print(interface_element.get_optical_element().get_surface_shape().get_conic_coefficients())
    #
    # trace
    #

    beam2, mirr2 = interface_element.trace_beam()
    #
    if False:
        from shadow4.tools.graphics import plotxy
        plotxy(beam2, 1, 3, nbins=100, title="FOCAL PLANE")
        plotxy(mirr2, 1, 3, nbins=100, title="LENS HEIGHT")
        # plotxy(mirr2, 4, 5, nbins=100, title="FOOT DIV")

    FX, FZ = (1e6 * beam2.get_standard_deviation(1), 1e6 * beam2.get_standard_deviation(3))
    print("Source dimensions: %f %f um" % (SX, SZ))
    print("Focal dimensions: %f %f um" % (FX, FZ))
    print("Demagnification: %g %g" % (SX / FX, SX / FZ))


    #
    # write script to file
    #
    if True:
        script = interface_element.to_python_code()
        f = open("tmp.py", 'w')
        f.write("""
from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical
src = SourceGeometrical()
src.set_energy_distribution_singleline(value=5000, unit='A')
src.set_spatial_type_rectangle(width=1e-3, height=1e-3)
src.set_angular_distribution_uniform(0, 0, 0, 0)

beam = src.get_beam()
""")
        f.write(script)

        f.write("""
from shadow4.tools.graphics import plotxy
plotxy(beam, 1, 3, nbins=100, title="FOCAL PLANE")
plotxy(mirr, 1, 3, nbins=100, title="LENS HEIGHT")
""")
        f.close()
        print("File tmp.py written to disk.")

