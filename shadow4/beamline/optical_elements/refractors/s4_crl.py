import numpy

from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.optical_elements.refractors.lens import CRL

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.s4_beamline_element import S4BeamlineElement
from shadow4.beamline.optical_elements.refractors.s4_conic_interface import S4ConicInterface, S4ConicInterfaceElement
from shadow4.beamline.s4_optical_element_decorators import S4LensOpticalElementDecorator

class S4CRL(CRL, S4LensOpticalElementDecorator):
    def __init__(self,
                 name="Undefined",
                 n_lens=1,
                 surface_shape=None,   # now: 0=plane, 1=sphere, 2=parabola, 3=conic coefficients
                                       # (in shadow3: 1=sphere 4=paraboloid, 5=plane)
                 boundary_shape=None,  # syned stuff, replaces "diameter" in the shadow3 append_lens
                 material="",          # syned stuff, not (yet) used
                 thickness=0.0,        # syned stuff, lens thickness [m] (distance between the two interfaces at the center of the lenses)
                 piling_thickness=0.0,  # syned stuff,
                 convex_to_the_beam=1,  # for surface_shape: convexity of the first interface exposed to the beam 0=No, 1=Yes
                                        # the second interface has opposite convexity
                 cylinder_angle=0,      # for surface_shape: 0=not cylindricaL, 1=meridional 2=sagittal
                 ri_calculation_mode=0,   # source of refraction indices and absorption coefficients
                                 # 0=User
                                 # 1=prerefl file
                                 # 2=direct calculation using xraylib
                                 # 3=direct calculation using dabax
                 prerefl_file=None,    # for ri_calculation_mode=0: file name (from prerefl) to get the refraction index.
                 refraction_index=1.0, # for ri_calculation_mode=1: n (real)
                 attenuation_coefficient=0.0, # for ri_calculation_mode=1: mu in cm^-1 (real)
                 radius=500e-6,        # for surface_shape=(1,2): lens radius [m] (for spherical, or radius at the tip for paraboloid)
                 conic_coefficients=numpy.zeros(10),   # for surface_shape = 3: the conic coefficients
                 ):
        S4LensOpticalElementDecorator.__init__(self,
                                               surface_shape,
                                               convex_to_the_beam,
                                               cylinder_angle,
                                               ri_calculation_mode,
                                               prerefl_file,
                                               refraction_index,
                                               attenuation_coefficient,
                                               radius,
                                               conic_coefficients)

        surface_shape1, surface_shape2 = self.get_surface_shape_instance()

        CRL.__init__(self,
                     name=name,
                     n_lens=n_lens,
                     surface_shape1=surface_shape1,
                     surface_shape2=surface_shape2,
                     boundary_shape=boundary_shape,
                     material=material,
                     thickness=thickness,
                     piling_thickness=piling_thickness)

        self.__inputs = {
            "name": name,
            "n_lens" : n_lens,
            "boundary_shape": boundary_shape,
            "material": material,
            "thickness": thickness,
            "piling_thickness" : piling_thickness,
            "surface_shape": surface_shape,
            "convex_to_the_beam": convex_to_the_beam,
            "cylinder_angle": cylinder_angle,
            "ri_calculation_mode": ri_calculation_mode,
            "prerefl_file": prerefl_file,
            "refraction_index": refraction_index,
            "attenuation_coefficient": attenuation_coefficient,
            "radius": radius,
            "conic_coefficients": repr(conic_coefficients),
        }

    def to_python_code(self, **kwargs):
        return "# ** not implemented python code for S4Lens() **"

    def apply_geometrical_model(self, beam):
        pass

    def get_lens_interfaces(self):
        interface_shape_1, interface_shape_2 = self.get_optical_surface_instance()

        interfaces = numpy.full((self._n_lens, 2), None)

        for lens_index in range(self._n_lens):
            interfaces[lens_index, 0] = \
                S4ConicInterface(
                    name="First half-lens - Lens #" + str(lens_index+1),
                    boundary_shape=self.get_boundary_shape(),
                    material_object=None,
                    material_image=None,
                    f_r_ind=2,  # source of optical constants, from constant value or PREREFL preprocessor (file):
                    #      (0) constant value in both object and image spaces
                    #      (1) file in object space, constant value in image space
                    #      (2) constant value in object space, file in image space
                    #      (3) file in both object and image space
                    r_ind_obj=1.0,  # (for f_r_ind=0,2): index of refraction in object space.
                    r_ind_ima=1.0,  # (for f_r_ind=0,1): index of refraction in image space.
                    r_attenuation_obj=0.0,  # (for f_r_ind=0,2): attenuation coefficient in object space. Units of UserUnitLength^(-1)
                    r_attenuation_ima=0.0,  # (for f_r_ind=0,1): attenuation coefficient in image space. Units of UserUnitLength^(-1)
                    file_r_ind_obj="",  # (for f_r_ind=1,3): file generated by PREREFL
                    file_r_ind_ima=self._prerefl_file,  # (for f_r_ind=2,3): file generated by PREREFL
                    conic_coefficients=interface_shape_1.get_coefficients())

            interfaces[lens_index, 1] = \
                S4ConicInterface(
                    name="Second half-lens - Lens #" + str(lens_index+1),
                    boundary_shape=self.get_boundary_shape(),
                    material_object=None,
                    material_image=None,
                    f_r_ind=1,  # source of optical constants, from constant value or PREREFL preprocessor (file):
                    #      (0) constant value in both object and image spaces
                    #      (1) file in object space, constant value in image space
                    #      (2) constant value in object space, file in image space
                    #      (3) file in both object and image space
                    r_ind_obj=1.0,  # (for f_r_ind=0,2): index of refraction in object space.
                    r_ind_ima=1.0,  # (for f_r_ind=0,1): index of refraction in image space.
                    r_attenuation_obj=0.0,  # (for f_r_ind=0,2): attenuation coefficient in object space. Units of UserUnitLength^(-1)
                    r_attenuation_ima=0.0,  # (for f_r_ind=0,1): attenuation coefficient in image space. Units of UserUnitLength^(-1)
                    file_r_ind_obj=self._prerefl_file,  # (for f_r_ind=1,3): file generated by PREREFL
                    file_r_ind_ima="",  # (for f_r_ind=2,3): file generated by PREREFL
                    conic_coefficients=interface_shape_2.get_coefficients()
                )

        return interfaces


class S4CRLElement(S4BeamlineElement):
    def __init__(self,
                 optical_element: S4CRL = None,
                 coordinates: ElementCoordinates = None,
                 input_beam: S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4CRL(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         input_beam=input_beam)

    def to_python_code(self, **kwargs):
        return "# ** not implemented python code for S4LensElement() **"

    # def trace_beam(self, **params):
    #     # raise NotImplementedError()
    def trace_beam(self, **params):
        # p = self.get_coordinates().p()
        # q = self.get_coordinates().q()
        # theta_grazing1 = numpy.pi / 2 - self.get_coordinates().angle_radial()
        # theta_grazing2 = numpy.pi / 2 - self.get_coordinates().angle_radial_out()
        # alpha1 = self.get_coordinates().angle_azimuthal()
        #
        # #
        input_beam = self.get_input_beam().duplicate()
        oe         = self.get_optical_element()

        interfaces = oe.get_lens_interfaces()
        p, q, angle_radial, angle_radial_out, angle_azimuthal = self.get_coordinates().get_positions()
        n_lens = oe.get_n_lens()
        for lens_index in range(n_lens):
            if lens_index==0: source_plane = p
            else:             source_plane = oe.get_piling_thickness()
            if lens_index==n_lens-1: image_plane = q
            else:                    image_plane = 0.0

            coordinates_1 = ElementCoordinates(p=source_plane, q=oe.get_thickness()*0.5, angle_radial=angle_radial, angle_radial_out=numpy.pi,          angle_azimuthal=angle_azimuthal)
            coordinates_2 = ElementCoordinates(p=oe.get_thickness()*0.5, q=image_plane,  angle_radial=0.0,           angle_radial_out=angle_radial_out, angle_azimuthal=0.0)

            beamline_element_1 = S4ConicInterfaceElement(optical_element=interfaces[lens_index, 0], coordinates=coordinates_1, input_beam=input_beam)
            if lens_index==0: beam1, footprint1 = beamline_element_1.trace_beam()
            else:             beam1, _          = beamline_element_1.trace_beam()

            beamline_element_2 = S4ConicInterfaceElement(optical_element=interfaces[lens_index, 1], coordinates=coordinates_2, input_beam=beam1)
            if lens_index==n_lens-1: beam2, footprint2 = beamline_element_2.trace_beam()
            else:                    beam2, _          = beamline_element_2.trace_beam()

            if lens_index < n_lens-1: input_beam = beam2.duplicate()

        return beam2, [footprint1, footprint2]


if __name__ == "__main__":
    import numpy
    from shadow4.physical_models.prerefl.prerefl import PreRefl

    #
    # collimated source
    #
    from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical

    src = SourceGeometrical()
    src.set_energy_distribution_singleline(value=5000, unit='eV')
    src.set_spatial_type_rectangle(width=1e-3, height=1e-3)
    src.set_angular_distribution_uniform(0, 0, 0, 0)

    beam = src.get_beam()

    filename = "/Users/lrebuffi/Documents/Workspace/OASYS/shadow4/TEST/Be.dat"

    PreRefl.prerefl(interactive=False, SYMBOL="Be", FILE=filename, DENSITY=1.848, E_MIN=4500, E_MAX=5500, E_STEP=1)

    #
    # lens
    #
    lens = S4CRL(n_lens=3, piling_thickness=2.5e-3, prerefl_file=filename)
    e = S4CRLElement(optical_element=lens,
                      coordinates=ElementCoordinates(p=10, q=20, angle_radial=0, angle_azimuthal=0, angle_radial_out=numpy.pi),
                      input_beam=beam)
    print(e.to_python_code())

    e.trace_beam()
