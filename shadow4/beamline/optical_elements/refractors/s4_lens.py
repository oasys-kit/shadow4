import numpy

from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.optical_elements.refractors.lens import Lens

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.s4_beamline_element import S4BeamlineElement
from shadow4.beamline.optical_elements.refractors.s4_conic_interface import S4ConicInterface, S4ConicInterfaceElement

class S4Lens(Lens):

    def __init__(self,
                 name="Undefined",
                 surface_shape1=None,  # syned stuff
                 surface_shape2=None,  # syned stuff, not (yet) used
                 boundary_shape=None,  # syned stuff, replaces "diameter" in the shadow3 append_lens
                 material="",          # syned stuff, not (yet) used
                 thickness=0.0,        # syned stuff, lens thickness [m] (distance between the two interfaces at the center of the lenses)
                 surface_shape=1,      # now: 0=plane, 1=sphere, 2=parabola, 3=conic coefficients
                                       # (in shadow3: 1=sphere 4=paraboloid, 5=plane)
                 convex_to_the_beam=1, # for surface_shape (1,2): convexity of the first interface exposed to the beam 0=No, 1=Yes
                                       # the second interface has opposite convexity
                 cylinder_angle=0,     # for surface_shape (1,2): 0=not cylindricaL, 1=meridional 2=sagittal
                 ri_calculation_mode=0,       # source of refraction indices and absorption coefficients
                                 # 0=User
                                 # 1=prerefl file
                                 # 2=direct calculation using xraylib
                                 # 3=direct calculation using dabax
                 prerefl_file=None,    # for ri_calculation_mode=0: file name (from prerefl) to get the refraction index.
                 refraction_index=1.0, # for ri_calculation_mode=1: n (real)
                 attenuation_coefficient=0.0, # for ri_calculation_mode=1: mu in cm^-1 (real)
                 radius=500e-6,        # for surface_shape=(1,2): lens radius [m] (for spherical, or radius at the tip for paraboloid)
                 conic_coefficients=[0]*10,   # for surface_shape = 3: the conic coefficients
                 ):




        Lens.__init__(self,
                      name=name,
                      surface_shape1=surface_shape1,
                      surface_shape2=surface_shape2,
                      boundary_shape=boundary_shape,
                      material=material,
                      thickness=thickness,
                        )

        self._boundary_shape = boundary_shape
        self._surface_shape = surface_shape
        self._convex_to_the_beam = convex_to_the_beam
        self._cylinder_angle = cylinder_angle
        self._ri_calculation_mode = ri_calculation_mode
        self._prerefl_file = prerefl_file
        self._refraction_index = refraction_index
        self._attenuation_coefficient = attenuation_coefficient
        self._radius = radius
        self._conic_coefficients = conic_coefficients

        self.__inputs = {
            "name": name,
            "surface_shape1":          surface_shape1,
            "surface_shape2":          surface_shape2,
            "boundary_shape":          boundary_shape,
            "material":                material,
            "thickness":               thickness,
            "surface_shape":           surface_shape,
            "convex_to_the_beam":      convex_to_the_beam,
            "cylinder_angle":          cylinder_angle,
            "ri_calculation_mode":     ri_calculation_mode,
            "prerefl_file":            prerefl_file,
            "refraction_index":        refraction_index,
            "attenuation_coefficient": attenuation_coefficient,
            "radius":                  radius,
            "conic_coefficients":      repr(conic_coefficients),
        }

    def to_python_code(self, **kwargs):
        return "# ** not implemented python code for S4Lens() **"


    def apply_geometrical_model(self, beam):
        pass


class S4LensElement(S4BeamlineElement):
    def __init__(self,
                 optical_element : S4Lens = None,
                 coordinates : ElementCoordinates = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4Lens(),
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
        oe = self.get_optical_element()

        #
        # first half lens
        #
        if oe._surface_shape == 0: # plane
            conic_coefficients = [0,0,0,0,0,0,0,0,-1,0]
        elif oe._surface_shape == 1: # sphere
            conic_coefficients = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0 * oe._radius, 0.0]
        elif oe._surface_shape == 2:  # parabola
            conic_coefficients = [1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0 * oe._radius, 0.0]
        elif oe._surface_shape == 3:  # conic coefficients
            conic_coefficients = oe._cocic_coefficients

        if oe._cylinder_angle == 0:
            pass
        elif oe._cylinder_angle == 1: # meridional: no dependence in x
            conic_coefficients[0] = 0
            conic_coefficients[3] = 0
            conic_coefficients[5] = 0
            conic_coefficients[6] = 0
        elif oe._cylinder_angle == 2: # sagittal: no dependence in y
            conic_coefficients[1] = 0
            conic_coefficients[3] = 0
            conic_coefficients[4] = 0
            conic_coefficients[7] = 0

        if oe._convex_to_the_beam: conic_coefficients[8] *= -1

        half_lens_1 = S4ConicInterface(
            name="First half-lens",
            boundary_shape=oe._boundary_shape,
            material_object=None,
            material_image=None,
            f_r_ind=2,               # source of optical constants, from constant value or PREREFL preprocessor (file):
                                     #      (0) constant value in both object and image spaces
                                     #      (1) file in object space, constant value in image space
                                     #      (2) constant value in object space, file in image space
                                     #      (3) file in both object and image space
            r_ind_obj=1.0,           # (for f_r_ind=0,2): index of refraction in object space.
            r_ind_ima=1.0,           # (for f_r_ind=0,1): index of refraction in image space.
            r_attenuation_obj=0.0,   # (for f_r_ind=0,2): attenuation coefficient in object space. Units of UserUnitLength^(-1)
            r_attenuation_ima=0.0,   # (for f_r_ind=0,1): attenuation coefficient in image space. Units of UserUnitLength^(-1)
            file_r_ind_obj="",       # (for f_r_ind=1,3): file generated by PREREFL
            file_r_ind_ima=oe._prerefl_file,       # (for f_r_ind=2,3): file generated by PREREFL
            conic_coefficients=conic_coefficients,
        )

        conic_coefficients2 = conic_coefficients.copy()
        conic_coefficients2[8] *= -1

        half_lens_2 = S4ConicInterface(
            name="Second half-lens",
            boundary_shape=oe._boundary_shape,
            material_object=None,
            material_image=None,
            f_r_ind=1,               # source of optical constants, from constant value or PREREFL preprocessor (file):
                                     #      (0) constant value in both object and image spaces
                                     #      (1) file in object space, constant value in image space
                                     #      (2) constant value in object space, file in image space
                                     #      (3) file in both object and image space
            r_ind_obj=1.0,           # (for f_r_ind=0,2): index of refraction in object space.
            r_ind_ima=1.0,           # (for f_r_ind=0,1): index of refraction in image space.
            r_attenuation_obj=0.0,   # (for f_r_ind=0,2): attenuation coefficient in object space. Units of UserUnitLength^(-1)
            r_attenuation_ima=0.0,   # (for f_r_ind=0,1): attenuation coefficient in image space. Units of UserUnitLength^(-1)
            file_r_ind_obj=oe._prerefl_file,       # (for f_r_ind=1,3): file generated by PREREFL
            file_r_ind_ima="",       # (for f_r_ind=2,3): file generated by PREREFL
            conic_coefficients=conic_coefficients2,
        )

        coordinates = self.get_coordinates()
        print(coordinates)
        p, q, angle_radial, angle_radial_out, angle_azimuthal = coordinates.get_positions()

        coordinates_1 = ElementCoordinates(p=p, q=oe._thickness/2,
                                angle_radial=angle_radial, angle_radial_out=numpy.pi, angle_azimuthal=angle_azimuthal)

        coordinates_2 = ElementCoordinates(p=oe._thickness/2, q=q,
                                angle_radial=0, angle_radial_out=angle_radial_out, angle_azimuthal=0)

        beamline_element_1 = S4ConicInterfaceElement(optical_element=half_lens_1,
                                                     coordinates=coordinates_1,
                                                     input_beam=beam)

        beam1, footprint1 = beamline_element_1.trace_beam()

        beamline_element_2 = S4ConicInterfaceElement(optical_element=half_lens_2,
                                                     coordinates=coordinates_2,
                                                     input_beam=beam)

        beam2, footprint2 = beamline_element_2.trace_beam()

        return beam2, footprint2


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


    #
    # lens
    #
    lens = S4Lens()
    e = S4LensElement(optical_element=lens,
                      coordinates=ElementCoordinates(p=10, q=20,
                                                    angle_radial=0, angle_azimuthal=0, angle_radial_out=numpy.pi),
                      input_beam=beam)
    print(e.to_python_code())

    e.trace_beam()
