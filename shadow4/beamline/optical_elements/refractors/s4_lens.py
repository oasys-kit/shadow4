from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.optical_elements.refractors.lens import Lens
from syned.beamline.shape import Rectangle, Ellipse, Circle

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.s4_beamline_element import S4BeamlineElement
from shadow4.beamline.optical_elements.refractors.s4_conic_interface import S4ConicInterface, S4ConicInterfaceElement
from shadow4.beamline.s4_optical_element_decorators import S4RefractiveLensOpticalElementDecorator

class S4Lens(Lens, S4RefractiveLensOpticalElementDecorator):
    def __init__(self,
                 name="Undefined",
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
                 conic_coefficients=[0.0]*10,   # for surface_shape = 3: the conic coefficients
                 ):
        S4RefractiveLensOpticalElementDecorator.__init__(self,
                                                         surface_shape,
                                                         convex_to_the_beam,
                                                         cylinder_angle,
                                                         ri_calculation_mode,
                                                         prerefl_file,
                                                         refraction_index,
                                                         attenuation_coefficient,
                                                         radius,
                                                         conic_coefficients)
        Lens.__init__(self,
                      name=name,
                      surface_shape1=self.get_surface_shape_instance()[0],
                      surface_shape2=self.get_surface_shape_instance()[1],
                      boundary_shape=boundary_shape,
                      material=material,
                      thickness=thickness)

        self.__inputs = {
            "name": name,
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

    def to_python_code_boundary_shape(self):
        txt = ""
        bs = self._boundary_shape
        if bs is None:
            txt += "\nboundary_shape = None"
        elif isinstance(bs, Rectangle):
            txt += "\nfrom syned.beamline.shape import Rectangle"
            txt += "\nboundary_shape = Rectangle(x_left=%g, x_right=%g, y_bottom=%g, y_top=%g)" % bs.get_boundaries()
        elif isinstance(bs, Circle):
            txt += "\nfrom syned.beamline.shape import Circle"
            txt += "\nboundary_shape = Circle(radius=%g, x_center=%g,y_center=%g)" % bs.get_boundaries()
        elif isinstance(bs, Ellipse):
            txt += "\nfrom syned.beamline.shape import Ellipse"
            txt += "\nboundary_shape = Ellipse(a_axis_min=%g, a_axis_max=%g, b_axis_min=%g, b_axis_max=%g)" % bs.get_boundaries()
        return txt

    def to_python_code(self, **kwargs):
        txt = self.to_python_code_boundary_shape()

        txt_pre = """
from shadow4.beamline.optical_elements.refractors.s4_lens import S4Lens

optical_element = S4Lens(name='{name:s}',
     boundary_shape=boundary_shape,         # syned stuff, replaces "diameter" in the shadow3 append_lens
     material="", # syned stuff, not (yet) used
     thickness={thickness}, # syned stuff, lens thickness [m] (distance between the two interfaces at the center of the lenses)
     surface_shape={surface_shape}, # now: 0=plane, 1=sphere, 2=parabola, 3=conic coefficients
                                    # (in shadow3: 1=sphere 4=paraboloid, 5=plane)
     convex_to_the_beam={convex_to_the_beam}, # for surface_shape: convexity of the first interface exposed to the beam 0=No, 1=Yes
     cylinder_angle={cylinder_angle}, # for surface_shape: 0=not cylindricaL, 1=meridional 2=sagittal
     ri_calculation_mode={ri_calculation_mode},   # source of refraction indices and absorption coefficients
                                     # 0=User, 1=prerefl file, 2=xraylib, 3=dabax
     prerefl_file='{prerefl_file:s}', # for ri_calculation_mode=0: file name (from prerefl) to get the refraction index.
     refraction_index={refraction_index:g}, # for ri_calculation_mode=1: n (real)
     attenuation_coefficient={attenuation_coefficient:g}, # for ri_calculation_mode=1: mu in cm^-1 (real)
     radius={radius:g}, # for surface_shape=(1,2): lens radius [m] (for spherical, or radius at the tip for paraboloid)
     conic_coefficients=None, # for surface_shape = 3: the conic coefficients [todo: noy yet implemented]
     )
    """
        txt += txt_pre.format(**self.__inputs)
        return txt

    def get_lens_interfaces(self):
        return _get_lens_interfaces(lens_optical_surfaces=self.get_optical_surface_instance(),
                                    boundary_shape=self.get_boundary_shape(),
                                    ri_calculation_mode=self._ri_calculation_mode,
                                    refraction_index=self._refraction_index,
                                    attenuation_coefficient=self._attenuation_coefficient,
                                    prerefl_file=self._prerefl_file)

def _get_lens_interfaces(lens_optical_surfaces, boundary_shape, ri_calculation_mode, refraction_index, attenuation_coefficient, prerefl_file):
    conic_coefficients1 = lens_optical_surfaces[0].get_coefficients()
    conic_coefficients2 = lens_optical_surfaces[1].get_coefficients()

    if ri_calculation_mode == 0:
        half_lens_1 = S4ConicInterface(
            name="First half-lens",
            boundary_shape=boundary_shape,
            material_object=None,
            material_image=None,
            f_r_ind=0,  # source of optical constants, from constant value or PREREFL preprocessor (file):
            #      (0) constant value in both object and image spaces
            #      (1) file in object space, constant value in image space
            #      (2) constant value in object space, file in image space
            #      (3) file in both object and image space
            r_ind_obj=1.0,  # (for f_r_ind=0,2): index of refraction in object space.
            r_ind_ima=refraction_index,  # (for f_r_ind=0,1): index of refraction in image space.
            r_attenuation_obj=0.0,  # (for f_r_ind=0,2): attenuation coefficient in object space. Units of UserUnitLength^(-1)
            r_attenuation_ima=attenuation_coefficient,  # (for f_r_ind=0,1): attenuation coefficient in image space. Units of UserUnitLength^(-1)
            file_r_ind_obj="",  # (for f_r_ind=1,3): file generated by PREREFL
            file_r_ind_ima="",  # (for f_r_ind=2,3): file generated by PREREFL
            conic_coefficients=conic_coefficients1,
        )

        half_lens_2 = S4ConicInterface(
            name="Second half-lens",
            boundary_shape=boundary_shape,
            material_object=None,
            material_image=None,
            f_r_ind=0,  # source of optical constants, from constant value or PREREFL preprocessor (file):
            #      (0) constant value in both object and image spaces
            #      (1) file in object space, constant value in image space
            #      (2) constant value in object space, file in image space
            #      (3) file in both object and image space
            r_ind_obj=refraction_index,  # (for f_r_ind=0,2): index of refraction in object space.
            r_ind_ima=1.0,  # (for f_r_ind=0,1): index of refraction in image space.
            r_attenuation_obj=attenuation_coefficient,  # (for f_r_ind=0,2): attenuation coefficient in object space. Units of UserUnitLength^(-1)
            r_attenuation_ima=0.0,  # (for f_r_ind=0,1): attenuation coefficient in image space. Units of UserUnitLength^(-1)
            file_r_ind_obj="",  # (for f_r_ind=1,3): file generated by PREREFL
            file_r_ind_ima="",  # (for f_r_ind=2,3): file generated by PREREFL
            conic_coefficients=conic_coefficients2,
        )
    else:
        half_lens_1 = S4ConicInterface(
            name="First half-lens",
            boundary_shape=boundary_shape,
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
            file_r_ind_ima=prerefl_file,  # (for f_r_ind=2,3): file generated by PREREFL
            conic_coefficients=conic_coefficients1,
        )

        half_lens_2 = S4ConicInterface(
            name="Second half-lens",
            boundary_shape=boundary_shape,
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
            file_r_ind_obj=prerefl_file,  # (for f_r_ind=1,3): file generated by PREREFL
            file_r_ind_ima="",  # (for f_r_ind=2,3): file generated by PREREFL
            conic_coefficients=conic_coefficients2,
        )

    return half_lens_1, half_lens_2


class S4LensElement(S4BeamlineElement):
    def __init__(self,
                 optical_element : S4Lens = None,
                 coordinates : ElementCoordinates = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4Lens(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         input_beam=input_beam)

    def to_python_code(self, **kwargs):
        txt = "\n\n# optical element number XX"
        txt += self.get_optical_element().to_python_code()
        coordinates = self.get_coordinates()
        txt += "\nfrom syned.beamline.element_coordinates import ElementCoordinates"
        txt += "\nimport numpy"
        txt += "\ncoordinates = ElementCoordinates(p=%g, q=%g, angle_radial=0, angle_azimuthal=0, angle_radial_out=numpy.pi)" % \
               (coordinates.p(), coordinates.q())
        txt += "\nfrom shadow4.beamline.optical_elements.refractors.s4_lens import S4LensElement"
        txt += "\nbeamline_element = S4LensElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)"
        txt += "\n\nbeam, mirr = beamline_element.trace_beam()"
        return txt

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

        half_lens_1, half_lens_2 = oe.get_lens_interfaces()
        print(half_lens_1.inputs)

        p, q, angle_radial, angle_radial_out, angle_azimuthal = self.get_coordinates().get_positions()

        coordinates_1 = ElementCoordinates(p=p, q=oe.get_thickness() * 0.5, angle_radial=angle_radial, angle_radial_out=numpy.pi, angle_azimuthal=angle_azimuthal)
        coordinates_2 = ElementCoordinates(p=oe.get_thickness() * 0.5, q=q, angle_radial=0, angle_radial_out=angle_radial_out, angle_azimuthal=0)

        beamline_element_1 = S4ConicInterfaceElement(optical_element=half_lens_1, coordinates=coordinates_1, input_beam=input_beam)
        beam1, footprint1  = beamline_element_1.trace_beam()

        # return beam1, footprint1

        beamline_element_2 = S4ConicInterfaceElement(optical_element=half_lens_2, coordinates=coordinates_2, input_beam=beam1)
        beam2, footprint2  = beamline_element_2.trace_beam()

        return beam2, [footprint1, footprint2]


if __name__ == "__main__":

    import numpy

    #
    # collimated source
    #
    from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical
    src = SourceGeometrical(nrays=10000)
    src.set_energy_distribution_singleline(value=15000, unit='eV')
    src.set_spatial_type_rectangle(width=1e-3, height=1e-3)
    src.set_angular_distribution_uniform(0, 0, 0, 0)

    beam = src.get_beam()


    #
    # lens (Be)
    #
    lens = S4Lens(name="Undefined",
                 boundary_shape=None,  # syned stuff, replaces "diameter" in the shadow3 append_lens
                 material="",          # syned stuff, not (yet) used
                 thickness=0.03,       # syned stuff, lens thickness [m] (distance between the two interfaces at the center of the lenses)
                 surface_shape=1,      # now: 0=plane, 1=sphere, 2=parabola, 3=conic coefficients
                                       # (in shadow3: 1=sphere 4=paraboloid, 5=plane)
                 convex_to_the_beam=0, # for surface_shape (1,2): convexity of the first interface exposed to the beam 0=No, 1=Yes
                                       # the second interface has opposite convexity
                 cylinder_angle=0,     # for surface_shape (1,2): 0=not cylindricaL, 1=meridional 2=sagittal
                 ri_calculation_mode=1,       # source of refraction indices and absorption coefficients
                                 # 0=User
                                 # 1=prerefl file
                                 # 2=direct calculation using xraylib
                                 # 3=direct calculation using dabax
                 prerefl_file="/users/srio/Oasys/reflec.dat",    # for ri_calculation_mode=0: file name (from prerefl) to get the refraction index.
                 refraction_index=1.5, # for ri_calculation_mode=1: n (real)
                 attenuation_coefficient=1e-3, # for ri_calculation_mode=1: mu in cm^-1 (real)
                 radius=0.1,        # for surface_shape=(1,2): lens radius [m] (for spherical, or radius at the tip for paraboloid)
                 conic_coefficients=[0.0]*10,   # for surface_shape = 3: the conic coefficients)
                  )

    e = S4LensElement(optical_element=lens,
                      coordinates=ElementCoordinates(p=10, q=20,
                                                    angle_radial=0, angle_azimuthal=0, angle_radial_out=numpy.pi),
                      input_beam=beam)
    print(e.to_python_code())

    beam, footprints = e.trace_beam()
    print("Intensity: ", beam.intensity(nolost=1))

    # test plot
    if True:
        from srxraylib.plot.gol import plot_scatter

        # plot_scatter(beam.get_photon_energy_eV(nolost=1), beam.get_column(23, nolost=1),
        #              title='(Intensity,Photon Energy)', plot_histograms=0)
        plot_scatter(1e6 * beam.get_column(1, nolost=1), 1e6 * beam.get_column(3, nolost=1), title="(X,Z) in microns")
        plot_scatter(1e6 * beam.get_column(4, nolost=1), 1e6 * beam.get_column(6, nolost=1), title="(X',Z') in microns")

    print(lens.info())

    print(e.to_python_code())
