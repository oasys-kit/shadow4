import numpy

from shadow4.syned.shape import Sphere, SphericalCylinder, Convexity, Direction

from shadow4.beamline.s4_optical_element import SurfaceCalculation, S4SphereOpticalElement
from shadow4.beamline.optical_elements.gratings.s4_grating import S4GratingElement, S4Grating, ElementCoordinates
from shadow4.optical_surfaces.s4_conic import S4Conic

class S4SphereGrating(S4Grating, S4SphereOpticalElement):
    def __init__(self,
                 name="Sphere Grating",
                 boundary_shape=None,
                 ruling=800e3,
                 ruling_coeff_linear=0.0,
                 ruling_coeff_quadratic=0.0,
                 ruling_coeff_cubic=0.0,
                 ruling_coeff_quartic=0.0,
                 coating=None,
                 coating_thickness=None,
                 f_central=False,
                 f_phot_cent=0,
                 phot_cent=8000.0,
                 f_reflec=0,
                 material_constants_library_flag=0,  # 0=xraylib, 1=dabax, 2=shadow preprocessor
                 file_refl="",
                 #
                 surface_calculation=SurfaceCalculation.EXTERNAL,
                 is_cylinder=False,
                 cylinder_direction=Direction.TANGENTIAL,
                 convexity=Convexity.UPWARD,
                 radius=1.0,
                 p_focus=0.0,
                 q_focus=0.0,
                 grazing_angle=0.0,
                 ):

        S4SphereOpticalElement.__init__(self, surface_calculation, is_cylinder, cylinder_direction, convexity,
                                        radius, p_focus, q_focus, grazing_angle)
        S4Grating.__init__(self,
                           name=name,
                           surface_shape=self.get_surface_shape_instance(),
                           boundary_shape=boundary_shape,
                           ruling=ruling,
                           ruling_coeff_linear=ruling_coeff_linear,
                           ruling_coeff_quadratic=ruling_coeff_quadratic,
                           ruling_coeff_cubic=ruling_coeff_cubic,
                           ruling_coeff_quartic=ruling_coeff_quartic,
                           coating=coating,
                           coating_thickness=coating_thickness,
                           f_central=False,
                           f_phot_cent=0,
                           phot_cent=8000.0,
                           f_reflec=0,
                           material_constants_library_flag=0,
                           file_refl="",
                           )

    # def get_optical_surface_instance(self):
    #     surface_shape = self.get_surface_shape()
    #
    #     switch_convexity = 0 if surface_shape.get_convexity() == Convexity.DOWNWARD else 1
    #
    #     if isinstance(surface_shape, SphericalCylinder):
    #         print(">>>>> SphericalCylinder mirror", surface_shape)
    #         cylindrical = 1
    #         cylangle = 0.0 if surface_shape.get_cylinder_direction() == Direction.TANGENTIAL else (0.5 * numpy.pi)
    #     elif isinstance(surface_shape, Sphere):
    #         print(">>>>> Sphere mirror", surface_shape)
    #         cylindrical = 0
    #         cylangle    = 0.0
    #
    #     return S4Conic.initialize_as_sphere_from_curvature_radius(surface_shape.get_radius(),
    #                                                              cylindrical=cylindrical, cylangle=cylangle, switch_convexity=switch_convexity)


    # def apply_grating_diffraction(self, beam):
    #     ccc = self.get_optical_surface_instance()
    #     print(">>>>>>> ccc: ", ccc.info())
    #     mirr, normal = ccc.apply_grating_diffraction_on_beam(beam)
    #     return mirr, normal

class S4SphereGratingElement(S4GratingElement):
    def __init__(self, optical_element=None, coordinates=None):
        super().__init__(optical_element if optical_element is not None else S4SphereMirror(),
                         coordinates if coordinates is not None else ElementCoordinates())
        if not (isinstance(self.get_optical_element().get_surface_shape(), SphericalCylinder) or
                isinstance(self.get_optical_element().get_surface_shape(), Sphere)):
            raise ValueError("Wrong Optical Element: only Sphere or Spherical Cylinder shape is accepted")

    # def apply_grating_diffraction(self, beam):
    #     return self.get_optical_element().apply_grating_diffraction(beam)

if __name__ == "__main__":

    from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical
    from shadow4.beam.beam import Beam

    #
    # source
    #
    src = SourceGeometrical(spatial_type="Point",
                    angular_distribution = "Flat",
                    energy_distribution = "Uniform",
                            )

    src.set_angular_distribution_flat(0,0,0,0)

    src.set_energy_distribution_uniform(value_min=999.8,value_max=1000.2,unit='eV')

    # print(src.info())

    beam = src.get_beam(N=5000, POL_DEG=1.0, POL_ANGLE=0.0, F_COHER=False)


    print(beam.info())

    # plotxy(Beam3.initialize_from_shadow4_beam(beam),1,3,nbins=100,title="SOURCE")

    #
    # grating
    #
    g = S4SphereGrating(
        name = "my_grating",
        boundary_shape = None, # BoundaryShape(),
        ruling = 800.0e3,
        ruling_coeff_linear = 0,
        ruling_coeff_quadratic = 0,
        ruling_coeff_cubic = 0,
        ruling_coeff_quartic = 0,
        coating = None,
        coating_thickness = None,
        f_central=False,
        f_phot_cent=0,
        phot_cent=8000.0,
        material_constants_library_flag=0,  # 0=xraylib, 1=dabax, 2=shadow preprocessor
        file_refl="",
        #
        surface_calculation=SurfaceCalculation.EXTERNAL,
        is_cylinder=False,
        cylinder_direction=Direction.TANGENTIAL,
        convexity=Convexity.UPWARD,
        radius=635757.0e-3,
        p_focus=0.0,
        q_focus=0.0,
        grazing_angle=0.0,
        )

    coordinates_syned = ElementCoordinates(p = 30.0,
                                           q = 9.93427,
                                           angle_radial = 87.29533343 * numpy.pi / 180,
                                           angle_radial_out= 89.10466657 * numpy.pi / 180,
                                           angle_azimuthal = 0.0)



    ge = S4SphereGratingElement(optical_element=g, coordinates=coordinates_syned)

    print(ge.info())

    beam_out = ge.trace_beam(beam)