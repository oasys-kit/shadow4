from shadow4.syned.shape import SurfaceData
from shadow4.beamline.optical_elements.mirrors.s4_mirror import S4MirrorElement, S4Mirror, ElementCoordinates

from shadow4.beamline.s4_optical_element import S4AdditiveSurfaceDataOpticalElement

class S4AdditiveSurfaceDataMirror(S4Mirror, S4AdditiveSurfaceDataOpticalElement):
    def __init__(self,
                 name="Surface Data Mirror",
                 boundary_shape=None,
                 xx=None,
                 yy=None,
                 zz=None,
                 surface_data_file=None,
                 # inputs related to mirror reflectivity
                 f_reflec=0,  # reflectivity of surface: 0=no reflectivity, 1=full polarization
                 f_refl=0,  # 0=prerefl file
                 # 1=electric susceptibility
                 # 2=user defined file (1D reflectivity vs angle)
                 # 3=user defined file (1D reflectivity vs energy)
                 # 4=user defined file (2D reflectivity vs energy and angle)
                 file_refl="",  # preprocessor file fir f_refl=0,2,3,4
                 refraction_index=1.0,  # refraction index (complex) for f_refl=1
                 base_surface_function=None,
                 ):
        S4AdditiveSurfaceDataOpticalElement.__init__(self, xx, yy, zz, surface_data_file, base_surface_function)
        S4Mirror.__init__(self, name, boundary_shape, self._curved_surface_shape,
                          f_reflec, f_refl, file_refl, refraction_index)

    def apply_geometrical_model(self, beam):
        num_mesh = self.get_optical_surface_instance()
        mirr, normal, _, _, _, _, _ = num_mesh.apply_specular_reflection_on_beam(beam)
        return mirr, normal

class S4AdditiveSurfaceDataMirrorElement(S4MirrorElement):
    def __init__(self, optical_element=None, coordinates=None):
        super().__init__(optical_element if optical_element is not None else S4AdditiveSurfaceDataMirror(),
                         coordinates if coordinates is not None else ElementCoordinates())
        if not isinstance(self.get_optical_element().get_surface_shape(), SurfaceData):
            raise ValueError("Wrong Optical Element: only Surface Data shape is accepted")



if __name__ == "__main__":
    import numpy

    from shadow4.sources.source_geometrical.source_gaussian import SourceGaussian

    from shadow4.beam.beam import Beam

    from shadow4.beamline.optical_elements.ideal_elements.s4_ideal_lens import S4IdealLens, S4IdealLensElement

    from shadow4.tools.graphics import plotxy

    from shadow4.syned.element_coordinates import ElementCoordinates

    from shadow4.beamline.optical_elements.mirrors.s4_toroidal_mirror import S4ToroidalMirror, S4ToroidalMirrorElement
    from shadow4.beamline.optical_elements.mirrors.s4_sphere_mirror import S4SphereMirror, S4SphereMirrorElement
    from shadow4.beamline.s4_optical_element import SurfaceCalculation
    from shadow4.beamline.optical_elements.mirrors.s4_surface_data_mirror import S4SurfaceDataMirror, \
        S4SurfaceDataMirrorElement

    from shadow4.syned.shape import Rectangle, Direction, Side  # TODO from syned.beamline.shape
    from srxraylib.plot.gol import set_qt

    set_qt()
    do_plot = True

    source = SourceGaussian.initialize_from_keywords(number_of_rays=10000,
                                                     sigmaX=0.0,
                                                     sigmaY=0.0,
                                                     sigmaZ=0.0,
                                                     sigmaXprime=1e-6,
                                                     sigmaZprime=1e-6, )
    beam0 = Beam()
    beam0.genSource(source)
    print(beam0.info())

    if do_plot:
        plotxy(beam0, 4, 6, title="Image 0", nbins=201)

    rlen1 = 0.6
    rlen2 = 0.6
    rwidx1 = 0.05
    rwidx2 = 0.05

    base_element = S4ToroidalMirror(name="M1b",
                                    surface_calculation=SurfaceCalculation.EXTERNAL,
                                    min_radius=0.157068,
                                    maj_radius=358.124803 - 0.157068,
                                    boundary_shape=Rectangle(x_left=-rwidx2, x_right=rwidx1, y_bottom=-rlen2,
                                                             y_top=rlen1))

    # base_element = S4SphereMirror(name="M1b",
    #                                 surface_calculation=SurfaceCalculation.EXTERNAL,
    #                                 radius=0.157068,
    #                                 boundary_shape=Rectangle(x_left=-rwidx2, x_right=rwidx1, y_bottom=-rlen2,
    #                                                          y_top=rlen1))

    base_surface_function = base_element.get_optical_surface_instance().surface_height

    mesh_element = S4AdditiveSurfaceDataMirror(name="M1",
                                       surface_data_file="../../../../oasys_workspaces/test_shadow4.hdf5",
                                       boundary_shape=Rectangle(x_left=-rwidx2, x_right=rwidx1, y_bottom=-rlen2,
                                                                y_top=rlen1),
                                       base_surface_function=base_surface_function)

    mirror1 = S4AdditiveSurfaceDataMirrorElement(optical_element=mesh_element,
                                                    coordinates=ElementCoordinates(p=10.0,
                                                                                   q=6.0,
                                                                                   angle_radial=numpy.radians(88.8)))

    print(mirror1.info())

    #
    # run
    #
    beam1, mirr1 = mirror1.trace_beam(beam0)
    print(mirr1.info())

    #
    # check
    #

    if do_plot:
        plotxy(beam1, 1, 3, title="Image 1", nbins=101, nolost=1)
        plotxy(mirr1, 2, 1, title="Footprint 1", nbins=101, nolost=1)


