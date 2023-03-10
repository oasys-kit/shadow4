from syned.beamline.shape import NumericalMesh
from syned.beamline.element_coordinates import ElementCoordinates
from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.optical_elements.mirrors.s4_mirror import S4MirrorElement, S4Mirror
from shadow4.beamline.optical_elements.mirrors.s4_numerical_mesh_mirror import S4NumericalMeshMirror

def add_mesh_to_ideal_surface(numerical_mesh, ideal_surface_ccc):
    return None

class S4AdditionalNumericalMeshMirror(S4NumericalMeshMirror):
    def __init__(self,
                 name="Mirror with Additional Numerical Mesh",
                 ideal_mirror : S4Mirror = None,
                 numerical_mesh_mirror : S4NumericalMeshMirror = None):
        S4NumericalMeshMirror.__init__(name)

        self.__ideal_mirror   = ideal_mirror
        self.__numerical_mesh_mirror = numerical_mesh_mirror

        # these attributes are necessary since they are called by the trace method.
        self._f_reflec         = self.__ideal_mirror.f_reflec
        self._f_refl           = self.__ideal_mirror.f_refl
        self._file_refl        = self.__ideal_mirror.file_refl
        self._refraction_index = self.__ideal_mirror.refraction_index

    def get_surface_shape(self):  return self.__numerical_mesh_mirror.get_surface_shape()

    def get_boundary_shape(self):
        # I think that the boundary shape should be dictated by the mirror and not by the error profile.
        # if the error profile is smaller, in the non-overlap area its value will be 0.
        return self.__ideal_mirror.get_boundary_shape()

    def set_boundaries_rectangle(self, x_left=-1e3, x_right=1e3, y_bottom=-1e3, y_top=1e3): # this method is for completeness
        self.__numerical_mesh_mirror.set_boundaries_rectangle(x_left=x_left, x_right=x_right, y_bottom=y_bottom, y_top=y_top)
        self.__ideal_mirror.set_boundaries_rectangle(x_left=x_left, x_right=x_right, y_bottom=y_bottom, y_top=y_top)

    def apply_geometrical_model(self, beam):
        numerical_mesh    = self.__numerical_mesh_mirror.get_optical_surface_instance()
        ideal_surface_ccc = self.__ideal_mirror.get_optical_surface_instance() # this mean that every S4Mirror must inherit from S4OpticalElementDecorator

        # here sum ideal surface to numerical mesh, and obtain a new numerical mesh:
        numerical_mesh = add_mesh_to_ideal_surface(numerical_mesh, ideal_surface_ccc)

        footprint, normal, _, _, _, _, _ = numerical_mesh.apply_specular_reflection_on_beam(beam)

        return footprint, normal

class S4AdditionalNumericalMeshMirrorElement(S4MirrorElement):
    def __init__(self,
                 optical_element: S4AdditionalNumericalMeshMirror = None,
                 coordinates: ElementCoordinates = None,
                 input_beam: S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4AdditionalNumericalMeshMirror(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         input_beam=input_beam)
        if not isinstance(self.get_optical_element().get_surface_shape(), NumericalMesh):
            raise ValueError("Wrong Optical Element: only Surface Data shape is accepted")


if __name__ == "__main__":
    import numpy

    from shadow4.sources.source_geometrical.source_gaussian import SourceGaussian
    from shadow4.beam.s4_beam import S4Beam
    from shadow4.tools.graphics import plotxy
    from shadow4.beamline.optical_elements.mirrors.s4_toroidal_mirror import S4ToroidalMirror
    from shadow4.beamline.s4_optical_element import SurfaceCalculation

    from syned.beamline.shape import Rectangle
    from srxraylib.plot.gol import set_qt

    set_qt()
    do_plot = True

    source = SourceGaussian.initialize_from_keywords(nrays=10000,
                                                     sigmaX=0.0,
                                                     sigmaY=0.0,
                                                     sigmaZ=0.0,
                                                     sigmaXprime=1e-6,
                                                     sigmaZprime=1e-6, )
    beam0 = S4Beam()
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

    mesh_element = S4AdditionalNumericalMeshMirror(name="M1",
                                                   ideal_mirror=base_element,
                                                   numerical_mesh_mirror=S4NumericalMeshMirror(surface_data_file="../../../../oasys_workspaces/test_shadow4.hdf5",
                                                                                               boundary_shape=Rectangle(x_left=-rwidx2,
                                                                                                                        x_right=rwidx1,
                                                                                                                        y_bottom=-rlen2,
                                                                                                                        y_top=rlen1)))

    mirror1 = S4AdditionalNumericalMeshMirrorElement(optical_element=mesh_element,
                                                     coordinates=ElementCoordinates(p=10.0,
                                                                                    q=6.0,
                                                                                    angle_radial=numpy.radians(88.8)),
                                                    input_beam=beam0)

    print(mirror1.info())

    #
    # run
    #
    beam1, mirr1 = mirror1.trace_beam()
    print(mirr1.info())

    #
    # check
    #

    if do_plot:
        plotxy(beam1, 1, 3, title="Image 1", nbins=101, nolost=1)
        plotxy(mirr1, 2, 1, title="Footprint 1", nbins=101, nolost=1)


