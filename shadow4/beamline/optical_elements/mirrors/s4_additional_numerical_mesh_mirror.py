from syned.beamline.shape import NumericalMesh
from syned.beamline.element_coordinates import ElementCoordinates
from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.optical_elements.mirrors.s4_mirror import S4MirrorElement, S4Mirror
from shadow4.beamline.optical_elements.mirrors.s4_numerical_mesh_mirror import S4NumericalMeshMirror


class S4AdditionalNumericalMeshMirror(S4NumericalMeshMirror):
    def __init__(self,
                 name="Mirror with Additional Numerical Mesh",
                 ideal_mirror : S4Mirror = None,
                 numerical_mesh_mirror : S4NumericalMeshMirror = None):
        S4NumericalMeshMirror.__init__(self, name=name,
                 boundary_shape=numerical_mesh_mirror.get_boundary_shape(),
                 xx=numerical_mesh_mirror._curved_surface_shape._xx,
                 yy=numerical_mesh_mirror._curved_surface_shape._yy,
                 zz=numerical_mesh_mirror._curved_surface_shape._zz,
                 surface_data_file=numerical_mesh_mirror._curved_surface_shape._surface_data_file,
                 # inputs related to mirror reflectivity
                 f_reflec=ideal_mirror._f_reflec,  # reflectivity of surface: 0=no reflectivity, 1=full polarization
                 f_refl=ideal_mirror._f_refl,  # 0=prerefl file
                 # 1=electric susceptibility
                 # 2=user defined file (1D reflectivity vs angle)
                 # 3=user defined file (1D reflectivity vs energy)
                 # 4=user defined file (2D reflectivity vs energy and angle)
                 file_refl=ideal_mirror._file_refl,  # preprocessor file fir f_refl=0,2,3,4
                 refraction_index=ideal_mirror._refraction_index)  # refraction index (complex) for f_refl=1)

        self.__ideal_mirror   = ideal_mirror
        self.__numerical_mesh_mirror = numerical_mesh_mirror

        # these attributes are necessary since they are called by the trace method.
        # self._f_reflec         = self.__ideal_mirror.f_reflec
        # self._f_refl           = self.__ideal_mirror.f_refl
        # self._file_refl        = self.__ideal_mirror.file_refl
        # self._refraction_index = self.__ideal_mirror.refraction_index

    # def get_surface_shape(self):  return self.__numerical_mesh_mirror.get_surface_shape()

    # def get_boundary_shape(self):
    #     # I think that the boundary shape should be dictated by the mirror and not by the error profile.
    #     # if the error profile is smaller, in the non-overlap area its value will be 0.
    #
    #     # srio: It is defined like that in __init__
    #     return self.__ideal_mirror.get_boundary_shape()

    # def set_boundaries_rectangle(self, x_left=-1e3, x_right=1e3, y_bottom=-1e3, y_top=1e3): # this method is for completeness
    #     self.__numerical_mesh_mirror.set_boundaries_rectangle(x_left=x_left, x_right=x_right, y_bottom=y_bottom, y_top=y_top)
    #     self.__ideal_mirror.set_boundaries_rectangle(x_left=x_left, x_right=x_right, y_bottom=y_bottom, y_top=y_top)

        self.__inputs = {
            "name": name,
            "ideal_mirror": ideal_mirror,
            "numerical_mesh_mirror": numerical_mesh_mirror,
        }

    def to_python_code(self, data=None):
        txt = "\nfrom shadow4.beamline.optical_elements.mirrors.s4_additional_numerical_mesh_mirror import S4AdditionalNumericalmeshMirror"
        txt_pre = """
optical_element = S4AdditionalNumericalMeshMirror(name='{name:s}',boundary_shape=None,
    ideal_mirror=ideal_mirror,numerical_mesh_mirror=numerical_mesh_mirror)
    """
        txt += txt_pre.format(**self.__inputs)
        return txt


    def apply_geometrical_model(self, beam):
        # numerical_mesh    = self.__numerical_mesh_mirror.get_optical_surface_instance()
        numerical_mesh = self.get_optical_surface_instance()
        ideal = self.__ideal_mirror.get_optical_surface_instance()

        # here sum ideal surface to numerical mesh, and obtain a new numerical mesh:
        # numerical_mesh = add_mesh_to_ideal_surface(numerical_mesh, ideal_surface_ccc)
        x, y = numerical_mesh.get_mesh_x_y()
        X = numpy.outer(x, numpy.ones_like(y))
        Y = numpy.outer(numpy.ones_like(x), y)
        Z = ideal.surface_height(X,Y)
        numerical_mesh.add_to_mesh(Z)
        # ideal_surface_ccc = self.__ideal_mirror.get_optical_surface_instance() # this mean that every S4Mirror must inherit from S4OpticalElementDecorator

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

    def to_python_code(self, data=None):
        txt = "\n\n# optical element number XX"
        txt += self.get_optical_element().to_python_code()
        coordinates = self.get_coordinates()
        txt += "\nfrom syned.beamline.element_coordinates import ElementCoordinates"
        txt += "\ncoordinates=ElementCoordinates(p=%g,q=%g,angle_radial=%g)" % \
               (coordinates.p(), coordinates.q(), coordinates.angle_radial())
        txt += "\nfrom shadow4.beamline.optical_elements.mirrors.s4_additional_numerical_mesh_mirror import S4AdditionalNumericalMeshMirrorElement"
        txt += "\nbeamline_element = S4AdditionalNumericalMeshMirrorElement(optical_element=optical_element,coordinates=coordinates,input_beam=beam)"
        txt += "\n\nbeam, mirr = beamline_element.trace_beam()"
        return txt

    def duplicate(self):
        return S4AdditionalNumericalMeshMirrorElement(optical_element=self.duplicate_coordinates(),
                                coordinates=self.duplicate_coordinates(),
                                input_beam=self.duplicate_input_beam())

if __name__ == "__main__":
    import numpy
    do_plot = True
    use_errors = True
    from srxraylib.plot.gol import plot_scatter

    #
    #
    #
    from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical

    light_source = SourceGeometrical(name='SourceGeometrical', nrays=10000, seed=5676561)
    light_source.set_spatial_type_gaussian(sigma_h=5e-06, sigma_v=0.000001)
    light_source.set_angular_distribution_gaussian(sigdix=0.000001, sigdiz=0.000001)
    light_source.set_energy_distribution_singleline(1.000000, unit='A')
    light_source.set_polarization(polarization_degree=1.000000, phase_diff=0.000000, coherent_beam=0)
    beam = light_source.get_beam()

    # optical element number XX
    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirror

    optical_element = S4EllipsoidMirror(name='Ellipsoid Mirror', boundary_shape=None,
                                        surface_calculation=0, is_cylinder=1, cylinder_direction=0,
                                        convexity=1, min_axis=0.000000, maj_axis=0.000000, p_focus=10.000000,
                                        q_focus=10.000000,
                                        grazing_angle=0.020944,
                                        f_reflec=0, f_refl=1, file_refl='<none>', refraction_index=0.99999 + 0.001j)

    from syned.beamline.element_coordinates import ElementCoordinates

    coordinates = ElementCoordinates(p=10, q=10, angle_radial=1.54985)


    if use_errors:
        from syned.beamline.shape import Rectangle
        mesh_element = S4AdditionalNumericalMeshMirror(name="M1",
                                                       ideal_mirror=optical_element,
                                                       numerical_mesh_mirror=S4NumericalMeshMirror(surface_data_file="/users/srio/Oasys/mirrors_branch3_mesh.hdf5",
                                                                                                   boundary_shape=Rectangle(x_left=-0.05,
                                                                                                                            x_right=0.05,
                                                                                                                            y_bottom=-0.5,
                                                                                                                            y_top=0.5)))

        beamline_element = S4AdditionalNumericalMeshMirrorElement(optical_element=mesh_element,
                                                         coordinates=ElementCoordinates(p=10.0,
                                                                                        q=6.0,
                                                                                        angle_radial=numpy.radians(88.8)),
                                                        input_beam=beam)

    else:
        from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirrorElement

        beamline_element = S4EllipsoidMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                                    input_beam=beam)
    #
    # run
    #
    beam1, mirr1 = beamline_element.trace_beam()
    #
    #
    #
    if do_plot:

        plot_scatter(1e6*beam1.get_column(1), 1e6*beam1.get_column(3), title="Cols 1,3 / um")



