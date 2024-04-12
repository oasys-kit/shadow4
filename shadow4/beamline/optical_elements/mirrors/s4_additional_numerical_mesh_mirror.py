import numpy
from syned.beamline.shape import NumericalMesh
from syned.beamline.element_coordinates import ElementCoordinates
from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.optical_elements.mirrors.s4_mirror import S4MirrorElement, S4Mirror
from shadow4.beamline.optical_elements.mirrors.s4_numerical_mesh_mirror import S4NumericalMeshMirror
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements


class S4AdditionalNumericalMeshMirror(S4NumericalMeshMirror):
    """
    Constructor.


    Parameters
    ----------
    ideal_mirror : instance of S4Mirror
        The mirror baseline.
    numerical_mesh_mirror : instance of S4NumericalMeshMirror
        The numerical mesh to be added to the ideal mirror.
    name : str, optional
        The name of the mirror.

    Returns
    -------
    instance of S4AdditionalNumericalMeshMirror.
    """
    def __init__(self,
                 ideal_mirror : S4Mirror = None,
                 numerical_mesh_mirror : S4NumericalMeshMirror = None,
                 name="Mirror with Additional Numerical Mesh"):
        """

        """
        S4NumericalMeshMirror.__init__(self, name=name,
                 boundary_shape=None if ideal_mirror is None else ideal_mirror.get_boundary_shape(),
                 xx=None if numerical_mesh_mirror is None else numerical_mesh_mirror._curved_surface_shape._xx,
                 yy=None if numerical_mesh_mirror is None else numerical_mesh_mirror._curved_surface_shape._yy,
                 zz=None if numerical_mesh_mirror is None else numerical_mesh_mirror._curved_surface_shape._zz,
                 surface_data_file="" if numerical_mesh_mirror is None else numerical_mesh_mirror._curved_surface_shape._surface_data_file,
                 # inputs related to mirror reflectivity
                 f_reflec=0 if ideal_mirror is None else ideal_mirror._f_reflec,  # reflectivity of surface: 0=no reflectivity, 1=full polarization
                 f_refl=0 if ideal_mirror is None else ideal_mirror._f_refl,  # 0=prerefl file
                     # 1=electric susceptibility
                     # 2=user defined file (1D reflectivity vs angle)
                     # 3=user defined file (1D reflectivity vs energy)
                     # 4=user defined file (2D reflectivity vs energy and angle)
                     # 5=direct calculation using xraylib
                     # 6=direct calculation using dabax
                file_refl="" if ideal_mirror is None else ideal_mirror._file_refl,  # preprocessor file fir f_refl=0,2,3,4
                refraction_index=1+0j if ideal_mirror is None else ideal_mirror._refraction_index,  # refraction index (complex) for f_refl=1)
                coating_material= ""  if ideal_mirror is None else ideal_mirror._coating,   # string with coating material formula for f_refl=5,6
                coating_density=  1.0 if ideal_mirror is None else ideal_mirror._coating_density,    # coating material density for f_refl=5,6
                coating_roughness=0.0 if ideal_mirror is None else ideal_mirror._coating_roughness,  # coating material roughness in A for f_refl=5,6
                )

        self.__ideal_mirror          = ideal_mirror
        self.__numerical_mesh_mirror = numerical_mesh_mirror

        self.__inputs = {
            "name": name,
            "ideal_mirror": ideal_mirror,
            "numerical_mesh_mirror": numerical_mesh_mirror,
        }

    def ideal_mirror(self):
        """
        get the ideal optical element.

        Returns
        -------
        instance of S4Mirror
        """
        return self.__ideal_mirror

    def get_ideal(self):
        """
        get the ideal optical element.

        Returns
        -------
        instance of S4Mirror
        """
        return self.__ideal_mirror

    def to_python_code(self, **kwargs):
        """
        Creates the python code for defining the element.

        Parameters
        ----------
        **kwargs

        Returns
        -------
        str
            Python code.
        """
        txt = self.__ideal_mirror.to_python_code()
        txt += "ideal_mirror = optical_element"
        txt += self.__numerical_mesh_mirror.to_python_code()
        txt += "numerical_mesh_mirror = optical_element"

        txt += self.to_python_code_boundary_shape()
        txt_pre = """

from shadow4.beamline.optical_elements.mirrors.s4_additional_numerical_mesh_mirror import S4AdditionalNumericalMeshMirror
optical_element = S4AdditionalNumericalMeshMirror(name='{name:s}', ideal_mirror=ideal_mirror, numerical_mesh_mirror=numerical_mesh_mirror)
    """
        txt += txt_pre.format(**self.__inputs)
        return txt

    #
    # overwrite this method combining ideal shape + error shape
    #
    def _apply_mirror_reflection(self, beam):
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
    """
    Constructor.

    Parameters
    ----------
    optical_element : instance of OpticalElement, optional
        The syned optical element.
    coordinates : instance of ElementCoordinates, optional
        The syned element coordinates.
    movements : instance of S4BeamlineElementMovements, optional
        The S4 element movements.
    input_beam : instance of S4Beam, optional
        The S4 incident beam.

    Returns
    -------
    instance of S4AdditionalNumericalMeshMirrorElement
    """
    def __init__(self,
                 optical_element: S4AdditionalNumericalMeshMirror = None,
                 coordinates: ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam: S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4AdditionalNumericalMeshMirror(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         movements=movements,
                         input_beam=input_beam)
        if not isinstance(self.get_optical_element().get_surface_shape(), NumericalMesh):
            raise ValueError("Wrong Optical Element: only Surface Data shape is accepted")

    def to_python_code(self, **kwargs):
        """
        Creates the python code for defining the element.

        Parameters
        ----------
        **kwargs

        Returns
        -------
        str
            Python code.
        """
        txt = "\n\n# optical element number XX"
        txt += self.get_optical_element().to_python_code()
        txt += self.to_python_code_coordinates()
        txt += self.to_python_code_movements()
        txt += "\nfrom shadow4.beamline.optical_elements.mirrors.s4_additional_numerical_mesh_mirror import S4AdditionalNumericalMeshMirrorElement"
        txt += "\nbeamline_element = S4AdditionalNumericalMeshMirrorElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)"
        txt += "\n\nbeam, mirr = beamline_element.trace_beam()"
        return txt

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


    # from shadow4.beamline.optical_elements.mirrors.s4_numerical_mesh_mirror import S4NumericalMeshMirror, S4NumericalMeshMirrorElement
    # m = S4NumericalMeshMirror()
    # e = S4NumericalMeshMirrorElement()
    # m = S4AdditionalNumericalMeshMirror()
    # e = S4AdditionalNumericalMeshMirrorElement(None, None, None)

