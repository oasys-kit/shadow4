import numpy
from syned.beamline.shape import NumericalMesh
from syned.beamline.element_coordinates import ElementCoordinates
from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.optical_elements.multilayers.s4_multilayer import S4MultilayerElement, S4Multilayer
from shadow4.beamline.optical_elements.multilayers.s4_numerical_mesh_multilayer import S4NumericalMeshMultilayer
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements


class S4AdditionalNumericalMeshMultilayer(S4NumericalMeshMultilayer):
    """
    Constructor.


    Parameters
    ----------
    ideal_multilayer : instance of S4Multilayer
        The multilayer baseline.
    numerical_mesh_multilayer : instance of S4NumericalMeshMultilayer
        The numerical mesh to be added to the ideal multilayer.
    name : str, optional
        The name of the multilayer.

    Returns
    -------
    instance of S4AdditionalNumericalMeshMultilayer.
    """
    def __init__(self,
                 ideal_multilayer : S4Multilayer = None,
                 numerical_mesh_multilayer : S4NumericalMeshMultilayer = None,
                 name="Multilayer with Additional Numerical Mesh"):
        """

        """
        S4NumericalMeshMultilayer.__init__(self, name=name,
                 boundary_shape=None if ideal_multilayer is None else ideal_multilayer.get_boundary_shape(),
                 xx=None if numerical_mesh_multilayer is None else numerical_mesh_multilayer._curved_surface_shape._xx,
                 yy=None if numerical_mesh_multilayer is None else numerical_mesh_multilayer._curved_surface_shape._yy,
                 zz=None if numerical_mesh_multilayer is None else numerical_mesh_multilayer._curved_surface_shape._zz,
                 surface_data_file="" if numerical_mesh_multilayer is None else numerical_mesh_multilayer._curved_surface_shape._surface_data_file,
                 # inputs related to multilayer reflectivity
                 f_refl=0 if ideal_multilayer is None else ideal_multilayer._f_refl,  # 0=pre_mlayer file
                             # 0=pre_mlayer file
                             # 1=user defined file (1D reflectivity vs angle)
                             # 2=user defined file (1D reflectivity vs energy)
                             # 3=user defined file (2D reflectivity vs energy and angle)
                             # 4=direct calculation using xraylib
                             # 5=direct calculation using dabax
                 file_refl="" if ideal_multilayer is None else ideal_multilayer._file_refl,  # reflectivity file for f_refl=0,1,2,3
                 structure='[B/W]x50+Si' if ideal_multilayer is None else ideal_multilayer._structure,
                 period=25.0 if ideal_multilayer is None else ideal_multilayer._period,
                 Gamma=0.5 if ideal_multilayer is None else ideal_multilayer._Gamma,
                )

        self.__ideal_multilayer          = ideal_multilayer
        self.__numerical_mesh_multilayer = numerical_mesh_multilayer

        self.__inputs = {
            "name": name,
            "ideal_multilayer": ideal_multilayer,
            "numerical_mesh_multilayer": numerical_mesh_multilayer,
        }

    def ideal_multilayer(self):
        """
        get the ideal optical element.

        Returns
        -------
        instance of S4Multilayer
        """
        return self.__ideal_multilayer

    def get_ideal(self):
        """
        get the ideal optical element.

        Returns
        -------
        instance of S4Multilayer
        """
        return self.__ideal_multilayer

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
        txt = self.__ideal_multilayer.to_python_code()
        txt += "ideal_multilayer = optical_element"
        txt += self.__numerical_mesh_multilayer.to_python_code()
        txt += "numerical_mesh_multilayer = optical_element"

        txt += self.to_python_code_boundary_shape()
        txt_pre = """

from shadow4.beamline.optical_elements.multilayers.s4_additional_numerical_mesh_multilayer import S4AdditionalNumericalMeshMultilayer
optical_element = S4AdditionalNumericalMeshMultilayer(name='{name:s}', ideal_multilayer=ideal_multilayer, numerical_mesh_multilayer=numerical_mesh_multilayer)
    """
        txt += txt_pre.format(**self.__inputs)
        return txt

    #
    # overwrite this method combining ideal shape + error shape
    #
    def _apply_multilayer_reflection(self, beam):
        # numerical_mesh    = self.__numerical_mesh_multilayer.get_optical_surface_instance()
        numerical_mesh = self.get_optical_surface_instance()
        ideal = self.__ideal_multilayer.get_optical_surface_instance()
        # here sum ideal surface to numerical mesh, and obtain a new numerical mesh:
        # numerical_mesh = add_mesh_to_ideal_surface(numerical_mesh, ideal_surface_ccc)
        x, y = numerical_mesh.get_mesh_x_y()
        X = numpy.outer(x, numpy.ones_like(y))
        Y = numpy.outer(numpy.ones_like(x), y)
        Z = ideal.surface_height(X,Y)
        numerical_mesh.add_to_mesh(Z)
        # ideal_surface_ccc = self.__ideal_multilayer.get_optical_surface_instance() # this mean that every S4Multilayer must inherit from S4OpticalElementDecorator
        footprint, normal, _, _, _, _, _ = numerical_mesh.apply_specular_reflection_on_beam(beam)

        return footprint, normal

class S4AdditionalNumericalMeshMultilayerElement(S4MultilayerElement):
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
    instance of S4AdditionalNumericalMeshMultilayerElement
    """
    def __init__(self,
                 optical_element: S4AdditionalNumericalMeshMultilayer = None,
                 coordinates: ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam: S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4AdditionalNumericalMeshMultilayer(),
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
        txt += "\nfrom shadow4.beamline.optical_elements.multilayers.s4_additional_numerical_mesh_multilayer import S4AdditionalNumericalMeshMultilayerElement"
        txt += "\nbeamline_element = S4AdditionalNumericalMeshMultilayerElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)"
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
    from shadow4.beamline.optical_elements.multilayers.s4_ellipsoid_multilayer import S4EllipsoidMultilayer

    optical_element = S4EllipsoidMultilayer(name='Ellipsoid Multilayer', boundary_shape=None,
                                        surface_calculation=0, is_cylinder=1, cylinder_direction=0,
                                        convexity=1, min_axis=0.000000, maj_axis=0.000000, p_focus=10.000000,
                                        q_focus=10.000000,
                                        grazing_angle=0.020944,
                                        f_refl=1, file_refl='<none>')

    from syned.beamline.element_coordinates import ElementCoordinates

    coordinates = ElementCoordinates(p=10, q=10, angle_radial=1.54985)


    if use_errors:
        from syned.beamline.shape import Rectangle
        mesh_element = S4AdditionalNumericalMeshMultilayer(name="M1",
                                                       ideal_multilayer=optical_element,
                                                       numerical_mesh_multilayer=S4NumericalMeshMultilayer(surface_data_file="/users/srio/Oasys/multilayers_branch3_mesh.hdf5",
                                                                                                   boundary_shape=Rectangle(x_left=-0.05,
                                                                                                                            x_right=0.05,
                                                                                                                            y_bottom=-0.5,
                                                                                                                            y_top=0.5)))

        beamline_element = S4AdditionalNumericalMeshMultilayerElement(optical_element=mesh_element,
                                                         coordinates=ElementCoordinates(p=10.0,
                                                                                        q=6.0,
                                                                                        angle_radial=numpy.radians(88.8)),
                                                        input_beam=beam)

    else:
        from shadow4.beamline.optical_elements.multilayers.s4_ellipsoid_multilayer import S4EllipsoidMultilayerElement

        beamline_element = S4EllipsoidMultilayerElement(optical_element=optical_element, coordinates=coordinates,
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


    # from shadow4.beamline.optical_elements.multilayers.s4_numerical_mesh_multilayer import S4NumericalMeshMultilayer, S4NumericalMeshMultilayerElement
    # m = S4NumericalMeshMultilayer()
    # e = S4NumericalMeshMultilayerElement()
    # m = S4AdditionalNumericalMeshMultilayer()
    # e = S4AdditionalNumericalMeshMultilayerElement(None, None, None)

