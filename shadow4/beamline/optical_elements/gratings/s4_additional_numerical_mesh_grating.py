import numpy
from syned.beamline.shape import NumericalMesh
from syned.beamline.element_coordinates import ElementCoordinates
from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.optical_elements.gratings.s4_grating import S4GratingElement, S4Grating
from shadow4.beamline.optical_elements.gratings.s4_numerical_mesh_grating import S4NumericalMeshGrating
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements


class S4AdditionalNumericalMeshGrating(S4NumericalMeshGrating):
    """
    Constructor.


    Parameters
    ----------
    ideal_grating : instance of S4Grating
        The grating baseline.
    numerical_mesh_grating : instance of S4NumericalMeshGrating
        The numerical mesh to be added to the ideal grating.
    name : str, optional
        The name of the grating.

    Returns
    -------
    instance of S4AdditionalNumericalMeshGrating.
    """
    def __init__(self,
                 ideal_grating : S4Grating = None,
                 numerical_mesh_grating : S4NumericalMeshGrating = None,
                 name="Grating with Additional Numerical Mesh"):

        if ideal_grating is not None:
            oe = ideal_grating
        else:
            oe = None

        S4NumericalMeshGrating.__init__(self, name=name,
                 boundary_shape=None if ideal_grating is None else ideal_grating.get_boundary_shape(),
                 xx=None if numerical_mesh_grating is None else numerical_mesh_grating._curved_surface_shape._xx,
                 yy=None if numerical_mesh_grating is None else numerical_mesh_grating._curved_surface_shape._yy,
                 zz=None if numerical_mesh_grating is None else numerical_mesh_grating._curved_surface_shape._zz,
                 surface_data_file="" if numerical_mesh_grating is None else numerical_mesh_grating._curved_surface_shape._surface_data_file,
                 # inputs related to grating
                 ruling=800e3 if   oe is None else oe._ruling,
                 ruling_coeff_linear=0.0 if   oe is None else oe._ruling_coeff_linear,
                 ruling_coeff_quadratic=0.0 if   oe is None else oe._ruling_coeff_quadratic,
                 ruling_coeff_cubic=0.0 if   oe is None else oe._ruling_coeff_cubic,
                 ruling_coeff_quartic=0.0 if   oe is None else oe._ruling_coeff_quartic,
                 order=0 if oe is None else oe._order,
                 f_ruling=0 if ideal_grating is None else ideal_grating._f_ruling,
                 )

        self.__ideal_grating         = ideal_grating
        self.__numerical_mesh_grating = numerical_mesh_grating

        self.__inputs = {
            "name": name,
            "ideal_grating": ideal_grating,
            "numerical_mesh_grating": numerical_mesh_grating,
        }

    def ideal_grating(self):
        """
        get the ideal optical element.

        Returns
        -------
        instance of S4Grating
        """
        return self.__ideal_grating

    def get_ideal(self):
        """
        get the ideal optical element.

        Returns
        -------
        instance of S4Grating
        """
        return self.__ideal_grating

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
        txt = self.__ideal_grating.to_python_code()
        txt += "ideal_grating = optical_element"
        txt += self.__numerical_mesh_grating.to_python_code()
        txt += "numerical_mesh_grating = optical_element"

        txt += self.to_python_code_boundary_shape()
        txt_pre = """

from shadow4.beamline.optical_elements.crystals.s4_additional_numerical_mesh_crystal import S4AdditionalNumericalMeshCrystal
optical_element = S4AdditionalNumericalMeshCrystal(name='{name:s}', ideal_crystal=ideal_crystal, numerical_mesh_crystal=numerical_mesh_crystal)
    """
        txt += txt_pre.format(**self.__inputs)
        return txt

class S4AdditionalNumericalMeshGratingElement(S4GratingElement):
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
    """
    def __init__(self,
                 optical_element: S4AdditionalNumericalMeshGrating = None,
                 coordinates: ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam: S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4AdditionalNumericalMeshGrating(),
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
        txt += "\nfrom shadow4.beamline.optical_elements.gratings.s4_additional_numerical_mesh_grating import S4AdditionalNumericalMeshGratingElement"
        txt += "\nbeamline_element = S4AdditionalNumericalMeshGratingElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)"
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
    from shadow4.beamline.optical_elements.gratings.s4_sphere_grating import S4SphereGrating

    optical_element = S4SphereGrating(name='xxx', boundary_shape=None,
                                        is_cylinder=1, cylinder_direction=0,
                                        convexity=1, radius=1.0,
                                        )

    from syned.beamline.element_coordinates import ElementCoordinates

    coordinates = ElementCoordinates(p=10, q=10, angle_radial=1.54985)


    if use_errors:
        from syned.beamline.shape import Rectangle
        mesh_element = S4AdditionalNumericalMeshGrating(name="C1",
                                                       ideal_grating=optical_element,
                                                       numerical_mesh_grating=S4NumericalMeshGrating(surface_data_file="/users/srio/Oasys/mirrors_branch3_mesh.hdf5",
                                                                                                   boundary_shape=Rectangle(x_left=-0.05,
                                                                                                                            x_right=0.05,
                                                                                                                            y_bottom=-0.5,
                                                                                                                            y_top=0.5)))

        beamline_element = S4AdditionalNumericalMeshGratingElement(optical_element=mesh_element,
                                                         coordinates=ElementCoordinates(p=10.0,
                                                                                        q=6.0,
                                                                                        angle_radial=numpy.radians(88.8)),
                                                        input_beam=beam)

    else:
        from shadow4.beamline.optical_elements.gratings.s4_sphere_grating import S4SphereGratingElement

        beamline_element = S4SphereGratingElement(optical_element=optical_element, coordinates=coordinates,
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

