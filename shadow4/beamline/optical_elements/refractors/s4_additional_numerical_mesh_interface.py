import numpy
from syned.beamline.shape import NumericalMesh
from syned.beamline.element_coordinates import ElementCoordinates
from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.optical_elements.refractors.s4_interface import S4InterfaceElement, S4Interface
from shadow4.beamline.optical_elements.refractors.s4_numerical_mesh_interface import S4NumericalMeshInterface
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements


class S4AdditionalNumericalMeshInterface(S4NumericalMeshInterface):
    """
    Constructor.


    Parameters
    ----------
    ideal_interface : instance of S4Interface
        The refractive interface baseline.
    numerical_mesh_interface : instance of S4NumericalMeshInterface
        The numerical mesh to be added to the ideal refractive interface.
    name : str, optional
        The name of the refractive interface.

    Returns
    -------
    instance of S4AdditionalNumericalMeshInterface.
    """
    def __init__(self,
                 ideal_interface : S4Interface = None,
                 numerical_mesh_interface : S4NumericalMeshInterface = None,
                 name="Refractive interface with Additional Numerical Mesh",
                 ):
        """

        """
        S4NumericalMeshInterface.__init__(self, name=name,
            xx=None if numerical_mesh_interface is None else numerical_mesh_interface._curved_surface_shape._xx,
            yy=None if numerical_mesh_interface is None else numerical_mesh_interface._curved_surface_shape._yy,
            zz=None if numerical_mesh_interface is None else numerical_mesh_interface._curved_surface_shape._zz,
            surface_data_file="" if numerical_mesh_interface is None else numerical_mesh_interface._curved_surface_shape._surface_data_file,
            # kwds for refraction index, copied from ideal
            f_r_ind            = 0    if ideal_interface is None else ideal_interface._f_r_ind,
            r_ind_obj          = 1.0  if ideal_interface is None else ideal_interface._r_ind_obj,
            r_ind_ima          = 1.0  if ideal_interface is None else ideal_interface._r_ind_ima,
            r_attenuation_obj  = 0.0  if ideal_interface is None else ideal_interface._r_attenuation_obj,
            r_attenuation_ima  = 0.0  if ideal_interface is None else ideal_interface._r_attenuation_ima,
            file_r_ind_obj     = ''   if ideal_interface is None else ideal_interface._file_r_ind_obj,
            file_r_ind_ima     = ''   if ideal_interface is None else ideal_interface._file_r_ind_ima,
            material_object    = ''   if ideal_interface is None else ideal_interface._material_object,
            material_image     = ''   if ideal_interface is None else ideal_interface._material_image,
            density_object     = 1.0  if ideal_interface is None else ideal_interface._density_object,
            density_image      = 1.0  if ideal_interface is None else ideal_interface._density_image,
            dabax              = None if ideal_interface is None else ideal_interface._dabax,
            )

        self.__ideal_interface          = ideal_interface
        self.__numerical_mesh_interface = numerical_mesh_interface

        self.__inputs = {
            "name": name,
            "ideal_interface": ideal_interface,
            "numerical_mesh_interface": numerical_mesh_interface,
        }

    def ideal_interface(self):
        """
        get the ideal optical element.

        Returns
        -------
        instance of S4Interface
        """
        return self.__ideal_interface

    def get_ideal(self):
        """
        get the ideal optical element.

        Returns
        -------
        instance of S4Interface
        """
        return self.__ideal_interface

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
        txt = self.__ideal_interface.to_python_code()
        txt += "ideal_interface = optical_element"
        txt += self.__numerical_mesh_interface.to_python_code()
        txt += "numerical_mesh_interface = optical_element"

        txt += self.to_python_code_boundary_shape()
        txt_pre = """

from shadow4.beamline.optical_elements.refractors.s4_additional_numerical_mesh_interface import S4AdditionalNumericalMeshInterface
optical_element = S4AdditionalNumericalMeshInterface(name='{name:s}', ideal_interface=ideal_interface, numerical_mesh_interface=numerical_mesh_interface)
    """
        txt += txt_pre.format(**self.__inputs)
        return txt

    #
    # overwrite this method combining ideal shape + error shape
    #
    def _apply_interface_refraction(self, beam,
                                    refraction_index_object, refraction_index_image, mu_object,
                                    apply_attenuation=1):

        numerical_mesh = self.get_optical_surface_instance()
        ideal = self.__ideal_interface.get_optical_surface_instance()
        # here sum ideal surface to numerical mesh, and obtain a new numerical mesh:
        # numerical_mesh = add_mesh_to_ideal_surface(numerical_mesh, ideal_surface_ccc)
        x, y = numerical_mesh.get_mesh_x_y()
        X = numpy.outer(x, numpy.ones_like(y))
        Y = numpy.outer(numpy.ones_like(x), y)
        Z = ideal.surface_height(X,Y)
        numerical_mesh.add_to_mesh(Z)

        footprint, normal = numerical_mesh.apply_refraction_on_beam(beam,
                                                         refraction_index_object, refraction_index_image,
                                                         apply_attenuation=apply_attenuation,
                                                         linear_attenuation_coefficient=mu_object)

        return footprint, normal


class S4AdditionalNumericalMeshInterfaceElement(S4InterfaceElement):
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
    instance of S4AdditionalNumericalMeshInterfaceElement
    """
    def __init__(self,
                 optical_element: S4AdditionalNumericalMeshInterface = None,
                 coordinates: ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam: S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4AdditionalNumericalMeshInterface(),
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
        txt += "\nfrom shadow4.beamline.optical_elements.refractors.s4_additional_numerical_mesh_interface import S4AdditionalNumericalMeshInterfaceElement"
        txt += "\nbeamline_element = S4AdditionalNumericalMeshInterfaceElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)"
        txt += "\n\nbeam, footprint = beamline_element.trace_beam()"
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
    from shadow4.beamline.optical_elements.refractors.s4_conic_interface import S4ConicInterface

    optical_element = S4ConicInterface(
        name="Conic Refractive Interface",
        boundary_shape=None,
        material_object="",
        material_image="",
        density_object=1.0,
        density_image=1.0,
        f_r_ind=0,  # source of optical constants, from constant value or PREREFL preprocessor (file):
        #      (0) constant value in both object and image spaces
        #      (1) file in object space, constant value in image space
        #      (2) constant value in object space, file in image space
        #      (3) file in both object and image space
        r_ind_obj=1.0,  # (for f_r_ind=0,2): index of refraction in object space.
        r_ind_ima=1.0,  # (for f_r_ind=0,1): index of refraction in image space.
        r_attenuation_obj=0.0,
        # (for f_r_ind=0,2): attenuation coefficient in object space. Units of UserUnitLength^(-1)
        r_attenuation_ima=0.0,
        # (for f_r_ind=0,1): attenuation coefficient in image space. Units of UserUnitLength^(-1)
        file_r_ind_obj="",  # (for f_r_ind=1,3): file generated by PREREFL
        file_r_ind_ima="",  # (for f_r_ind=2,3): file generated by PREREFL
        dabax=None,
        conic_coefficients=numpy.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
        )

    from syned.beamline.element_coordinates import ElementCoordinates

    coordinates = ElementCoordinates(p=10, q=10, angle_radial=1.54985)


    if use_errors:
        from syned.beamline.shape import Rectangle
        mesh_element = S4AdditionalNumericalMeshInterface(name="M1",
                                                       ideal_interface=optical_element,
                                                       numerical_mesh_interface=S4NumericalMeshInterface(surface_data_file="/home/srio/Oasys/bump.h5",
                                                                                                   boundary_shape=Rectangle(x_left=-0.05,
                                                                                                                            x_right=0.05,
                                                                                                                            y_bottom=-0.5,
                                                                                                                            y_top=0.5)))

        beamline_element = S4AdditionalNumericalMeshInterfaceElement(optical_element=mesh_element,
                                                         coordinates=ElementCoordinates(p=10.0,
                                                                                        q=6.0,
                                                                                        angle_radial=0,
                                                                                        angle_radial_out=numpy.pi),
                                                        input_beam=beam)

    else:
        from shadow4.beamline.optical_elements.refractors.s4_conic_interface import S4ConicInterfaceElement

        beamline_element = S4ConicInterfaceElement(optical_element=optical_element, coordinates=coordinates,
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

    print(beamline_element.to_python_code())

