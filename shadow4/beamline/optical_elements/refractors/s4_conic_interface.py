import numpy
from syned.beamline.shape import Conic
from syned.beamline.element_coordinates import ElementCoordinates

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.optical_elements.refractors.s4_interface import S4InterfaceElement, S4Interface
from shadow4.optical_surfaces.s4_conic import S4Conic
from shadow4.beamline.s4_optical_element_decorators import S4ConicOpticalElementDecorator
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements

class S4ConicInterface(S4Interface, S4ConicOpticalElementDecorator):
    """
    Constructor.

    Parameters
    ----------
    name : str, optional
        The name of the mirror.
    boundary_shape : instance of BoundaryShape, optional
        The boundary shape of the mirror.
    conic_coefficients : list, ndarray, optional
        The list of the 10 conic coefficients.
    material_object : str, optional
        string with material formula (just a nema, not used, passed to syned).
    material_image : str, optional
        string with material formula (just a nema, not used, passed to syned).
    density_object : float, optional
        density for material_object (used when f_r_ind>3).
    density_image : float, optional
        density for material_image (used when f_r_ind>3).
    f_r_ind : int, optional
        source of optical constants, from constant value or PREREFL preprocessor (file):
            - (0) constant value in both object and image spaces,
            - (1) file in object space, constant value in image space,
            - (2) constant value in object space, file in image space,
            - (3) file in both object and image space.
            - (4) xraylib in object space, constant value in image space,
            - (5) constant value in object space, xraylib in image space,
            - (6) xraylib in both object and image space.
            - (7) dabax in object space, constant value in image space,
            - (8) constant value in object space, dabax in image space,
            - (9) dabax in both object and image space.
    r_ind_obj : float or numpy array
        (for f_r_ind=0,2): index of refraction (real) in object space.
    r_ind_ima : float or numpy array
        (for f_r_ind=0,1): index of refraction (real) in image space.
    r_attenuation_obj : float or numpy array
        (for f_r_ind=0,2): attenuation coefficient in object space. Units of m^(-1)
    r_attenuation_ima : float or numpy array
        (for f_r_ind=0,1): attenuation coefficient in image space. Units of m^(-1)
    file_r_ind_obj : str, optional
        (for f_r_ind=1,3): file generated by PREREFL preprocessor.
    file_r_ind_ima : str, optional
        (for f_r_ind=2,3): file generated by PREREFL preprocessor.
    dabax : None or instance of DabaxXraylib,
        The pointer to the dabax library  (used for f_r_ind > 6).

    Returns
    -------
    instance of S4ConicInterface.
    """
    def __init__(self,
                 name="Conic Refractive Interface",
                 boundary_shape=None,
                 material_object="",
                 material_image="",
                 density_object=1.0,
                 density_image=1.0,
                 f_r_ind=0,             # source of optical constants, from constant value or PREREFL preprocessor (file):
                                        #      (0) constant value in both object and image spaces
                                        #      (1) file in object space, constant value in image space
                                        #      (2) constant value in object space, file in image space
                                        #      (3) file in both object and image space
                 r_ind_obj=1.0,         # (for f_r_ind=0,2): index of refraction in object space.
                 r_ind_ima=1.0,         # (for f_r_ind=0,1): index of refraction in image space.
                 r_attenuation_obj=0.0, # (for f_r_ind=0,2): attenuation coefficient in object space. Units of UserUnitLength^(-1)
                 r_attenuation_ima=0.0, # (for f_r_ind=0,1): attenuation coefficient in image space. Units of UserUnitLength^(-1)
                 file_r_ind_obj="",     # (for f_r_ind=1,3): file generated by PREREFL
                 file_r_ind_ima="",     # (for f_r_ind=2,3): file generated by PREREFL
                 dabax=None,
                 conic_coefficients=numpy.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                 ):

        S4ConicOpticalElementDecorator.__init__(self, conic_coefficients)
        S4Interface.__init__(self,
                             name=name,
                             boundary_shape=boundary_shape,
                             surface_shape=self._curved_surface_shape,
                             material_object=material_object,
                             material_image=material_image,
                             density_object=density_object,
                             density_image =density_image,
                             f_r_ind=f_r_ind,
                             r_ind_obj=r_ind_obj,
                             r_ind_ima=r_ind_ima,
                             r_attenuation_obj=r_attenuation_obj,
                             r_attenuation_ima=r_attenuation_ima,
                             file_r_ind_obj=file_r_ind_obj,
                             file_r_ind_ima=file_r_ind_ima,
                             dabax=dabax,
                          )

        self.__inputs = {
             "name": name,
             "boundary_shape": boundary_shape,
             "material_object": material_object,
             "material_image": material_image,
             "density_object": density_object,
             "density_image":  density_image,
             "f_r_ind": f_r_ind,
             "r_ind_obj": r_ind_obj,
             "r_ind_ima": r_ind_ima,
             "r_attenuation_obj": r_attenuation_obj,
             "r_attenuation_ima": r_attenuation_ima,
             "file_r_ind_obj": file_r_ind_obj,
             "file_r_ind_ima": file_r_ind_ima,
             "conic_coefficients": repr(conic_coefficients),
             "dabax": repr(dabax),
        }

        self.inputs = self.__inputs

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
        txt = self.to_python_code_boundary_shape()
        txt_pre = """

from shadow4.beamline.optical_elements.refractors.s4_conic_interface import S4ConicInterface
optical_element = S4ConicInterface(name='{name:s}', boundary_shape=boundary_shape,
    f_r_ind={f_r_ind:g}, # source of optical constants:
               # (0) cte in both object (O) and image (I) spaces,
               # (1) file in O, cte in I, (2) cte in O, file in I, (3) file in O and I
               # (4) xraylib in O, cte in I, (5) cte O, xraylib in I, (6) xraylib in O and I
               # (7) dabax O, cte in I, (8) cte value in O, dabax in I, (9) dabax in O and I
    material_object='{material_object:s}', material_image='{material_image:s}',
    density_object={density_object:g}, density_image={density_image:g},
    r_ind_obj={r_ind_obj:g}, r_ind_ima={r_ind_ima:g},
    r_attenuation_obj={r_attenuation_obj:g}, r_attenuation_ima={r_attenuation_ima:g},
    file_r_ind_obj='{file_r_ind_obj:s}', file_r_ind_ima='{file_r_ind_ima:s}',
    dabax={dabax},
    conic_coefficients={conic_coefficients:s})
"""
        txt += txt_pre.format(**self.__inputs)
        return txt


class S4ConicInterfaceElement(S4InterfaceElement):
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
                 optical_element : S4ConicInterface = None,
                 coordinates : ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4ConicInterface(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         movements=movements,
                         input_beam=input_beam)
        if not isinstance(self.get_optical_element().get_surface_shape(), Conic):
            raise ValueError("Wrong Optical Element: only Conic shape is accepted")


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
        txt += "\nfrom shadow4.beamline.optical_elements.refractors.s4_conic_interface import S4ConicInterfaceElement"
        txt += "\nbeamline_element = S4ConicInterfaceElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)"
        txt += "\n\nbeam, footprint = beamline_element.trace_beam()"
        return txt

if __name__ == "__main__":

    if False: # using external refraction index
        from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical

        light_source = SourceGeometrical(name='SourceGeometrical', nrays=5000, seed=5676561)
        light_source.set_spatial_type_rectangle(width=0.001000, height=0.001000)
        light_source.set_angular_distribution_flat(hdiv1=0.000000, hdiv2=0.000000, vdiv1=0.000000, vdiv2=0.000000)
        light_source.set_energy_distribution_singleline(5000.000000, unit='eV')
        light_source.set_polarization(polarization_degree=1.000000, phase_diff=0.000000, coherent_beam=0)
        beam = light_source.get_beam()

        # test plot
        from srxraylib.plot.gol import plot_scatter

        # rays = beam.get_rays()
        # plot_scatter(1e6 * rays[:, 0], 1e6 * rays[:, 2], title='(X,Z) in microns')

        i = S4ConicInterface(name="Conic Refractive Interface",
                     boundary_shape=None,
                     material_object="Be",
                     material_image="Be",
                     density_object=1.848,
                     density_image=1.848,
                     f_r_ind=4,              # source of optical constants, from constant value or PREREFL preprocessor (file):
                                             #      (0) constant value in both object and image spaces
                                             #      (1) file in object space, constant value in image space
                                             #      (2) constant value in object space, file in image space
                                             #      (3) file in both object and image space
                     r_ind_obj=1.0,          # (for f_r_ind=0,2): index of refraction in object space.
                     r_ind_ima=1.5,          # (for f_r_ind=0,1): index of refraction in image space.
                     r_attenuation_obj=0.0,  # (for f_r_ind=0,2): attenuation coefficient in object space. Units of UserUnitLength^(-1)
                     r_attenuation_ima=1e-3, # (for f_r_ind=0,1): attenuation coefficient in image space. Units of UserUnitLength^(-1)
                     file_r_ind_obj="",      # (for f_r_ind=1,3): file generated by PREREFL
                     file_r_ind_ima="",      # (for f_r_ind=2,3): file generated by PREREFL
                     conic_coefficients=[1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.35, 0.0],)

        ie = S4ConicInterfaceElement(optical_element=i,
                                     coordinates=ElementCoordinates(p=0,q=5,angle_radial=0,angle_radial_out=numpy.pi,angle_azimuthal=0),
                                     input_beam=beam)
        beam, footprint = ie.trace_beam()

        print("Intensity: ", beam.intensity(nolost=1))
        print("Estimated intensity (**buggy in shadow3??**): ", numpy.exp(-0.001 * 5) * 5000)
        # test plot
        if True:
            from srxraylib.plot.gol import plot_scatter

            # plot_scatter(beam.get_photon_energy_eV(nolost=1), beam.get_column(23, nolost=1),
            #              title='(Intensity,Photon Energy)', plot_histograms=0)
            plot_scatter(1e6 * beam.get_column(1, nolost=1), 1e6 * beam.get_column(3, nolost=1), title='(X,Z) in microns')
            plot_scatter(1e6 * beam.get_column(4, nolost=1), 1e6 * beam.get_column(6, nolost=1), title="(X',Z') in microrads")


    if True:  # using prerefl in image space

        from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical

        light_source = SourceGeometrical(name='SourceGeometrical', nrays=5000, seed=5676561)
        light_source.set_spatial_type_rectangle(width=0.001000, height=0.001000)
        light_source.set_angular_distribution_flat(hdiv1=0.000000, hdiv2=0.000000, vdiv1=0.000000, vdiv2=0.000000)
        light_source.set_energy_distribution_singleline(5000.000000, unit='eV')
        light_source.set_polarization(polarization_degree=1.000000, phase_diff=0.000000, coherent_beam=0)
        beam = light_source.get_beam()

        # test plot
        from srxraylib.plot.gol import plot_scatter

        # rays = beam.get_rays()
        # plot_scatter(1e6 * rays[:, 0], 1e6 * rays[:, 2], title='(X,Z) in microns')

        i = S4ConicInterface(name="Conic Refractive Interface",
                     boundary_shape=None,
                     material_object="Be",
                     material_image="Be",
                     density_object=1.848,
                     density_image=1.848,
                     f_r_ind=2,              # source of optical constants, from constant value or PREREFL preprocessor (file):
                                             #      (0) constant value in both object and image spaces
                                             #      (1) file in object space, constant value in image space
                                             #      (2) constant value in object space, file in image space
                                             #      (3) file in both object and image space
                     r_ind_obj=1.0,          # (for f_r_ind=0,2): index of refraction in object space.
                     r_ind_ima=1.0,          # (for f_r_ind=0,1): index of refraction in image space.
                     r_attenuation_obj=0.0,  # (for f_r_ind=0,2): attenuation coefficient in object space. Units of UserUnitLength^(-1)
                     r_attenuation_ima=0.0,  # (for f_r_ind=0,1): attenuation coefficient in image space. Units of UserUnitLength^(-1)
                     file_r_ind_obj="",      # (for f_r_ind=1,3): file generated by PREREFL
                     file_r_ind_ima="/nobackup/gurb1/srio/Oasys/Be.dat",    # (for f_r_ind=2,3): file generated by PREREFL
                     conic_coefficients=[1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.35, 0.0],)

        ie = S4ConicInterfaceElement(optical_element=i,
                                     coordinates=ElementCoordinates(p=0,q=5e-3,angle_radial=0,angle_radial_out=numpy.pi,angle_azimuthal=0),
                                     input_beam=beam)
        beam, footprint = ie.trace_beam()

        print("Intensity: ", beam.intensity(nolost=1))
        import xraylib
        mu = xraylib.CS_Total_CP('Be', 5.000000)  * 1.848
        mu *= 100 # m^-1
        print("Estimated intensity: ", numpy.exp(-mu * 5e-3) * 5000)

        # test plot
        if False:
            from srxraylib.plot.gol import plot_scatter

            # plot_scatter(beam.get_photon_energy_eV(nolost=1), beam.get_column(23, nolost=1),
            #              title='(Intensity,Photon Energy)', plot_histograms=0)
            plot_scatter(1e6 * beam.get_column(1, nolost=1), 1e6 * beam.get_column(3, nolost=1), title='(X,Z) in microns')
            plot_scatter(1e6 * beam.get_column(4, nolost=1), 1e6 * beam.get_column(6, nolost=1), title="(X',Z') in microrads")



        print(ie.to_python_code())