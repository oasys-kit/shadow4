"""
The s4 plane mosaic crystal (optical element and beamline element).
"""
from syned.beamline.element_coordinates import ElementCoordinates

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.optical_elements.mosaic_crystals.s4_mosaic_crystal import S4MosaicCrystalElement, S4MosaicCrystal
from shadow4.beamline.s4_optical_element_decorators import S4PlaneOpticalElementDecorator
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements


class S4PlaneMosaicCrystal(S4MosaicCrystal, S4PlaneOpticalElementDecorator):
    """
    Shadow4 Plane Crystal Class
    This is a plane mosaic crystal in reflection geometry (Bragg), using the diffracted beam.

    Constructor.

    Parameters
    ----------
    name :  str, optional
        A name for the crystal
    boundary_shape : instance of BoundaryShape, optional
        The information on the crystal boundaries.
    material : str, optional
        The crystal material name (a name accepted by crystalpy).
    miller_index_h : int, optional
        The Miller index H.
    miller_index_k : int, optional
        The Miller index K.
    miller_index_l : int, optional
        The Miller index L.
    thickness : float, optional
        For is_thick=0, the crystal thickness in m.
    f_central : int, optional
        Flag for autosetting the crystal to the corrected Bragg angle.
    f_phot_cent : int, optional
        0: setting photon energy in eV, 1:setting photon wavelength in m.
    phot_cent : float, optional
        for f_central=1, the value of the photon energy (f_phot_cent=0) or photon wavelength (f_phot_cent=1).
    mosaicity_fwhm_deg : float, optional
        The crystal mosaicity (FWHM) in degrees.
    mosaicity_profile_flag : int, optional
        The distribution function of the crystallites:
        0: Gaussian,
        1: External.
    material_constants_library_flag : int, optional
        Flag for indicating the origin of the crystal data:
        0: xraylib, 1: dabax, 2: preprocessor file v1, 3: preprocessor file v2.
    file_refl : str, optional
        for material_constants_library_flag=2,3, the name of the file containing the crystal parameters.
    dabax : None or instance of DabaxXraylib,
        The pointer to the dabax library  (used for material_constants_library_flag=1).
    calculation_method : int, optional
        The model used for the diffraction in the mosaic crystal:
        0: macroscopic model (1992 paper), 1: Monte Carlo model (crystallite by crystallite).
        See S4MosaicCrystal.
    mc_crystallite_thickness_flag : int, optional
        For calculation_method=1, how the crystallite thickness t0 is set:
        0: automatic (mc_crystallite_factor times the primary extinction depth),
        1: user-defined (mc_crystallite_thickness).
    mc_crystallite_factor : float, optional
        For calculation_method=1 and mc_crystallite_thickness_flag=0, the ratio of the
        crystallite thickness to the primary extinction depth (amplitude, normal to the
        surface, sigma polarization).
    mc_crystallite_thickness : float, optional
        For calculation_method=1 and mc_crystallite_thickness_flag=1, the crystallite
        thickness in m.
    mc_max_energies : int, optional
        For calculation_method=1, the maximum number of photon energies where the crystallite
        rocking curve is computed (otherwise interpolated). See S4MosaicCrystal.

    Returns
    -------
    instance of S4PlaneMosaicCrystal.
    """
    def __init__(self,
                 name="Undefined",
                 boundary_shape=None,
                 material="graphite",
                 miller_index_h=1,
                 miller_index_k=1,
                 miller_index_l=1,
                 thickness=0.010,
                 f_central=0,
                 f_phot_cent=0,
                 phot_cent=8000.0,
                 material_constants_library_flag=0, # 0=xraylib, 1=dabax
                                                    # 2=shadow preprocessor file v1
                                                    # 3=shadow preprocessor file v2
                 file_refl="",
                 dabax=None,
                 mosaicity_fwhm_deg=0.4,
                 mosaicity_profile_flag=0,  # 0=Gaussian, 1=External
                 calculation_method=0,            # 0=macroscopic (1992 paper), 1=Monte Carlo
                 mc_crystallite_thickness_flag=0, # for calculation_method=1: 0=automatic, 1=user-defined
                 mc_crystallite_factor=0.3,       # for mc_crystallite_thickness_flag=0: t0 / extinction depth
                 mc_crystallite_thickness=1e-6,   # for mc_crystallite_thickness_flag=1: t0 in m
                 mc_max_energies=21,              # for calculation_method=1: max number of crystallite curves
                 ):
        S4PlaneOpticalElementDecorator.__init__(self)
        S4MosaicCrystal.__init__(self,
                        name=name,
                        surface_shape=self.get_surface_shape_instance(),
                        boundary_shape=boundary_shape,
                        material=material,
                        miller_index_h=miller_index_h,
                        miller_index_k=miller_index_k,
                        miller_index_l=miller_index_l,
                        thickness=thickness,
                        f_central=f_central,
                        f_phot_cent=f_phot_cent,
                        phot_cent=phot_cent,
                        material_constants_library_flag=material_constants_library_flag, # 0=xraylib, 1=dabax
                        file_refl=file_refl,
                        dabax=dabax,
                        mosaicity_fwhm_deg=mosaicity_fwhm_deg,
                        mosaicity_profile_flag=mosaicity_profile_flag,  # 0=Gaussian, 1=External
                        calculation_method=calculation_method,
                        mc_crystallite_thickness_flag=mc_crystallite_thickness_flag,
                        mc_crystallite_factor=mc_crystallite_factor,
                        mc_crystallite_thickness=mc_crystallite_thickness,
                        mc_max_energies=mc_max_energies,
                        )

        self.__inputs = {
            "name": name,
            "boundary_shape": boundary_shape,
            "material": material,
            # "diffraction_geometry": diffraction_geometry,
            "miller_index_h": miller_index_h,
            "miller_index_k": miller_index_k,
            "miller_index_l": miller_index_l,
            "thickness": thickness,
            "f_central": f_central,
            "f_phot_cent": f_phot_cent,
            "phot_cent": phot_cent,
            "file_refl": file_refl,
            "material_constants_library_flag": material_constants_library_flag,
            "dabax": self._get_dabax_txt(),
            "mosaicity_fwhm_deg": mosaicity_fwhm_deg,
            "mosaicity_profile_flag": mosaicity_profile_flag,  # 0=Gaussian, 1=External
            "calculation_method": calculation_method,
            "mc_crystallite_thickness_flag": mc_crystallite_thickness_flag,
            "mc_crystallite_factor": mc_crystallite_factor,
            "mc_crystallite_thickness": mc_crystallite_thickness,
            "mc_max_energies": mc_max_energies,
            }

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

from shadow4.beamline.optical_elements.mosaic_crystals.s4_plane_mosaic_crystal import S4PlaneMosaicCrystal
optical_element = S4PlaneMosaicCrystal(name='{name}',
    boundary_shape=boundary_shape, material='{material}',
    miller_index_h={miller_index_h}, miller_index_k={miller_index_k}, miller_index_l={miller_index_l},
    thickness={thickness},
    f_central={f_central}, f_phot_cent={f_phot_cent}, phot_cent={phot_cent},
    file_refl='{file_refl}',
    material_constants_library_flag={material_constants_library_flag}, # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
    dabax={dabax}, # used when material_constants_library_flag=1,
    mosaicity_fwhm_deg={mosaicity_fwhm_deg},
    mosaicity_profile_flag={mosaicity_profile_flag},  # 0=Gaussian, 1=External
    calculation_method={calculation_method},  # 0=macroscopic (1992 paper), 1=Monte Carlo
    mc_crystallite_thickness_flag={mc_crystallite_thickness_flag},  # for Monte Carlo: 0=automatic, 1=user-defined
    mc_crystallite_factor={mc_crystallite_factor},  # for automatic: crystallite thickness / extinction depth
    mc_crystallite_thickness={mc_crystallite_thickness},  # for user-defined: crystallite thickness in m
    mc_max_energies={mc_max_energies},  # for Monte Carlo: max number of energies for the crystallite curves
    )"""
        txt += txt_pre.format(**self.__inputs)

        return txt

class S4PlaneMosaicCrystalElement(S4MosaicCrystalElement):
    """
    The Shadow4 plane mosaic crystal element.
    It is made of a S4PlaneMosaicCrystal and an ElementCoordinates instance. It also includes the input beam.

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
                 optical_element : S4PlaneMosaicCrystal = None,
                 coordinates : ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4PlaneMosaicCrystal(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         movements=movements,
                         input_beam=input_beam)

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
        txt += "\nfrom shadow4.beamline.optical_elements.mosaic_crystals.s4_plane_mosaic_crystal import S4PlaneMosaicCrystalElement"
        txt += "\nbeamline_element = S4PlaneMosaicCrystalElement(optical_element=optical_element,coordinates=coordinates, movements=movements, input_beam=beam)"
        txt += "\n\nbeam, footprint = beamline_element.trace_beam()"
        if self.get_optical_element()._calculation_method == 1: txt += self.to_python_code_mc_crystallite_plot()
        return txt


if __name__ == "__main__":
    from dabax.dabax_xraylib import DabaxXraylib
    from shadow4.beamline.s4_beamline import S4Beamline

    beamline = S4Beamline()

    #
    #
    #
    from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical

    light_source = SourceGeometrical(name='Geometrical Source', nrays=15000, seed=5676561)
    light_source.set_spatial_type_point()
    light_source.set_depth_distribution_off()
    light_source.set_angular_distribution_gaussian(sigdix=1e-06, sigdiz=1e-06)
    light_source.set_energy_distribution_uniform(value_min=9000, value_max=11000, unit='eV')
    light_source.set_polarization(polarization_degree=1, phase_diff=0, coherent_beam=0)
    beam = light_source.get_beam()

    beamline.set_light_source(light_source)

    # optical element number XX
    boundary_shape = None

    optical_element = S4PlaneMosaicCrystal(name='Generic Crystal',
                                     boundary_shape=boundary_shape, material='Si',
                                     miller_index_h=1, miller_index_k=1, miller_index_l=1,
                                     thickness=0.001,
                                     f_central=1, f_phot_cent=0, phot_cent=10000.0,
                                     file_refl='bragg.dat',
                                     material_constants_library_flag=1,
                                     # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                     # method_efields_management=0,  # 0=new in S4; 1=like in S3
                                     dabax=DabaxXraylib(file_f0="f0_InterTables.dat", file_f1f2="f1f2_Windt.dat"),
                                     # used when material_constants_library_flag=1,
                                     calculation_method=1,
                                     )
    from syned.beamline.element_coordinates import ElementCoordinates

    coordinates = ElementCoordinates(p=30, q=0, angle_radial=1.371743969, angle_azimuthal=0,
                                     angle_radial_out=1.371743969)
    movements = None

    beamline_element = S4PlaneMosaicCrystalElement(optical_element=optical_element, coordinates=coordinates,
                                             movements=movements, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # test plot
    if True:
        from srxraylib.plot.gol import plot_scatter, plot_show

        plot_scatter(beam.get_photon_energy_eV(nolost=1), beam.get_column(23, nolost=1),
                     title='(Intensity,Photon Energy)', plot_histograms=0, show=0)#, yrange=[0,1.1])
        plot_scatter(1e6 * beam.get_column(1, nolost=1), 1e6 * beam.get_column(3, nolost=1), title='(X,Z) in microns',
                     show=0)
        plot_show()


    print(beamline_element.info())

    print(beam.get_intensity(nolost=1))

    print(beamline.to_python_code())
