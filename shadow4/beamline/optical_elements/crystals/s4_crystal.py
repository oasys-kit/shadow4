import numpy

from syned.beamline.shape import Plane
from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.optical_elements.crystals.crystal import Crystal, DiffractionGeometry

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.s4_beamline_element import S4BeamlineElement

from crystalpy.diffraction.DiffractionSetup import DiffractionSetup
from crystalpy.diffraction.DiffractionSetupDabax import DiffractionSetupDabax
from crystalpy.diffraction.DiffractionSetupShadowPreprocessorV1 import DiffractionSetupShadowPreprocessorV1
from crystalpy.diffraction.DiffractionSetupShadowPreprocessorV2 import DiffractionSetupShadowPreprocessorV2

from crystalpy.diffraction.GeometryType import BraggDiffraction
from crystalpy.diffraction.Diffraction import Diffraction
from crystalpy.util.Vector import Vector
from crystalpy.util.Photon import Photon
from crystalpy.util.ComplexAmplitudePhoton import ComplexAmplitidePhoton
from crystalpy.util.ComplexAmplitudePhotonBunch import ComplexAmplitudePhotonBunch

import scipy.constants as codata

class S4Crystal(Crystal):

    def __init__(self,
                 name="Undefined",
                 boundary_shape=None,
                 surface_shape=None,
                 material=None,
                 diffraction_geometry=DiffractionGeometry.BRAGG, #?? not supposed to be in syned...
                 miller_index_h=1,
                 miller_index_k=1,
                 miller_index_l=1,
                 asymmetry_angle=0.0,
                 thickness=0.010, ###########################
                 f_central=False,
                 f_phot_cent=0,
                 phot_cent=8000.0,
                 file_refl="",
                 f_bragg_a=False,
                 # a_bragg=0.0,
                 f_johansson=False,
                 r_johansson=1.0,
                 f_mosaic=False,
                 spread_mos=0.4*numpy.pi/180,
                 f_ext=0,
                 material_constants_library_flag=0, # 0=xraylib, 1=dabax
                                                    # 2=shadow preprocessor file v1
                                                    # 3=shadow preprocessor file v1
                 ):

        """
         f_crystal	=             1 - flag: crystal -- yes (1), no (0).
         f_mosaic	=	  1 - if f_crystal=1; flag: mosaic crystal - yes (1), no (0).

         f_central	 =    1 - flag: autotuning of grating or crystal - yes (1), no (0).
         f_phot_cent =    0 - for f_central=1: tune to eV(0) or Angstroms (1).
         phot_cent	=  11160.0 - for f_phot_cent=1: photon energ
         file_refl	= 'GAAS.SHA         - for f_crystal=1: file containing the crystal parameters.

         f_bragg_a	=     0 - flag: is the crystal asymmetric - yes (1), no (0).
         f_johansson =	  0 - if f_crystal=1; flag: johansson geometry - yes (1), no (0).
         a_bragg =  0.0 - f_bragg_a=1: angle between crystal planes and surface.

         spread_mos	=  0.4 - f_mosaic=1: mosaic spread FWHM (degrees).
         thickness	=  0.1 - crystal thickness in m.
         f_ext	=  0 - flag for internal/calculated (0) parameters vs. external/user defined parameters (1).
         r_johansson =  0.0 - f_ext=1: johansson radius.

        """

        Crystal.__init__(self,
                         name=name,
                         surface_shape=surface_shape,
                         boundary_shape=boundary_shape,
                         material=material,
                         diffraction_geometry=diffraction_geometry,
                         miller_index_h=miller_index_h,
                         miller_index_k=miller_index_k,
                         miller_index_l=miller_index_l,
                         asymmetry_angle=asymmetry_angle,
                         thickness=thickness,
                        )

        self._f_mosaic = f_mosaic
        self._f_central = f_central
        self._f_phot_cent = f_phot_cent
        self._phot_cent = phot_cent
        self._file_refl = file_refl
        self._f_bragg_a = f_bragg_a
        self._f_johansson = f_johansson
        self._spread_mos = spread_mos
        self._f_ext = f_ext
        self._r_johansson = r_johansson
        self._material_constants_library_flag = material_constants_library_flag

        self.congruence()

    def congruence(self):
        print(self._material)
        if self._f_mosaic or \
            self._f_bragg_a or \
            self._f_johansson:
            raise Exception("Not implemented")

class S4CrystalElement(S4BeamlineElement):
    def __init__(self,
                 optical_element : S4Crystal = None,
                 coordinates : ElementCoordinates = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4Crystal(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         input_beam=input_beam)

        self._crystalpy_diffraction_setup = None


    def set_crystalpy_diffraction_setup(self):
        oe = self.get_optical_element()
        coor = self.get_coordinates()

        if oe._material_constants_library_flag == 0:
            print("\nCreating a diffraction setup (XRAYLIB) for material:", oe._material)
            diffraction_setup = DiffractionSetup(geometry_type=BraggDiffraction(),  # todo: use oe._diffraction_geometry
                                                 crystal_name=oe._material,  # string
                                                 thickness=oe._thickness,  # meters
                                                 miller_h=oe._miller_index_h,  # int
                                                 miller_k=oe._miller_index_k,  # int
                                                 miller_l=oe._miller_index_l,  # int
                                                 asymmetry_angle=oe._asymmetry_angle,                            # radians
                                                 azimuthal_angle=0.0)
        elif oe._material_constants_library_flag == 1:
            print("\nCreating a diffraction setup (DABAX) for material:", oe._material)
            diffraction_setup = DiffractionSetupDabax(geometry_type=BraggDiffraction(),  # todo: use oe._diffraction_geometry
                                                 crystal_name=oe._material,  # string
                                                 thickness=oe._thickness,  # meters
                                                 miller_h=oe._miller_index_h,  # int
                                                 miller_k=oe._miller_index_k,  # int
                                                 miller_l=oe._miller_index_l,  # int
                                                 asymmetry_angle=oe._asymmetry_angle,  # radians
                                                 azimuthal_angle=0.0)
        elif oe._material_constants_library_flag == 2:
            print("\nCreating a diffraction setup (shadow preprocessor file V1)...")
            diffraction_setup = DiffractionSetupShadowPreprocessorV1(geometry_type=BraggDiffraction(),  # todo: use oe._diffraction_geometry
                                                 crystal_name=oe._material,            # string
                                                 thickness=oe._thickness,              # meters
                                                 miller_h=oe._miller_index_h,          # int
                                                 miller_k=oe._miller_index_k,          # int
                                                 miller_l=oe._miller_index_l,          # int
                                                 asymmetry_angle=oe._asymmetry_angle,  # radians
                                                 azimuthal_angle=0.0,
                                                 preprocessor_file=oe._file_refl)
        elif oe._material_constants_library_flag == 3:
            print("\nCreating a diffraction setup (shadow preprocessor file V2)...")
            diffraction_setup = DiffractionSetupShadowPreprocessorV2(geometry_type=BraggDiffraction(),  # todo: use oe._diffraction_geometry
                                                 crystal_name=oe._material,            # string
                                                 thickness=oe._thickness,              # meters
                                                 miller_h=oe._miller_index_h,          # int
                                                 miller_k=oe._miller_index_k,          # int
                                                 miller_l=oe._miller_index_l,          # int
                                                 asymmetry_angle=oe._asymmetry_angle,  # radians
                                                 azimuthal_angle=0.0,
                                                 preprocessor_file=oe._file_refl)
        else:
            raise NotImplementedError

        self._crystalpy_diffraction_setup = diffraction_setup

    def align_crystal(self):
        oe = self.get_optical_element()
        coor = self.get_coordinates()

        if oe is None:
            raise Exception("Undefined optical element")

        if oe._f_central:
            if oe._f_phot_cent == 0:
                energy = oe._phot_cent
            else:
                energy = codata.h * codata.c / codata.e * 1e2 / (oe._phot_cent * 1e-8)

            setting_angle = self._crystalpy_diffraction_setup.angleBraggCorrected(energy)

            print("Bragg angle for E=%f eV is %f deg" % (energy, setting_angle * 180.0 / numpy.pi))

            coor.set_angles(angle_radial=numpy.pi/2-setting_angle,
                            angle_radial_out=numpy.pi/2-setting_angle,
                            angle_azimuthal=0.0)
        else:
            print("align_crystal: nothing to align: f_central=0")

        print(coor.info())

    def trace_beam(self, **params):
        flag_lost_value = params.get("flag_lost_value", -1)

        if self._crystalpy_diffraction_setup is None:  # todo: supress if?
            self.set_crystalpy_diffraction_setup()
            self.align_crystal()

        p = self.get_coordinates().p()
        q = self.get_coordinates().q()
        theta_grazing1 = numpy.pi / 2 - self.get_coordinates().angle_radial()
        theta_grazing2 = numpy.pi / 2 - self.get_coordinates().angle_radial_out()
        alpha1 = self.get_coordinates().angle_azimuthal()

        #
        input_beam = self.get_input_beam().duplicate()
        #
        # put beam in mirror reference system
        #
        input_beam.rotate(alpha1, axis=2)
        input_beam.rotate(theta_grazing1, axis=1)
        input_beam.translation([0.0, -p * numpy.cos(theta_grazing1), p * numpy.sin(theta_grazing1)])

        #
        # reflect beam in the mirror surface
        #
        soe = self.get_optical_element()

        beam_in_crystal_frame_before_reflection = input_beam.duplicate()
        if not isinstance(soe, Crystal): # undefined
            raise Exception("Undefined Crystal")
        else:
            footprint, normal = self.apply_crystal_diffraction(input_beam) # warning, beam is also changed!!

        #
        # apply mirror boundaries
        #
        footprint.apply_boundaries_syned(soe.get_boundary_shape(), flag_lost_value=flag_lost_value)


        ########################################################################################
        #
        # TODO" apply crystal reflectivity
        #
        nrays = footprint.get_number_of_rays()
        # energy = 8000.0  # eV

        # Create a Diffraction object (the calculator)
        diffraction = Diffraction()


        scan_type = 1 # 0=scan, 1=loop on rays, 2=bunch of photons (not functional)  # TODO: delete 0,2
        if scan_type == 0: # scan
            energy = 8000.0  # eV
            # setting_angle = self._crystalpy_diffraction_setup.angleBragg(energy)
            setting_angle = self._crystalpy_diffraction_setup.angleBraggCorrected(energy)

            angle_deviation_points = nrays
            # initialize arrays for storing outputs
            intensityS = numpy.zeros(nrays)
            intensityP = numpy.zeros(nrays)

            angle_deviation_min = -100e-6  # radians
            angle_deviation_max = 100e-6  # radians
            angle_step = (angle_deviation_max - angle_deviation_min) / angle_deviation_points
            deviations = numpy.zeros(angle_deviation_points)
            for ia in range(angle_deviation_points):
                deviation = angle_deviation_min + ia * angle_step
                angle = deviation + setting_angle

                # calculate the components of the unitary vector of the incident photon scan
                # Note that diffraction plane is YZ
                yy = numpy.cos(angle)
                zz = - numpy.abs(numpy.sin(angle))
                photon = Photon(energy_in_ev=energy, direction_vector=Vector(0.0, yy, zz))
                # if ia < 10: print(ia, 0.0, yy, zz)

                # perform the calculation
                coeffs = diffraction.calculateDiffractedComplexAmplitudes(self._crystalpy_diffraction_setup, photon)

                # store results
                deviations[ia] = deviation
                intensityS[ia] = coeffs['S'].intensity()
                intensityP[ia] = coeffs['P'].intensity()
        elif scan_type == 1: # from beam, loop
            # initialize arrays for storing outputs
            complex_reflectivity_S = numpy.zeros(nrays, dtype=complex)
            complex_reflectivity_P = numpy.zeros(nrays, dtype=complex)

            # we retrieve data from "beam" meaning the beam before reflection, in the crystal frame (incident beam...)
            xp = beam_in_crystal_frame_before_reflection.get_column(4)
            yp = beam_in_crystal_frame_before_reflection.get_column(5)
            zp = beam_in_crystal_frame_before_reflection.get_column(6)
            energies = beam_in_crystal_frame_before_reflection.get_photon_energy_eV()
            for ia in range(nrays):
                photon = Photon(energy_in_ev=energies[ia], direction_vector=Vector(xp[ia], yp[ia], zp[ia]))
                # if ia < 10: print(ia, xp[ia], yp[ia], zp[ia])
                # perform the calculation
                coeffs = diffraction.calculateDiffractedComplexAmplitudes(self._crystalpy_diffraction_setup, photon)
                # store results
                complex_reflectivity_S[ia] = coeffs['S'].complexAmplitude()
                complex_reflectivity_P[ia] = coeffs['P'].complexAmplitude()

            footprint.apply_complex_reflectivities(complex_reflectivity_S, complex_reflectivity_P)
        elif scan_type == 2: # from beam, bunch
            # this is complicated... and not faster...
            # todo: accelerate crystalpy create calculateDiffractedComplexAmplitudes for a PhotonBunch

            # we retrieve data from "beam" meaning the beam before reflection, in the crystal frame (incident beam...)
            xp = beam_in_crystal_frame_before_reflection.get_column(4)
            yp = beam_in_crystal_frame_before_reflection.get_column(5)
            zp = beam_in_crystal_frame_before_reflection.get_column(6)
            energies = beam_in_crystal_frame_before_reflection.get_photon_energy_eV()

            Esigma = numpy.sqrt(beam_in_crystal_frame_before_reflection.get_column(24)) * \
                numpy.exp(1j * beam_in_crystal_frame_before_reflection.get_column(14))
            Epi = numpy.sqrt(beam_in_crystal_frame_before_reflection.get_column(25)) * \
                numpy.exp(1j * beam_in_crystal_frame_before_reflection.get_column(15))

            photons = ComplexAmplitudePhotonBunch()
            for ia in range(nrays):
                photons.addPhoton(
                    ComplexAmplitidePhoton(energy_in_ev=energies[ia],
                                    direction_vector=Vector(xp[ia], yp[ia], zp[ia]),
                                    Esigma= 1.0, # Esigma[ia],
                                    Epi   = 1.0, # [ia],
                                           )
                    )
            bunch_out = diffraction.calculateDiffractedComplexAmplitudePhotonBunch(self._crystalpy_diffraction_setup, photons)
            bunch_out_dict = bunch_out.toDictionary()
            reflectivity_S = numpy.sqrt(numpy.array(bunch_out_dict["intensityS"]))
            reflectivity_P = numpy.sqrt(numpy.array(bunch_out_dict["intensityP"]))

            footprint.apply_reflectivities(reflectivity_S, reflectivity_P)
            footprint.add_phases(numpy.array(bunch_out_dict["intensityS"]),
                                 numpy.array(bunch_out_dict["intensityP"]))



########################################################################################
        #
        # from element reference system to image plane
        #

        output_beam = footprint.duplicate()
        output_beam.change_to_image_reference_system(theta_grazing2, q)

        # plot results
        if False:
            if scan_type == 0:
                pass
            else:
                deviations = output_beam.get_column(6)
                intensityS = output_beam.get_column(24)
                intensityP = output_beam.get_column(25)

            from srxraylib.plot.gol import plot
            plot(1e6 * deviations, intensityS,
                 1e6 * deviations, intensityP,
                 xtitle="deviation angle [urad]",
                 ytitle="Reflectivity",
                 legend=["Sigma-polarization", "Pi-polarization"],
                 linestyle=['',''],
                 marker=['+','.'])

        return output_beam, footprint

    def apply_crystal_diffraction(self, beam):

        oe = self.get_optical_element()
        ssi = oe.get_surface_shape_instance()
        print(">>>>>>surface_shape_instalce: ", ssi)
        ccc = oe.get_optical_surface_instance()

        if isinstance(ssi, Plane):
            if oe._diffraction_geometry == DiffractionGeometry.BRAGG and oe._asymmetry_angle == 0.0:
                beam_mirr, normal = ccc.apply_crystal_diffraction_bragg_symmetric_on_beam(beam)
            else:
                raise Exception(NotImplementedError)
        else:
            raise NotImplementedError

        return beam_mirr, normal



if __name__ == "__main__":
    c = S4Crystal(
            name="Undefined",
            boundary_shape=None,
            surface_shape=None,
            material="Si",
            diffraction_geometry=DiffractionGeometry.BRAGG, #?? not supposed to be in syned...
            miller_index_h=1,
            miller_index_k=1,
            miller_index_l=1,
            asymmetry_angle=0.0,
            thickness=0.010, ###########################
            f_central=False,
            f_phot_cent=0,
            phot_cent=8000.0,
            file_refl="",
            f_bragg_a=False,
            # a_bragg=0.0,
            f_johansson=False,
            r_johansson=1.0,
            f_mosaic=False,
            spread_mos=0.4*numpy.pi/180,
            f_ext=0,)
    # print(c.info())

    # ce = S4CrystalElement(optical_element=c)
    # print(ce.info())

    ce = S4CrystalElement()
    ce.set_optical_element(c)
