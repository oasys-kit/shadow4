import numpy

from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.optical_elements.crystals.crystal import Crystal, DiffractionGeometry
from syned.beamline.shape import Rectangle, Ellipse

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.s4_beamline_element import S4BeamlineElement
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements

from crystalpy.diffraction.DiffractionSetupXraylib import DiffractionSetupXraylib
from crystalpy.diffraction.DiffractionSetupDabax import DiffractionSetupDabax
from crystalpy.diffraction.DiffractionSetupShadowPreprocessorV1 import DiffractionSetupShadowPreprocessorV1
from crystalpy.diffraction.DiffractionSetupShadowPreprocessorV2 import DiffractionSetupShadowPreprocessorV2

from crystalpy.diffraction.GeometryType import BraggDiffraction
from crystalpy.diffraction.Diffraction import Diffraction
from crystalpy.util.Vector import Vector
from crystalpy.util.Photon import Photon
from crystalpy.util.ComplexAmplitudePhoton import ComplexAmplitudePhoton
from crystalpy.util.ComplexAmplitudePhotonBunch import ComplexAmplitudePhotonBunch
from crystalpy.diffraction.PerfectCrystalDiffraction import PerfectCrystalDiffraction

from shadow4.optical_surfaces.s4_mesh import S4Mesh
from shadow4.optical_surfaces.s4_toroid import S4Toroid

from shadow4.tools.logger import is_verbose, is_debug

import scipy.constants as codata

class S4Crystal(Crystal):
    """
    Shadow4 Crystal Class
    This is a base class for perfect crystal in reflection geometry (Bragg), using the diffracted beam.

    Use derived classes for plane or other curved crystal surfaces.

    Use other classes for (to be developed):
        * S4TransmissionCrystal : Perfect crystal in transmission (Bragg-transmitted beam, Laue-diffracted and Laue-transmited)
        * S4JohanssonCrystal : Johanssong curved perfect crystals (in Bragg reflection).
        * S4MosaicCrystal : Mosaic crystals (in Bragg reflection).

    Constructor.

    Parameters
    ----------
    name :  str, optional
        A name for the crystal
    boundary_shape : instance of BoundaryShape, optional
        The information on the crystal boundaries.
    surface_shape : instance of SurfaceShape, optional
        The information on crystal surface.
    material : str, optional
        The crystal material name (a name accepted by crystalpy).
    miller_index_h : int, optional
        The Miller index H.
    miller_index_k : int, optional
        The Miller index K.
    miller_index_l : int, optional
        The Miller index L.
    f_bragg_a : int, optional
        Asymmetric crystal 0:No, 1:Yes.
    asymmetry_angle : float, optional
        For f_bragg_a=1, the asymmetry angle (angle between crystal planes and surface) in rads.
    is_thick : int, optional
        Use thick crystal approximation.
    thickness : float, optional
        For is_thick=0, the crystal thickness in m.
    f_central : int, optional
        Flag for autosetting the crystal to the corrected Bragg angle.
    f_phot_cent : int, optional
        0: setting photon energy in eV, 1:setting photon wavelength in A.
    phot_cent : float, optional
        for f_central=1, the value of the photon energy (f_phot_cent=0) or photon wavelength (f_phot_cent=1).
    f_ext : inf, optional
        Flag for autosetting the crystal surface parameters.
        0: internal/calculated parameters, 1:external/user defined parameters. TODO: delete?
    material_constants_library_flag : int, optional
        Flag for indicating the origin of the crystal data:
        0: xraylib, 1: dabax, 2: preprocessor file v1, 3: preprocessor file v2.
    file_refl : str, optional
        for material_constants_library_flag=2,3, the name of the file containing the crystal parameters.

    Returns
    -------
    instance of S4Crystal.
    """
    def __init__(self,
                 name="Undefined",
                 boundary_shape=None,
                 surface_shape=None,
                 material=None,
                 # diffraction_geometry=DiffractionGeometry.BRAGG,
                 miller_index_h=1,
                 miller_index_k=1,
                 miller_index_l=1,
                 f_bragg_a=False,
                 asymmetry_angle=0.0,
                 is_thick=0,          # 1=Use thick crystal approximation
                 thickness=0.010,
                 f_central=False,
                 f_phot_cent=0,
                 phot_cent=8000.0,
                 # f_johansson=False,
                 # r_johansson=1.0,
                 # f_mosaic=False,
                 # spread_mos=0.4*numpy.pi/180,
                 f_ext=0,
                 material_constants_library_flag=0, # 0=xraylib, 1=dabax
                                                    # 2=shadow preprocessor file v1
                                                    # 3=shadow preprocessor file v2
                 file_refl="",
                 ):


        Crystal.__init__(self,
                         name=name,
                         surface_shape=surface_shape,
                         boundary_shape=boundary_shape,
                         material=material,
                         diffraction_geometry=DiffractionGeometry.BRAGG,
                         miller_index_h=miller_index_h,
                         miller_index_k=miller_index_k,
                         miller_index_l=miller_index_l,
                         asymmetry_angle=asymmetry_angle,
                         thickness=thickness,
                        )


        self._f_central = f_central
        self._f_phot_cent = f_phot_cent
        self._phot_cent = phot_cent
        self._f_bragg_a = f_bragg_a
        self._is_thick = is_thick
        self._f_ext = f_ext
        self._material_constants_library_flag = material_constants_library_flag
        self._file_refl = file_refl

        # self._f_mosaic = f_mosaic
        # self._r_johansson = r_johansson
        # self._f_johansson = f_johansson
        # self._spread_mos = spread_mos

        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._add_support_text([
                    ("f_central",           "S4: autotuning",                              ""),
                    ("f_phot_cent",         "S4: for f_central=1: tune to eV(0) or A (1)", ""),
                    ("phot_cent",           "S4: for f_central=1: value in eV or A",       ""),
                    ("f_bragg_a",           "S4: use asymmetruc cut",                      ""),
                    ("is_thick",            "S4: use thick crystal approximation", ""),
                    ("f_ext",               "S4: autosetting curved surface parms.",       ""),
                    ("material_constants_library_flag", "S4: crystal data from: 0=xraylib, 1=dabax, 2=file v1, 3=file v1", ""),
                    ("file_refl",           "S4: preprocessor file name",                  ""),
            ] )


    def get_info(self):
        """
        Returns the specific information of the S4 crystal optical element.

        Returns
        -------
        str
        """
        txt = "\n\n"
        txt += "CRYSTAL\n"
        if self._material_constants_library_flag == 0:
            txt += "Crystal data using xraylib for %s %d%d%d\n" % (self._material,
                                                                   self._miller_index_h,
                                                                   self._miller_index_k,
                                                                   self._miller_index_l)
        elif self._material_constants_library_flag == 1:
            txt += "Crystal data using dabax for %s %d%d%d\n" % (self._material,
                                                                   self._miller_index_h,
                                                                   self._miller_index_k,
                                                                   self._miller_index_l)
        elif self._material_constants_library_flag == 2:
           txt += "Crystal data using preprocessor (bragg V1) file: %s \n" % self._file_refl
        elif self._material_constants_library_flag == 3:
           txt += "Crystal data using preprocessor (bragg V2) file: %s \n" % self._file_refl

        if self._f_central == 0:
            txt += "Using EXTERNAL incidence and reflection angles.\n"
        else:
            txt += "Using INTERNAL or calculated incidence and reflection angles for "
            if self._f_phot_cent == 0:
                txt += "photon energy %.6f eV\n" % self._phot_cent
            else:
                txt += "photon wavelength %f A\n" % (self._phot_cent)


        txt += "\n"
        ss = self.get_surface_shape()
        if ss is None:
            txt += "Surface shape is: Plane (** UNDEFINED?? **)\n"
        else:
            txt += "Surface shape is: %s\n" % ss.__class__.__name__

        #
        if ss is not None: txt += "\nParameters:\n %s\n" % ss.info()

        txt += self.get_optical_surface_instance().info() + "\n"

        boundary = self.get_boundary_shape()
        if boundary is None:
            txt += "Surface boundaries not considered (infinite)"
        else:
            txt += "Surface boundaries are: %s\n" % boundary.__class__.__name__
            txt += "    Limits: " + repr( boundary.get_boundaries()) + "\n"
            txt += boundary.info()

        return txt

    def to_python_code_boundary_shape(self):
        """
        Creates a code block with information of boundary shape.

        Returns
        -------
        str
            The text with the code.
        """
        txt = "" # "\nfrom shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirror"
        bs = self._boundary_shape
        if bs is None:
            txt += "\nboundary_shape = None"
        elif isinstance(bs, Rectangle):
            txt += "\nfrom syned.beamline.shape import Rectangle"
            txt += "\nboundary_shape = Rectangle(x_left=%g, x_right=%g, y_bottom=%g, y_top=%g)" % bs.get_boundaries()
        elif isinstance(bs, Ellipse):
            txt += "\nfrom syned.beamline.shape import Ellipse"
            txt += "\nboundary_shape = Ellipse(a_axis_min=%g, a_axis_max=%g, b_axis_min=%g, b_axis_max=%g)" % bs.get_boundaries()
        return txt

class S4CrystalElement(S4BeamlineElement):
    """
    The base class for Shadow4 crystal element.
    It is made of a S4Crystal and an ElementCoordinates instance. It also includes the input beam.

    Use derived classes for plane or other curved crystal surfaces.

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
    instance of S4CrystalElement.
    """

    def __init__(self,
                 optical_element : S4Crystal = None,
                 coordinates : ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4Crystal(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         movements=movements,
                         input_beam=input_beam)

        self._crystalpy_diffraction_setup = None


    def set_crystalpy_diffraction_setup(self):
        """
        Returns the crystalpy DiffractionSetup.

        Returns
        -------
        instance of crystalpy DiffractionSetupAbstract
        """
        oe = self.get_optical_element()
        coor = self.get_coordinates()

        if oe._material_constants_library_flag == 0:
            if is_verbose(): print("\nCreating a diffraction setup (XRAYLIB) for material:", oe._material)
            diffraction_setup = DiffractionSetupXraylib(geometry_type=BraggDiffraction(),  # todo: use oe._diffraction_geometry
                                                 crystal_name=oe._material,  # string
                                                 thickness=oe._thickness,  # meters
                                                 miller_h=oe._miller_index_h,  # int
                                                 miller_k=oe._miller_index_k,  # int
                                                 miller_l=oe._miller_index_l,  # int
                                                 asymmetry_angle=oe._asymmetry_angle,                            # radians
                                                 azimuthal_angle=0.0)
        elif oe._material_constants_library_flag == 1:
            if is_verbose(): print("\nCreating a diffraction setup (DABAX) for material:", oe._material)
            diffraction_setup = DiffractionSetupDabax(geometry_type=BraggDiffraction(),  # todo: use oe._diffraction_geometry
                                                 crystal_name=oe._material,  # string
                                                 thickness=oe._thickness,  # meters
                                                 miller_h=oe._miller_index_h,  # int
                                                 miller_k=oe._miller_index_k,  # int
                                                 miller_l=oe._miller_index_l,  # int
                                                 asymmetry_angle=oe._asymmetry_angle,  # radians
                                                 azimuthal_angle=0.0)
        elif oe._material_constants_library_flag == 2:
            if is_verbose(): print("\nCreating a diffraction setup (shadow preprocessor file V1)...")
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
            if is_verbose(): print("\nCreating a diffraction setup (shadow preprocessor file V2)...")
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
        """
        Sets the adequate incident and reflection angles to match the tuning energy.
        """
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
            theta_in_grazing  = setting_angle + oe._asymmetry_angle

            if is_verbose():
                print("    align_crystal: dSpacingSI: " , (self._crystalpy_diffraction_setup.dSpacingSI()))
                print("    align_crystal: Bragg angle (uncorrected) for E=%f eV is %f deg" % (energy, numpy.degrees(self._crystalpy_diffraction_setup.angleBragg(energy))))
                print("    align_crystal: Bragg angle (corrected) for E=%f eV is %f deg" % (energy, numpy.degrees(setting_angle)))
                print("    align_crystal: (normal) Incident   angle",  numpy.degrees(numpy.pi/2 - (theta_in_grazing ) ))
                print("    align_crystal: grazing incident angle: ", numpy.degrees(theta_in_grazing ))

                theta_out_grazing = setting_angle - oe._asymmetry_angle # wrong because this just applies the Laue equation
                print("    align_crystal: (normal) Reflection angle [LAUE EQUATION]",  numpy.degrees(numpy.pi/2 - (theta_out_grazing) ))
                print("    align_crystal: grazing output angle [LAUE EQUATION]: ", numpy.degrees(theta_out_grazing))


            KIN = self._crystalpy_diffraction_setup.vectorKscattered(energy=energy)
            theta_out = KIN.angle(self._crystalpy_diffraction_setup.vectorNormalSurface())
            if is_verbose(): print("    align_crystal: (normal) Reflection angle [SCATTERING EQUATION]: ", numpy.degrees(theta_out))
            _,_,angle_azimuthal = coor.get_angles()
            coor.set_angles(angle_radial=    numpy.pi/2 - theta_in_grazing ,
                            angle_radial_out=theta_out,
                            angle_azimuthal=angle_azimuthal)
        else:
            if is_verbose(): print("align_crystal: nothing to align: f_central=0")

        if is_verbose(): print(coor.info())

    def trace_beam(self, **params):
        """
        Runs (ray tracing) the input beam through the element.

        Parameters
        ----------
        **params

        Returns
        -------
        tuple
            (output_beam, footprint) instances of S4Beam.
        """
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

        # mirror movement:
        movements = self.get_movements()
        if movements is not None:
            if movements.f_move:
                input_beam.rot_for(OFFX=movements.offset_x,
                                   OFFY=movements.offset_y,
                                   OFFZ=movements.offset_z,
                                   X_ROT=movements.rotation_x,
                                   Y_ROT=movements.rotation_y,
                                   Z_ROT=movements.rotation_z)

        #
        # reflect beam in the crystal surface and apply crystal reflectivity
        #
        soe = self.get_optical_element()

        if not isinstance(soe, Crystal): raise Exception("Undefined Crystal")

        # if 0:
        #     # two steps (diffraction delegated to optical surface, reflectivity with crystalpy)
        #     footprint, normal = self.apply_crystal_diffraction(input_beam) # warning, beam is also changed!!
        #     beam_in_crystal_frame_before_reflection = input_beam.duplicate()
        #     footprint = self.apply_crystal_reflectivities(beam_in_crystal_frame_before_reflection, footprint)
        # else:
        #     # one steps (diffraction and reflectivity with crystalpy)
        #     if 1:
        #         footprint, normal = self.apply_crystal_diffraction_and_reflectivities(input_beam)  # warning, beam is also changed!!
        #     else:
        #         #
        #         # intercept calculation
        #         #
        #         ccc = soe.get_optical_surface_instance()
        #         footprint = input_beam.duplicate()
        #
        #         x1 = footprint.get_columns([1, 2, 3])  # numpy.array(a3.getshcol([1,2,3]))
        #         v1 = footprint.get_columns([4, 5, 6])  # numpy.array(a3.getshcol([4,5,6]))
        #         flag = footprint.get_column(10)  # numpy.array(a3.getshonecol(10))
        #         optical_path = footprint.get_column(13)
        #
        #         t1, t2 = ccc.calculate_intercept(x1, v1)
        #         reference_distance = -footprint.get_column(2).mean() + footprint.get_column(3).mean()
        #         t, iflag = ccc.choose_solution(t1, t2, reference_distance=reference_distance)
        #
        #         x2 = x1 + v1 * t
        #         for i in range(flag.size):
        #             if iflag[i] < 0: flag[i] = -100
        #
        #         normal = ccc.get_normal(x2)
        #
        #         footprint.set_column(1, x2[0])
        #         footprint.set_column(2, x2[1])
        #         footprint.set_column(3, x2[2])
        #         # footprint.set_column(4, v2[0])
        #         # footprint.set_column(5, v2[1])
        #         # footprint.set_column(6, v2[2])
        #         footprint.set_column(10, flag)
        #         footprint.set_column(13, optical_path + t)
        #
        #         #
        #         # direction and reflectivity calculation using crystalpy
        #         #
        #
        #         # incident  crystalpy photon stack
        #         # we retrieve data from "beam" meaning the beam before reflection, in the crystal frame (incident beam...)
        #         xp = v1[0] # beam_in_crystal_frame_before_reflection.get_column(4)
        #         yp = v1[1] # beam_in_crystal_frame_before_reflection.get_column(5)
        #         zp = v1[2] # beam_in_crystal_frame_before_reflection.get_column(6)
        #         energies = footprint.get_photon_energy_eV()
        #
        #         Esigma = numpy.sqrt(footprint.get_column(24)) * \
        #                  numpy.exp(1j * footprint.get_column(14))
        #         Epi = numpy.sqrt(footprint.get_column(25)) * \
        #               numpy.exp(1j * footprint.get_column(15))
        #
        #         photons_in = ComplexAmplitudePhoton(energies, Vector(xp, yp, zp), Esigma=Esigma, Epi=Epi)
        #
        #
        #         # create crystalpy PerfectCrystalDiffraction instance
        #
        #         # photons_out = Diffraction.calculateDiffractedComplexAmplitudePhoton(self._crystalpy_diffraction_setup,
        #         #                                                                     photons_in) # this is only for flat crystals
        #
        #         #
        #         surface_normal = Vector(normal[0], normal[1], normal[2]).scalarMultiplication(-1.0) # normal is inwards!
        #
        #         # calculate vector H
        #         # Geometrical convention from M.Sanchez del Rio et al., J.Appl.Cryst.(2015). 48, 477-491.
        #
        #         g_modulus = 2.0 * numpy.pi / (self._crystalpy_diffraction_setup.dSpacingSI())
        #         # Let's start from a vector parallel to the surface normal (z axis).
        #         temp_normal_bragg = surface_normal.scalarMultiplication(g_modulus)
        #
        #         # Let's now rotate this vector of an angle alphaX around the y axis (according to the right-hand-rule).
        #         alpha_x = self._crystalpy_diffraction_setup.asymmetryAngle()
        #         axis = self._crystalpy_diffraction_setup.vectorParallelSurface().crossProduct(surface_normal)  # should be ~(1, 0, 0)
        #         temp_normal_bragg = temp_normal_bragg.rotateAroundAxis(axis, -alpha_x)
        #
        #         # Let's now rotate this vector of an angle phi around the z axis (following the ISO standard 80000-2:2009).
        #         phi = self._crystalpy_diffraction_setup.azimuthalAngle()
        #         bragg_normal = temp_normal_bragg.rotateAroundAxis(temp_normal_bragg, phi)
        #
        #         perfect_crystal = PerfectCrystalDiffraction.initializeFromDiffractionSetupAndEnergy(
        #             self._crystalpy_diffraction_setup,
        #             energies,
        #             geometry_type=None,
        #             bragg_normal=bragg_normal,
        #             surface_normal=surface_normal,
        #             # bragg_angle=None,
        #             # psi_0=None,
        #             # psi_H=None,
        #             # psi_H_bar=None,
        #             thickness=None,
        #             d_spacing=None,
        #         )
        #
        #         # Calculate outgoing Photon.
        #         apply_reflectivity = True #todo set always  True
        #         outgoing_complex_amplitude_photon = perfect_crystal._calculatePhotonOut(photons_in,
        #                 apply_reflectivity=apply_reflectivity,
        #                 calculation_method=1,
        #                 is_thick=0,
        #                 use_transfer_matrix=0
        #                 )
        #
        #         if not apply_reflectivity: #todo delete
        #             coeffs = perfect_crystal.calculateDiffraction(photons_in,
        #                                                           calculation_method=1,
        #                                                           is_thick=0,
        #                                                           use_transfer_matrix=0)
        #             outgoing_complex_amplitude_photon.rescaleEsigma(coeffs["S"])
        #             outgoing_complex_amplitude_photon.rescaleEpi(coeffs["P"])
        #
        #
        #
        #         # copy values from crystalpy photon stack to shadow4 beam
        #
        #         footprint.apply_reflectivities(
        #             numpy.sqrt(outgoing_complex_amplitude_photon.getIntensityS()),
        #             numpy.sqrt(outgoing_complex_amplitude_photon.getIntensityP()))
        #
        #         footprint.add_phases(outgoing_complex_amplitude_photon.getPhaseS(),
        #                            outgoing_complex_amplitude_photon.getPhaseP())
        #
        #         footprint.set_column(4, outgoing_complex_amplitude_photon.unitDirectionVector().components()[0])
        #         footprint.set_column(5, outgoing_complex_amplitude_photon.unitDirectionVector().components()[1])
        #         footprint.set_column(6, outgoing_complex_amplitude_photon.unitDirectionVector().components()[2])

        footprint, normal = self.apply_crystal_diffraction_and_reflectivities(input_beam)  # warning, beam is also changed!!


        if movements is not None:
            if movements.f_move:
                footprint.rot_back(OFFX=movements.offset_x,
                                   OFFY=movements.offset_y,
                                   OFFZ=movements.offset_z,
                                   X_ROT=movements.rotation_x,
                                   Y_ROT=movements.rotation_y,
                                   Z_ROT=movements.rotation_z)

        #
        # apply mirror boundaries
        #
        footprint.apply_boundaries_syned(soe.get_boundary_shape(), flag_lost_value=flag_lost_value)


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

    # def apply_crystal_diffraction(self, beam):
    #
    #     oe = self.get_optical_element()
    #     ccc = oe.get_optical_surface_instance()
    #
    #     dSpacingSI = self._crystalpy_diffraction_setup.dSpacingSI()
    #     alphaX = oe._asymmetry_angle
    #
    #     if oe._diffraction_geometry == DiffractionGeometry.BRAGG and oe._asymmetry_angle == 0.0:
    #         print(">>>>>> Using non-dispersive reflection (BRAGG-SYMMETRIC)")
    #         beam_mirr, normal = ccc.apply_crystal_diffraction_bragg_symmetric_on_beam(beam)
    #     else:
    #         print(">>>>>> Using dispersive reflection (BRAGG or LAUE), dSpacingSI, alphaX[deg]", dSpacingSI, numpy.degrees(oe._asymmetry_angle))
    #         beam_mirr, normal = ccc.apply_crystal_diffraction_dispersive_on_beam(beam, dSpacingSI=dSpacingSI, alphaX=alphaX)
    #
    #     return beam_mirr, normal
    #
    # def apply_crystal_reflectivities(self, beam_in_crystal_frame_before_reflection, footprint):
    #     ########################################################################################
    #     #
    #     # apply crystal reflectivity
    #     #
    #     nrays = footprint.get_number_of_rays()
    #     # energy = 8000.0  # eV
    #
    #     # Create a Diffraction object (the calculator)
    #     diffraction = Diffraction()
    #
    #
    #     scan_type = 2 # 0=scan, 1=loop on rays, 2=bunch of photons, 3=bunch in, c.a. out # TODO: delete 0,1,2
    #     if scan_type == 0: # scan
    #         energy = 8000.0  # eV
    #         # setting_angle = self._crystalpy_diffraction_setup.angleBragg(energy)
    #         setting_angle = self._crystalpy_diffraction_setup.angleBraggCorrected(energy)
    #         theta_in_grazing = setting_angle + self.get_optical_element()._asymmetry_angle
    #
    #         angle_deviation_points = nrays
    #         # initialize arrays for storing outputs
    #         intensityS = numpy.zeros(nrays)
    #         intensityP = numpy.zeros(nrays)
    #
    #         angle_deviation_min = -100e-6  # radians
    #         angle_deviation_max = 100e-6  # radians
    #         angle_step = (angle_deviation_max - angle_deviation_min) / angle_deviation_points
    #         deviations = numpy.zeros(angle_deviation_points)
    #         for ia in range(angle_deviation_points):
    #             deviation = angle_deviation_min + ia * angle_step
    #             angle = deviation + theta_in_grazing
    #
    #             # calculate the components of the unitary vector of the incident photon scan
    #             # Note that diffraction plane is YZ
    #             yy = numpy.cos(angle)
    #             zz = - numpy.abs(numpy.sin(angle))
    #             photon = Photon(energy_in_ev=energy, direction_vector=Vector(0.0, yy, zz))
    #             # if ia < 10: print(ia, 0.0, yy, zz)
    #
    #             # perform the calculation
    #             coeffs = diffraction.calculateDiffractedComplexAmplitudes(self._crystalpy_diffraction_setup, photon)
    #
    #             # store results
    #             deviations[ia] = deviation
    #             intensityS[ia] = numpy.abs(coeffs['S'])**2
    #             intensityP[ia] = numpy.abs(coeffs['P'])**2
    #         from srxraylib.plot.gol import plot
    #         plot(deviations, intensityS)
    #
    #     elif scan_type == 1: # from beam, loop
    #         # initialize arrays for storing outputs
    #         complex_reflectivity_S = numpy.zeros(nrays, dtype=complex)
    #         complex_reflectivity_P = numpy.zeros(nrays, dtype=complex)
    #
    #         # we retrieve data from "beam" meaning the beam before reflection, in the crystal frame (incident beam...)
    #         xp = beam_in_crystal_frame_before_reflection.get_column(4)
    #         yp = beam_in_crystal_frame_before_reflection.get_column(5)
    #         zp = beam_in_crystal_frame_before_reflection.get_column(6)
    #         energies = beam_in_crystal_frame_before_reflection.get_photon_energy_eV()
    #         for ia in range(nrays):
    #             photon = Photon(energy_in_ev=energies[ia], direction_vector=Vector(xp[ia], yp[ia], zp[ia]))
    #             # if ia < 10: print(ia, xp[ia], yp[ia], zp[ia])
    #             # perform the calculation
    #             coeffs = diffraction.calculateDiffractedComplexAmplitudes(self._crystalpy_diffraction_setup, photon)
    #             # store results
    #             complex_reflectivity_S[ia] = coeffs['S']
    #             complex_reflectivity_P[ia] = coeffs['P']
    #
    #         footprint.apply_complex_reflectivities(complex_reflectivity_S, complex_reflectivity_P)
    #
    #     elif scan_type == 2: # from beam, bunch
    #         # we retrieve data from "beam" meaning the beam before reflection, in the crystal frame (incident beam...)
    #         xp = beam_in_crystal_frame_before_reflection.get_column(4)
    #         yp = beam_in_crystal_frame_before_reflection.get_column(5)
    #         zp = beam_in_crystal_frame_before_reflection.get_column(6)
    #         energies = beam_in_crystal_frame_before_reflection.get_photon_energy_eV()
    #
    #         Esigma = numpy.sqrt(beam_in_crystal_frame_before_reflection.get_column(24)) * \
    #             numpy.exp(1j * beam_in_crystal_frame_before_reflection.get_column(14))
    #         Epi = numpy.sqrt(beam_in_crystal_frame_before_reflection.get_column(25)) * \
    #             numpy.exp(1j * beam_in_crystal_frame_before_reflection.get_column(15))
    #
    #
    #         photons_in = ComplexAmplitudePhoton(energies, Vector(xp,yp,zp), Esigma=Esigma, Epi=Epi)
    #         photons_out = diffraction.calculateDiffractedComplexAmplitudePhoton(self._crystalpy_diffraction_setup, photons_in)
    #
    #         footprint.apply_reflectivities(
    #             numpy.sqrt(photons_out.getIntensityS()),
    #             numpy.sqrt(photons_out.getIntensityP()))
    #
    #         footprint.add_phases(photons_out.getPhaseS(),
    #                              photons_out.getPhaseP())
    #
    #
    #     elif scan_type == 3:  # from beam, complex amplitides
    #
    #         # we retrieve data from "beam" meaning the beam before reflection, in the crystal frame (incident beam...)
    #         xp = beam_in_crystal_frame_before_reflection.get_column(4)
    #         yp = beam_in_crystal_frame_before_reflection.get_column(5)
    #         zp = beam_in_crystal_frame_before_reflection.get_column(6)
    #         energies = beam_in_crystal_frame_before_reflection.get_photon_energy_eV()
    #
    #         Esigma = numpy.sqrt(beam_in_crystal_frame_before_reflection.get_column(24)) * \
    #                  numpy.exp(1j * beam_in_crystal_frame_before_reflection.get_column(14))
    #         Epi = numpy.sqrt(beam_in_crystal_frame_before_reflection.get_column(25)) * \
    #               numpy.exp(1j * beam_in_crystal_frame_before_reflection.get_column(15))
    #
    #         photons_in = ComplexAmplitudePhoton(energies, Vector(xp, yp, zp), Esigma=Esigma, Epi=Epi)
    #
    #         # photons_out = diffraction.calculateDiffractedComplexAmplitudePhoton(self._crystalpy_diffraction_setup,
    #         #                                                                     photons_in)
    #
    #         coeff = diffraction.calculateDiffractedComplexAmplitudes(self._crystalpy_diffraction_setup,
    #                                                                             photons_in)
    #
    #         footprint.apply_reflectivities(
    #             numpy.sqrt(numpy.abs(coeff['S'])**2),
    #             numpy.sqrt(numpy.abs(coeff['P'])**2))
    #
    #         footprint.add_phases(
    #             numpy.angle(coeff['S']),
    #             numpy.angle(coeff['P']))
    #
    #         # just for check they are equal (for plane crystals only....)
    #         # print(">>>> shadow4 vz: ", footprint.rays[:,5][0:20])
    #         # print(">>>> crystalpy vz: ", photons_out.unitDirectionVector().components()[2][0:20])
    #
    #     ########################################################################################
    #
    #     return footprint

    def apply_crystal_diffraction_and_reflectivities(self, beam):
        """
        Applies the changes in direction and in reflectivity in a beam due to crystal diffraction.

        Parameters
        ----------
        beam : instance of S4Beam
            the beam (already transformed to the local crystal reference system).

        Returns
        -------
        tuple
            (footprint, normal), with footprint an instance of S4beam and the vector normal as numpy arrat (3,:)

        """
        #
        # intercept calculation
        #
        soe = self.get_optical_element()
        ccc = soe.get_optical_surface_instance()
        footprint = beam.duplicate()

        x1 = footprint.get_columns([1, 2, 3])  # numpy.array(a3.getshcol([1,2,3]))
        v1 = footprint.get_columns([4, 5, 6])  # numpy.array(a3.getshcol([4,5,6]))
        flag = footprint.get_column(10)  # numpy.array(a3.getshonecol(10))
        optical_path = footprint.get_column(13)

        # t1, t2 = ccc.calculate_intercept(x1, v1)
        # reference_distance = -footprint.get_column(2).mean() + footprint.get_column(3).mean()
        # t, iflag = ccc.choose_solution(t1, t2, reference_distance=reference_distance)

        reference_distance = -footprint.get_column(2).mean() + footprint.get_column(3).mean()
        t, iflag = ccc.calculate_intercept_and_choose_solution(x1, v1, reference_distance=reference_distance)

        x2 = x1 + v1 * t
        for i in range(flag.size):
            if iflag[i] < 0: flag[i] = -100

        normal = ccc.get_normal(x2)

        footprint.set_column(1, x2[0])
        footprint.set_column(2, x2[1])
        footprint.set_column(3, x2[2])
        # footprint.set_column(4, v2[0])
        # footprint.set_column(5, v2[1])
        # footprint.set_column(6, v2[2])
        footprint.set_column(10, flag)
        footprint.set_column(13, optical_path + t)

        #
        # direction and reflectivity calculation using crystalpy
        #

        energies = footprint.get_photon_energy_eV()

        # incident  crystalpy photon stack
        photons_in = ComplexAmplitudePhoton(
            energies,
            Vector(v1[0], v1[1], v1[2]),
            Esigma=numpy.ones(energies.size, dtype=complex), #numpy.sqrt(footprint.get_column(24)) * numpy.exp(1j * footprint.get_column(14)),
            Epi   =numpy.ones(energies.size, dtype=complex), #numpy.sqrt(footprint.get_column(25)) * numpy.exp(1j * footprint.get_column(15)),
            )

        # create crystalpy PerfectCrystalDiffraction instance
        # Warning:
        # S4Conic, S4Toroid give normal_z < 0 for concave surfaces (and >0 for convex)
        # S4mesh give always normal_z > 0
        # We need for crystalpy the upwards normal
        if isinstance(ccc, S4Mesh):
            surface_normal = Vector(normal[0], normal[1], normal[2])  # normal is outwards!
        elif isinstance(ccc, S4Toroid):
            if ccc.f_torus == 0 or ccc.f_torus == 2:
                surface_normal = Vector(normal[0], normal[1], normal[2]).scalarMultiplication(-1.0) # normal is inwards!
            else:
                surface_normal = Vector(normal[0], normal[1], normal[2])  # normal is outwards!
        else: # todo: check with convex surfaces
            surface_normal = Vector(normal[0], normal[1], normal[2]).scalarMultiplication(-1.0)  # normal is inwards!

        # calculate vector H
        # Geometrical convention from M.Sanchez del Rio et al., J.Appl.Cryst.(2015). 48, 477-491.
        bragg_normal = surface_normal.getVectorH(
            surface_normal,
            self._crystalpy_diffraction_setup.dSpacingSI(),
            asymmetry_angle=self._crystalpy_diffraction_setup.asymmetryAngle(),
            azimuthal_angle=self._crystalpy_diffraction_setup.azimuthalAngle())

        # g_modulus = 2.0 * numpy.pi / (self._crystalpy_diffraction_setup.dSpacingSI())
        # # Let's start from a vector parallel to the surface normal (z axis).
        # temp_normal_bragg = surface_normal.scalarMultiplication(g_modulus)
        #
        # # Let's now rotate this vector of an angle alphaX around the y axis (according to the right-hand-rule).
        # alpha_x = self._crystalpy_diffraction_setup.asymmetryAngle()
        # axis = self._crystalpy_diffraction_setup.vectorParallelSurface().crossProduct(
        #     surface_normal)  # should be ~(1, 0, 0)
        # temp_normal_bragg = temp_normal_bragg.rotateAroundAxis(axis, -alpha_x)
        #
        # # Let's now rotate this vector of an angle phi around the z axis (following the ISO standard 80000-2:2009).
        # phi = self._crystalpy_diffraction_setup.azimuthalAngle()
        # bragg_normal = temp_normal_bragg.rotateAroundAxis(temp_normal_bragg, phi)

        perfect_crystal = PerfectCrystalDiffraction.initializeFromDiffractionSetupAndEnergy(
            self._crystalpy_diffraction_setup,
            energies,
            geometry_type=None,
            bragg_normal=bragg_normal,
            surface_normal=surface_normal,
            # bragg_angle=None,
            # psi_0=None,
            # psi_H=None,
            # psi_H_bar=None,
            thickness=None,
            d_spacing=None,
        )

        # Calculate outgoing Photon.
        apply_reflectivity = 1 # todo: set always  True
        if apply_reflectivity:
            outgoing_complex_amplitude_photon = perfect_crystal.calculatePhotonOut(photons_in,
                                                                                    apply_reflectivity=True,
                                                                                    calculation_method=1,
                                                                                    is_thick=soe._is_thick,
                                                                                    use_transfer_matrix=0
                                                                                    )
        else: # in two steps: todo delete
            outgoing_complex_amplitude_photon = perfect_crystal.calculatePhotonOut(photons_in,
                                                                                    apply_reflectivity=False,
                                                                                    calculation_method=1,
                                                                                    is_thick=soe._is_thick,
                                                                                    use_transfer_matrix=0
                                                                                    )
            coeffs = perfect_crystal.calculateDiffraction(photons_in,
                                                          calculation_method=1,
                                                          is_thick=1, #soe._is_thick,
                                                          use_transfer_matrix=0)
            print(coeffs["S"])
            outgoing_complex_amplitude_photon.rescaleEsigma(coeffs["S"])
            outgoing_complex_amplitude_photon.rescaleEpi(coeffs["P"])


        # copy/apply values from crystalpy photon stack to shadow4 beam
        footprint.apply_reflectivities(
            numpy.sqrt(outgoing_complex_amplitude_photon.getIntensityS()),
            numpy.sqrt(outgoing_complex_amplitude_photon.getIntensityP()))

        footprint.add_phases(outgoing_complex_amplitude_photon.getPhaseS(),
                             outgoing_complex_amplitude_photon.getPhaseP())

        footprint.set_column(4, outgoing_complex_amplitude_photon.unitDirectionVector().components()[0])
        footprint.set_column(5, outgoing_complex_amplitude_photon.unitDirectionVector().components()[1])
        footprint.set_column(6, outgoing_complex_amplitude_photon.unitDirectionVector().components()[2])

        return footprint, normal

if __name__ == "__main__":
    c = S4Crystal(
            name="Undefined",
            boundary_shape=None,
            surface_shape=None,
            material="Si",
            # diffraction_geometry=DiffractionGeometry.BRAGG, #?? not supposed to be in syned...
            miller_index_h=1,
            miller_index_k=1,
            miller_index_l=1,
            asymmetry_angle=0.0,
            is_thick=0,
            thickness=0.010,
            f_central=False,
            f_phot_cent=0,
            phot_cent=8000.0,
            file_refl="",
            f_bragg_a=False,
            # f_johansson=False,
            # r_johansson=1.0,
            # f_mosaic=False,
            # spread_mos=0.4*numpy.pi/180,
            f_ext=0,)

    print(c.info())


    ce = S4CrystalElement(optical_element=c)
    print(ce.info())

    ce = S4CrystalElement()
    ce.set_optical_element(c)

    print(ce.info())
