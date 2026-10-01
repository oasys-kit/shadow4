"""
The s4 mosaic crystal base class (optical element and beamline element).
"""
import numpy

from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.optical_elements.crystals.crystal import Crystal, DiffractionGeometry
from syned.beamline.shape import Rectangle, Ellipse

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.s4_beamline_element import S4BeamlineElement
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
from shadow4.tools.arrayofvectors import vector_modulus, vector_dot, vector_cross, vector_norm, vector_rotate_around_axis, vector_reflection
from shadow4.tools.logger import is_verbose, is_debug

from crystalpy.diffraction.DiffractionSetupXraylib import DiffractionSetupXraylib
from crystalpy.diffraction.DiffractionSetupDabax import DiffractionSetupDabax
from crystalpy.diffraction.DiffractionSetupShadowPreprocessorV1 import DiffractionSetupShadowPreprocessorV1
from crystalpy.diffraction.DiffractionSetupShadowPreprocessorV2 import DiffractionSetupShadowPreprocessorV2
from crystalpy.diffraction.GeometryType import BraggDiffraction

from dabax.dabax_xraylib import DabaxXraylib

import scipy.constants as codata


class S4MosaicCrystal(Crystal):
    """
    Shadow4 Mosaic Crystal Class
    This is a base class for mosaic crystal in reflection geometry (Bragg), using the diffracted beam.

    Use derived classes for plane or other curved crystal surfaces.

    Use other classes for (to be developed):
        * S4TransmissionCrystal : Perfect crystal in transmission (Bragg-transmitted beam, Laue-diffracted and Laue-transmited)
        * S4JohanssonCrystal : Johanssong curved mosaic crystals (in Bragg reflection).

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
    thickness : float, optional
        For is_thick=0, the crystal thickness in m.
    f_central : int, optional
        Flag for autosetting the crystal to the corrected Bragg angle.
    f_phot_cent : int, optional
        0: setting photon energy in eV, 1:setting photon wavelength in A.
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
        0: macroscopic model (Sanchez del Rio et al., Rev. Sci. Instrum. 63, 932 (1992)):
        Zachariasen/Sears reflectivity as the ray weight, Bragg exit direction, and a
        sampled penetration point.
        1: Monte Carlo model (module s4_mosaic_crystal_monte_carlo): each ray is traced crystallite by
        crystallite through the crystal, as in the code RTbent (L. Alianelli, 2002).
        Multiple reflections, primary extinction (through the crystallite thickness) and
        the exit point on the crystal surface are then included.
    mc_crystallite_thickness_flag : int, optional
        For calculation_method=1, how the thickness t0 of the crystallites is set:
        0: automatic, t0 = mc_crystallite_factor times the primary extinction depth,
        1: user-defined, t0 = mc_crystallite_thickness.
    mc_crystallite_factor : float, optional
        For calculation_method=1 and mc_crystallite_thickness_flag=0, the ratio of the
        crystallite thickness to the primary extinction depth. The depth is the amplitude
        extinction depth, normal to the surface, for sigma polarization. Values of 0.1-0.3
        reproduce the Zachariasen/Sears reflectivity; larger values include primary
        extinction.
    mc_crystallite_thickness : float, optional
        For calculation_method=1 and mc_crystallite_thickness_flag=1, the crystallite
        thickness t0 (normal to the surface) in m.
    mc_max_energies : int, optional
        For calculation_method=1, the maximum number of photon energies where the crystallite
        rocking curve is computed with crystalpy. If the beam has more distinct energies, the
        curves are computed on a regular grid of mc_max_energies energies between the lowest
        and the highest energy, and the crystallite parameters are interpolated linearly.

    Returns
    -------
    instance of S4MosaicCrystal.
    """
    def __init__(self,
                 name="Undefined",
                 boundary_shape=None,
                 surface_shape=None,
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


        Crystal.__init__(self,
                         name=name,
                         surface_shape=surface_shape,
                         boundary_shape=boundary_shape,
                         material=material,
                         diffraction_geometry=DiffractionGeometry.BRAGG,
                         miller_index_h=miller_index_h,
                         miller_index_k=miller_index_k,
                         miller_index_l=miller_index_l,
                         thickness=thickness,
                        )


        self._f_central = f_central
        self._f_phot_cent = f_phot_cent
        self._phot_cent = phot_cent
        self._material_constants_library_flag = material_constants_library_flag
        self._file_refl = file_refl
        self._calculation_method = calculation_method
        self._mc_crystallite_thickness_flag = mc_crystallite_thickness_flag
        self._mc_crystallite_factor = mc_crystallite_factor
        self._mc_crystallite_thickness = mc_crystallite_thickness
        self._mc_max_energies = mc_max_energies

        self._dabax = dabax

        if mosaicity_profile_flag != 0:
            raise NotImplementedError("Only the Gaussian mosaicity profile (mosaicity_profile_flag=0) is implemented.")
        if calculation_method not in (0, 1):
            raise ValueError("calculation_method must be 0 (macroscopic) or 1 (Monte Carlo).")
        if calculation_method == 1:
            if mc_crystallite_thickness_flag not in (0, 1):
                raise ValueError("mc_crystallite_thickness_flag must be 0 (automatic) or 1 (user-defined).")
            if mc_crystallite_thickness_flag == 0 and not (numpy.isfinite(mc_crystallite_factor) and mc_crystallite_factor > 0):
                raise ValueError("mc_crystallite_factor must be finite and positive.")
            if mc_crystallite_thickness_flag == 1 and not (numpy.isfinite(mc_crystallite_thickness) and mc_crystallite_thickness > 0):
                raise ValueError("mc_crystallite_thickness must be finite and positive (in m).")
            if int(mc_max_energies) != mc_max_energies or mc_max_energies < 1:
                raise ValueError("mc_max_energies must be an integer >= 1.")

        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._mosaicity_fwhm_deg = mosaicity_fwhm_deg
        self._mosaicity_profile_flag = mosaicity_profile_flag

        self._add_support_text([
                    ("f_central",           "S4: autotuning",                              ""),
                    ("f_phot_cent",         "S4: for f_central=1: tune to eV(0) or A (1)", ""),
                    ("phot_cent",           "S4: for f_central=1: value in eV or A",       ""),
                    ("material_constants_library_flag", "S4: crystal data from: 0=xraylib, 1=dabax, 2=file v1, 3=file v1", ""),
                    ("file_refl",           "S4: preprocessor file name",                  ""),
                    ("mosaicity_fwhm_deg", "Mosaicity fwhm", "deg"),
                    ("mosaicity_profile_flag", "Mosaic distribution profile 0=Gaussian, 1=External", ""),
                    ("calculation_method",     "Calculation method: 0=macroscopic, 1=Monte Carlo", ""),
                    ("mc_crystallite_thickness_flag", "MC Crystallite thickness 0=automatic, 1=user-defined", ""),
                    ("mc_crystallite_factor",  "MC Crystallite ratio (for automatic)", ""),
                    ("mc_crystallite_thickness", "MC Crystallite thickness (for user-defined)", "m"),
                    ("mc_max_energies", "MC Max number of energies for the crystallite curves", ""),
            ] )


    def get_info(self):
        """
        Returns the specific information of the S4 crystal optical element.

        Returns
        -------
        str
        """
        txt = "\n\n"
        txt += "MOSAICCRYSTAL\n"
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

        txt += "\n"
        if self._calculation_method == 0:
            txt += "Calculation method: macroscopic ray-tracing (1992 model)\n"
        elif self._calculation_method == 1:
            txt += "Calculation method: Monte Carlo (crystallite by crystallite)\n"
            if self._mc_crystallite_thickness_flag == 0:
                txt += "    Crystallite thickness equal to %.3f times the primary extinction depth\n" % self._mc_crystallite_factor
            elif self._mc_crystallite_thickness_flag == 1:
                txt += "    Crystallite thickness = %.4f um\n" % (1e6 * self._mc_crystallite_thickness)
            txt += "    Crystallite rocking curves computed at up to %d energies\n" % self._mc_max_energies

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

    def _get_dabax_txt(self):
        if self._material_constants_library_flag == 1:
            if isinstance(self._dabax, DabaxXraylib):
                dabax_txt = 'DabaxXraylib(file_f0="%s", file_f1f2="%s")' % (self._dabax.get_file_f0(), self._dabax.get_file_f1f2())
            else:
                dabax_txt = "DabaxXraylib()"
        else:
            dabax_txt = "None"

        return dabax_txt

class S4MosaicCrystalElement(S4BeamlineElement):
    """
    The base class for Shadow4 crystal element.
    It is made of a S4MosaicCrystal and an ElementCoordinates instance. It also includes the input beam.

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
    instance of S4MosaicCrystalElement.
    """

    def __init__(self,
                 optical_element : S4MosaicCrystal = None,
                 coordinates : ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4MosaicCrystal(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         movements=movements,
                         input_beam=input_beam)

        self._crystalpy_diffraction_setup = None
        self._mc_crystallite_profiles = None   # set by trace_beam with calculation_method=1


    def get_mc_crystallite_profiles(self):
        """
        Returns the crystallite rocking curves used in the last Monte Carlo trace
        (calculation_method=1).

        The curves are computed with crystalpy at the distinct photon energies of the beam, or at
        mc_max_energies energies on a regular grid between its lowest and highest energy.

        Returns
        -------
        list or None
            None if no Monte Carlo trace has been done. Otherwise, a list with one dict per
            computed energy, with keys:

            * energy [eV], theta_B [rad] (Bragg angle), depth [m] (amplitude primary extinction
              depth, sigma polarization, normal to the surface), t0 [m] (crystallite thickness),
              deviation [rad] (numpy array, angle from theta_B of the scan points);
            * "s" and "p": dicts for each polarization, with reflectivity (numpy array, the
              rocking curve of a perfect t0 slab at the deviation points), gaussian (numpy
              array, the Gaussian used by the Monte Carlo, same integral and FWHM), width [rad]
              (FWHM), peak (of the Gaussian), maximum (of the rocking curve), integrated [rad]
              (integrated reflectivity, including the tails beyond the scan) and f_primary
              (integrated / kinematic integrated reflectivity Q t0 / sin(theta_B); 1 means no
              primary extinction).

        Examples
        --------
        >>> beam, footprint = beamline_element.trace_beam()
        >>> for prof in beamline_element.get_mc_crystallite_profiles():
        ...     print(prof["energy"], prof["s"]["width"], prof["s"]["integrated"])
        ...     # plot(1e6 * prof["deviation"], prof["s"]["reflectivity"], 1e6 * prof["deviation"], prof["s"]["gaussian"])
        """
        return self._mc_crystallite_profiles

    def set_crystalpy_diffraction_setup(self):
        """
        Returns the crystalpy DiffractionSetup.

        Returns
        -------
        instance of crystalpy DiffractionSetupAbstract
        """
        oe = self.get_optical_element()

        if oe._material_constants_library_flag == 0:
            if is_verbose(): print("\nCreating a diffraction setup (XRAYLIB) for material:", oe._material)
            diffraction_setup = DiffractionSetupXraylib(geometry_type=BraggDiffraction(),
                                                 crystal_name=oe._material,  # string
                                                 thickness=oe._thickness,  # meters
                                                 miller_h=oe._miller_index_h,  # int
                                                 miller_k=oe._miller_index_k,  # int
                                                 miller_l=oe._miller_index_l,  # int
                                                 asymmetry_angle=0.0,  # radians (mosaic crystals are symmetric)
                                                 azimuthal_angle=0.0)
        elif oe._material_constants_library_flag == 1:
            if is_verbose(): print("\nCreating a diffraction setup (DABAX) for material:", oe._material)
            diffraction_setup = DiffractionSetupDabax(geometry_type=BraggDiffraction(),
                                                 crystal_name=oe._material,  # string
                                                 thickness=oe._thickness,  # meters
                                                 miller_h=oe._miller_index_h,  # int
                                                 miller_k=oe._miller_index_k,  # int
                                                 miller_l=oe._miller_index_l,  # int
                                                 asymmetry_angle=0.0,  # radians (mosaic crystals are symmetric)
                                                 azimuthal_angle=0.0,
                                                 dabax=oe._dabax)
        elif oe._material_constants_library_flag == 2:
            if is_verbose(): print("\nCreating a diffraction setup (shadow preprocessor file V1)...")
            diffraction_setup = DiffractionSetupShadowPreprocessorV1(geometry_type=BraggDiffraction(),
                                                 crystal_name=oe._material,            # string
                                                 thickness=oe._thickness,              # meters
                                                 miller_h=oe._miller_index_h,          # int
                                                 miller_k=oe._miller_index_k,          # int
                                                 miller_l=oe._miller_index_l,          # int
                                                 asymmetry_angle=0.0,  # radians (mosaic crystals are symmetric)
                                                 azimuthal_angle=0.0,
                                                 preprocessor_file=oe._file_refl)
        elif oe._material_constants_library_flag == 3:
            if is_verbose(): print("\nCreating a diffraction setup (shadow preprocessor file V2)...")
            diffraction_setup = DiffractionSetupShadowPreprocessorV2(geometry_type=BraggDiffraction(),
                                                 crystal_name=oe._material,            # string
                                                 thickness=oe._thickness,              # meters
                                                 miller_h=oe._miller_index_h,          # int
                                                 miller_k=oe._miller_index_k,          # int
                                                 miller_l=oe._miller_index_l,          # int
                                                 asymmetry_angle=0.0,  # radians (mosaic crystals are symmetric)
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


            setting_angle = self._crystalpy_diffraction_setup.angleBragg(energy)
            if isinstance(setting_angle, (list, tuple, numpy.ndarray)): setting_angle = setting_angle[0]

            theta_in_grazing  = setting_angle

            if is_verbose():
                print("    align_crystal: dSpacingSI [m]: " , (self._crystalpy_diffraction_setup.dSpacingSI()))
                print("    align_crystal: Bragg angle (uncorrected) for E=%f eV is %f deg" % (energy, numpy.degrees(self._crystalpy_diffraction_setup.angleBragg(energy))))
                print("    align_crystal: angle set at %f deg" % (numpy.degrees(setting_angle)))
                print("    align_crystal: (normal) Incident   angle [deg]",  numpy.degrees(numpy.pi/2 - (theta_in_grazing ) ))
                print("    align_crystal: grazing incident angle [deg]: ", numpy.degrees(theta_in_grazing ))

                theta_out_grazing = setting_angle
                print("    align_crystal: (normal) Reflection angle [LAUE EQUATION] [deg]",  numpy.degrees(numpy.pi/2 - (theta_out_grazing) ))
                print("    align_crystal: grazing output angle [LAUE EQUATION] [deg]: ", numpy.degrees(theta_out_grazing))


            _, _, angle_azimuthal = coor.get_angles()

            coor.set_angles(angle_radial     = numpy.pi/2 - theta_in_grazing,
                            angle_radial_out = numpy.pi/2 - theta_in_grazing,
                            angle_azimuthal  = angle_azimuthal)
        else:
            if is_verbose(): print("align_crystal: nothing to align: f_central=0")

        if is_verbose(): print(coor.info())

    def to_python_code_mc_crystallite_plot(self):
        """
        Creates a (commented) code block that plots the crystallite rocking curves used in the
        Monte Carlo calculation (calculation_method=1), with their Gaussian approximations.
        To be appended to to_python_code() after the call to trace_beam().

        Returns
        -------
        str
            The text with the code.
        """
        return r"""
# Uncomment to plot the crystallite diffraction profiles (one plot per computed energy)
# import numpy
# from srxraylib.plot.gol import plot
# for prof in beamline_element.get_mc_crystallite_profiles():
#     plot(1e6 * prof["deviation"], prof["s"]["reflectivity"],
#          1e6 * prof["deviation"], prof["s"]["gaussian"],
#          1e6 * prof["deviation"], prof["p"]["reflectivity"],
#          1e6 * prof["deviation"], prof["p"]["gaussian"],
#          legend=[r"Crystallite profile S (FWHM=%.1f $\mu$rad, int-ref=%.1f $\mu$rad)" % (1e6 * prof["s"]["width"], 1e6 * prof["s"]["integrated"]),
#                  "Approximated Gaussian profile S",
#                  r"Crystallite profile P (FWHM=%.1f $\mu$rad, int-ref=%.1f $\mu$rad)" % (1e6 * prof["p"]["width"], 1e6 * prof["p"]["integrated"]),
#                  "Approximated Gaussian profile P"],
#          title=r"Crystallite E=%.3f eV; $\theta_B$=%.2f deg; ext-depth(ampl)=%.3f $\mu$m; thickness=%.3f $\mu$m" % (
#              prof["energy"], numpy.degrees(prof["theta_B"]), 1e6 * prof["depth"], 1e6 * prof["t0"]),
#          xtitle=r"$\theta-\theta_B$ [$\mu$rad]", ytitle="Reflectivity", grid=1)
"""

    def trace_beam(self, **params):
        """
        Runs (ray tracing) the input beam through the element.

        Parameters
        ----------
        **params : accepted parameters, in particular:

        flag_lost_value: float
            numeric value to set in the flag column when ray is lost.

        Returns
        -------
        tuple
            (output_beam, footprint) instances of S4Beam.
        """

        if not isinstance(self.get_optical_element(), Crystal): raise Exception("Undefined Crystal")
        flag_lost_value = params.get("flag_lost_value", -1)
        change_reference_system_in = params.get("change_reference_system_in", True)
        change_reference_system_out = params.get("change_reference_system_out", True)

        if is_verbose():
            print(">>>>>> change_reference_system: ", change_reference_system_in, change_reference_system_out)
            if not change_reference_system_in:
                print("change_reference_system_in = False: skipping reference change to o.e.")
            if not change_reference_system_out:
                print("change_reference_system_out = False: skipping reference change from o.e. to image")

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

        soe = self.get_optical_element()

        if is_verbose():
            b_S, b_P = input_beam.get_efield_directions()
            print("\n\n")
            print(">>> input beam e_S, mod e_s", b_S[0], vector_modulus(b_S)[0])
            print(">>> input beam e_P, mod e_P, e_S.e_P: ", b_P[0], vector_modulus(b_P)[0], vector_dot(b_S, b_P)[0])

        #
        # put input_beam in crystal reference system
        #
        if change_reference_system_in:
            input_beam.rotate(alpha1,         axis=2)
            input_beam.rotate(theta_grazing1, axis=1)

            if is_verbose():
                b_S, b_P = input_beam.get_efield_directions()
                print("")
                print(">>> local beam e_S, mod e_s", b_S[0], vector_modulus(b_S)[0])
                print(">>> local beam e_P, mod e_P, e_S.e_P: ", b_P[0], vector_modulus(b_P)[0], vector_dot(b_S, b_P)[0])

            input_beam.translation([0.0, -p * numpy.cos(theta_grazing1), p * numpy.sin(theta_grazing1)])

        # crystal movement (forward):
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
        # crystal diffraction
        #
        if soe._calculation_method == 0:
            footprint, normal = self._apply_crystal_diffraction(input_beam)
        elif soe._calculation_method == 1:
            footprint, normal = self._apply_monte_carlo(input_beam, flag_lost_value=flag_lost_value)
        else:
            raise ValueError("calculation_method must be 0 (macroscopic) or 1 (Monte Carlo).")
        #
        # apply crystal movements (backwards) and boundaries
        #
        if movements is not None:
            if movements.f_move:
                footprint.rot_back(OFFX=movements.offset_x,
                                   OFFY=movements.offset_y,
                                   OFFZ=movements.offset_z,
                                   X_ROT=movements.rotation_x,
                                   Y_ROT=movements.rotation_y,
                                   Z_ROT=movements.rotation_z)

        footprint.apply_boundaries_syned(soe.get_boundary_shape(), flag_lost_value=flag_lost_value)

        #
        # from element reference system to image plane
        #
        output_beam = footprint.duplicate()
        if change_reference_system_out:
            output_beam.change_to_image_reference_system(theta_grazing2, q)

            if is_verbose():
                b_S, b_P = output_beam.get_efield_directions()
                print("")
                print(">>> image e_S, mod e_s", b_S[0], vector_modulus(b_S)[0])
                print(">>> image e_P, mod e_P, e_S.e_P: ", b_P[0], vector_modulus(b_P)[0], vector_dot(b_S, b_P)[0])

        return output_beam, footprint

    def _apply_crystal_diffraction(self, input_beam):
        """
        Applies mosaic crystal diffraction to the input beam.

        Calculates the surface intercepts, mosaic reflectivity, outgoing directions,
        and internal diffraction positions. Updates the beam positions, directions,
        Jones components, and electric field directions.

        Parameters
        ----------
        input_beam : instance of S4Beam
            The incident beam in the optical element reference system.

        Returns
        -------
        tuple
            (footprint, normal), where footprint is the updated S4Beam and normal
            is a numpy array of shape (3, nrays) containing the surface normals
            at the intercept points.
        """

        footprint, normal = self.get_optical_element().get_optical_surface_instance().calculate_intercept_on_beam(input_beam)

        if is_debug():
            print("    >>>>>> intercept: ", footprint.get_columns([1, 2, 3])[:, 0])
            print("    >>>>>> vout: ", footprint.get_columns([4, 5, 6])[:, 0])
            print("    >>>>>> normal: ", normal.shape, normal[:, 0])

        r_SS, r_PP = self._calculate_mosaic_reflectivity(footprint, normal)
        vIn, vOut = self._calculate_mosaic_reflection(footprint, normal)
        rIn, rOut, path_in = self._sample_mosaic_penetration(footprint, normal)

        jv_out_0, jv_out_1, ee_S, ee_P = self._calculate_jones_and_efield_directions(footprint, normal,
                                                                                        vIn, vOut, r_SS, r_PP)
        # update beam array with the new position
        footprint.set_column(1, rOut[:, 0])
        footprint.set_column(2, rOut[:, 1])
        footprint.set_column(3, rOut[:, 2])
        # update beam array with the new direction
        footprint.set_column(4, vOut[:, 0])
        footprint.set_column(5, vOut[:, 1])
        footprint.set_column(6, vOut[:, 2])
        # update beam array with the new electric fields
        footprint.set_jones_components(jv_out_0, jv_out_1, e_S=ee_S, e_P=ee_P)
        # update optical path
        # note that
        # 1) we add here the path going down into the crystal. The (x,y,z) for the outgoing ray
        #    will be inside the crystal, so the outgoing path will be implemented.
        # 2) we are adding the path not the optical path as we are not considering the refraction index with
        #    the crystal. This is an approximation!
        footprint.set_column(13, footprint.get_column(13) + path_in)

        if is_verbose():
            print(">>> Orthogonal footprint: ", footprint.efields_orthogonal(),
                vector_dot(ee_S, ee_P)[0],
                vector_dot(ee_S, vOut)[0],
                vector_dot(ee_P, vOut)[0])

            b_S, b_P = footprint.get_efield_directions()
            print("")
            print(">>> reflected beam e_S, mod e_s", b_S[0], vector_modulus(b_S)[0])
            print(">>> reflected beam e_P, mod e_P, e_S.e_P: ", b_P[0], vector_modulus(b_P)[0], vector_dot(b_S, b_P)[0])


            print(">>> Intensity foot s, beam in s, foot p,  beam in p:",
                    footprint.get_column(24)[0], input_beam.get_column(24)[0],
                    footprint.get_column(25)[0], input_beam.get_column(25)[0],)


        return footprint, normal

    @staticmethod
    def _incident_facing_normal(footprint, normal):
        """Orient surface normals toward the incident beam, independently of shape.

        Multiplying conic coefficients by -1 reverses their normals without
        changing the surface. Hyperboloids and toroid branches also differ in
        this convention. All mosaic calculations require an incident-facing normal.
        """
        v_in = footprint.get_columns([4, 5, 6]).T
        reverse = vector_dot(v_in, normal.T) > 0
        return normal * numpy.where(reverse, -1.0, 1.0)[numpy.newaxis, :]

    def _calculate_mosaic_reflectivity(self, footprint, normal):
        """
        Calculates the amplitude reflectivity for a symmetric mosaic crystal.

        Uses a Gaussian crystallite orientation profile and the finite-thickness
        reflectivity of the 1992 model, with corrected Q coefficients.
        No orientation or penetration sampling is performed. The footprint is
        not modified.

        Parameters
        ----------
        footprint : instance of S4Beam
            The incident beam at the surface intercepts, in the optical element
            reference system.
        normal : numpy array shape (3, nrays)
            The surface normals at the intercept points.

        Returns
        -------
        tuple
            (r_SS, r_PP), two complex numpy arrays of shape (nrays,) containing
            the amplitude factors for the S and P polarizations. Their squared
            moduli give the corresponding intensity reflectivities.

        Raises
        ------
        ValueError
            If the thickness is not finite and non-negative, the Gaussian
            mosaicity FWHM is not finite and positive, or the absorption
            coefficients are not finite and non-negative.
        """
        if self._crystalpy_diffraction_setup is None:
            self.set_crystalpy_diffraction_setup()
            
        setup = self._crystalpy_diffraction_setup
        soe = self.get_optical_element()
        if soe._thickness < 0 or not numpy.isfinite(soe._thickness):
            raise ValueError("Crystal thickness must be finite and non-negative.")

        surface_normal = self._incident_facing_normal(footprint, normal)

        energies = footprint.get_photon_energy_eV()
        v_in = footprint.get_columns([4, 5, 6]).T
        sin_theta = -vector_dot(v_in, surface_normal.T)
        theta = numpy.arcsin(numpy.clip(sin_theta, -1.0, 1.0))
        theta_bragg = setup.angleBragg(energies)
        theta_diff = theta - theta_bragg

        lambda_cm = codata.h * codata.c / (codata.e * energies) * 100
        mu = -2 * numpy.pi / lambda_cm * numpy.imag(setup.psi0(energies))
        Q_s = (numpy.pi**2 * numpy.abs(setup.psiH(energies) * setup.psiH_bar(energies))
               / (lambda_cm * numpy.sin(2 * theta_bragg)))
        Q_p = Q_s * numpy.cos(2 * theta_bragg)**2

        # Gaussian profile is intentionally fixed; its FWHM comes from the element.
        kappa = numpy.radians(soe._mosaicity_fwhm_deg) / numpy.sqrt(8 * numpy.log(2))
        if not numpy.isfinite(kappa) or kappa <= 0:
            raise ValueError("Gaussian mosaicity FWHM must be finite and positive.")
        w = numpy.exp(-0.5 * (theta_diff / kappa)**2) / (kappa * numpy.sqrt(2 * numpy.pi))
        if numpy.any(~numpy.isfinite(mu)) or numpy.any(mu < 0):
            raise ValueError("Absorption coefficients must be finite and non-negative.")

        path_cm = soe._thickness * 100 / numpy.sin(theta_bragg)

        def amplitude(Q):
            eta = w * Q
            # Algebraically equivalent to Eq. (8), without division by mu or
            # tanh(0). Includes the zero-thickness and zero-absorption limits.
            x = path_cm * numpy.sqrt(mu * (mu + 2 * eta))
            tanhc = numpy.ones_like(x)
            numpy.divide(numpy.tanh(x), x, out=tanhc, where=x != 0)
            effective_path = path_cm * tanhc
            reflectivity = eta * effective_path / (1 + (mu + eta) * effective_path)
            reflectivity = numpy.where(sin_theta > 0, reflectivity, 0.0)
            return numpy.sqrt(reflectivity).astype(complex)

        return amplitude(Q_s), amplitude(Q_p)

    def _calculate_mosaic_reflection(self, footprint, normal):
        """
        Samples the outgoing directions for a symmetric mosaic crystal.

        Rotates the surface normal to satisfy the Bragg condition, then samples
        an additional rotation around the incident direction using the Gaussian
        approximation of the 2013 model. Calculates the outgoing directions by
        reflection about the sampled crystallite normals. The footprint is
        not modified.

        Parameters
        ----------
        footprint : instance of S4Beam
            The incident beam at the surface intercepts, in the optical element
            reference system.
        normal : numpy array shape (3, nrays)
            The surface normals at the intercept points.

        Returns
        -------
        tuple
            (vIn, vOut), two numpy arrays of shape (nrays, 3) containing the
            incident and sampled outgoing unit directions.

        Raises
        ------
        ValueError
            If the thickness is not finite and non-negative, the Gaussian
            mosaicity FWHM is not finite and positive, or the beta sampling
            coefficient s1 is not finite and positive.
        """
        if self._crystalpy_diffraction_setup is None:
            self.set_crystalpy_diffraction_setup()
            
        setup = self._crystalpy_diffraction_setup
        soe = self.get_optical_element()
        if soe._thickness < 0 or not numpy.isfinite(soe._thickness):
            raise ValueError("Crystal thickness must be finite and non-negative.")

        surface_normal = self._incident_facing_normal(footprint, normal)

        energies = footprint.get_photon_energy_eV()
        vIn = footprint.get_columns([4, 5, 6]).T

        # 1. Calculate Delta and a1
        sin_theta = -vector_dot(vIn, surface_normal.T)
        theta = numpy.arcsin(numpy.clip(sin_theta, -1.0, 1.0))
        theta_bragg = setup.angleBragg(energies)
        delta = theta_bragg - theta
        a1 = vector_cross(vIn, surface_normal.T)

        # 2. Calculate the reflected direction using the Bragg condition
        n1 = vector_rotate_around_axis(surface_normal.T, a1, delta)

        # 3. Sample beta using the Gaussian approximation of the 2013 model
        cos_alpha = -vector_dot(vIn, surface_normal.T)
        alpha = numpy.arccos(numpy.clip(cos_alpha, -1.0, 1.0))
        theta_D = numpy.pi / 2 - theta_bragg

        kappa = (
            numpy.radians(soe._mosaicity_fwhm_deg)
            / numpy.sqrt(8 * numpy.log(2))
        )
        if not numpy.isfinite(kappa) or kappa <= 0:
            raise ValueError("Gaussian mosaicity FWHM must be finite and positive.")

        # sinc(u/pi) = sin(u)/u, with the correct limit at u = 0
        u = alpha - theta_D
        s1 = numpy.sin(alpha) * numpy.sin(theta_D) / numpy.sinc(u / numpy.pi)

        if numpy.any(~numpy.isfinite(s1)) or numpy.any(s1 <= 0):
            raise ValueError("Gaussian beta sampling requires finite, positive s1.")

        sigma_beta = kappa / numpy.sqrt(s1)
        beta = numpy.random.normal(loc=0.0, scale=sigma_beta)

        # 4. Calculate n2
        n2 = vector_rotate_around_axis(n1, vIn, beta)

        # 5. Calculcate the reflected direction vOut
        vOut = vector_reflection(vIn, n2)     

        return vIn, vOut

    def _sample_mosaic_penetration(self, footprint, normal):
        """
        Samples the internal diffraction positions for a symmetric mosaic crystal.

        Samples the distance along each incident ray from a truncated exponential
        distribution, limited by the crystal thickness. Uses the Gaussian
        orientation profile and the intensity-weighted S/P scattering coefficient.
        The footprint is not modified.

        Parameters
        ----------
        footprint : instance of S4Beam
            The incident beam at the surface intercepts, in the optical element
            reference system.
        normal : numpy array shape (3, nrays)
            The surface normals at the intercept points.

        Returns
        -------
        tuple
            (rin, rout, path_in). rin and rout are numpy arrays of shape
            (nrays, 3) containing the surface entry positions and sampled
            internal diffraction positions in meters, both expressed in the
            optical element reference system. path_in is a numpy array of
            shape (nrays,) with the sampled path length (in meters) travelled
            inside the crystal from rin to rout.

        Raises
        ------
        ValueError
            If the thickness is not finite and non-negative or the Gaussian
            mosaicity FWHM is not finite and positive.
        """
        if self._crystalpy_diffraction_setup is None:
            self.set_crystalpy_diffraction_setup()
            
        setup = self._crystalpy_diffraction_setup
        soe = self.get_optical_element()
        if soe._thickness < 0 or not numpy.isfinite(soe._thickness):
            raise ValueError("Crystal thickness must be finite and non-negative.")

        surface_normal = self._incident_facing_normal(footprint, normal)

        energies = footprint.get_photon_energy_eV()
        vIn = footprint.get_columns([4, 5, 6]).T
        sin_theta = -vector_dot(vIn, surface_normal.T)
        theta = numpy.arcsin(numpy.clip(sin_theta, -1.0, 1.0))
        theta_bragg = setup.angleBragg(energies)
        theta_diff = theta - theta_bragg

        lambda_cm = codata.h * codata.c / (codata.e * energies) * 100
        # Only s polarization for now!
        Q_s = (numpy.pi**2 * numpy.abs(setup.psiH(energies) * setup.psiH_bar(energies))
               / (lambda_cm * numpy.sin(2 * theta_bragg)))
        Q_p = Q_s * numpy.cos(2 * theta_bragg) ** 2

        kappa = numpy.radians(soe._mosaicity_fwhm_deg) / numpy.sqrt(8 * numpy.log(2))
        if not numpy.isfinite(kappa) or kappa <= 0:
            raise ValueError("Gaussian mosaicity FWHM must be finite and positive.")
        w = numpy.exp(-0.5 * (theta_diff / kappa)**2) / (kappa * numpy.sqrt(2 * numpy.pi))


        Itot = footprint.get_column(23)
        Is = footprint.get_column(24)
        Ip = footprint.get_column(25)

        # Itot = Is + Ip exactly but can be zero for degenerate zero-amplitude rays.
        # Guard the division so a single such ray does not raise below.
        eta_cm = w * Q_s # start with sigma polarization
        numpy.divide(w * (Q_s * Is + Q_p * Ip), Itot, out=eta_cm, where=Itot > 0) # if Itot>0 use average instead

        if numpy.any(~numpy.isfinite(eta_cm)) or numpy.any(eta_cm < 0):
            raise ValueError("Scattering coefficients must be finite and non-negative.")

        # Absolute cosine of the incidence angle relative to the surface normal
        cos_incidence = numpy.abs(vector_dot(vIn, surface_normal.T))
        if numpy.any(~numpy.isfinite(cos_incidence)) or numpy.any(cos_incidence == 0):
            raise ValueError("Penetration sampling requires finite, non-zero incidence cosines.")

        # Maximum path length in cm: _thickness is in meters
        L_cm = 100.0 * soe._thickness / cos_incidence
        if numpy.any(~numpy.isfinite(L_cm)) or numpy.any(L_cm < 0):
            raise ValueError("Maximum path lengths must be finite and non-negative.")

        # Uniform random numbers in [0, 1)
        u = numpy.random.random(size=vIn.shape[0])

        # Uniform limit for eta_cm == 0
        s_cm = u * L_cm

        # Truncated exponential for eta_cm > 0
        numpy.divide(
            -numpy.log1p(u * numpy.expm1(-eta_cm * L_cm)),
            eta_cm,
            out=s_cm,
            where=eta_cm > 0,
        )

        # Entry positions in meters, shape (N, 3)
        rin = footprint.get_columns([1, 2, 3]).T

        # Internal diffraction positions in meters
        rout = rin + (s_cm / 100.0)[:, None] * vIn

        return rin, rout, s_cm / 100.0
                

    def _calculate_jones_and_efield_directions(self, footprint, normal, vIn, vOut, r_SS, r_PP):
        """
        Calculates the Jones vector after crystal diffraction. It also returns the directions of the
        S and P polarized components of the electric field.


        Parameters
        ----------
        footprint : instance of S4Beam
            The input beam
        normal : numpy array shape (nrays, 3)
            The normal to the surface at the intercept points.
        vIn :  numpy array shape (nrays, 3)
            The incident directions
        vOut :  numpy array shape (nrays, 3)
            The incident directions
        r_SS : numpy array complex shape (nrays)
            The crystal reflectivity for the S polarization
        r_PP : numpy array complex shape (nrays)
            The crystal reflectivity for the P polarization

        Returns
        -------
        tuple
            (jv_out_0, jv_out_1, ee_S, ee_P) the two components of the Jones vector and the two vectors of
            shape(nrays, 3) with the electric vectors for the S and P polarizations.

        """
        #
        # get versors with the sigma and pi directions for:
        #     e_S, e_P: the incident beam (as it is)
        #     es_S, es_P: the scattering plane spanned by vIn (incident)
        #     ee_S, ee_P: the scattering plane spanned by vOut (incident)
        #
        # Note that vector_norm() is not needed (for the vectors that should be unitary),
        # but renormalizing them improves accuracy in the calculation of c, s
        #
        e_S, e_P = footprint.get_efield_directions()  # these are \hat{u}_{\sigma,\pi} in Eq. 3

        axis = vector_norm(vector_cross(vIn, vOut))

        es_S = axis  # \hat{u}_{\sigma,i} in Eq. 12
        es_P = vector_norm(vector_cross(es_S, vIn))  # \hat{u}_{\pi,i} in Eq. 12

        ee_S = axis  # \hat{u}_{\sigma,f} in Eq. 13
        ee_P = vector_norm(vector_cross(ee_S, vOut))  # \hat{u}_{\pi,f} in Eq. 13

        if is_verbose():
            print(">>>>> e_S, perp vIn: ", e_S[0], vector_dot(e_S, vIn)[0])
            print(">>>>> e_P, perp vIn: ", e_P[0], vector_dot(e_P, vIn)[0])

            print(">>>>> axis, mod, perp vIn: ", axis[0], vector_modulus(axis)[0], vector_dot(axis, vIn)[0])
            print(">>>>> final ee_S, perp vOut: ", ee_S[0], vector_dot(ee_S, vOut)[0])
            print(">>>>> final ee_P, perp vOut: ", ee_P[0], vector_dot(ee_P, vOut)[0])

        #
        # Jones calculus of refletivity
        #

        # Jones matrix (local)
        J00 = r_SS
        J01 = 0
        J10 = 0
        J11 = r_PP

        # rotation matrix R
        R00 = vector_dot(e_S, es_S)
        R01 = vector_dot(e_S, es_P)
        R10 = vector_dot(e_P, es_S)
        R11 = vector_dot(e_P, es_P)

        # J x R(alpha), the Jones matrix to apply to the Jones vector of the incident rays
        Jrotated_00 = J00 * R00 + J01 * R10  # r_SS * c
        Jrotated_01 = J00 * R01 + J01 * R11  # -r_SS * s
        Jrotated_10 = J10 * R00 + J11 * R10  # r_PP * s
        Jrotated_11 = J10 * R01 + J11 * R11  # r_PP * c

        if is_verbose():
            print(">>> R dotprd: ", R00[0], R01[0], R10[0], R11[0])
            print(">>> J: ", Jrotated_00[0], Jrotated_01[0], Jrotated_10[0], Jrotated_11[0])
            print(">>> |J|: ", numpy.abs(Jrotated_00[0]), numpy.abs(Jrotated_01[0]), numpy.abs(Jrotated_10[0]),
                  numpy.abs(Jrotated_11[0]))

        # Jones vector of incident rays
        jv_in_0, jv_in_1 = footprint.get_jones_components()
        # Jones vector or reflected rays
        jv_out_0 = Jrotated_00 * jv_in_0 + Jrotated_01 * jv_in_1
        jv_out_1 = Jrotated_10 * jv_in_0 + Jrotated_11 * jv_in_1

        return jv_out_0, jv_out_1, ee_S, ee_P

    #
    # Monte Carlo routines
    #
    def _apply_monte_carlo(self, input_beam, flag_lost_value=-1):
        """
        Applies mosaic crystal diffraction to the input beam using the Monte Carlo model
        (calculation_method=1).

        Calculates the surface intercepts and traces each ray crystallite by crystallite
        through the crystal (_apply_monte_carlo_to_rays). Then updates the beam positions
        (exit points on the crystal surface), directions, Jones components, electric field
        directions and optical paths. The steps after the Monte Carlo are the same as in
        _apply_crystal_diffraction.

        Parameters
        ----------
        input_beam : instance of S4Beam
            The incident beam in the optical element reference system.
        flag_lost_value : float, optional
            The value set in the flag column for the rays that are not reflected
            (transmitted through the crystal or stopped at the safety cap).

        Returns
        -------
        tuple
            (footprint, normal), where footprint is the updated S4Beam and normal
            is a numpy array of shape (3, nrays) containing the surface normals
            at the intercept points.
        """

        footprint, normal = self.get_optical_element().get_optical_surface_instance().calculate_intercept_on_beam(input_beam)

        if is_debug():
            print("    >>>>>> intercept: ", footprint.get_columns([1, 2, 3])[:, 0])
            print("    >>>>>> vout: ", footprint.get_columns([4, 5, 6])[:, 0])
            print("    >>>>>> normal: ", normal.shape, normal[:, 0])


        # Monte Carlo: exit points, directions, amplitudes and paths inside the crystal
        rIn, rOut, vIn, vOut, r_SS, r_PP, path_in, lost = self._apply_monte_carlo_to_rays(footprint, normal)
        flag = footprint.get_column(10)
        flag[lost] = flag_lost_value
        footprint.set_column(10, flag)

        jv_out_0, jv_out_1, ee_S, ee_P = self._calculate_jones_and_efield_directions(footprint, normal,
                                                                                        vIn, vOut, r_SS, r_PP)
        # update beam array with the new position
        footprint.set_column(1, rOut[:, 0])
        footprint.set_column(2, rOut[:, 1])
        footprint.set_column(3, rOut[:, 2])
        # update beam array with the new direction
        footprint.set_column(4, vOut[:, 0])
        footprint.set_column(5, vOut[:, 1])
        footprint.set_column(6, vOut[:, 2])
        # update beam array with the new electric fields
        footprint.set_jones_components(jv_out_0, jv_out_1, e_S=ee_S, e_P=ee_P)
        # update optical path
        # note that
        # 1) unlike the macroscopic model, rOut is the exit point on the crystal surface, so path_in
        #    is the whole path travelled inside the crystal (all hops between entry and exit).
        # 2) we are adding the path not the optical path as we are not considering the refraction index with
        #    the crystal. This is an approximation!
        footprint.set_column(13, footprint.get_column(13) + path_in)

        if is_verbose():
            print(">>> Orthogonal footprint: ", footprint.efields_orthogonal(),
                vector_dot(ee_S, ee_P)[0],
                vector_dot(ee_S, vOut)[0],
                vector_dot(ee_P, vOut)[0])

            b_S, b_P = footprint.get_efield_directions()
            print("")
            print(">>> reflected beam e_S, mod e_s", b_S[0], vector_modulus(b_S)[0])
            print(">>> reflected beam e_P, mod e_P, e_S.e_P: ", b_P[0], vector_modulus(b_P)[0], vector_dot(b_S, b_P)[0])


            print(">>> Intensity foot s, beam in s, foot p,  beam in p:",
                    footprint.get_column(24)[0], input_beam.get_column(24)[0],
                    footprint.get_column(25)[0], input_beam.get_column(25)[0],)


        return footprint, normal

    def _apply_monte_carlo_to_rays(self, footprint, normal):
        """
        Traces the rays crystallite by crystallite through the mosaic crystal (calculation_method=1).

        The model is that of RTbent (L. Alianelli, 2002), as validated in the project
        monte-carlo-mosaics; the algorithm and its approximations are described in the module
        s4_mosaic_crystal_monte_carlo:

        * The crystal is a slab of thickness self.get_optical_element()._thickness under the
          tangent plane at the entry point of each ray. It is filled with crystallites of
          thickness t0 (mc_crystallite_thickness_flag, mc_crystallite_factor,
          mc_crystallite_thickness), whose normals follow the Gaussian mosaic distribution
          (mosaicity_fwhm_deg) about the surface normal.
        * Each crystallite reflects with the rocking curve of a perfect t0 slab, computed with
          crystalpy for the photon energy of the ray, approximated by a Gaussian of the same
          FWHM and integrated reflectivity. Otherwise the ray crosses it. Reflections take the
          exact Bragg direction.
        * The ray is traced until it leaves the crystal. It is reflected if it has been
          reflected an odd number of times, and transmitted otherwise. Only reflected rays are
          kept; their intensity is multiplied by exp(-mu * path_in).

        Polarization: the crystallite curve differs for S and P. Each ray is traced once, with
        the S (P) curve with probability Is/I (Ip/I), and its S (P) amplitude is divided by the
        square root of that probability, while the other amplitude is set to zero. The expected
        S and P reflected intensities are then exact (incoherent sum, as usual for mosaic
        crystals).

        Parameters
        ----------
        footprint : instance of S4Beam
            The incident beam at the surface intercepts, in the optical element reference system.
        normal : numpy array of shape (3, nrays)
            The surface normals at the intercept points (any orientation).

        Returns
        -------
        tuple
            (rIn, rOut, vIn, vOut, r_SS, r_PP, path_in, lost):
            rIn : numpy array (nrays, 3), the entry points on the surface [m].
            rOut : numpy array (nrays, 3), the exit points on the surface [m].
            vIn : numpy array (nrays, 3), the incident directions.
            vOut : numpy array (nrays, 3), the exit directions (the specular direction on the
            surface for rays that are not reflected).
            r_SS, r_PP : complex numpy arrays (nrays,), the amplitude factors for S and P
            polarizations (0 for rays that are not reflected).
            path_in : numpy array (nrays,), the path length inside the crystal [m].
            lost : boolean numpy array (nrays,), True for good incident rays that are not
            reflected (transmitted, or stopped at the safety cap).
        """
        from shadow4.beamline.optical_elements.mosaic_crystals.s4_mosaic_crystal_monte_carlo import \
            crystallite_parameters, trace_mosaic_monte_carlo

        if self._crystalpy_diffraction_setup is None:
            self.set_crystalpy_diffraction_setup()

        setup = self._crystalpy_diffraction_setup
        soe = self.get_optical_element()
        if soe._thickness <= 0 or not numpy.isfinite(soe._thickness):
            raise ValueError("Crystal thickness must be finite and positive for the Monte Carlo model.")
        kappa = numpy.radians(soe._mosaicity_fwhm_deg) / numpy.sqrt(8 * numpy.log(2))
        if not numpy.isfinite(kappa) or kappa <= 0:
            raise ValueError("Gaussian mosaicity FWHM must be finite and positive.")

        surface_normal = self._incident_facing_normal(footprint, normal).T   # (nrays, 3), toward the beam
        rIn = footprint.get_columns([1, 2, 3]).T
        vIn = footprint.get_columns([4, 5, 6]).T
        nrays = rIn.shape[0]
        good = footprint.get_column(10) > 0

        # defaults for the rays that are not reflected
        rOut = rIn.copy()
        vOut = vector_reflection(vIn, surface_normal)
        r_SS = numpy.zeros(nrays, dtype=complex)
        r_PP = numpy.zeros(nrays, dtype=complex)
        path_in = numpy.zeros(nrays)
        lost = numpy.zeros(nrays, dtype=bool)
        self._mc_crystallite_profiles = None
        if not good.any():
            return rIn, rOut, vIn, vOut, r_SS, r_PP, path_in, lost

        rng = numpy.random.default_rng(numpy.random.randint(0, 2 ** 31 - 1))  # reproducible with numpy.random.seed

        energies = footprint.get_photon_energy_eV()[good]
        theta_B = setup.angleBragg(energies)
        wavelength = codata.h * codata.c / (codata.e * energies)                  # m
        mu = -2 * numpy.pi / wavelength * numpy.imag(setup.psi0(energies))       # m^-1

        cp, self._mc_crystallite_profiles = crystallite_parameters(setup, energies,
                                    soe._mc_crystallite_thickness_flag,
                                    soe._mc_crystallite_factor,
                                    soe._mc_crystallite_thickness,
                                    max_energies=int(soe._mc_max_energies),
                                    return_profiles=True)

        # polarization channel of each ray
        Itot = footprint.get_column(23)[good]
        Is = footprint.get_column(24)[good]
        p_s = numpy.ones_like(Itot)
        numpy.divide(Is, Itot, out=p_s, where=Itot > 0)
        channel_s = rng.random(energies.size) < p_s
        peak = numpy.where(channel_s, cp["peak_s"], cp["peak_p"])
        width = numpy.where(channel_s, cp["width_s"], cp["width_p"])

        res = trace_mosaic_monte_carlo(rIn[good], vIn[good], surface_normal[good],
                                       theta_B, mu, soe._thickness, kappa,
                                       cp["t0"] / numpy.tan(theta_B), peak, width,
                                       flag_direction=1, rng=rng)

        reflected = res["flag"] == 1
        weight = numpy.exp(-mu * res["path"])
        amp_s = numpy.zeros(energies.size)
        amp_p = numpy.zeros(energies.size)
        numpy.divide(weight, p_s, out=amp_s, where=reflected & channel_s)
        numpy.divide(weight, 1 - p_s, out=amp_p, where=reflected & ~channel_s)

        idx = numpy.flatnonzero(good)
        refl_idx = idx[reflected]
        rOut[refl_idx] = res["X_out"][reflected]
        vOut[refl_idx] = res["V_out"][reflected]
        path_in[refl_idx] = res["path"][reflected]
        r_SS[idx] = numpy.sqrt(amp_s)
        r_PP[idx] = numpy.sqrt(amp_p)
        lost[idx[~reflected]] = True

        st = res["stats"]
        if st["violations"] > 0:
            print("Warning: mosaic Monte Carlo thinning bound exceeded %d times (max ratio %.3f)" %
                  (st["violations"], st["max_ratio"]))
        if is_verbose():
            def stat(label, values, weights=None):
                if values.size == 0:
                    print(f"   {label:>36s}: (no rays)")
                    return
                mean = numpy.average(values, weights=weights)
                std = numpy.sqrt(numpy.average((values - mean) ** 2, weights=weights))
                print(f"   {label:>36s}: mean: {mean:10.5g}, stdev: {round(std, 10):10.5g}, "
                      f"min: {values.min():10.5g}, max: {values.max():10.5g}")

            n_good = energies.size
            n_refl = int(reflected.sum())
            w_refl = weight[reflected]
            Ip = footprint.get_column(25)[good]
            I_out = amp_s * Is + amp_p * Ip
            nref_refl = res["n_reflections"][reflected]
            offset = numpy.linalg.norm(res["X_out"][reflected] - rIn[good][reflected], axis=1)
            print("\nMonte Carlo mosaic crystal: tracing")
            print(f"   crystal thickness: {1e3 * soe._thickness:g} mm, mosaicity FWHM: {soe._mosaicity_fwhm_deg:g} deg "
                  f"(rms tilt per component {1e6 * kappa:.5g} urad), exit direction: exact Bragg")
            print(f"   polarization channels: S {int(channel_s.sum())}, P {int((~channel_s).sum())}")
            print(f"   good incident rays: {n_good}, reflected: {n_refl} ({n_refl / n_good:.4f}), "
                  f"transmitted: {int(numpy.sum(res['flag'] == 0))}, "
                  f"stopped at the safety cap: {int(numpy.sum(res['flag'] == -1))}")
            print(f"   reflected intensity / incident intensity (good rays): {I_out.sum() / Itot.sum():.5g}")
            stat("Bragg angle [deg]", numpy.degrees(theta_B))
            stat("absorption mu [cm-1]", 1e-2 * mu)
            stat("tthi = t0/tan(thetaB) [um]", 1e6 * cp["t0"] / numpy.tan(theta_B))
            stat("weight exp(-mu path) (refl.)", w_refl)
            stat("path inside [um] (refl., weighted)", 1e6 * res["path"][reflected], w_refl)
            stat("exit offset [um] (refl., weighted)", 1e6 * offset, w_refl)
            stat("reflections (refl., weighted)", nref_refl.astype(float), w_refl)
            stat("crystallites crossed (all rays)", res["n_hops"].astype(float))
            if n_refl > 0:
                share = [w_refl[nref_refl == k].sum() / w_refl.sum() for k in (1, 3, 5)]
                print(f"   {'share of 1/3/5 reflections':>36s}: {share[0]:.4f} / {share[1]:.4f} / {share[2]:.4f} "
                      f"(of the reflected intensity)")
            print(f"   thinning: {st['iterations']} iterations, {st['candidates']} candidate crystallites, "
                  f"{st['reflections']} reflections, bound violations: {st['violations']} "
                  f"(max ratio {st['max_ratio']:.3f}, must be <= 1)")

        return rIn, rOut, vIn, vOut, r_SS, r_PP, path_in, lost


if __name__ == "__main__":
    c = S4MosaicCrystal(
        name="Undefined",
        boundary_shape=None,
        surface_shape=None,
        material="graphite",
        miller_index_h=1,
        miller_index_k=1,
        miller_index_l=1,
        thickness=0.010,
        f_central=0,
        f_phot_cent=0,
        phot_cent=8000.0,
        material_constants_library_flag=0,  # 0=xraylib, 1=dabax
        # 2=shadow preprocessor file v1
        # 3=shadow preprocessor file v2
        file_refl="",
        dabax=None,
        mosaicity_fwhm_deg=0.4,
        mosaicity_profile_flag=0,  # 0=Gaussian, 1=External
    )

    print(c.info())


    ce = S4MosaicCrystalElement(optical_element=c)
    print(ce.info())
    print(ce.to_python_code_mc_crystallite_plot())

