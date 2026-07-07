import numpy

from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.optical_elements.crystals.crystal import Crystal, DiffractionGeometry
from syned.beamline.shape import Rectangle, Ellipse

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.s4_beamline_element import S4BeamlineElement
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
from shadow4.tools.arrayofvectors import vector_modulus, vector_dot

from crystalpy.diffraction.DiffractionSetupXraylib import DiffractionSetupXraylib
from crystalpy.diffraction.DiffractionSetupDabax import DiffractionSetupDabax
from crystalpy.diffraction.DiffractionSetupShadowPreprocessorV1 import DiffractionSetupShadowPreprocessorV1
from crystalpy.diffraction.DiffractionSetupShadowPreprocessorV2 import DiffractionSetupShadowPreprocessorV2

from crystalpy.diffraction.GeometryType import LaueDiffraction

from dabax.dabax_xraylib import DabaxXraylib

from shadow4.tools.logger import is_verbose
from shadow4.tools.arrayofvectors import vector_modulus, vector_dot, vector_cross, vector_norm
from shadow4.tools.arrayofvectors import vector_multiply_scalar, vector_diff

import scipy.constants as codata

from shadow4.tools.logger import is_verbose, is_debug
from shadow4.tools.arrayofvectors import vector_modulus_square, vector_modulus, vector_norm, vector_rotate_around_axis
from shadow4.tools.logger import is_verbose, is_debug

from crystalpy.diffraction.GeometryType import BraggDiffraction
from crystalpy.util.Vector import Vector
from crystalpy.util.ComplexAmplitudePhoton import ComplexAmplitudePhoton
from crystalpy.diffraction.PerfectCrystalDiffraction import PerfectCrystalDiffraction

from shadow4.optical_surfaces.s4_mesh import S4Mesh
from shadow4.optical_surfaces.s4_toroid import S4Toroid

class S4GeometricLaueCrystal(Crystal):
    """
    Shadow4 Crystal Class
    This is a base class for a perfect crystal in Laue (transmission) diffraction geometry, treated with a
    purely geometric (ray-tracing) formalism.

    Use derived classes for plane or other curved crystal surfaces.

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
    poisson_ratio : float, optional
        The (isotropic) Poisson-type elastic constant rho = nu = -s23/s22 of the crystal cut for the
        given bending axis (free-plate, meridional cylindrical bending -- see Guigay & Sanchez del
        Rio, Appendix B). Used only by the geometric (bent-crystal) diffraction model; for Si,
        typically 0.22-0.28 depending on the cut.
    dabax : None or instance of DabaxXraylib,
        The pointer to the dabax library  (used for material_constants_library_flag=1).

    Returns
    -------
    instance of S4GeometricLaueCrystal.
    """
    def __init__(self,
                 name="Undefined",
                 boundary_shape=None,
                 surface_shape=None,
                 material=None,
                 miller_index_h=1,
                 miller_index_k=1,
                 miller_index_l=1,
                 f_bragg_a=False,
                 asymmetry_angle=0.0,
                 is_thick=0,          # 1=Use thick crystal approximation
                 thickness=0.010,
                 f_central=0,
                 f_phot_cent=0,
                 phot_cent=8000.0,
                 f_ext=0,
                 material_constants_library_flag=0, # 0=xraylib, 1=dabax
                                                    # 2=shadow preprocessor file v1
                                                    # 3=shadow preprocessor file v2
                 file_refl="",
                 poisson_ratio=0.22,
                 dabax=None,
                 ):


        Crystal.__init__(self,
                         name=name,
                         surface_shape=surface_shape,
                         boundary_shape=boundary_shape,
                         material=material,
                         diffraction_geometry=DiffractionGeometry.LAUE,
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
        self._poisson_ratio = poisson_ratio

        self._dabax = dabax

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
                    ("poisson_ratio",       "S4: elastic constant rho=nu for bent-crystal geometric diffraction", ""),
            ] )


    def get_info(self):
        """
        Returns the specific information of the S4 crystal optical element.

        Returns
        -------
        str
        """
        txt = "\n\n"
        txt += "GEOMETRIC LAUE CRYSTAL\n"
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
        txt = ""
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

class S4GeometricLaueCrystalElement(S4BeamlineElement):
    """
    The base class for Shadow4 geometric Laue crystal element.
    It is made of a S4GeometricLaueCrystal and an ElementCoordinates instance. It also includes the input beam.

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
    instance of S4GeometricLaueCrystalElement.
    """

    def __init__(self,
                 optical_element : S4GeometricLaueCrystal = None,
                 coordinates : ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4GeometricLaueCrystal(),
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
            diffraction_setup = DiffractionSetupXraylib(geometry_type=LaueDiffraction(),
                                                 crystal_name=oe._material,  # string
                                                 thickness=oe._thickness,  # meters
                                                 miller_h=oe._miller_index_h,  # int
                                                 miller_k=oe._miller_index_k,  # int
                                                 miller_l=oe._miller_index_l,  # int
                                                 asymmetry_angle=oe._asymmetry_angle,                            # radians
                                                 azimuthal_angle=0.0)
        elif oe._material_constants_library_flag == 1:
            if is_verbose(): print("\nCreating a diffraction setup (DABAX) for material:", oe._material)
            diffraction_setup = DiffractionSetupDabax(geometry_type=LaueDiffraction(),
                                                 crystal_name=oe._material,  # string
                                                 thickness=oe._thickness,  # meters
                                                 miller_h=oe._miller_index_h,  # int
                                                 miller_k=oe._miller_index_k,  # int
                                                 miller_l=oe._miller_index_l,  # int
                                                 asymmetry_angle=oe._asymmetry_angle,  # radians
                                                 azimuthal_angle=0.0,
                                                 dabax=oe._dabax)
        elif oe._material_constants_library_flag == 2:
            if is_verbose(): print("\nCreating a diffraction setup (shadow preprocessor file V1)...")
            diffraction_setup = DiffractionSetupShadowPreprocessorV1(geometry_type=LaueDiffraction(),
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
            diffraction_setup = DiffractionSetupShadowPreprocessorV2(geometry_type=LaueDiffraction(),
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
            if isinstance(setting_angle, (list, tuple, numpy.ndarray)): setting_angle = setting_angle[0]

            theta_in_grazing  = setting_angle + oe._asymmetry_angle

            if is_verbose():
                print("    align_crystal: dSpacingSI [m]: " , (self._crystalpy_diffraction_setup.dSpacingSI()))
                print("    align_crystal: Bragg angle (uncorrected) for E=%f eV is %f deg" % (energy, numpy.degrees(self._crystalpy_diffraction_setup.angleBragg(energy))))
                print("    align_crystal: Bragg angle (corrected) for E=%f eV is %f deg" % (energy, numpy.degrees(setting_angle)))
                print("    align_crystal: (normal) Incident   angle [deg]",  numpy.degrees(numpy.pi/2 - (theta_in_grazing ) ))
                print("    align_crystal: grazing incident angle [deg]: ", numpy.degrees(theta_in_grazing ))

                theta_out_grazing = setting_angle - oe._asymmetry_angle # wrong because this just applies the Laue equation
                print("    align_crystal: (normal) Reflection angle [LAUE EQUATION] [deg]",  numpy.degrees(numpy.pi/2 - (theta_out_grazing) ))
                print("    align_crystal: grazing output angle [LAUE EQUATION] [deg]: ", numpy.degrees(theta_out_grazing))

            KIN = self._crystalpy_diffraction_setup.vectorKscattered(energy=energy)
            theta_out = KIN.angle(self._crystalpy_diffraction_setup.vectorNormalSurface())
            if isinstance(theta_out, (list, tuple, numpy.ndarray)): theta_out = theta_out[0]

            if is_verbose(): print("    align_crystal: (normal) Reflection angle [SCATTERING EQUATION] [deg]: ", numpy.degrees(theta_out))
            _, _, angle_azimuthal = coor.get_angles()

            coor.set_angles(angle_radial     = numpy.pi/2 - theta_in_grazing,
                            angle_radial_out = theta_out,
                            angle_azimuthal  = angle_azimuthal)
        else:
            if is_verbose(): print("align_crystal: nothing to align: f_central=0")

        if is_verbose(): print(coor.info())

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

        footprint, normal = self._apply_crystal_diffraction(input_beam, flag_lost_value=flag_lost_value)

        #
        # for i in range(10):
        #     print(">>>> ray:", i)
        #     print("     footprint: ", footprint.rays.shape, footprint.rays[i, 0:3])
        #     print("     normal: ", normal.shape, normal[:, i])

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

    def _apply_crystal_diffraction(self, input_beam, flag_lost_value=-1):
        """
        Geometric (purely kinematic) ray tracing through a bent Laue crystal, following the
        closed-form algorithm of Section 5 ("Ray-tracing methodology for bent Laue crystals") of
        Guigay & Sanchez del Rio, "Geometrical focusing in Laue crystals".

        For each ray: (1) find the entry point on the (possibly curved) entrance surface; (2) build
        the local, strain-corrected reciprocal-lattice vector h'(t) along the ray path [Eq. (hpt)];
        (3) solve in closed form for the diffraction depth t_d at which the ray meets the exact
        local Bragg condition [Eq. (td)] -- rays for which t_d falls outside [0, T] do not diffract
        anywhere in the crystal and are flagged lost; (4) compute the exit direction from the local
        diffraction vector at t_d [Eq. (ih)]; (5) advance the ray to the diffraction point.

        The reflectivity model itself is only a placeholder (unit reflectivity for both
        polarizations, r_S = r_P = 1, applied through the same Jones-matrix rotation machinery as
        the Bragg-reflection crystal classes) -- the reference paper is purely geometric/kinematic
        and gives no reflectivity/dynamical-theory model; that is left as a later TODO (e.g. a
        Penning-Polder-type model), see _calculate_penetration_and_scattering().

        Coordinate mapping to the local shadow4 beamline-element frame (see trace_beam()): the
        local (zeta, eta) meridional-plane basis used in _calculate_penetration_and_scattering() is
        built from the existing (surface-type-aware) surface_normal -- not hardcoded local axes --
        so it self-consistently matches whichever inward/outward sign convention that block already
        applies for S4Mesh / S4Toroid / S4Conic; the sagittal direction is local X, per the
        S4SphereOpticalElementDecorator / S4Beam.rotate() convention used throughout shadow4, and is
        not affected by the (purely meridional) cylindrical-bending model used here.

        NOTE: the sign of R (bending radius, positive when the incident beam impinges on the
        concave side of the bent crystal, matching the free-plate displacement field of Eq. (displ))
        against syned's Convexity/Sphere.get_radius() convention has not been independently verified
        against a worked numerical example; check against a known magic-condition case (e.g. Qi
        et al., 2019/2021) before relying on the sign of the focusing effects.

        NOTE: align_crystal() is inherited unchanged from the Bragg-reflection polariser skeleton
        and still encodes the Bragg angle relation (phi_0 + phi_h = 2*alpha + pi); the Laue relation
        used throughout this paper is (phi_0 + phi_h = 2*alpha). align_crystal() most likely needs
        its own Laue-specific rewrite; this method does not rely on it beyond using the diffraction
        setup's Bragg-angle machinery (angleBraggCorrected), which is geometry-independent.

        Parameters
        ----------
        input_beam : instance of S4Beam
            The beam (already transformed to the local crystal reference system).
        flag_lost_value : float
            Numeric value to set in the flag column (10) for rays that do not satisfy the local
            Bragg condition anywhere inside the crystal thickness (transmitted, undiffracted, and
            therefore not modeled by this diffracted-beam element).

        Returns
        -------
        tuple
            (footprint, normal), with footprint an instance of S4Beam (direction and position
            updated to the diffracted ray and its diffraction point) and normal the (3, nrays)
            inward unit surface normal at the entry point.
        """


        # geometric and physics for the scattering process:
        # reflect beam in the crystal and apply crystal reflectivity

        footprint, normal = self.get_optical_element().get_optical_surface_instance().calculate_intercept_on_beam(input_beam)

        if is_debug():
            print("    >>>>>> intercept: ", footprint.get_columns([1, 2, 3])[:, 0])
            print("    >>>>>> vout: ", footprint.get_columns([4, 5, 6])[:, 0])
            print("    >>>>>> normal: ", normal.shape, normal[:, 0])

        penetration_distance, x2, vIn, vOut, r_SS, r_PP = self._calculate_penetration_and_scattering(footprint, normal)
        jv_out_0, jv_out_1, ee_S, ee_P = self._calculate_jones_and_efield_directions(footprint, normal,
                                                                                        vIn, vOut, r_SS, r_PP)

        # update beam array with the new direction
        footprint.set_column(4, vOut[:, 0])
        footprint.set_column(5, vOut[:, 1])
        footprint.set_column(6, vOut[:, 2])
        # update beam array with the new electric fields
        footprint.set_jones_components(jv_out_0, jv_out_1, e_S=ee_S, e_P=ee_P)

        # advance the ray to the actual diffraction point (x2, at depth penetration_distance along
        # the incident direction) and update the optical path accordingly; rays for which
        # penetration_distance falls outside [0, thickness] do not satisfy the local Bragg
        # condition anywhere inside the crystal (transmitted, undiffracted -- not modeled by this
        # diffracted-beam element) and are flagged lost.
        footprint.set_column(1, x2[0])
        footprint.set_column(2, x2[1])
        footprint.set_column(3, x2[2])
        footprint.set_column(13, footprint.get_column(13) + penetration_distance)

        thickness = self.get_optical_element()._thickness
        lost = ~numpy.isfinite(penetration_distance) | (penetration_distance < 0.0) | (penetration_distance > thickness)
        footprint.set_column(10, numpy.where(lost, flag_lost_value, footprint.get_column(10)))

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

    def _calculate_penetration_and_scattering(self, footprint1, normal):
        """
        Computes, for each ray, the diffraction (penetration) depth inside a bent Laue crystal and
        the resulting exit direction, following the closed-form algorithm of Section 5
        ("Ray-tracing methodology for bent Laue crystals") of Guigay & Sanchez del Rio,
        "Geometrical focusing in Laue crystals" (main23.tex in this directory):

        - the local, strain-corrected reciprocal-lattice vector h'(t) along the ray path [Eqs.
          (hpP), (gradphi), (hpt)];
        - the diffraction depth t_d at which the ray meets the exact local Bragg condition, in
          closed form [Eq. (td), using the crystal-level constant of Eq. (phi-sum)];
        - the exit direction from h'(t_d) [Eq. (ih)].

        Reflectivity is only a placeholder (r_S = r_P = 1): the reference paper is purely
        geometric/kinematic and gives no reflectivity model; TODO: introduce one (e.g.
        Penning-Polder) in a second phase.

        Parameters
        ----------
        footprint1 : instance of S4Beam
            The beam at the entrance-surface intercept (position = entry point, direction =
            incident direction), in the local crystal reference system.
        normal : numpy array shape (3, nrays)
            The surface normal at the intercept, as returned by calculate_intercept_on_beam().

        Returns
        -------
        tuple
            (penetration_distance, x2, vIn, vOut, r_S, r_P): penetration_distance and x2 (shape
            (3, nrays)) are the diffraction depth (along the incident direction) and diffraction
            point for each ray; vIn, vOut (shape (nrays, 3)) are the incident and exit directions;
            r_S, r_P (shape (nrays,), complex) are the (placeholder) sigma/pi reflectivities.
        """
        footprint = footprint1.duplicate()
        v1 = footprint.get_columns([4, 5, 6])  # shape (3,nrays), incident direction i0
        x1 = footprint.get_columns([1, 2, 3])  # shape (3,nrays), entry point P
        nrays = v1.shape[1]

        energies = footprint.get_photon_energy_eV()

        # using crystalpy Vector
        direction_in = Vector(v1[0], v1[1], v1[2])

        # Warning:
        # S4Conic, S4Toroid give normal_z < 0 for concave surfaces (and >0 for convex)
        # S4mesh give always normal_z > 0
        # We need for crystalpy the upwards normal
        soe = self.get_optical_element()
        ccc = soe.get_optical_surface_instance()
        if isinstance(ccc, S4Mesh):
            surface_normal = Vector(normal[0], normal[1], normal[2])  # normal is outwards!
        elif isinstance(ccc, S4Toroid):
            if ccc.f_torus == 0 or ccc.f_torus == 2:
                surface_normal = Vector(normal[0], normal[1], normal[2]).scalarMultiplication(-1.0) # normal is inwards!
            else:
                surface_normal = Vector(normal[0], normal[1], normal[2])  # normal is outwards!
        else: # todo: check with convex surfaces
            surface_normal = Vector(normal[0], normal[1], normal[2]).scalarMultiplication(-1.0)  # normal is inwards!

        if is_debug():
            print(">>> direction_in: ", direction_in.componentsStack()[:, 0])
            print(">>> surface_normal: ", surface_normal.componentsStack()[:, 0])

        r_S = numpy.ones(nrays, dtype=complex) # TODO: introduce model, Penning-Polder?
        r_P = numpy.ones(nrays, dtype=complex) # TODO: introduce model, Penning-Polder?

        if is_verbose():
            print(">> r_S: ", r_S[0], numpy.abs(r_S[0]) ** 2)
            print(">> r_P: ", r_P[0], numpy.abs(r_P[0]) ** 2)

        #
        # local (zeta, eta) meridional-plane basis. zeta_hat is the paper's inward normal n --
        # reusing the surface_normal computed above (already handling the S4Mesh/S4Toroid/S4Conic
        # sign differences) rather than a hardcoded local axis. eta_hat is tied to the crystal's
        # own (fixed) local frame -- not to the incident ray's direction, to avoid a sign ambiguity
        # for divergent/misaligned beams -- via a cross product with the sagittal direction, which
        # is local X throughout shadow4 (S4SphereOpticalElementDecorator / S4Beam.rotate()
        # convention): eta_hat = zeta_hat x sagittal_hat.
        #
        # zeta_hat is re-oriented, if needed, so that the incident ray always has a positive
        # component along it (i0 . zeta_hat > 0): a ray reaching the entrance surface must be
        # travelling INTO the crystal bulk, i.e. along +zeta by the paper's own definition of zeta
        # (Section 2.1). This self-corrects for the surface_normal sign, which the "todo: check with
        # convex surfaces" comment above already flags as unverified for S4Conic, and which was
        # found (empirically) to flip between convexity=Convexity.UPWARD and
        # convexity=Convexity.DOWNWARD for the same incident direction.
        #
        # The zeta_hat x sagittal_hat order (rather than sagittal_hat x zeta_hat) was fixed against
        # an experimental reference diagram (Si 111, 10 keV, symmetric Laue, theta_B = 11 deg): with
        # outward normal n = +Z and H = +Y (in the surface, chi = 0, symmetric), h(O) =
        # h*(sin(alpha0), -cos(alpha0)) in (zeta, eta) only reproduces H = +h*Y for alpha0 = 0 with
        # this cross-product order (eta_hat = -Y here); the opposite order gives eta_hat = +Y, which
        # would put h(O) along -Y -- the wrong side.
        #
        zeta_hat = surface_normal.componentsStack()  # (3, nrays), already unit
        zeta_sign = numpy.sign(numpy.sum(v1 * zeta_hat, axis=0))
        zeta_hat = zeta_hat * zeta_sign[numpy.newaxis, :]

        sagittal_hat = numpy.zeros_like(zeta_hat)
        sagittal_hat[0, :] = 1.0
        eta_hat = numpy.cross(zeta_hat, sagittal_hat, axis=0)
        eta_hat = eta_hat / numpy.linalg.norm(eta_hat, axis=0, keepdims=True)

        i0_zeta = numpy.sum(v1 * zeta_hat, axis=0)
        i0_eta  = numpy.sum(v1 * eta_hat,  axis=0)

        # entry point in (zeta, eta) coordinates, relative to the crystal-design reference point O
        # (the local o.e. origin)
        zeta_P = numpy.sum(x1 * zeta_hat, axis=0)
        eta_P  = numpy.sum(x1 * eta_hat,  axis=0)

        #
        # crystal / elastic constants (crystal-level, not per-ray)
        #
        # soe._asymmetry_angle follows the SHADOW convention: the angle between the Bragg planes
        # and the surface (90 deg for symmetric Laue). The paper's alpha (main23.tex, "the rotation
        # of h from its direction in the symmetric case") is instead the angle between H and the
        # surface -- zero for symmetric Laue -- i.e. the complement: alpha = pi/2 - asymmetry_angle.
        alpha0 = numpy.pi / 2.0 - soe._asymmetry_angle
        rho = soe._poisson_ratio
        thickness = soe._thickness  # normal thickness, m

        shape = soe.get_surface_shape()
        R = shape.get_radius() if hasattr(shape, "get_radius") else numpy.inf

        d_spacing = self._crystalpy_diffraction_setup.dSpacingSI()  # m
        h_mod = 2.0 * numpy.pi / d_spacing                          # 1/m

        # design (tuning) Bragg angle -- used, as in the reference paper, as a fixed crystal-level
        # property in the closed-form elastic terms below (Eqs. phi-sum, td). Uses the UNCORRECTED
        # (pure vacuum/kinematic) Bragg angle, angleBragg -- not angleBraggCorrected -- since the
        # refractive-index correction is a dynamical-theory effect outside the scope of this purely
        # geometric/kinematic derivation.
        if soe._f_phot_cent == 0:
            energy_design = soe._phot_cent
        else:
            energy_design = codata.h * codata.c / codata.e * 1e2 / (soe._phot_cent * 1e-8)
        theta_B = self._crystalpy_diffraction_setup.angleBragg(energy_design)
        if isinstance(theta_B, (list, tuple, numpy.ndarray)): theta_B = theta_B[0]

        # per-ray photon energy -> wavenumber
        wavelength = codata.h * codata.c / codata.e / energies     # m
        k = 2.0 * numpy.pi / wavelength                             # 1/m

        #
        # local (strained) diffraction vector at the entry point P, Eq. (hpP) + (gradphi)
        #
        dphi_dzeta_P = (h_mod / R) * (rho * zeta_P * numpy.sin(alpha0) + eta_P * numpy.cos(alpha0))
        dphi_deta_P  = (h_mod / R) * (eta_P * numpy.sin(alpha0) + zeta_P * numpy.cos(alpha0))

        hP_zeta =  h_mod * numpy.sin(alpha0) - dphi_dzeta_P
        hP_eta  = -h_mod * numpy.cos(alpha0) - dphi_deta_P
        hP_mod = numpy.sqrt(hP_zeta ** 2 + hP_eta ** 2)

        #
        # variation of the diffraction vector along the ray, Eq. (hpt)
        #
        phi_0zeta = (h_mod / R) * (rho * numpy.sin(alpha0) * i0_zeta + numpy.cos(alpha0) * i0_eta)
        phi_0eta  = (h_mod / R) * (numpy.sin(alpha0) * i0_eta + numpy.cos(alpha0) * i0_zeta)
        ht_zeta = -phi_0zeta
        ht_eta  = -phi_0eta

        #
        # diffraction (penetration) depth t_d, Eq. (td), using the closed form of (phi_00+phi_0h)
        # for cylindrical bending, Eq. (phi-sum) -- a crystal-level constant (not per-ray).
        #
        denom = -(h_mod * numpy.cos(theta_B) / R) * (
                    2.0 * numpy.sin(theta_B - alpha0)
                    - (1.0 + rho) * numpy.sin(2.0 * alpha0) * numpy.cos(theta_B - alpha0))

        sin_theta_inc = -(hP_zeta * i0_zeta + hP_eta * i0_eta) / hP_mod   # geometric glancing angle at P
        sin_theta_local_bragg = hP_mod / (2.0 * k)                        # Bragg's law, local spacing

        with numpy.errstate(divide="ignore", invalid="ignore"):
            penetration_distance = -2.0 * h_mod * (sin_theta_inc - sin_theta_local_bragg) / denom #  !!!! CHANGED SIGN  !!!!!



        # if is_debug():
        print(">>> min, max, size penetration_distance [m]: ", penetration_distance.min(), penetration_distance.max(),
              penetration_distance.size)
        count = numpy.sum((penetration_distance >= 0) & (penetration_distance <= thickness))
        print("thickness, Number of elements in [0, thickness]:", thickness, count)
        count = numpy.sum((penetration_distance >= -thickness) & (penetration_distance <= 0))
        print("thickness, Number of elements in [-thickness, 0]:", thickness, count)
        count = numpy.sum((penetration_distance >= -10*thickness) & (penetration_distance <= 10*thickness))
        print("thickness, Number of elements in [-thickness, thichness]:", thickness, count)

        from srxraylib.plot.gol import plot
        plot(energies, penetration_distance, linestyle='', marker='.',
             xtitle="Energy [eV]", ytitle='Penetration depth [m]', show=0)

        #
        # interaction (diffraction) point x2 = x1 + t_d * v1
        #
        x2 = x1 + penetration_distance[numpy.newaxis, :] * v1  # shape (3, nrays)

        #
        # local diffraction vector at the interaction point, Eq. (hp-at-t), and exit direction,
        # Eq. (ih): i_h = (k0 + h'(t_d)) / k. The (zeta, eta) correction is embedded back into 3D
        # via zeta_hat, eta_hat; the incident ray's own sagittal component (already in v1/k0_vec)
        # is carried through unchanged, consistent with the purely meridional bending model.
        #
        h_zeta_td = hP_zeta + penetration_distance * ht_zeta
        h_eta_td  = hP_eta  + penetration_distance * ht_eta

        h_vec = h_zeta_td[numpy.newaxis, :] * zeta_hat + h_eta_td[numpy.newaxis, :] * eta_hat

        k0_vec = k[numpy.newaxis, :] * v1
        kh_vec = k0_vec + h_vec
        kh_mod = numpy.sqrt(kh_vec[0] ** 2 + kh_vec[1] ** 2 + kh_vec[2] ** 2)
        vOut_arr = (kh_vec / kh_mod).T  # shape (nrays, 3)

        if is_debug():
            scattering_angle = numpy.arccos(numpy.clip(numpy.sum(v1 * (kh_vec / kh_mod), axis=0), -1.0, 1.0))
            print(">>> scattering angle [deg]: ", numpy.degrees(scattering_angle[0]))
            print(">>> direction_out: ", vOut_arr[0])

        scattering_angle = numpy.arccos(numpy.clip(numpy.sum(v1 * (kh_vec / kh_mod), axis=0), -1.0, 1.0))
        from srxraylib.plot.gol import plot
        plot(energies, 0.5 * numpy.degrees(scattering_angle), linestyle='', marker='.',
             xtitle="Energy [eV]", ytitle='Half of Scattering angle [deg]', show=0)

        vIn = v1.T  # shape (npoints, 3)

        return penetration_distance, x2, vIn, vOut_arr, r_S, r_P

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
        if True:  # todo delete, only for test
            # c = e_S[:, 0] # cos of angle between e_S and the x axis
            c = vector_dot(e_S, ee_S)
            s = numpy.sqrt(1 - c ** 2)  # sin
            if is_verbose(): print(">>> s, c, angle: ", s[0], c[0], numpy.degrees(numpy.arctan2(s[0], c[0])))
            R00 = c
            R01 = -s
            R10 = s
            R11 = c
            if is_verbose(): print(">>> R angles: ", R00[0], R01[0], R10[0], R11[0])

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

if __name__ == "__main__":
    c = S4GeometricLaueCrystal(
            name="Undefined",
            boundary_shape=None,
            surface_shape=None,
            material="Si",
            miller_index_h=1,
            miller_index_k=1,
            miller_index_l=1,
            asymmetry_angle=numpy.pi/2,
            is_thick=0,
            thickness=15.5e-6,
            f_central=0,
            f_phot_cent=0,
            phot_cent=10000.0,
            file_refl="",
            f_bragg_a=True,
            f_ext=0,)

    print(c.info())


    ce = S4GeometricLaueCrystalElement(optical_element=c)
    print(ce.info())

    ce = S4GeometricLaueCrystalElement()
    ce.set_optical_element(c)

    print(ce.info())
