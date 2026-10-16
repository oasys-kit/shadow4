"""
The s4 phase plate crystal base class (optical element and beamline element).

A phase plate (or phase retarder, phase shifter) is a thin perfect crystal plate set close to (a few Darwin widths
off) a Bragg reflection, and used in transmission: the transmitted (or forward-diffracted) beam has different
complex amplitudes t_sigma, t_pi for the sigma and pi polarizations (referred to the diffraction plane), and
therefore a phase retardation that converts, e.g., linear into circular polarization when the diffraction plane is
at 45 deg from the incident polarization.
"""
import numpy

from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.optical_elements.crystals.crystal import Crystal, DiffractionGeometry
from syned.beamline.shape import Rectangle, Ellipse

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.s4_beamline_element import S4BeamlineElement
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
from shadow4.tools.arrayofvectors import vector_modulus, vector_dot, vector_cross

from crystalpy.diffraction.DiffractionSetupXraylib import DiffractionSetupXraylib
from crystalpy.diffraction.DiffractionSetupDabax import DiffractionSetupDabax
from crystalpy.diffraction.DiffractionSetupShadowPreprocessorV1 import DiffractionSetupShadowPreprocessorV1
from crystalpy.diffraction.DiffractionSetupShadowPreprocessorV2 import DiffractionSetupShadowPreprocessorV2
from crystalpy.diffraction.GeometryType import BraggTransmission
from crystalpy.diffraction.PerfectCrystalDiffraction import PerfectCrystalDiffraction
from crystalpy.util.Vector import Vector
from crystalpy.util.ComplexAmplitudePhoton import ComplexAmplitudePhoton

from dabax.dabax_xraylib import DabaxXraylib

from shadow4.optical_surfaces.s4_mesh import S4Mesh
from shadow4.optical_surfaces.s4_toroid import S4Toroid
from shadow4.tools.logger import is_verbose

import scipy.constants as codata


class S4PhasePlateCrystal(Crystal):
    """
    Shadow4 Phase Plate Crystal Class
    This is a base class for perfect crystal phase plates, consisting in a thin crystal plate set to diffract in
    reflection geometry (Bragg), but using the transmitted (or forward-diffracted) beam.

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
    thickness : float, optional
        The crystal thickness in m.
    f_central : int, optional
        Flag for autosetting the crystal to the corrected Bragg angle (plus the detuning).
    f_phot_cent : int, optional
        0: setting photon energy in eV, 1:setting photon wavelength in A.
    phot_cent : float, optional
        for f_central=1, the value of the photon energy (f_phot_cent=0) or photon wavelength (f_phot_cent=1).
    f_detune : int, optional
        For f_central=1, detuning of the crystal angle from the corrected Bragg angle:
        0: no detuning, 1: detuning angle detune_angle, 2: detuning of detune_half_darwin_widths Darwin half widths.
    detune_angle : float, optional
        For f_detune=1, the detuning angle in rad (positive: larger glancing angle).
    detune_half_darwin_widths : float, optional
        For f_detune=2, the detuning in units of the Darwin half width (sigma polarization, symmetric reflection,
        as given by crystalpy darwinHalfwidthS), as in Hirano et al. (1993) NIM A 336, 343.
    material_constants_library_flag : int, optional
        Flag for indicating the origin of the crystal data:
        0: xraylib, 1: dabax, 2: preprocessor file v1, 3: preprocessor file v2.
    file_refl : str, optional
        for material_constants_library_flag=2,3, the name of the file containing the crystal parameters.
    dabax : None or instance of DabaxXraylib,
        The pointer to the dabax library  (used for material_constants_library_flag=1).

    Returns
    -------
    instance of S4PhasePlateCrystal.
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
                 thickness=100e-6,
                 f_central=0,
                 f_phot_cent=0,
                 phot_cent=8000.0,
                 f_detune=0,
                 detune_angle=0.0,
                 detune_half_darwin_widths=0.0,
                 material_constants_library_flag=0, # 0=xraylib, 1=dabax
                                                    # 2=shadow preprocessor file v1
                                                    # 3=shadow preprocessor file v2
                 file_refl="",
                 dabax=None,
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
        self._f_detune = f_detune
        self._detune_angle = detune_angle
        self._detune_half_darwin_widths = detune_half_darwin_widths
        self._f_bragg_a = f_bragg_a
        self._material_constants_library_flag = material_constants_library_flag
        self._file_refl = file_refl
        self._dabax = dabax

        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._add_support_text([
                    ("f_central",           "S4: autotuning",                              ""),
                    ("f_phot_cent",         "S4: for f_central=1: tune to eV(0) or A (1)", ""),
                    ("phot_cent",           "S4: for f_central=1: value in eV or A",       ""),
                    ("f_detune",            "S4: for f_central=1: detuning 0=No, 1=angle, 2=Darwin half widths", ""),
                    ("detune_angle",        "S4: for f_detune=1: detuning angle",          "rad"),
                    ("detune_half_darwin_widths", "S4: for f_detune=2: detuning in Darwin half widths", ""),
                    ("f_bragg_a",           "S4: use asymmetric cut",                      ""),
                    ("material_constants_library_flag", "S4: crystal data from: 0=xraylib, 1=dabax, 2=file v1, 3=file v1", ""),
                    ("file_refl",           "S4: preprocessor file name",                  ""),
            ] )

    def get_info(self):
        """
        Returns the specific information of the S4 phase plate crystal optical element.

        Returns
        -------
        str
        """
        txt = "\n\n"
        txt += "PHASE PLATE CRYSTAL (Bragg geometry, transmitted beam)\n"
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

        txt += "Crystal thickness: %g m\n" % self._thickness

        if self._f_central == 0:
            txt += "Using EXTERNAL incidence angle.\n"
        else:
            txt += "Using INTERNAL or calculated incidence angle for "
            if self._f_phot_cent == 0:
                txt += "photon energy %.6f eV\n" % self._phot_cent
            else:
                txt += "photon wavelength %f A\n" % (self._phot_cent)
            if self._f_detune == 1:
                txt += "Detuning angle: %g rad\n" % self._detune_angle
            elif self._f_detune == 2:
                txt += "Detuning: %g Darwin half widths\n" % self._detune_half_darwin_widths

        txt += "\n"
        ss = self.get_surface_shape()
        if ss is None:
            txt += "Surface shape is: Plane (** UNDEFINED?? **)\n"
        else:
            txt += "Surface shape is: %s\n" % ss.__class__.__name__

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


class S4PhasePlateCrystalElement(S4BeamlineElement):
    """
    The base class for Shadow4 phase plate crystal element.
    It is made of a S4PhasePlateCrystal and an ElementCoordinates instance. It also includes the input beam.

    The rays are transmitted without change of direction (the refraction is neglected), from the intercept with the
    crystal surface. Their electric fields are multiplied by the transmission amplitudes t_sigma, t_pi calculated by
    crystalpy (Bragg transmission geometry, Takagi-Taupin solution), with the Jones calculus in the basis of the
    diffraction plane (sigma direction v_in x H, with H the diffraction vector). The image plane is perpendicular to
    the transmitted beam (angle_radial_out = pi - angle_radial).

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
    instance of S4PhasePlateCrystalElement.
    """

    def __init__(self,
                 optical_element : S4PhasePlateCrystal = None,
                 coordinates : ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4PhasePlateCrystal(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         movements=movements,
                         input_beam=input_beam)

        self._crystalpy_diffraction_setup = None

    def set_crystalpy_diffraction_setup(self):
        """
        Sets the crystalpy DiffractionSetup (Bragg transmission geometry).
        """
        oe = self.get_optical_element()

        kwargs = dict(geometry_type=BraggTransmission(),
                      crystal_name=oe._material,
                      thickness=oe._thickness,
                      miller_h=oe._miller_index_h,
                      miller_k=oe._miller_index_k,
                      miller_l=oe._miller_index_l,
                      asymmetry_angle=oe._asymmetry_angle,
                      azimuthal_angle=0.0)

        if oe._material_constants_library_flag == 0:
            if is_verbose(): print("\nCreating a diffraction setup (XRAYLIB) for material:", oe._material)
            diffraction_setup = DiffractionSetupXraylib(**kwargs)
        elif oe._material_constants_library_flag == 1:
            if is_verbose(): print("\nCreating a diffraction setup (DABAX) for material:", oe._material)
            diffraction_setup = DiffractionSetupDabax(dabax=oe._dabax, **kwargs)
        elif oe._material_constants_library_flag == 2:
            if is_verbose(): print("\nCreating a diffraction setup (shadow preprocessor file V1)...")
            diffraction_setup = DiffractionSetupShadowPreprocessorV1(preprocessor_file=oe._file_refl, **kwargs)
        elif oe._material_constants_library_flag == 3:
            if is_verbose(): print("\nCreating a diffraction setup (shadow preprocessor file V2)...")
            diffraction_setup = DiffractionSetupShadowPreprocessorV2(preprocessor_file=oe._file_refl, **kwargs)
        else:
            raise NotImplementedError

        self._crystalpy_diffraction_setup = diffraction_setup

    def get_detuning_angle(self, energy):
        """
        Returns the detuning angle (rad) from the corrected Bragg angle.

        Parameters
        ----------
        energy : float
            The photon energy in eV.

        Returns
        -------
        float
        """
        oe = self.get_optical_element()
        if oe._f_detune == 0:
            return 0.0
        elif oe._f_detune == 1:
            return oe._detune_angle
        elif oe._f_detune == 2:
            darwin_half_width = self._crystalpy_diffraction_setup.darwinHalfwidthS(energy)
            return oe._detune_half_darwin_widths * float(numpy.ravel(darwin_half_width)[0])
        else:
            raise ValueError("f_detune must be 0, 1 or 2")

    def align_crystal(self):
        """
        Sets the incident angle to the corrected Bragg angle plus the detuning (for f_central=1), and the exit angle
        for the transmitted (undeviated) beam: angle_radial_out = pi - angle_radial.
        """
        oe = self.get_optical_element()
        coor = self.get_coordinates()

        if oe is None:
            raise Exception("Undefined optical element")

        angle_radial, _, angle_azimuthal = coor.get_angles()

        if oe._f_central:
            if oe._f_phot_cent == 0:
                energy = oe._phot_cent
            else:
                energy = codata.h * codata.c / codata.e * 1e2 / (oe._phot_cent * 1e-8)

            setting_angle = self._crystalpy_diffraction_setup.angleBraggCorrected(energy)
            if isinstance(setting_angle, (list, tuple, numpy.ndarray)): setting_angle = setting_angle[0]

            detuning = self.get_detuning_angle(energy)
            theta_in_grazing = setting_angle + oe._asymmetry_angle + detuning
            angle_radial = numpy.pi / 2 - theta_in_grazing

            if is_verbose():
                print("    align_crystal: Bragg angle (corrected) for E=%f eV is %f deg" % (energy, numpy.degrees(setting_angle)))
                print("    align_crystal: detuning %g urad" % (1e6 * detuning))
                print("    align_crystal: grazing incident angle [deg]: ", numpy.degrees(theta_in_grazing))
        else:
            if is_verbose(): print("align_crystal: f_central=0: using the given incident angle")

        # transmitted beam: the image plane is perpendicular to the (undeviated) incident beam
        coor.set_angles(angle_radial=angle_radial,
                        angle_radial_out=numpy.pi - angle_radial,
                        angle_azimuthal=angle_azimuthal)

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

        if self._crystalpy_diffraction_setup is None:
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

        #
        # put input_beam in crystal reference system
        #
        if change_reference_system_in:
            input_beam.rotate(alpha1,         axis=2)
            input_beam.rotate(theta_grazing1, axis=1)
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
        # crystal transmission
        #
        footprint, normal = self._apply_crystal_transmission(input_beam)

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

        return output_beam, footprint

    def _apply_crystal_transmission(self, input_beam):
        """
        Calculates the intercepts with the crystal surface, and the electric fields of the transmitted rays.

        Parameters
        ----------
        input_beam : instance of S4Beam
            The incident beam in the crystal reference frame.

        Returns
        -------
        tuple
            (footprint, normal): the beam at the crystal surface with the transmitted fields, and the normals
            (3, nrays) at the intercept points.
        """
        footprint, normal = self.get_optical_element().get_optical_surface_instance().calculate_intercept_on_beam(input_beam)

        v_in, t_s, t_p, vector_h = self._calculate_transmission_amplitudes(footprint, normal)

        # Jones calculus: transmission amplitudes in the basis of the diffraction plane (k0, kh)
        sigma_axis = vector_cross(v_in, vector_h)
        j0, j1, e_S, e_P = footprint.calculate_jones_scattering(v_in, t_s, t_p, sigma_axis=sigma_axis)
        footprint.set_jones_components(j0, j1, e_S=e_S, e_P=e_P)

        return footprint, normal

    def _calculate_transmission_amplitudes(self, footprint, normal):
        """
        Calculates with crystalpy the complex transmission amplitudes for the sigma and pi polarizations.

        Parameters
        ----------
        footprint : instance of S4Beam
            The beam at the crystal surface (directions are the incident directions).
        normal : numpy array
            The normals to the surface at the intercept points, shape (3, nrays).

        Returns
        -------
        tuple
            (v_in, t_s, t_p, vector_h): incident (= transmitted) directions (nrays, 3), complex transmission
            amplitudes (nrays), and diffraction vectors H (nrays, 3).
        """
        v1 = footprint.get_columns([4, 5, 6])
        energies = footprint.get_photon_energy_eV()

        photons_in = ComplexAmplitudePhoton(energies,
                                            Vector(v1[0], v1[1], v1[2]),
                                            Esigma=numpy.ones(energies.size, dtype=complex),
                                            Epi=numpy.ones(energies.size, dtype=complex))

        # crystalpy needs the normal pointing into the crystal (same convention as in S4Crystal)
        ccc = self.get_optical_element().get_optical_surface_instance()
        if isinstance(ccc, S4Mesh):
            surface_normal = Vector(normal[0], normal[1], normal[2])
        elif isinstance(ccc, S4Toroid):
            if ccc.f_torus == 0 or ccc.f_torus == 2:
                surface_normal = Vector(normal[0], normal[1], normal[2]).scalarMultiplication(-1.0)
            else:
                surface_normal = Vector(normal[0], normal[1], normal[2])
        else:
            surface_normal = Vector(normal[0], normal[1], normal[2]).scalarMultiplication(-1.0)

        bragg_normal = surface_normal.getVectorH(
            surface_normal,
            self._crystalpy_diffraction_setup.dSpacingSI(),
            asymmetry_angle=self._crystalpy_diffraction_setup.asymmetryAngle(),
            azimuthal_angle=self._crystalpy_diffraction_setup.azimuthalAngle())

        perfect_crystal = PerfectCrystalDiffraction.initializeFromDiffractionSetupAndEnergy(
            self._crystalpy_diffraction_setup,
            energies,
            geometry_type=None,          # from the setup: BraggTransmission
            bragg_normal=bragg_normal,
            surface_normal=surface_normal,
            thickness=None,
            d_spacing=None,
        )

        coeffs = perfect_crystal.calculateDiffraction(photons_in, calculation_method=1, is_thick=0)
        t_s = numpy.asarray(coeffs["S"], dtype=complex) * numpy.ones(energies.size)
        t_p = numpy.asarray(coeffs["P"], dtype=complex) * numpy.ones(energies.size)

        if is_verbose():
            print(">> t_s: ", t_s[0], numpy.abs(t_s[0]) ** 2)
            print(">> t_p: ", t_p[0], numpy.abs(t_p[0]) ** 2)

        return v1.T, t_s, t_p, bragg_normal.components().T


if __name__ == "__main__":
    c = S4PhasePlateCrystal(
            name="Undefined",
            boundary_shape=None,
            surface_shape=None,
            material="Si",
            miller_index_h=1,
            miller_index_k=1,
            miller_index_l=1,
            asymmetry_angle=0.0,
            thickness=100e-6,
            f_central=0,
            f_phot_cent=0,
            phot_cent=8000.0,
            file_refl="",
            f_bragg_a=False,)

    print(c.info())
    ce = S4PhasePlateCrystalElement(optical_element=c)
    print(ce.info())
