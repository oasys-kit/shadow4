import numpy

from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.optical_elements.crystals.crystal import Crystal, DiffractionGeometry
from syned.beamline.shape import Rectangle, Ellipse

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.s4_beamline_element import S4BeamlineElement
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
from shadow4.tools.arrayofvectors import vector_modulus, vector_dot, vector_cross, vector_norm
from shadow4.tools.arrayofvectors import vector_multiply_scalar, vector_diff

from crystalpy.diffraction.DiffractionSetupXraylib import DiffractionSetupXraylib
from crystalpy.diffraction.DiffractionSetupDabax import DiffractionSetupDabax
from crystalpy.diffraction.DiffractionSetupShadowPreprocessorV1 import DiffractionSetupShadowPreprocessorV1
from crystalpy.diffraction.DiffractionSetupShadowPreprocessorV2 import DiffractionSetupShadowPreprocessorV2

from crystalpy.diffraction.GeometryType import BraggDiffraction
from crystalpy.util.Vector import Vector
from crystalpy.util.ComplexAmplitudePhoton import ComplexAmplitudePhoton
from crystalpy.diffraction.PerfectCrystalDiffraction import PerfectCrystalDiffraction

from shadow4.optical_surfaces.s4_mesh import S4Mesh
from shadow4.optical_surfaces.s4_toroid import S4Toroid

from shadow4.tools.logger import is_verbose, is_debug
from shadow4.tools.arrayofvectors import vector_modulus_square, vector_modulus, vector_norm, vector_rotate_around_axis

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
                 method_efields_management=0, # 0=S4, 1=S3
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

        self._method_efields_management = method_efields_management

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
                    ("method_efields_management", "flag 0:new in S4; 1=like S3",           ""),
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
            if isinstance(setting_angle, (list, tuple, numpy.ndarray)): setting_angle = setting_angle[0]

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
            if isinstance(theta_out, (list, tuple, numpy.ndarray)): theta_out = theta_out[0]

            if is_verbose(): print("    align_crystal: (normal) Reflection angle [SCATTERING EQUATION]: ", numpy.degrees(theta_out))
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
        footprint, normal = self._apply_crystal_diffraction(input_beam)

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
        output_beam.change_to_image_reference_system(theta_grazing2, q)

        if is_verbose():
            b_S, b_P = output_beam.get_efield_directions()
            print("")
            print(">>> image e_S, mod e_s", b_S[0], vector_modulus(b_S)[0])
            print(">>> image e_P, mod e_P, e_S.e_P: ", b_P[0], vector_modulus(b_P)[0], vector_dot(b_S, b_P)[0])

        return output_beam, footprint

    def _apply_crystal_diffraction(self, input_beam):
        #
        # geometric and physics for the scattering process:
        # reflect beam in the crystal and apply crystal reflectivity
        #

        if self.get_optical_element()._method_efields_management == 0: # S4
            footprint, normal = self.get_optical_element().get_optical_surface_instance().calculate_intercept_on_beam(input_beam)

            vIn, vOut, r_SS, r_PP = self._calculate_perfect_crystal_scattering(footprint, normal)
            jv_out_0, jv_out_1, ee_S, ee_P = self._calculate_jones_and_efield_directions(footprint, normal,
                                                                                            vIn, vOut, r_SS, r_PP)
            # update beam array with the new direction
            footprint.set_column(4, vOut[:, 0])
            footprint.set_column(5, vOut[:, 1])
            footprint.set_column(6, vOut[:, 2])
            # update beam array with the new electric fields
            footprint.set_jones_components(jv_out_0, jv_out_1, e_S=ee_S, e_P=ee_P)

            if is_verbose():
                print(">>>> Orthogonal footprint: ", footprint.efields_orthogonal(),
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

        else: # S3
            footprint, normal = self._apply_crystal_diffraction_and_reflectivities_S3(input_beam)

        return footprint, normal

    def _calculate_perfect_crystal_scattering(self, footprint1, normal):
        """
        Compute the scattering by the crystal beamline element.

        It calls crystalpy to compute the crystal parameters:

        - output direction vOut

        - reflectivities (complex amplitudes) r__S and r_P

        Parameters
        ----------
        footprint : S4Beam instance
            The beam in the local reference of the beamline element with only the intercepts calculated (directions
            and electric fields will be later updated with the results of this method).

        normal : numpy array
            The array with the normal to the surface for all rays.

        Returns
        -------
        tuple
            (vIn, vOut, r_S, r_P) : input and output direction with shape (npoints, 3) and complex amplitude
            reflectivity for sigma and pi.

        """

        footprint = footprint1.duplicate()
        v1 = footprint.get_columns([4, 5, 6])  # numpy.array(a3.getshcol([4,5,6]))

        #
        # direction and reflectivity calculation using crystalpy
        #
        energies = footprint.get_photon_energy_eV()

        # incident  crystalpy photon stack
        photons_in = ComplexAmplitudePhoton(
            energies,
            Vector(v1[0], v1[1], v1[2]),
            Esigma=numpy.ones(energies.size, dtype=complex),
            Epi   =numpy.ones(energies.size, dtype=complex),
            )

        # create crystalpy PerfectCrystalDiffraction instance
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
            thickness=None,
            d_spacing=None,
        )

        # Calculate outgoing Photon.

        # NEW USING Jones matrix JR
        outgoing_complex_amplitude_photon = perfect_crystal.calculatePhotonOut(photons_in,
                                                                                apply_reflectivity=True,
                                                                                calculation_method=1,
                                                                                is_thick=soe._is_thick,
                                                                                use_transfer_matrix=0
                                                                                )
        r_S = outgoing_complex_amplitude_photon.getComplexAmplitudeS()
        r_P = outgoing_complex_amplitude_photon.getComplexAmplitudeP()

        if is_verbose():
            print(">> r_S: ", r_S[0], numpy.abs(r_S[0]) ** 2)
            print(">> r_P: ", r_P[0], numpy.abs(r_P[0]) ** 2)
        vv = outgoing_complex_amplitude_photon.unitDirectionVector().components()  # (3, npoints)
        vOut = vv.T  # shape (npoints, 3)
        vIn = v1.T
        return vIn, vOut, r_S, r_P

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

    def _apply_crystal_diffraction_and_reflectivities_S3(self, beam): #TODO delete
        """
        Applies the changes in direction and in reflectivity in a beam due to crystal diffraction.
        It mimics the procedure in SHADOW3.
        Note that it is a new way in SHADOW4 using Jones algebra.

        Parameters
        ----------
        beam : instance of S4Beam
            the beam (already transformed to the local crystal reference system).

        Returns
        -------
        tuple
            (footprint, normal), with footprint an instance of S4Beam and the vector normal as numpy arrat (3,:)

        """
        #
        soe = self.get_optical_element()
        ccc = soe.get_optical_surface_instance()
        v1 = beam.get_columns([4, 5, 6])

        # intercept calculation
        footprint, normal = ccc.calculate_intercept_on_beam(beam)


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

        if is_verbose(): print(">> method_efields_management=1: traditional S3")
        outgoing_complex_amplitude_photon = perfect_crystal.calculatePhotonOut(photons_in,
                                                                                apply_reflectivity=True,
                                                                                calculation_method=1,
                                                                                is_thick=soe._is_thick,
                                                                                use_transfer_matrix=0
                                                                                )


        # copy data to vectors
        E_S = footprint.get_columns([7,8,9]).T
        E_P = footprint.get_columns([16, 17, 18]).T
        PhiS = footprint.get_column(14)
        PhiP = footprint.get_column(15)
        vNormal = normal.T
        vIn = v1.T
        vH = vNormal # todo: set for asymmetry!!!

        if is_verbose():
            print(">>>>> incident E_S, norm: ", E_S[0], vector_modulus(E_S)[0])
            print(">>>>> incident E_P, norm: ", E_P[0], vector_modulus(E_P)[0])

        #
        # STEP 1: get local sigma and pi directions v_S and v_P (uS and uP unitary vectors)
        #

        uS, uP = footprint.get_local_directions_sigma_pi(vH, vIn=vIn)

        #
        # calculate the electric field in the local coordinate system
        #

        # CALL	DOT	(AS_VEC,AS_TEMP,A11)	! matrix element of rotation
        # CALL	DOT	(AP_VEC,AS_TEMP,A12)	! matrix element of rotation
        # CALL	DOT	(AS_VEC,AP_TEMP,A21)	! matrix element of rotation
        # CALL	DOT	(AP_VEC,AP_TEMP,A22)	! matrix element of rotation
        # ! ** Now recompute the ampltitude and phase of the local S- and P- component.
        # AS_NEW	= SQRT(ABS(A11**2 + A12**2 +  &
        #         2.0D0*A11*A12*COS(PHS-PHP)))
        # AP_NEW	= SQRT(ABS(A21**2 + A22**2 +  &
        #         2.0D0*A21*A22*COS(PHS-PHP)))
        # CALL	SCALAR	(AS_TEMP,AS_NEW,AS_VEC)	! Local As vector
        # CALL	SCALAR	(AP_TEMP,AP_NEW,AP_VEC)	! Local Ap vector
        # PHTS	= A11*SIN(PHS) + A12*SIN(PHP)
        # PHBS	= A11*COS(PHS) + A12*COS(PHP)
        # PHTP	= A21*SIN(PHS) + A22*SIN(PHP)
        # PHBP	= A21*COS(PHS) + A22*COS(PHP)
        # CALL	ATAN_2	(PHTS,PHBS,PHS)		! Phase of local As vector
        # CALL	ATAN_2	(PHTP,PHBP,PHP)		! Phase of local Ap vector
        # ! C
        # ! C
        # CALL	DOT	(VVIN,VNOR,SIN_VAL)	! sin(graz. ang)
        # CALL	DOT	(Q_OUT,VNOR,SIN_REF)	! sin(graz.ref.ang)

        a11 = vector_dot(E_S, uS)
        a12 = vector_dot(E_P, uS)
        a21 = vector_dot(E_S, uP)
        a22 = vector_dot(E_P, uP)

        if is_verbose():
            print(">>>>    Matrix a11: ", a11[0])
            print(">>>>    Matrix a12: ", a12[0])
            print(">>>>    Matrix a21: ", a21[0])
            print(">>>>    Matrix a22: ", a22[0])

        M2_S = a11 ** 2 + a12 ** 2 + 2 * a11 * a12 * numpy.cos(PhiS - PhiP)
        M2_P = a21 ** 2 + a22 ** 2 + 2 * a21 * a22 * numpy.cos(PhiS - PhiP)

        E_local_S = vector_multiply_scalar(uS, numpy.sqrt(M2_S))
        E_local_P = vector_multiply_scalar(uP, numpy.sqrt(M2_P))


        if is_verbose():
            print(">>>>> E_local_S, norm: ", E_local_S[0], vector_modulus(E_local_S)[0])
            print(">>>>> E_local_P, norm: ", E_local_P[0], vector_modulus(E_local_P)[0])

        local_PhiS = numpy.arctan2((a11 * numpy.sin(PhiS) + a12 * numpy.sin(PhiP)) , (a11 * numpy.cos(PhiS) + a12 * numpy.cos(PhiP)))
        local_PhiP = numpy.arctan2((a21 * numpy.sin(PhiS) + a22 * numpy.sin(PhiP)) , (a21 * numpy.cos(PhiS) + a22 * numpy.cos(PhiP)))
        # END OF STEP 1

        #
        # STEP 2: now apply crystal reflectivity
        #
        crystal_reflectivity_S = numpy.sqrt(outgoing_complex_amplitude_photon.getIntensityS())
        crystal_reflectivity_P = numpy.sqrt(outgoing_complex_amplitude_photon.getIntensityP())
        crystal_phase_S = outgoing_complex_amplitude_photon.getPhaseS()
        crystal_phase_P = outgoing_complex_amplitude_photon.getPhaseP()

        if is_verbose():
            print(">>>> crystal_reflectivity_S: ", crystal_reflectivity_S[0])
            print(">>>> crystal_reflectivity_P: ", crystal_reflectivity_P[0])


        E_diffracted_S = numpy.zeros_like(E_local_S)
        E_diffracted_P = numpy.zeros_like(E_local_S)

        E_diffracted_S[:, 0] = E_local_S[:, 0] * crystal_reflectivity_S
        E_diffracted_S[:, 1] = E_local_S[:, 1] * crystal_reflectivity_S
        E_diffracted_S[:, 2] = E_local_S[:, 2] * crystal_reflectivity_S

        E_diffracted_P[:, 0] = E_local_P[:, 0] * crystal_reflectivity_P
        E_diffracted_P[:, 1] = E_local_P[:, 1] * crystal_reflectivity_P
        E_diffracted_P[:, 2] = E_local_P[:, 2] * crystal_reflectivity_P

        Phi_diffracted_S = local_PhiS + crystal_phase_S
        Phi_diffracted_P = local_PhiP + crystal_phase_P

        #
        # update footprint with the new values (direction and phases)
        # (electric vector are not updated as they will change in STEP 3)
        #

        footprint.set_column(14, Phi_diffracted_S)
        footprint.set_column(15, Phi_diffracted_P)

        # print("orthogonal INCIDENT: ", footprint.efields_orthogonal(verbose=1))

        vv = outgoing_complex_amplitude_photon.unitDirectionVector().components() # (3, npoints)
        footprint.set_column(4, vv[0])
        footprint.set_column(5, vv[1])
        footprint.set_column(6, vv[2])

        # print("orthogonal BEFORE: ", footprint.efields_orthogonal(verbose=1))

        vOut = vv.T # (npoints, 3)

        #
        # STEP 3
        #

        # ! Electric vectors are changed to assure orthogonality with the new direction VVOUT
        # ! To conserve intensity, the moduli of Es and Ep must not change
        # ! AS_VEC and VVOUT are not orthogonal so a projection of S and P coordinates into the
        # ! new ones do not work as it may be a component of the electric field along VVOUT
        #
        # CALL CROSS_M_FLAG  (VVOUT,VNOR,AS_TEMP,M_FLAG) ! vector pp. to inc.pl.
        # CALL DOT (AS_VEC,AS_VEC,AS2)
        # CALL DOT (AP_VEC,AP_VEC,AP2)
        #
        # IF (M_FLAG.EQ.1) THEN
        #  IF (AS2.NE.0) THEN
        #    DO I=1,3
        #      AS_TEMP(I) = AS_VEC(I)
        #    END DO
        #  ELSE
        #   DO I=1,3
        #    AS_TEMP(I) = AP_VEC(I)
        #   END DO
        #  END IF
        # END IF
        #
        # CALL NORM   (AS_TEMP,AS_TEMP) ! Local unit As vector perp to vvout
        # CALL CROSS (AS_TEMP,VVOUT,AP_TEMP)
        # CALL NORM (AP_TEMP,AP_TEMP) ! Local unit Ap vector perp to vvout
        #
        # do i=1,3
        #   as_vec(i) = as_temp(i) * sqrt(as2)
        #   ap_vec(i) = ap_temp(i) * sqrt(ap2)
        # end do

        aS = vector_modulus(E_diffracted_S)
        aP = vector_modulus(E_diffracted_P)

        if is_verbose():
            print(">>>>>>  Modulus of initial E vectors s, p",   vector_modulus(E_S)[0], vector_modulus(E_P)[0])
            print(">>>>>>  Modulus of diffracted E vectors s, p", aS[0], aP[0])

        v_S = vector_cross(vOut, vH) # AS_TEMP
        v_Smod = vector_modulus(v_S)
        mask = (v_Smod == 0.0)
        if mask.any():
            print(">>>>>>>>>> FOUND A ZERO!!!!!")
            if v_Smod.sum() > 0:
                v_S[mask, 0] = E_diffracted_S[mask, 0]
                v_S[mask, 1] = E_diffracted_S[mask, 1]
                v_S[mask, 2] = E_diffracted_S[mask, 2]
            else:
                v_S[mask, 0] = E_diffracted_P[mask, 0]
                v_S[mask, 1] = E_diffracted_P[mask, 1]
                v_S[mask, 2] = E_diffracted_P[mask, 2]

        v_P = vector_cross(v_S, vOut) # AP_TEMP

        uS = vector_norm(v_S)
        uP = vector_norm(v_P)

        E_diffracted_S[:, 0] = uS[:, 0] * aS
        E_diffracted_S[:, 1] = uS[:, 1] * aS
        E_diffracted_S[:, 2] = uS[:, 2] * aS

        E_diffracted_P[:, 0] = uP[:, 0] * aP
        E_diffracted_P[:, 1] = uP[:, 1] * aP
        E_diffracted_P[:, 2] = uP[:, 2] * aP

        # END STEP 3

        #
        # replace electric fields and phases in footprint.
        # Note that reflectivity_S affects the X components, and reflectivity_P the Z ones.
        #
        footprint.set_column(7, E_diffracted_S[:, 0])
        footprint.set_column(8, E_diffracted_S[:, 1])
        footprint.set_column(9, E_diffracted_S[:, 2])
        footprint.set_column(16, E_diffracted_P[:, 0])
        footprint.set_column(17, E_diffracted_P[:, 1])
        footprint.set_column(18, E_diffracted_P[:, 2])

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
