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

from crystalpy.diffraction.GeometryType import BraggDiffraction

from dabax.dabax_xraylib import DabaxXraylib

from shadow4.tools.logger import is_verbose

import scipy.constants as codata

class S4PolariserCrystal(Crystal):
    """
    Shadow4 Crystal Class
    This is a base class for perfect crystal Polariser in transmission geometry (Laue), using the transmitted beam.

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
    dabax : None or instance of DabaxXraylib,
        The pointer to the dabax library  (used for material_constants_library_flag=1).

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
                 f_central=0,
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
        self._f_bragg_a = f_bragg_a
        self._is_thick = is_thick
        self._f_ext = f_ext
        self._material_constants_library_flag = material_constants_library_flag
        self._file_refl = file_refl

        self._dabax = dabax
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

    def _get_dabax_txt(self):
        if self._material_constants_library_flag == 1:
            if isinstance(self._dabax, DabaxXraylib):
                dabax_txt = 'DabaxXraylib(file_f0="%s", file_f1f2="%s")' % (self._dabax.get_file_f0(), self._dabax.get_file_f1f2())
            else:
                dabax_txt = "DabaxXraylib()"
        else:
            dabax_txt = "None"

        return dabax_txt

class S4PolariserCrystalElement(S4BeamlineElement):
    """
    The base class for Shadow4 polariser crystal element.
    It is made of a S4PolariserCrystal and an ElementCoordinates instance. It also includes the input beam.

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
    instance of S4PolariserCrystalElement.
    """

    def __init__(self,
                 optical_element : S4PolariserCrystal = None,
                 coordinates : ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4PolariserCrystal(),
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
                                                 azimuthal_angle=0.0,
                                                 dabax=oe._dabax)
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
        print(">>>>>> change_reference_system: ", change_reference_system_in, change_reference_system_out)

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
        if change_reference_system_out:
            output_beam.change_to_image_reference_system(theta_grazing2, q)

            if is_verbose():
                b_S, b_P = output_beam.get_efield_directions()
                print("")
                print(">>> image e_S, mod e_s", b_S[0], vector_modulus(b_S)[0])
                print(">>> image e_P, mod e_P, e_S.e_P: ", b_P[0], vector_modulus(b_P)[0], vector_dot(b_S, b_P)[0])

        return output_beam, footprint

    def _apply_crystal_diffraction(self, input_beam):
        #
        # geometric and physics for the scattering process (transmission/Laue polariser):
        # to be completely rewritten.
        #
        pass

if __name__ == "__main__":
    c = S4PolariserCrystal(
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
            f_central=0,
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


    ce = S4PolariserCrystalElement(optical_element=c)
    print(ce.info())

    ce = S4PolariserCrystalElement()
    ce.set_optical_element(c)

    print(ce.info())
