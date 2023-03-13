import numpy
from syned.beamline.element_coordinates import ElementCoordinates

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.optical_elements.crystals.s4_crystal import S4CrystalElement, S4Crystal
from shadow4.beamline.s4_optical_element import S4PlaneOpticalElementDecorator

from syned.beamline.optical_elements.crystals.crystal import DiffractionGeometry

class S4PlaneCrystal(S4Crystal, S4PlaneOpticalElementDecorator):
    def __init__(self,
                 name="Plane crystal",
                 boundary_shape=None,
                 # surface_shape=None,
                 material=None,
                 diffraction_geometry=DiffractionGeometry.BRAGG,  # ?? not supposed to be in syned...
                 miller_index_h=1,
                 miller_index_k=1,
                 miller_index_l=1,
                 asymmetry_angle=0.0,
                 thickness=0.010,  ###########################
                 f_central=False,
                 f_phot_cent=0,
                 phot_cent=8000.0,
                 file_refl="",
                 f_bragg_a=False,
                 # a_bragg=0.0,
                 f_johansson=False,
                 r_johansson=1.0,
                 f_mosaic=False,
                 spread_mos=0.4 * numpy.pi / 180,
                 f_ext=0,
                 material_constants_library_flag=0,  # 0=xraylib, 1=dabax
                                                     # 2=shadow preprocessor file v1
                                                     # 3=shadow preprocessor file v2
                 ):

        S4PlaneOpticalElementDecorator.__init__(self)
        S4Crystal.__init__(self,
                           name=name,
                           boundary_shape=boundary_shape,
                           surface_shape=self.get_surface_shape_instance(),
                           material=material,
                           diffraction_geometry=diffraction_geometry,  # ?? not supposed to be in syned...
                           miller_index_h=miller_index_h,
                           miller_index_k=miller_index_k,
                           miller_index_l=miller_index_l,
                           asymmetry_angle=asymmetry_angle,
                           thickness=thickness,  ###########################
                           f_central=f_central,
                           f_phot_cent=f_phot_cent,
                           phot_cent=phot_cent,
                           file_refl=file_refl,
                           f_bragg_a=f_bragg_a,
                           # a_bragg=0.0,
                           f_johansson=f_johansson,
                           r_johansson=r_johansson,
                           f_mosaic=f_mosaic,
                           spread_mos=spread_mos,
                           f_ext=f_ext,
                           material_constants_library_flag=material_constants_library_flag,
                           )

    def get_input_dictionary(self):
        fmt = {
            "name": self.get_name(),
            "material": self._material                       ,
            "diffraction_geometry": self._diffraction_geometry           ,
            "miller_index_h": self._miller_index_h                 ,
            "miller_index_k": self._miller_index_k                 ,
            "miller_index_l": self._miller_index_l                 ,
            "asymmetry_angle": self._asymmetry_angle                ,
            "thickness": self._thickness                      ,
            "f_mosaic": self._f_mosaic                       ,
            "f_central": self._f_central                      ,
            "f_phot_cent": self._f_phot_cent                    ,
            "phot_cent": self._phot_cent                      ,
            "file_refl": self._file_refl                      ,
            "f_bragg_a": self._f_bragg_a                      ,
            "f_johansson": self._f_johansson                    ,
            "spread_mos": self._spread_mos                     ,
            "f_ext": self._f_ext                          ,
            "r_johansson": self._r_johansson                    ,
            "material_constants_library_flag": self._material_constants_library_flag,
            }
        return fmt

    def to_python_code(self, data=None):

        txt = "\nfrom shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystal"

        fmt = self.get_input_dictionary()
        txt_pre = """\noptical_element = S4PlaneCrystal(name='{name}',
    boundary_shape=None,material='{material}',diffraction_geometry={diffraction_geometry},
    miller_index_h={miller_index_h},miller_index_k={miller_index_k},miller_index_l={miller_index_l},
    asymmetry_angle={asymmetry_angle},
    thickness={thickness},
    f_central={f_central},f_phot_cent={f_phot_cent},phot_cent={phot_cent},
    file_refl='{file_refl}',
    f_bragg_a={f_bragg_a},
    f_johansson={f_johansson},r_johansson = {r_johansson},
    f_mosaic={f_mosaic},spread_mos={spread_mos},
    f_ext={f_ext},
    material_constants_library_flag={material_constants_library_flag},
    )"""
        txt += txt_pre.format(**fmt)

        return txt

class S4PlaneCrystalElement(S4CrystalElement):
    def __init__(self,
                 optical_element : S4PlaneCrystal = None,
                 coordinates : ElementCoordinates = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4PlaneCrystal(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         input_beam=input_beam)

    def to_python_code(self, data=None):
        txt = "\n\n# optical element number XX"
        txt += self.get_optical_element().to_python_code()
        coordinates = self.get_coordinates()
        txt += "\nfrom syned.beamline.element_coordinates import ElementCoordinates"
        txt += "\ncoordinates=ElementCoordinates(p=%g,q=%g,angle_radial=%g)" % \
               (coordinates.p(), coordinates.q(), coordinates.angle_radial())
        txt += "\nfrom shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystalElement"
        txt += "\nbeamline_element = S4PlaneCrystalElement(optical_element=optical_element,coordinates=coordinates,input_beam=beam)"
        txt += "\n\nbeam, mirr = beamline_element.trace_beam()"
        return txt

if __name__ == "__main__":
    # c = S4PlaneCrystal(
    #         name="Undefined",
    #         boundary_shape=None,
    #         material="Si",
    #         diffraction_geometry=DiffractionGeometry.BRAGG, #?? not supposed to be in syned...
    #         miller_index_h=1,
    #         miller_index_k=1,
    #         miller_index_l=1,
    #         asymmetry_angle=0.0,
    #         thickness=0.010, ###########################
    #         f_central=False,
    #         f_phot_cent=0,
    #         phot_cent=8000.0,
    #         file_refl="",
    #         f_bragg_a=False,
    #         # a_bragg=0.0,
    #         f_johansson=False,
    #         r_johansson=1.0,
    #         f_mosaic=False,
    #         spread_mos=0.4*numpy.pi/180,
    #         f_ext=0,)
    # # print(c.info())
    # # print(c.to_python_code())
    #
    #
    # ce = S4PlaneCrystalElement(optical_element=c)
    # print(ce.info())
    # print(ce.to_python_code())

    cc = S4PlaneCrystalElement()
