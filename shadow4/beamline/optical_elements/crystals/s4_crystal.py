import numpy

from shadow4.syned.shape import Rectangle

from shadow4.syned.element_coordinates import ElementCoordinates
from syned.beamline.optical_elements.crystals.crystal import Crystal, DiffractionGeometry

from shadow4.physical_models.prerefl.prerefl import PreRefl

from shadow4.beamline.s4_beamline_element import S4BeamlineElement



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

        if self._f_central:
            self.align_crystal()

        if self._f_johansson and self._f_ext:
            self.set_default_johansson_radius()

    def load_material_file(self):
        raise Exception(NotImplementedError)

    def align_crystal(self):
        raise Exception(NotImplementedError)

    def set_default_johansson_radius(self):
        raise Exception(NotImplementedError)



class S4CrystalElement(S4BeamlineElement):
    
    def __init__(self, optical_element=None, coordinates=None):
        super().__init__(optical_element if optical_element is not None else S4Crystal(),
                         coordinates if coordinates is not None else ElementCoordinates())
    
    def trace_beam(self, beam_in, flag_lost_value=-1):

        p = self.get_coordinates().p()
        q = self.get_coordinates().q()
        theta_grazing1 = numpy.pi / 2 - self.get_coordinates().angle_radial()
        theta_grazing2 = numpy.pi / 2 - self.get_coordinates().angle_radial_out()
        alpha1 = self.get_coordinates().angle_azimuthal()

        #
        beam = beam_in.duplicate()
        #
        # put beam in mirror reference system
        #
        beam.rotate(alpha1, axis=2)
        beam.rotate(theta_grazing1, axis=1)
        beam.translation([0.0, -p * numpy.cos(theta_grazing1), p * numpy.sin(theta_grazing1)])

        #
        # reflect beam in the mirror surface
        #
        soe = self.get_optical_element()


        if not isinstance(soe, Crystal): # undefined
            raise Exception("Undefined Crystal")
        else:
            beam_mirr, normal = self.apply_crystal_diffraction(beam)


        #
        # apply mirror boundaries
        #
        beam_mirr.apply_boundaries_syned(soe.get_boundary_shape(), flag_lost_value=flag_lost_value)
        #
        # TODO" apply lens absorption
        #

        #
        # from element reference system to image plane
        #

        beam_out = beam_mirr.duplicate()
        beam_out.change_to_image_reference_system(theta_grazing2, q)

        return beam_out, beam_mirr

    def apply_crystal_diffraction(self, beam):
        raise NotImplementedError()


if __name__ == "__main__":
    c = S4Crystal()
    print(c.info())
