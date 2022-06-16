import numpy
import scipy.constants as codata

from syned.beamline.optical_elements.gratings.grating import GratingVLS
from shadow4.syned.shape import Plane, Sphere
from shadow4.syned.element_coordinates import ElementCoordinates # TODO from shadow4.syned.element_coordinates

from shadow4.beamline.s4_optical_element import S4OpticalElement
from shadow4.beamline.s4_beamline_element import S4BeamlineElement


class S4Grating(GratingVLS, S4OpticalElement):
    def __init__(self,
                 # inputs related tosyned
                 name="Undefined",
                 surface_shape=None,
                 boundary_shape=None,
                 ruling=800e3,
                 ruling_coeff_linear=0.0,
                 ruling_coeff_quadratic=0.0,
                 ruling_coeff_cubic=0.0,
                 ruling_coeff_quartic=0.0,
                 coating=None,
                 coating_thickness=None,
                 # inputs related autosetting
                 f_central=False,
                 f_phot_cent=0,
                 phot_cent=8000.0,
                 # inputs related to mirror reflectivity
                 f_reflec=0,  # reflectivity of surface: 0=no reflectivity, 1=full polarization
                 material_constants_library_flag=0,  # 0=xraylib, 1=dabax, 2=shadow preprocessor
                 file_refl="",
                 # inputs related to observation direction
                 order=0,
                 f_ruling=0,    # - for f_grating=1 - ruling type: (0) constant on X-Y plane (0)
                                # constant on mirror surface (1),
                                # holographic (2),
                                # fan type	(3),
                                # reserved (4),
                                # polynomial line density (5).
                 # inputs NOT USED ANYMORE....
                 # f_mono=0,      #- f_grating, f_central=1 - monochromator type:
                 #                # TGM / Seya(0)
                 #                # ERG(1),
                 #                # Constant Incidence Angle(2),
                 #                # Constant diffraction angle(3),
                 #                # Hunter(4)
                 # f_hunt=1,      #- for f_mono = 4: first(1) or second(2) grating.
                 # dist_fan=0.0,  # - f_ruling= 3: distance from grating center (cm).
                 ):
        GratingVLS.__init__(self,
                 name=name,
                 surface_shape=surface_shape,
                 boundary_shape=boundary_shape,
                 ruling=ruling,
                 ruling_coeff_linear=ruling_coeff_linear,
                 ruling_coeff_quadratic=ruling_coeff_quadratic,
                 ruling_coeff_cubic=ruling_coeff_cubic,
                 ruling_coeff_quartic=ruling_coeff_quartic,
                 coating=coating,
                 coating_thickness=coating_thickness,)


        self._f_central = f_central
        self._f_phot_cent = f_phot_cent
        self._phot_cent = phot_cent
        self._f_reflec = f_reflec
        self._material_constants_library_flag = material_constants_library_flag
        self._file_refl = file_refl
        self._order = order
        self._f_ruling = f_ruling

        self.congruence()

    def congruence(self):
        if not self._f_ruling in [0,1,5]:
            raise Exception("Not implemented grating with f_ruling=%d" % self._fruling)


class S4GratingElement(S4BeamlineElement):

    def __init__(self, optical_element=None, coordinates=None):
        super().__init__(optical_element if optical_element is not None else S4Grating(),
                         coordinates if coordinates is not None else ElementCoordinates())

        self.align_grating()

    def align_grating(self):

        oe = self.get_optical_element()
        coor = self.get_coordinates()


        if oe._f_central:
            if oe._f_phot_cent == 0:
                energy = oe._phot_cent
            else:
                energy = codata.h * codata.c / codata.e * 1e2 / (oe._phot_cent * 1e-8)

            raise Exception(NotImplementedError)

            angle_radial = 0.0
            angle_radial_out = 0.0

            coor.set_angles(angle_radial=angle_radial,
                            angle_radial_out=angle_radial_out,
                            angle_azimuthal=0.0)
        else:
            print("Info: nothing to align: f_central=0")

        print(coor.info())

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

        beam_in_crystal_frame_before_reflection = beam.duplicate()
        if not isinstance(soe, GratingVLS):  # undefined
            raise Exception("Undefined Grating")
        else:
            beam_mirr, normal = self.apply_grating_diffraction(beam)  # warning, beam is also changed!!

        #
        # apply mirror boundaries
        #
        beam_mirr.apply_boundaries_syned(soe.get_boundary_shape(), flag_lost_value=flag_lost_value)

        ########################################################################################
        #
        # TODO" apply grating reflectivity/efficiency
        #
        # beam_mirr.apply_complex_reflectivities(complex_reflectivity_S, complex_reflectivity_P)
        ########################################################################################


        #
        # from element reference system to image plane
        #

        beam_out = beam_mirr.duplicate()
        beam_out.change_to_image_reference_system(theta_grazing2, q)

        # plot results
        if False:
            if scan_type == 0:
                pass
            else:
                deviations = beam_out.get_column(6)
                intensityS = beam_out.get_column(24)
                intensityP = beam_out.get_column(25)

            from srxraylib.plot.gol import plot
            plot(1e6 * deviations, intensityS,
                 1e6 * deviations, intensityP,
                 xtitle="deviation angle [urad]",
                 ytitle="Reflectivity",
                 legend=["Sigma-polarization", "Pi-polarization"],
                 linestyle=['', ''],
                 marker=['+', '.'])

        return beam_out, beam_mirr

    def apply_grating_diffraction(self, beam):  # to be implemented in the children classes

        oe = self.get_optical_element()
        ssi = oe.get_surface_shape_instance()
        ccc = oe.get_optical_surface_instance()

        if oe._f_ruling == 5:
            ruling = [oe._ruling,
                        oe._ruling_coeff_linear,
                        oe._ruling_coeff_quadratic,
                        oe._ruling_coeff_cubic,
                        oe._ruling_coeff_quartic]
        else:
            ruling = [oe._ruling]

        if isinstance(ssi, Plane):
            beam_mirr, normal = ccc.apply_grating_diffraction_on_beam(
                beam,
                ruling=ruling,
                order=oe._order,
                f_ruling=oe._f_ruling)
        elif isinstance(ssi, Sphere):
            beam_mirr, normal = ccc.apply_grating_diffraction_on_beam(
                beam,
                ruling=ruling,
                order=oe._order,
                f_ruling=oe._f_ruling)
        else:
            raise NotImplementedError

        return beam_mirr, normal


if __name__ == "__main__":
    from shadow4.sources.source_geometrical.source_gaussian import SourceGaussian
    from shadow4.beam.beam import Beam

    from syned.beamline.shape import Plane
    #
    # source
    #
    src = SourceGaussian.initialize_from_keywords(
                 number_of_rays=10000,
                 sigmaX=1.0e-6,
                 sigmaY=0.0,
                 sigmaZ=1.0e-6,
                 sigmaXprime=0.0002,
                 sigmaZprime=0.0002,
                 real_space_center=[0.0,0.0,0.0],
                 direction_space_center=[0.0,0.0]
                                 )
    beam = Beam()

    beam.genSource(src)
    beam.set_photon_energy_eV(1000.0)

    print(beam.info())

    # plotxy(Beam3.initialize_from_shadow4_beam(beam),1,3,nbins=100,title="SOURCE")

    #
    # grating
    #
    g = S4Grating(
        name = "my_grating",
        surface_shape =  Plane(), # SurfaceShape(),
        boundary_shape = None, # BoundaryShape(),
        ruling = 600000.0,
        ruling_coeff_linear = 260818.35944225,
        ruling_coeff_quadratic = 260818.35944225,
        ruling_coeff_cubic = 13648.21037618,
        ruling_coeff_quartic = 0.0,
        coating = None,
        coating_thickness = None,
        f_central=False,
        f_phot_cent=0,
        phot_cent=8000.0,
        material_constants_library_flag=0,  # 0=xraylib, 1=dabax, 2=shadow preprocessor
        file_refl="",
        order=0,
        )

    coordinates_syned = ElementCoordinates(p = 10.0,
                                           q = 6.0,
                                           angle_radial = 88.840655 * numpy.pi / 180,
                                           angle_radial_out= 87.588577 * numpy.pi / 180,
                                           angle_azimuthal = 0.0)



    ge = S4GratingElement(optical_element=g, coordinates=coordinates_syned)

    print(ge.info())

    beam_out = ge.trace_beam(beam)
