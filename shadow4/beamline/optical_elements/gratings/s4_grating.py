import numpy
import scipy.constants as codata

from syned.beamline.optical_elements.gratings.grating import GratingVLS
from syned.beamline.shape import Plane, Sphere, Conic, Toroid, Paraboloid, Hyperboloid, Ellipsoid, NumericalMesh
from syned.beamline.element_coordinates import ElementCoordinates

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.s4_optical_element_decorators import S4OpticalElementDecorator
from shadow4.beamline.s4_beamline_element import S4BeamlineElement
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements


class S4Grating(GratingVLS, S4OpticalElementDecorator):
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
                 coating=None,            # not used, passed to syned
                 coating_thickness=None,  # not used, passed to syned
                 # inputs related to observation direction
                 order=0,
                 f_ruling=0,    # (0) constant on X-Y plane (0)
                                # (1) polynomial line density (5 in shadow5).
                 # inputs NOT USED ANYMORE....
                 # # inputs related autosetting
                 # f_central=False,
                 # f_phot_cent=0,
                 # phot_cent=8000.0,
                 # inputs related to mirror reflectivity
                 # f_reflec=0,  # reflectivity of surface: 0=no reflectivity, 1=full polarization
                 # material_constants_library_flag=0,  # 0=xraylib, 1=dabax, 2=shadow preprocessor
                 # file_refl="",
                 # f_mono=0,      #- f_grating, f_central=1 - monochromator type:
                 #                # TGM / Seya(0)
                 #                # ERG(1),
                 #                # Constant Incidence Angle(2),
                 #                # Constant diffraction angle(3),
                 #                # Hunter(4)
                 # f_hunt=1,      #- for f_mono = 4: first(1) or second(2) grating.
                 # dist_fan=0.0,  # - f_ruling= 3: distance from grating center (cm).
                 ):
        if f_ruling == 0:
            ruling_coeff_linear = 0
            ruling_coeff_quadratic = 0
            ruling_coeff_cubic = 0
            ruling_coeff_quartic = 0

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


        # self._f_central = f_central
        # self._f_phot_cent = f_phot_cent
        # self._phot_cent = phot_cent
        # self._f_reflec = f_reflec
        # self._material_constants_library_flag = material_constants_library_flag
        # self._file_refl = file_refl
        self._order = order
        self._f_ruling = f_ruling

        self.congruence()

    def congruence(self):
        if not self._f_ruling in [0,1,5]:
            raise Exception("Not implemented grating with f_ruling=%d" % self._fruling)

class S4GratingElement(S4BeamlineElement):
    def __init__(self,
                 optical_element : S4Grating = None,
                 coordinates : ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4Grating(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         movements=movements,
                         input_beam=input_beam)

    def trace_beam(self, **params):
        flag_lost_value = params.get("flag_lost_value", -1)

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
        # reflect beam in the mirror surface
        #
        soe = self.get_optical_element()

        beam_in_crystal_frame_before_reflection = input_beam.duplicate()
        footprint, normal = self.apply_grating_diffraction(input_beam)  # warning, beam is also changed!!


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

        ########################################################################################
        #
        # TODO" apply grating reflectivity/efficiency
        #
        # beam_mirr.apply_complex_reflectivities(complex_reflectivity_S, complex_reflectivity_P)
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
                 linestyle=['', ''],
                 marker=['+', '.'])

        return output_beam, footprint

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

        # todo: remove if's the interface is the same for Conic, Toroid and NumericalMesh
        if isinstance(ssi, Plane) or isinstance(ssi, Sphere) or isinstance(ssi, Paraboloid) \
                or isinstance(ssi, Ellipsoid) or isinstance(ssi, Hyperboloid) or isinstance(ssi, Conic):
            beam_mirr, normal = ccc.apply_grating_diffraction_on_beam(
                beam,
                ruling=ruling,
                order=oe._order,
                f_ruling=oe._f_ruling)
        # elif isinstance(ssi, Sphere):
        #     beam_mirr, normal = ccc.apply_grating_diffraction_on_beam(
        #         beam,
        #         ruling=ruling,
        #         order=oe._order,
        #         f_ruling=oe._f_ruling)
        # elif isinstance(ssi, Conic):
        #     beam_mirr, normal = ccc.apply_grating_diffraction_on_beam(
        #         beam,
        #         ruling=ruling,
        #         order=oe._order,
        #         f_ruling=oe._f_ruling)
        elif isinstance(ssi, Toroid):
            beam_mirr, normal = ccc.apply_grating_diffraction_on_beam(
                beam,
                ruling=ruling,
                order=oe._order,
                f_ruling=oe._f_ruling)
        elif isinstance(ssi, NumericalMesh):
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
    from shadow4.beam.s4_beam import S4Beam

    from syned.beamline.shape import Plane
    #
    # source
    #
    src = SourceGaussian(
                 nrays=10000,
                 sigmaX=1.0e-6,
                 sigmaY=0.0,
                 sigmaZ=1.0e-6,
                 sigmaXprime=0.0002,
                 sigmaZprime=0.0002,
                 real_space_center=[0.0,0.0,0.0],
                 direction_space_center=[0.0,0.0]
                                 )
    beam = S4Beam()

    beam.generate_source(src)
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
        order=0,
        )

    coordinates_syned = ElementCoordinates(p = 10.0,
                                           q = 6.0,
                                           angle_radial = 88.840655 * numpy.pi / 180,
                                           angle_radial_out= 87.588577 * numpy.pi / 180,
                                           angle_azimuthal = 0.0)



    ge = S4GratingElement(optical_element=g, coordinates=coordinates_syned, input_beam=beam)

    print(ge.info())

    beam_out = ge.trace_beam()
