import numpy
from shadow4.sources.wiggler.s4_wiggler_light_source import S4WigglerLightSource
from shadow4.beam.s4_beam import S4Beam

from syned.beamline.shape import Rectangle
from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
from syned.beamline.element_coordinates import ElementCoordinates
from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement


class S4WigglerOptimizedLightSource(S4WigglerLightSource):
    """
    Defines an "optimized" wiggler light source. "Optimized" means that a reject mechanist is in place to avoid
    storing the rays that will not go through a slit.

    Parameters
    ----------
    name : str, optional
        The name of the light source.
    electron_beam : instance of ElectronBeam
        The electron beam parameters.
    magnetic_structure : instance of S4BendingMagnet
        The shadow4 bending magnet magnetic structure.
    nrays : int, optional
        The number of rays.
    seed : int, optional
        The Monte Carlo seed.
    optim_method: int, optional
        0=No optimization, 1=Optimization by defined slit, 2=Optimization by entered S4BeamlineElement.
    optim_max_iterations: int, optional
        The maximum number of iterations (source creation).
    optim_slit_distance: float, optional
        The position of the slit in m (for optim_method=1).
    optim_slit_center_x: float, optional
        The X coordinate of the slit center in m (for optim_method=1).
    optim_slit_center_z: float, optional
        The Z coordinate of the slit center in m (for optim_method=1).
    optim_slit_gap_x: float, optional
        The slit gap in X direction, in m (for optim_method=1).
    optim_slit_gap_z: float, optional
        The slit gap in Z direction, in m (for optim_method=1).
    optim_beamline_element: None or instance of S4BeamlineElement, optional
        For optim_method=2, the beamline element.
    """
    def __init__(self,
                 name="Undefined",
                 electron_beam=None,
                 magnetic_structure=None,
                 nrays=5000,
                 seed=12345,
                 optim_method=0,    # 0=None, 1=Rect slit, 2=beamline_element
                 optim_slit_d=1.0,
                 optim_slit_center_x=0.0,
                 optim_slit_center_z=0.0,
                 optim_slit_gap_x=1.0,
                 optim_slit_gap_z=1.0,
                 optim_max_iterations=10,
                 optim_beamline_element=None,
                 ):
        super().__init__(name=name,
                         electron_beam=electron_beam,
                         magnetic_structure=magnetic_structure,
                         nrays=nrays,
                         seed=seed,
                         )

        self._optim_method           = optim_method
        self._optim_slit_d           = optim_slit_d
        self._optim_slit_center_x    = optim_slit_center_x
        self._optim_slit_center_z    = optim_slit_center_z
        self._optim_slit_gap_x       = optim_slit_gap_x
        self._optim_slit_gap_z       = optim_slit_gap_z
        self._optim_max_iterations   = optim_max_iterations
        self._optim_beamline_element = optim_beamline_element

    ############################################################################
    #
    ############################################################################
    def get_beam(self, F_COHER=0, psi_interval_in_units_one_over_gamma=None):
        """
        Creates the beam as emitted by the wiggler.

        Parameters
        ----------
        F_COHER : int, optional
            A flag to indicate that the phase for the s-component is set to zero (coherent_beam=1) or is random for incoherent.
        psi_interval_in_units_one_over_gamma : None or float, optional
            The interval of psi*gamma for sampling rays.

        Returns
        -------
        instance of S4Beam
        """

        if self.get_seed() != 0:
            numpy.random.seed(self.get_seed())

        beam = S4Beam.initialize_from_array(
            self._S4WigglerLightSource__calculate_rays( # TODO: find a nicer solution?
                user_unit_to_m                       = 1.0,
                F_COHER                              = F_COHER,
                psi_interval_in_units_one_over_gamma = psi_interval_in_units_one_over_gamma,
                psi_interval_number_of_points        = self.get_magnetic_structure()._psi_interval_number_of_points,
                )
            )

        if self._optim_method == 0:
            return beam
        elif self._optim_method == 1:
            boundary_shape = Rectangle(x_left=self._optim_slit_center_x - 0.5 * self._optim_slit_gap_x,
                                       x_right=self._optim_slit_center_x + 0.5 * self._optim_slit_gap_x,
                                       y_bottom=self._optim_slit_center_z - 0.5 * self._optim_slit_gap_z,
                                       y_top=self._optim_slit_center_z + 0.5 * self._optim_slit_gap_z,
                                       )
            optical_element = S4Screen(name='Generic Beam Screen/Slit/Stopper/Attenuator',
                                       boundary_shape=boundary_shape,
                                       i_abs=0,  # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
                                       i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)
            coordinates = ElementCoordinates(p=self._optim_slit_d, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
            beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

        else:
            beamline_element = self._optim_beamline_element

        N0 = beam.N

        for iter in range(self._optim_max_iterations):
            print("======================================== iteration %d ========================================" % iter)

            beamline_element.set_input_beam(beam)

            beam_on_slit, footprint = beamline_element.trace_beam()

            print("=============================== on slit: N, Ngood: ", beam_on_slit.N, beam_on_slit.Ngood)
            #
            # if beam_on_slit.Ngood == 0:
            #     raise ("Error: optimization slit does nor receive any rays")

            beam.rays[:, 9] = beam_on_slit.rays[:, 9]

            beam.clean_lost_rays()

            if iter == 0:
                beam_accumulated = beam.duplicate()
            else:
                beam_accumulated.append_beam(beam)

            print("=============================== on accumulated beam: N, Ngood: ", beam_accumulated.N, beam_accumulated.Ngood)

            if beam_accumulated.Ngood >= N0:
                print("\n=============================== Reached Ngood %d >= %d" % (beam_accumulated.Ngood, N0))
                break

            beam = S4Beam.initialize_from_array(
                self._S4WigglerLightSource__calculate_rays( # TODO: find a nicer solution?
                    user_unit_to_m                       = 1.0,
                    F_COHER                              = F_COHER,
                    psi_interval_in_units_one_over_gamma = psi_interval_in_units_one_over_gamma,
                    psi_interval_number_of_points        = self.get_magnetic_structure()._psi_interval_number_of_points,
                    )
                )

        if iter == (self._optim_max_iterations - 1):
            print("\n=============================== Reached maximum number of iterations %d:  Ngood=%d < %d" % \
                  (self._optim_max_iterations, beam_accumulated.Ngood, N0))

        return beam_accumulated

    def to_python_code(self, **kwargs):
        """
        Returns the python code for calculating the wiggler source.

        Returns
        -------
        str
            The python code.
        """
        script = ''
        try:
            script += self.get_electron_beam().to_python_code()
        except:
            script += "\n\n#Error retrieving electron_beam code"

        try:
            script += self.get_magnetic_structure().to_python_code()
        except:
            script += "\n\n#Error retrieving magnetic structure code"

        script += "\n\n\n#light source\nfrom shadow4.sources.wiggler.s4_wiggler_optimized_light_source import S4WigglerOptimizedLightSource"
        script += "\nlight_source = S4WigglerOptimizedLightSource(name='%s', electron_beam=electron_beam, magnetic_structure=source," % self.get_name()
        script += "\n    nrays=%d, seed=%s," % (self.get_nrays(), self.get_seed())
        script += "\n    optim_method=1,  # 0=None, 1=Rect slit, 2=beamline_element"
        script += "\n    optim_max_iterations=%d," % self._optim_max_iterations
        script += "\n    optim_slit_d=%g, optim_slit_center_x=%g, optim_slit_center_z=%g, " % (self._optim_slit_d, self._optim_slit_center_x, self._optim_slit_center_z)
        script += "\n    optim_slit_gap_x=%g, optim_slit_gap_z=%g)" % (self._optim_slit_gap_x, self._optim_slit_gap_z)
        #
        # script += "\n\n\n#beamline\nfrom shadow4.beamline.s4_beamline import S4Beamline"
        # script += "\nbeamline = S4Beamline(light_source=light_source)"
        script += "\nbeam = light_source.get_beam()"
        return script


if __name__ == "__main__":

    from srxraylib.plot.gol import plot_scatter

    # electron beam
    from shadow4.sources.s4_electron_beam import S4ElectronBeam

    electron_beam = S4ElectronBeam(energy_in_GeV=6, energy_spread=0.001, current=0.2)
    electron_beam.set_sigmas_all(sigma_x=2.377e-05, sigma_y=2.472e-05, sigma_xp=3.58e-06, sigma_yp=3.04e-06)

    # magnetic structure
    from shadow4.sources.wiggler.s4_wiggler import S4Wiggler

    source = S4Wiggler(
        magnetic_field_periodic=0,  # 0=external, 1=periodic
        file_with_magnetic_field="/nobackup/gurb1/srio/Oasys/SW_BM18_Joel.txt",  # used only if magnetic_field_periodic=0
        K_vertical=10.0,       # syned Wiggler pars: used only if magnetic_field_periodic=1
        period_length=0.1,     # syned Wiggler pars: used only if magnetic_field_periodic=1
        number_of_periods=10,  # syned Wiggler pars: used only if magnetic_field_periodic=1
        emin=100.0,     # Photon energy scan from energy (in eV)
        emax=200000.0,  # Photon energy scan to energy (in eV)
        ng_e=51,        # Photon energy scan number of points
        ng_j=501,       # Number of points in electron trajectory (per period) for internal calculation only
        flag_interpolation=2,  #
        psi_interval_number_of_points=101,
        flag_emittance=1,       # Use emittance (0=No, 1=Yes)
        shift_x_flag=4,         # 0="No shift", 1="Half excursion", 2="Minimum", 3="Maximum", 4="Value at zero", 5="User value"
        shift_x_value=0.001,    # used only if shift_x_flag=5
        shift_betax_flag=4,     # 0="No shift", 1="Half excursion", 2="Minimum", 3="Maximum", 4="Value at zero", 5="User value"
        shift_betax_value=0.0,  # used only if shift_betax_flag=5
    )

    # light source
    from shadow4.sources.wiggler.s4_wiggler_light_source import S4WigglerLightSource

    light_source = S4WigglerOptimizedLightSource(name='wiggler',
                                        electron_beam=electron_beam,
                                        magnetic_structure=source,
                                        nrays=20000,
                                        seed=0, # 5676561,
                                        optim_method=1,  # 0=None, 1=Rect slit, 2=beamline_element
                                        optim_slit_d=30.0,
                                        optim_slit_center_x=0.0,
                                        optim_slit_center_z=0.0,
                                        optim_slit_gap_x=0.01,
                                        optim_slit_gap_z=0.01,
                                        optim_beamline_element=None,
                                        )


    beam = light_source.get_beam()
    print("beam N, Ngood, Nbad: ", beam.N, beam.Ngood, beam.Nbad)

    # # test plot
    # from srxraylib.plot.gol import plot_scatter
    # rays = beam.get_rays()
    # plot_scatter(1e6 * rays[:, 0], 1e6 * rays[:, 2], title='(X,Z) in microns',   xrange=[-1800,100], yrange=[-100,100], show=0)
    # plot_scatter(1e6 * rays[:, 3], 1e6 * rays[:, 5], title='(Xp,Zp) in microns', xrange=[-10000,8000], yrange=[-300,300])
    #
    print(light_source.to_python_code())
