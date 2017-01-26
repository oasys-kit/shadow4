#
# examples for SourceUndulator to be used in ShadowOui
#

def Setting(hello):
    return hello

if __name__ == "__main__":

    import Shadow
    import numpy
    import os
    from SourceUndulatorInputOutput import plot_undul_cdf,plot_undul_phot
    from SourceUndulator import SourceUndulator



    NONE_SPECIFIED = "NONE SPECIFIED"

    plot_graph = Setting(0)

    number_of_rays=Setting(5000)
    seed=Setting(5676561)
    set_to_resonance_combo=Setting(1)

    e_min=Setting(7990)
    e_max=Setting(8010)
    harmonic_number=Setting(1)


    # optimize_source_combo=Setting(0)
    # max_number_of_rejected_rays = Setting(10000000)
    # slit_distance = Setting(1000.0)
    # min_x = Setting(-1.0)
    # max_x = Setting(1.0)
    # min_z = Setting(-1.0)
    # max_z = Setting(1.0)
    # file_with_phase_space_volume = Setting(NONE_SPECIFIED)

    energy=Setting(6.04)
    electron_current = Setting(200)
    use_emittances_combo=Setting(1)
    sigma_x=Setting(0.04)
    sigma_z=Setting(0.001)
    # emittance_x=Setting(3.8E-7)
    # emittance_z=Setting(3.8E-9)
    # New
    sigma_xprime = Setting(10e-6)
    sigma_zprime = Setting(4e-6)

    number_of_periods=Setting(50)
    k_value=Setting(0.25)
    id_period=Setting(0.032)

    # New precalculations

    number_of_points_slit_theta = Setting(51)
    number_of_points_slit_phi = Setting(11)
    number_of_points_photon_energy = Setting(11)
    number_of_points_per_trajectory_period = Setting(20)
    maximum_radiation_semiaperture_mrad = Setting(0.1)

            #
            # (traj, pars) = srfunc.wiggler_trajectory(b_from=self.type_combo,
            #                                          inData=inData,
            #                                          nPer=self.number_of_periods,
            #                                          nTrajPoints=501,
            #                                          ener_gev=self.energy,
            #                                          per=self.id_period,
            #                                          kValue=self.k_value,
            #                                          trajFile=congruence.checkFileName("tmp.traj"),
            #                                          shift_x_flag=self.shift_x_flag,
            #                                          shift_x_value=self.shift_x_value,
            #                                          shift_betax_flag=self.shift_betax_flag,
            #                                          shift_betax_value=self.shift_betax_value)






    u = SourceUndulator()

    u.set_from_keywords(E_ENERGY=energy,INTENSITY=electron_current,FLAG_EMITTANCE=use_emittances_combo,
                        SX=sigma_x,SZ=sigma_z,SXP=sigma_xprime,SZP=sigma_xprime,
                        LAMBDAU=id_period,NPERIODS=number_of_periods,K=k_value,
                        EMIN=e_min,EMAX=e_max,
                        NG_E=number_of_points_photon_energy,MAXANGLE=maximum_radiation_semiaperture_mrad,
                        NG_T=number_of_points_slit_theta,NG_P=number_of_points_slit_phi,N_J=number_of_points_per_trajectory_period,
                        SEED=seed,NRAYS=number_of_rays)


    if harmonic_number != 0:
        u.set_energy_monochromatic_at_resonance(harmonic_number=harmonic_number)
        angle = u.get_resonance_central_cone(harmonic_number)

    print(u.info(debug=True))

    os.system("rm begin.dat start.00 uphot.dat")

    method = 1 # 0=direct, 1=get intermediate radiation
    if method == 0:
        # direct calculation
        beam = u.calculate_shadow3_beam(code_undul_phot='internal',dump_undul_phot_file=True,dump_start_files=True)
        plot_undul_phot("uphot.dat",do_plot_intensity=True,do_plot_polarization=False)
    else:
        # via intermediate radiation (without writing uphot.dat)

        dict_radiation = u.calculate_radiation(code_undul_phot='internal')
        # plot_undul_phot(dict_radiation,do_plot_intensity=True,do_plot_polarization=False)

        tmp = u.calculate_cdf(code_undul_phot='intrernal',
                                 use_existing_undul_phot_output=dict_radiation,
                                 dump_undul_phot_file=False)

        # initialize shadow3 source (oe0) and beam
        oe0 = u.get_shadow3_source_object()

        oe0.write("start.00")

        beam = Shadow.Beam()
        beam.genSource(oe0)

        oe0.write("end.00")


    beam.write("begin.dat")


    tkt = Shadow.ShadowTools.plotxy("begin.dat",4,6,nbins=150)
    print("FWHM: h histo: %f, v histo: %f urad"%(1e6*tkt['fwhm_h'],1e6*tkt['fwhm_h']))

    if harmonic_number!= 0:
        print("FWHM: calculated: %f (or %f) urad"%(1e6*2.35*0.69*u.get_resonance_central_cone(harmonic_number),
              1e6*2.35*0.69*numpy.sqrt(u.get_resonance_wavelength(harmonic_number)/(u.NPERIODS*u.LAMBDAU))))


    # plot_undul_cdf("xshundul.sha")
