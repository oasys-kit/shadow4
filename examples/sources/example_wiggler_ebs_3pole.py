# electron beam
from shadow4.sources.s4_electron_beam import S4ElectronBeam

electron_beam = S4ElectronBeam(energy_in_GeV=6, energy_spread=0.001, current=0.2)
electron_beam.set_sigmas_all(sigma_x=2.377e-05, sigma_y=2.472e-05, sigma_xp=3.58e-06, sigma_yp=3.05e-06)

# magnetic structure
from shadow4.sources.wiggler.s4_wiggler import S4Wiggler

source = S4Wiggler(
    magnetic_field_periodic=0,  # 0=external, 1=periodic
    file_with_magnetic_field="https://raw.githubusercontent.com/oasys-kit/ShadowOui-Tutorial/master/SOS-WORKSHOP/EBS-WIGGLERS/SW_3Pcut.txt",
    # used only if magnetic_field_periodic=0
    K_vertical=10.0,  # syned Wiggler pars: used only if magnetic_field_periodic=1
    period_length=0.1,  # syned Wiggler pars: used only if magnetic_field_periodic=1
    number_of_periods=10,  # syned Wiggler pars: used only if magnetic_field_periodic=1
    emin=20000.0,  # Photon energy scan from energy (in eV)
    emax=20040.0,  # Photon energy scan to energy (in eV)
    ng_e=101,  # Photon energy scan number of points for spectrum and internal calculation
    ng_j=501,  # Number of points in electron trajectory (per period) for internal calculation only
    epsi_dx=-0.894,  # position y of waist X [m]
    epsi_dz=1.048,  # position y of waist Z [m]
    psi_interval_number_of_points=101,  # the number psi (vertical angle) points for internal calculation only
    flag_interpolation=1,  # Use interpolation to sample psi (0=No, 1=Yes)
    flag_emittance=1,  # Use emittance (0=No, 1=Yes)
    shift_x_flag=1,  # 0="No shift", 1="Half excursion", 2="Minimum", 3="Maximum", 4="Value at zero", 5="User value"
    shift_x_value=0.0,  # used only if shift_x_flag=5
    shift_betax_flag=0,  # 0="No shift", 1="Half excursion", 2="Minimum", 3="Maximum", 4="Value at zero", 5="User value"
    shift_betax_value=0.0,  # used only if shift_betax_flag=5
)

# light source
from shadow4.sources.wiggler.s4_wiggler_light_source import S4WigglerLightSource

light_source = S4WigglerLightSource(name='wiggler', electron_beam=electron_beam, magnetic_structure=source, nrays=20000,
                                    seed=5676561)
beam = light_source.get_beam()

# test plot
from srxraylib.plot.gol import plot_scatter

rays = beam.get_rays()
plot_scatter(1e6 * rays[:, 0], 1e6 * rays[:, 2], title='(X,Z) in microns')