import numpy as np
from dabax.dabax_xraylib import DabaxXraylib
from shadow4.beamline.s4_beamline import S4Beamline

beamline = S4Beamline()

# electron beam
from shadow4.sources.s4_electron_beam import S4ElectronBeam

electron_beam = S4ElectronBeam(energy_in_GeV=6, energy_spread=0.00094, current=0.2)
electron_beam.set_moments_horizontal(9.16257e-10, -7.2086e-13, 1.93063e-11)
electron_beam.set_moments_vertical(2.66211e-11, 1.298e-14, 3.75642e-12)
electron_beam.set_dispersion_all(-3.5e-06, 0, 0, -0)

# magnetic structure
from shadow4.sources.undulator.s4_undulator_gaussian import S4UndulatorGaussian

source = S4UndulatorGaussian(
    period_length=0.027,  # syned Undulator parameter (length in m)
    number_of_periods=52.0,  # syned Undulator parameter
    photon_energy=9000.0,  # Photon energy (in eV)
    delta_e=0.0,  # Photon energy width (in eV)
    ng_e=1,  # Photon energy scan number of points
    flag_emittance=1,  # when sampling rays: Use emittance (0=No, 1=Yes)
    flag_energy_spread=1,  # when sampling rays: Use e- energy spread (0=No, 1=Yes)
    harmonic_number=3,  # harmonic number
    flag_autoset_flux_central_cone=1,  # value to set the flux peak
    flux_central_cone=773567252092211.8,  # value to set the flux peak
)

# light source
from shadow4.sources.undulator.s4_undulator_gaussian_light_source import S4UndulatorGaussianLightSource

light_source = S4UndulatorGaussianLightSource(name='Undulator Gaussian', electron_beam=electron_beam,
                                              magnetic_structure=source, nrays=500000, seed=5676561)
beam = light_source.get_beam()

beamline.set_light_source(light_source)

# optical element number XX
from syned.beamline.shape import Rectangle

boundary_shape = Rectangle(x_left=-0.000585, x_right=0.000585, y_bottom=-0.000585, y_top=0.000585)

from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen

optical_element = S4Screen(name='10 µrad (z=117 m)', boundary_shape=boundary_shape,
                           i_abs=0,  # attenuation: 0=No, 1=prerefl file, 2=xraylib, 3=dabax
                           i_stop=0,  # 0=slit or aperture, 1=beam stop
                           thick=0,  # for i_abs>0
                           file_abs='<specify file name>',  # for i_abs=1
                           material='Au', density=19.3,  # for i_abs=2,3
                           dabax=None,  # if using dabax (i_abs=3), instance of DabaxXraylib() (use None for default)
                           )

from syned.beamline.element_coordinates import ElementCoordinates

coordinates = ElementCoordinates(p=117, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement

beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

beam, footprint = beamline_element.trace_beam()

beamline.append_beamline_element(beamline_element)

# optical element number XX
from syned.beamline.shape import Ellipse

boundary_shape = Ellipse(a_axis_min=-0.00022, a_axis_max=0.00022, b_axis_min=-0.00022, b_axis_max=0.00022)

from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen

optical_element = S4Screen(name='pinhole (D=440µm)', boundary_shape=boundary_shape,
                           i_abs=0,  # attenuation: 0=No, 1=prerefl file, 2=xraylib, 3=dabax
                           i_stop=0,  # 0=slit or aperture, 1=beam stop
                           thick=0,  # for i_abs>0
                           file_abs='<specify file name>',  # for i_abs=1
                           material='Au', density=19.3,  # for i_abs=2,3
                           dabax=None,  # if using dabax (i_abs=3), instance of DabaxXraylib() (use None for default)
                           )

from syned.beamline.element_coordinates import ElementCoordinates

coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement

beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

beam, footprint = beamline_element.trace_beam()

beamline.append_beamline_element(beamline_element)

# optical element number XX
from syned.beamline.shape import Circle

boundary_shape = Circle(radius=0.00022, x_center=0, y_center=0)
from shadow4.beamline.optical_elements.refractors.s4_lens import S4Lens

optical_element = S4Lens(name='Be 50µm (x1 lens)',
                         boundary_shape=boundary_shape,  # syned stuff, replaces "diameter" in the shadow3 append_lens
                         material='Be',  # the material for ri_calculation_mode > 1
                         density=1.484,  # the density for ri_calculation_mode > 1
                         thickness=2.9999999999999997e-05,
                         # syned stuff, lens thickness [m] (distance between the two interfaces at the center of the lenses)
                         surface_shape=2,  # 0=plane, 1=sphere, 2=parabola, 3=conic coefficients
                         convex_to_the_beam=0,
                         # for surface_shape: convexity of the first interface exposed to the beam 0=No, 1=Yes
                         cylinder_angle=0,  # for surface_shape: 0=not cylindricaL, 1=meridional 2=sagittal
                         ri_calculation_mode=1,
                         # source of refr indices and absorp coeff 0=User, 1=prerefl file, 2=xraylib, 3=dabax
                         prerefl_file='Be.dat',
                         # for ri_calculation_mode=0: file name (from prerefl) to get the refraction index.
                         refraction_index=1,  # for ri_calculation_mode=1: n (real)
                         attenuation_coefficient=0,  # for ri_calculation_mode=1: mu in m^-1 (real)
                         dabax=None,
                         # if using dabax (ri_calculation_mode=3), instance of DabaxXraylib() (use None for default)
                         radius=5e-05,
                         # for surface_shape=(1,2): lens radius [m] (for spherical, or radius at the tip for paraboloid)
                         conic_coefficients1=None,  # for surface_shape = 3: the conic coefficients for interface 1
                         conic_coefficients2=None,  # for surface_shape = 3: the conic coefficients for interface 2
                         flag_add_mesh_surface_entrance=0,
                         # 0=No, 1=Yes: add a numerical mesh on top of the entrance interface native shape
                         flag_add_mesh_surface_exit=0,
                         # 0=No, 1=Yes: add a numerical mesh on top of the exit interface native shape
                         mesh_surface_entrance_h5file='<none>.hdf5',
                         # for flag_add_mesh_surface_entrance=1: h5 file with the entrance mesh
                         mesh_surface_exit_h5file='<none>.hdf5',
                         # for flag_add_mesh_surface_exit=1: h5 file with the exit mesh
                         )

import numpy
from syned.beamline.element_coordinates import ElementCoordinates

coordinates = ElementCoordinates(p=0, q=0.002, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
movements = None
from shadow4.beamline.optical_elements.refractors.s4_lens import S4LensElement

beamline_element = S4LensElement(optical_element=optical_element, coordinates=coordinates, movements=movements,
                                 input_beam=beam)

beam, footprint = beamline_element.trace_beam()

beamline.append_beamline_element(beamline_element)

# test plot
if True:
    from srxraylib.plot.gol import plot_scatter

    # plot_scatter(beam.get_photon_energy_eV(nolost=1), beam.get_column(23, nolost=1), title='(Intensity,Photon Energy)',
    #              plot_histograms=0)
    plot_scatter(1e6 * beam.get_column(1, nolost=1), 1e6 * beam.get_column(3, nolost=1),
                 title='(X,Z) in microns I= %d' % (beam.get_intensity(nolost=1)))