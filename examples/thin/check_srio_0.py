import numpy as np
from dabax.dabax_xraylib import DabaxXraylib
from shadow4.beamline.s4_beamline import S4Beamline

beamline = S4Beamline()

#
#
#
from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical

light_source = SourceGeometrical(name='Geometrical Source', nrays=5000, seed=5676561)
light_source.set_spatial_type_point()
light_source.set_depth_distribution_off()
light_source.set_angular_distribution_gaussian(sigdix=1e-05, sigdiz=1e-05)
light_source.set_energy_distribution_uniform(value_min=9950, value_max=10050, unit='eV')
light_source.set_polarization(polarization_degree=1, phase_diff=0, coherent_beam=0)
beam = light_source.get_beam()

beamline.set_light_source(light_source)

# optical element number XX
from syned.beamline.shape import Circle

boundary_shape = Circle(radius=0.0002, x_center=0, y_center=0)
from shadow4.beamline.optical_elements.refractors.s4_lens import S4Lens

optical_element = S4Lens(name='Refractive Lens (mesh)',
                         boundary_shape=boundary_shape,  # syned stuff, replaces "diameter" in the shadow3 append_lens
                         material='Be',  # the material for ri_calculation_mode > 1
                         density=1.85,  # the density for ri_calculation_mode > 1
                         thickness=0.0,
                         # syned stuff, lens thickness [m] (distance between the two interfaces at the center of the lenses)
                         surface_shape=0,  # 0=plane, 1=sphere, 2=parabola, 3=conic coefficients
                         convex_to_the_beam=0,
                         # for surface_shape: convexity of the first interface exposed to the beam 0=No, 1=Yes
                         cylinder_angle=0,  # for surface_shape: 0=not cylindricaL, 1=meridional 2=sagittal
                         ri_calculation_mode=3,
                         # source of refr indices and absorp coeff 0=User, 1=prerefl file, 2=xraylib, 3=dabax
                         prerefl_file='lens_si_O2.dat',
                         # for ri_calculation_mode=0: file name (from prerefl) to get the refraction index.
                         refraction_index=1.4682,  # for ri_calculation_mode=1: n (real)
                         attenuation_coefficient=0,  # for ri_calculation_mode=1: mu in m^-1 (real)
                         dabax=DabaxXraylib(file_f1f2="f1f2_Windt.dat", file_CrossSec="CrossSec_EPDL97.dat"),
                         # if using dabax (ri_calculation_mode=3), instance of DabaxXraylib() (use None for default)
                         radius=1.538e-06,
                         # for surface_shape=(1,2): lens radius [m] (for spherical, or radius at the tip for paraboloid)
                         conic_coefficients1=None,  # for surface_shape = 3: the conic coefficients for interface 1
                         conic_coefficients2=None,  # for surface_shape = 3: the conic coefficients for interface 2
                         flag_add_mesh_surface_entrance=1,
                         # 0=No, 1=Yes: add a numerical mesh on top of the entrance interface native shape
                         flag_add_mesh_surface_exit=0,
                         # 0=No, 1=Yes: add a numerical mesh on top of the exit interface native shape
                         mesh_surface_entrance_h5file='/home/srio/Oasys2/lens_interface_1b.h5',
                         # for flag_add_mesh_surface_entrance=1: h5 file with the entrance mesh
                         mesh_surface_exit_h5file='/home/srio/Oasys2/lens_interface_2.h5',
                         # for flag_add_mesh_surface_exit=1: h5 file with the exit mesh
                         )

import numpy
from syned.beamline.element_coordinates import ElementCoordinates

coordinates = ElementCoordinates(p=10, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
movements = None
from shadow4.beamline.optical_elements.refractors.s4_lens import S4LensElement

beamline_element = S4LensElement(optical_element=optical_element, coordinates=coordinates, movements=movements,
                                 input_beam=beam)

beam, footprint = beamline_element.trace_beam()

beamline.append_beamline_element(beamline_element)

# test plot
if True:
    from srxraylib.plot.gol import plot_scatter

    beam.retrace(27.4327)
    plot_scatter(1e6 * beam.get_column(1, nolost=1), 1e6 * beam.get_column(3, nolost=1), title='(X,Z) in microns I= %d' % (beam.get_intensity(nolost=1)))