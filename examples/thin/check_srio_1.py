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
light_source.set_angular_distribution_gaussian(sigdix=1e-05,sigdiz=1e-05)
light_source.set_energy_distribution_uniform(value_min=9950, value_max=10050, unit='eV')
light_source.set_polarization(polarization_degree=1, phase_diff=0, coherent_beam=0)
beam = light_source.get_beam()

beamline.set_light_source(light_source)

# optical element number XX
from syned.beamline.shape import Ellipse
boundary_shape = Ellipse(a_axis_min=-0.0002, a_axis_max=0.0002, b_axis_min=-0.0002, b_axis_max=0.0002)
from shadow4.beamline.optical_elements.phase_deflectors.s4_numerical_mesh_phase_deflector import S4NumericalMeshPhaseDeflector
optical_element = S4NumericalMeshPhaseDeflector(name='Phase Deflector', boundary_shape=boundary_shape,
    surface_data_file='/home/srio/Oasys2/lens_interface_1b.h5',
    material='Be', density=1.85,
    f_r_ind=3, # 0=constant, 1=file, 2=xraylib, 3=dabax
    file_r_ind='<none>',
    r_ind=1,
    r_attenuation=0,
    dabax=DabaxXraylib(file_f1f2="f1f2_Windt.dat", file_CrossSec="CrossSec_EPDL97.dat"),
    apply_to_lost=True,
    shift_thickness_to_zero=False,
    thickness_scaling=1,
    coordinate_scaling=1,
    invert_surface=False,
    )
from syned.beamline.element_coordinates import ElementCoordinates
coordinates = ElementCoordinates(p=10, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
movements = None
from shadow4.beamline.optical_elements.phase_deflectors.s4_numerical_mesh_phase_deflector import S4NumericalMeshPhaseDeflectorElement
beamline_element = S4NumericalMeshPhaseDeflectorElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)

beam, footprint = beamline_element.trace_beam()

beamline.append_beamline_element(beamline_element)


# test plot
if True:
   from srxraylib.plot.gol import plot_scatter

   beam.retrace(27.4327)
   plot_scatter(1e6 * beam.get_column(1, nolost=1), 1e6 * beam.get_column(3, nolost=1),
                title='(X,Z) in microns I= %d' % (beam.get_intensity(nolost=1)))