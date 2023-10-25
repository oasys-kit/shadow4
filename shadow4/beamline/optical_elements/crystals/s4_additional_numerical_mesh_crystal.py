import numpy
from syned.beamline.shape import NumericalMesh
from syned.beamline.element_coordinates import ElementCoordinates
from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.optical_elements.crystals.s4_crystal import S4CrystalElement, S4Crystal
from shadow4.beamline.optical_elements.crystals.s4_numerical_mesh_crystal import S4NumericalMeshCrystal
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements


class S4AdditionalNumericalMeshCrystal(S4NumericalMeshCrystal):
    def __init__(self,
                 ideal_crystal : S4Crystal = None,
                 numerical_mesh_crystal : S4NumericalMeshCrystal = None,
                 name="Crystal with Additional Numerical Mesh"):
        S4NumericalMeshCrystal.__init__(self, name=name,
                 boundary_shape=None if ideal_crystal is None else ideal_crystal.get_boundary_shape(),
                 xx=None if numerical_mesh_crystal is None else numerical_mesh_crystal._curved_surface_shape._xx,
                 yy=None if numerical_mesh_crystal is None else numerical_mesh_crystal._curved_surface_shape._yy,
                 zz=None if numerical_mesh_crystal is None else numerical_mesh_crystal._curved_surface_shape._zz,
                 surface_data_file="" if numerical_mesh_crystal is None else numerical_mesh_crystal._curved_surface_shape._surface_data_file,
                 # inputs related to crystal
                 material                       =None if ideal_crystal is None else ideal_crystal._material,
                 miller_index_h                 =1 if ideal_crystal is None else ideal_crystal._miller_index_h,
                 miller_index_k                 =1 if ideal_crystal is None else ideal_crystal._miller_index_k,
                 miller_index_l                 =1 if ideal_crystal is None else ideal_crystal._miller_index_l,
                 f_bragg_a                      =0 if ideal_crystal is None else ideal_crystal._f_bragg_a,
                 asymmetry_angle                =0.0 if ideal_crystal is None else ideal_crystal._asymmetry_angle,
                 is_thick                       =0   if ideal_crystal is None else ideal_crystal._is_thick,
                 thickness                      =0.010 if ideal_crystal is None else ideal_crystal._thickness,
                 f_central                      =0 if ideal_crystal is None else ideal_crystal._f_central,
                 f_phot_cent                    =0 if ideal_crystal is None else ideal_crystal._f_phot_cent,
                 phot_cent                      =8000.0 if ideal_crystal is None else ideal_crystal._phot_cent,
                 f_ext                          =0 if ideal_crystal is None else ideal_crystal._f_ext,
                 material_constants_library_flag=0 if ideal_crystal is None else ideal_crystal._material_constants_library_flag,
                 file_refl                      ="" if ideal_crystal is None else ideal_crystal._file_refl,
                 )

        self.__ideal_crystal         = ideal_crystal
        self.__numerical_mesh_crystal = numerical_mesh_crystal

        self.__inputs = {
            "name": name,
            "ideal_crystal": ideal_crystal,
            "numerical_mesh_crystal": numerical_mesh_crystal,
        }

    def to_python_code(self, **kwargs):


        txt = self.__ideal_crystal.to_python_code()
        txt += "ideal_crystal = optical_element"
        txt += self.__numerical_mesh_crystal.to_python_code()
        txt += "numerical_mesh_crystal = optical_element"

        txt += self.to_python_code_boundary_shape()
        txt_pre = """

from shadow4.beamline.optical_elements.crystals.s4_additional_numerical_mesh_crystal import S4AdditionalNumericalMeshCrystal
optical_element = S4AdditionalNumericalMeshCrystal(name='{name:s}', ideal_crystal=ideal_crystal, numerical_mesh_crystal=numerical_mesh_crystal)
    """
        txt += txt_pre.format(**self.__inputs)
        return txt

class S4AdditionalNumericalMeshCrystalElement(S4CrystalElement):
    def __init__(self,
                 optical_element: S4AdditionalNumericalMeshCrystal = None,
                 coordinates: ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam: S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4AdditionalNumericalMeshCrystal(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         movements=movements,
                         input_beam=input_beam)
        if not isinstance(self.get_optical_element().get_surface_shape(), NumericalMesh):
            raise ValueError("Wrong Optical Element: only Surface Data shape is accepted")

    def to_python_code(self, **kwargs):
        txt = "\n\n# optical element number XX"
        txt += self.get_optical_element().to_python_code()
        coordinates = self.get_coordinates()
        txt += "\nfrom syned.beamline.element_coordinates import ElementCoordinates"
        txt += "\ncoordinates = ElementCoordinates(p=%g, q=%g, angle_radial=%g, angle_azimuthal=%g, angle_radial_out=%g)" % \
               (coordinates.p(), coordinates.q(), coordinates.angle_radial(), coordinates.angle_azimuthal(), coordinates.angle_radial_out())

        txt += self.to_python_code_movements()

        txt += "\nfrom shadow4.beamline.optical_elements.crystals.s4_additional_numerical_mesh_crystal import S4AdditionalNumericalMeshCrystalElement"
        txt += "\nbeamline_element = S4AdditionalNumericalMeshCrystalElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)"
        txt += "\n\nbeam, mirr = beamline_element.trace_beam()"
        return txt

if __name__ == "__main__":
    import numpy
    do_plot = True
    use_errors = True
    from srxraylib.plot.gol import plot_scatter

    #
    #
    #
    from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical

    light_source = SourceGeometrical(name='SourceGeometrical', nrays=10000, seed=5676561)
    light_source.set_spatial_type_gaussian(sigma_h=5e-06, sigma_v=0.000001)
    light_source.set_angular_distribution_gaussian(sigdix=0.000001, sigdiz=0.000001)
    light_source.set_energy_distribution_singleline(1.000000, unit='A')
    light_source.set_polarization(polarization_degree=1.000000, phase_diff=0.000000, coherent_beam=0)
    beam = light_source.get_beam()

    # optical element number XX
    from shadow4.beamline.optical_elements.crystals.s4_sphere_crystal import S4SphereCrystal

    optical_element = S4SphereCrystal(name='xxx', material="Si", boundary_shape=None,
                                        is_cylinder=1, cylinder_direction=0,
                                        convexity=1, radius=1.0,
                                        )

    from syned.beamline.element_coordinates import ElementCoordinates

    coordinates = ElementCoordinates(p=10, q=10, angle_radial=1.54985)


    if use_errors:
        from syned.beamline.shape import Rectangle
        mesh_element = S4AdditionalNumericalMeshCrystal(name="C1",
                                                       ideal_crystal=optical_element,
                                                       numerical_mesh_crystal=S4NumericalMeshCrystal(surface_data_file="/users/srio/Oasys/mirrors_branch3_mesh.hdf5",
                                                                                                   boundary_shape=Rectangle(x_left=-0.05,
                                                                                                                            x_right=0.05,
                                                                                                                            y_bottom=-0.5,
                                                                                                                            y_top=0.5)))

        beamline_element = S4AdditionalNumericalMeshCrystalElement(optical_element=mesh_element,
                                                         coordinates=ElementCoordinates(p=10.0,
                                                                                        q=6.0,
                                                                                        angle_radial=numpy.radians(88.8)),
                                                        input_beam=beam)

    else:
        from shadow4.beamline.optical_elements.crystals.s4_sphere_crystal import S4SphereCrystalElement

        beamline_element = S4SphereCrystalElement(optical_element=optical_element, coordinates=coordinates,
                                                    input_beam=beam)
    #
    # run
    #
    beam1, mirr1 = beamline_element.trace_beam()
    #
    #
    #
    if do_plot:

        plot_scatter(1e6*beam1.get_column(1), 1e6*beam1.get_column(3), title="Cols 1,3 / um")


    # from shadow4.beamline.optical_elements.mirrors.s4_numerical_mesh_mirror import S4NumericalMeshMirror, S4NumericalMeshMirrorElement
    # m = S4NumericalMeshMirror()
    # e = S4NumericalMeshMirrorElement()
    # m = S4AdditionalNumericalMeshMirror()
    # e = S4AdditionalNumericalMeshMirrorElement(None, None, None)

