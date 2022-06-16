import numpy

from shadow4.syned.element_coordinates import ElementCoordinates
from shadow4.syned.shape import Plane

from shadow4.beamline.optical_elements.gratings.s4_grating import S4GratingElement, S4Grating
from shadow4.beamline.s4_optical_element import S4PlaneOpticalElement


class S4PlaneGrating(S4Grating, S4PlaneOpticalElement):
    def __init__(self,
                 name="Undefined",
                 boundary_shape=None,
                 ruling=800e3,
                 ruling_coeff_linear=0.0,
                 ruling_coeff_quadratic=0.0,
                 ruling_coeff_cubic=0.0,
                 ruling_coeff_quartic=0.0,
                 coating=None,
                 coating_thickness=None,
                 f_central=False,
                 f_phot_cent=0,
                 phot_cent=8000.0,
                 f_reflec=0,
                 material_constants_library_flag=0,  # 0=xraylib, 1=dabax, 2=shadow preprocessor
                 file_refl="",
                 order=0,
                 f_ruling=0,
                 ):

        S4PlaneOpticalElement.__init__(self)
        S4Grating.__init__(self,
                           name=name,
                           surface_shape=self.get_surface_shape_instance(),
                           boundary_shape=boundary_shape,
                           ruling=ruling,
                           ruling_coeff_linear=ruling_coeff_linear,
                           ruling_coeff_quadratic=ruling_coeff_quadratic,
                           ruling_coeff_cubic=ruling_coeff_cubic,
                           ruling_coeff_quartic=ruling_coeff_quartic,
                           coating=coating,
                           coating_thickness=coating_thickness,
                           f_central=f_central,
                           f_phot_cent=f_phot_cent,
                           f_reflec=f_reflec,
                           phot_cent=phot_cent,
                           material_constants_library_flag=material_constants_library_flag,  # 0=xraylib, 1=dabax, 2=shadow preprocessor
                           file_refl=file_refl,
                           order=order,
                           f_ruling=f_ruling,
                           )

    # def get_optical_surface_instance(self):
    #     return S4Conic.initialize_as_plane()

class S4PlaneGratingElement(S4GratingElement):
    def __init__(self, optical_element=None, coordinates=None):
        super().__init__(optical_element if optical_element is not None else S4Grating(),
                         coordinates if coordinates is not None else ElementCoordinates())
        if not isinstance(self.get_optical_element().get_surface_shape(), Plane):
            raise ValueError("Wrong Optical Element: only Plane shape is accepted")


if __name__ == "__main__":

    from shadow4.sources.source_geometrical.source_gaussian import SourceGaussian
    from shadow4.beam.beam import Beam
    from shadow4.tools.graphics import plotxy

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
    g = S4PlaneGrating(
        name = "my_grating",
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
        f_ruling=0,
        )

    coordinates_syned = ElementCoordinates(p = 10.0,
                                           q = 6.0,
                                           angle_radial = 88.840655 * numpy.pi / 180,
                                           angle_radial_out= 87.588577 * numpy.pi / 180,
                                           angle_azimuthal = 0.0)



    ge = S4PlaneGratingElement(optical_element=g, coordinates=coordinates_syned)

    print(ge.info())

    beam_out = ge.trace_beam(beam)
    plotxy(beam_out[0], 1, 3, title="Image 0", nbins=201)

