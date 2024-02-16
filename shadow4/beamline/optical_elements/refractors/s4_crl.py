import numpy

from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.optical_elements.refractors.crl import CRL
from syned.beamline.shape import Rectangle, Ellipse, Circle

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.s4_beamline_element import S4BeamlineElement
from shadow4.beamline.optical_elements.refractors.s4_conic_interface import S4ConicInterfaceElement
from shadow4.beamline.s4_optical_element_decorators import S4RefractiveLensOpticalElementDecorator
from shadow4.beamline.optical_elements.refractors.s4_lens import _get_lens_interfaces

class S4CRL(CRL, S4RefractiveLensOpticalElementDecorator):
    def __init__(self,
                 name="Undefined",
                 n_lens=1,
                 piling_thickness=0.0,  # syned stuff,
                 boundary_shape=None,  # syned stuff, replaces "diameter" in the shadow3 append_lens
                 material="",          # syned stuff, not (yet) used
                 thickness=0.0,        # syned stuff, lens thickness [m] (distance between the two interfaces at the center of the lenses)
                 surface_shape=1,      # now: 0=plane, 1=sphere, 2=parabola, 3=conic coefficients
                                       # (in shadow3: 1=sphere 4=paraboloid, 5=plane)
                 convex_to_the_beam=1,  # for surface_shape: convexity of the first interface exposed to the beam 0=No, 1=Yes
                                        # the second interface has opposite convexity
                 cylinder_angle=0,      # for surface_shape: 0=not cylindricaL, 1=meridional 2=sagittal
                 ri_calculation_mode=0,   # source of refraction indices and absorption coefficients
                                 # 0=User
                                 # 1=prerefl file
                                 # 2=direct calculation using xraylib
                                 # 3=direct calculation using dabax
                 prerefl_file=None,    # for ri_calculation_mode=0: file name (from prerefl) to get the refraction index.
                 refraction_index=1.0, # for ri_calculation_mode=1: n (real)
                 attenuation_coefficient=0.0, # for ri_calculation_mode=1: mu in cm^-1 (real)
                 radius=500e-6,        # for surface_shape=(1,2): lens radius [m] (for spherical, or radius at the tip for paraboloid)
                 conic_coefficients=[0.0]*10,   # for surface_shape = 3: the conic coefficients
                 ):
        S4RefractiveLensOpticalElementDecorator.__init__(self,
                                                         surface_shape,
                                                         convex_to_the_beam,
                                                         cylinder_angle,
                                                         ri_calculation_mode,
                                                         prerefl_file,
                                                         refraction_index,
                                                         attenuation_coefficient,
                                                         radius,
                                                         conic_coefficients)
        CRL.__init__(self,
                     name=name,
                     n_lens=n_lens,
                     surface_shape1=self.get_surface_shape_instance()[0],
                     surface_shape2=self.get_surface_shape_instance()[1],
                     boundary_shape=boundary_shape,
                     material=material,
                     thickness=thickness,
                     piling_thickness=piling_thickness)

        self.__inputs = {
            "name": name,
            "n_lens" : n_lens,
            "boundary_shape": boundary_shape,
            "material": material,
            "thickness": thickness,
            "piling_thickness" : piling_thickness,
            "surface_shape": surface_shape,
            "convex_to_the_beam": convex_to_the_beam,
            "cylinder_angle": cylinder_angle,
            "ri_calculation_mode": ri_calculation_mode,
            "prerefl_file": prerefl_file,
            "refraction_index": refraction_index,
            "attenuation_coefficient": attenuation_coefficient,
            "radius": radius,
            "conic_coefficients": repr(conic_coefficients),
        }

    def to_python_code_boundary_shape(self):
        txt = "" # "\nfrom shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirror"
        bs = self._boundary_shape
        if bs is None:
            txt += "\nboundary_shape = None"
        elif isinstance(bs, Rectangle):
            txt += "\nfrom syned.beamline.shape import Rectangle"
            txt += "\nboundary_shape = Rectangle(x_left=%g, x_right=%g, y_bottom=%g, y_top=%g)" % bs.get_boundaries()
        elif isinstance(bs, Circle):
            txt += "\nfrom syned.beamline.shape import Circle"
            txt += "\nboundary_shape = Circle(radius=%g, x_center=%g,y_center=%g)" % bs.get_boundaries()
        elif isinstance(bs, Ellipse):
            txt += "\nfrom syned.beamline.shape import Ellipse"
            txt += "\nboundary_shape = Ellipse(a_axis_min=%g, a_axis_max=%g, b_axis_min=%g, b_axis_max=%g)" % bs.get_boundaries()
        return txt

    def to_python_code(self, **kwargs):
        txt = self.to_python_code_boundary_shape()

        txt_pre = """
from shadow4.beamline.optical_elements.refractors.s4_crl import S4CRL

optical_element = S4CRL(name='{name:s}',
     n_lens={n_lens:d},
     piling_thickness={piling_thickness:g}, # syned stuff
     boundary_shape=boundary_shape,         # syned stuff, replaces "diameter" in the shadow3 append_lens
     material="", # syned stuff, not (yet) used
     thickness={thickness}, # syned stuff, lens thickness [m] (distance between the two interfaces at the center of the lenses)
     surface_shape={surface_shape}, # now: 0=plane, 1=sphere, 2=parabola, 3=conic coefficients
                                    # (in shadow3: 1=sphere 4=paraboloid, 5=plane)
     convex_to_the_beam={convex_to_the_beam}, # for surface_shape: convexity of the first interface exposed to the beam 0=No, 1=Yes
     cylinder_angle={cylinder_angle}, # for surface_shape: 0=not cylindricaL, 1=meridional 2=sagittal
     ri_calculation_mode={ri_calculation_mode},   # source of refraction indices and absorption coefficients
                                     # 0=User, 1=prerefl file, 2=xraylib, 3=dabax
     prerefl_file='{prerefl_file:s}', # for ri_calculation_mode=0: file name (from prerefl) to get the refraction index.
     refraction_index={refraction_index:g}, # for ri_calculation_mode=1: n (real)
     attenuation_coefficient={attenuation_coefficient:g}, # for ri_calculation_mode=1: mu in cm^-1 (real)
     radius={radius:g}, # for surface_shape=(1,2): lens radius [m] (for spherical, or radius at the tip for paraboloid)
     conic_coefficients=None, # for surface_shape = 3: the conic coefficients [todo: noy yet implemented]
     )
    """
        txt += txt_pre.format(**self.__inputs)
        return txt

    def get_lens_interfaces(self):
        single_lens_optical_surfaces = self.get_optical_surface_instance()
        lens_interfaces              = numpy.full((self._n_lens, 2), None)
        boundary_shape               = self.get_boundary_shape()

        for lens_index in range(self._n_lens):
            lens_interfaces[lens_index, 0], \
            lens_interfaces[lens_index, 1] = \
                _get_lens_interfaces(lens_optical_surfaces=single_lens_optical_surfaces,
                                     boundary_shape=boundary_shape,
                                     ri_calculation_mode=self._ri_calculation_mode,
                                     refraction_index=self._refraction_index,
                                     attenuation_coefficient=self._attenuation_coefficient,
                                     prerefl_file=self._prerefl_file)

        return lens_interfaces


class S4CRLElement(S4BeamlineElement):
    def __init__(self,
                 optical_element: S4CRL = None,
                 coordinates: ElementCoordinates = None,
                 input_beam: S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4CRL(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         input_beam=input_beam)

    def to_python_code(self, **kwargs):
        txt = "\n\n# optical element number XX"
        txt += self.get_optical_element().to_python_code()
        coordinates = self.get_coordinates()
        txt += "\nfrom syned.beamline.element_coordinates import ElementCoordinates"
        txt += "\nimport numpy"
        txt += "\ncoordinates = ElementCoordinates(p=%g, q=%g, angle_radial=0, angle_azimuthal=0, angle_radial_out=numpy.pi)" % \
               (coordinates.p(), coordinates.q())
        txt += "\nfrom shadow4.beamline.optical_elements.refractors.s4_crl import S4CRLElement"
        txt += "\nbeamline_element = S4CRLElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)"
        txt += "\n\nbeam, mirr = beamline_element.trace_beam()"
        return txt

    def trace_beam(self, **params):
        input_beam = self.get_input_beam().duplicate()
        oe         = self.get_optical_element()

        optical_surfaces = oe.get_lens_interfaces()
        p, q, angle_radial, angle_radial_out, angle_azimuthal = self.get_coordinates().get_positions()
        n_lens = oe.get_n_lens()
        for lens_index in range(n_lens):
            if lens_index==0: source_plane = p
            else:             source_plane = oe.get_piling_thickness()
            if lens_index==n_lens-1: image_plane = q
            else:                    image_plane = 0.0

            coordinates_1 = ElementCoordinates(p=source_plane, q=oe.get_thickness()*0.5, angle_radial=angle_radial, angle_radial_out=numpy.pi,          angle_azimuthal=angle_azimuthal)
            coordinates_2 = ElementCoordinates(p=oe.get_thickness()*0.5, q=image_plane,  angle_radial=0.0,           angle_radial_out=angle_radial_out, angle_azimuthal=0.0)

            beamline_element_1 = S4ConicInterfaceElement(optical_element=optical_surfaces[lens_index, 0], coordinates=coordinates_1, input_beam=input_beam)
            if lens_index==0: beam1, footprint1 = beamline_element_1.trace_beam()
            else:             beam1, _          = beamline_element_1.trace_beam()

            beamline_element_2 = S4ConicInterfaceElement(optical_element=optical_surfaces[lens_index, 1], coordinates=coordinates_2, input_beam=beam1)
            if lens_index==n_lens-1: beam2, footprint2 = beamline_element_2.trace_beam()
            else:                    beam2, _          = beamline_element_2.trace_beam()

            if lens_index < n_lens-1: input_beam = beam2.duplicate()

            print(beamline_element_2.info())

        return beam2, [footprint1, footprint2]


if __name__ == "__main__":
    from shadow4.physical_models.prerefl.prerefl import PreRefl

    #
    # collimated source
    #
    from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical

    src = SourceGeometrical()
    src.set_energy_distribution_singleline(value=5000, unit='eV')
    src.set_spatial_type_rectangle(width=1e-3, height=1e-3)
    src.set_angular_distribution_uniform(0, 0, 0, 0)

    beam = src.get_beam()

    filename = "Be.dat"

    PreRefl.prerefl(interactive=False, SYMBOL="Be", FILE=filename, DENSITY=1.848, E_MIN=4500, E_MAX=5500, E_STEP=1)

    #
    # lens
    #
    lens = S4CRL(n_lens=3, piling_thickness=2.5e-3, prerefl_file=filename)
    e = S4CRLElement(optical_element=lens,
                      coordinates=ElementCoordinates(p=10, q=20, angle_radial=0, angle_azimuthal=0, angle_radial_out=numpy.pi),
                      input_beam=beam)

    e.trace_beam()

    print(e.to_python_code())