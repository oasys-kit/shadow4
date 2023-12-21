#
# the screen optical element:
#       deals with screens, slits, beam-stoppers and absorbers (as in shadow3)
#       it is a stand-alone optical element (contrary to shadow4)
#
#
#
import numpy

from syned.beamline.optical_elements.absorbers.absorber import Absorber

from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.shape import Rectangle, Ellipse, MultiplePatch

from shadow4.physical_models.prerefl.prerefl import PreRefl

from shadow4.beamline.s4_optical_element_decorators import S4OpticalElementDecorator
from shadow4.beamline.s4_beamline_element import S4BeamlineElement
from shadow4.beam.s4_beam import S4Beam

class S4Screen(Absorber, S4OpticalElementDecorator):
    def __init__(self,
                 name="Undefined", boundary_shape=None,  # for syned absorber
                 i_abs=0,      # include absorption: 0=No, 1=using preprocessor file
                               # 2=direct calculation using xraylib
                               # 3=direct calculation using dabax
                 i_stop=False, # aperture/stop
                 thick=0.0,    # thickness of the absorber (in SI)
                 file_abs="",  # if i_abs=1, the material file (from prerefl)
                 material="",  # if i_abs=2,3, the material name
                 density=1.0,  # if i_abs=2,3, the material density in g/cm^3
                ):
        super().__init__(name=name, boundary_shape=boundary_shape)
        self._i_abs = i_abs
        self._i_stop = i_stop
        self._thick = thick
        self._file_abs = file_abs
        self._material = material
        self._density = density

        self.__inputs = {
            "name": name,
            "boundary_shape": boundary_shape,
            "i_abs":    i_abs,
            "i_stop":   i_stop,
            "thick":    thick ,
            "material": material,
            "file_abs": file_abs,
            "density": density,
        }

    def to_python_code_boundary_shape(self):
        txt = ""
        bs = self._boundary_shape
        if bs is None:
            txt += "\nboundary_shape = None"
        if isinstance(bs, Rectangle):
            txt += "\nfrom syned.beamline.shape import Rectangle"
            txt += "\nboundary_shape = Rectangle(x_left=%g, x_right=%g, y_bottom=%g, y_top=%g)" % bs.get_boundaries()
        elif isinstance(bs, Ellipse):
            txt += "\nfrom syned.beamline.shape import Ellipse"
            txt += "\nboundary_shape = Ellipse(a_axis_min=%g, a_axis_max=%g, b_axis_min=%g, b_axis_max=%g)" % bs.get_boundaries()
        else:
            txt += "\nboundary_shape = None"
        return txt

    def to_python_code(self, **kwargs):
        txt = self.to_python_code_boundary_shape()
        txt_pre = """

from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
optical_element = S4Screen(name='{name:s}', boundary_shape=boundary_shape,
    i_abs={i_abs:d}, # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
    i_stop={i_stop:d}, thick={thick:g}, file_abs='{file_abs:s}', material='{material:s}', density={density:g})
"""
        txt += txt_pre.format(**self.__inputs)
        return txt

class S4ScreenElement(S4BeamlineElement):

    def __init__(self, optical_element : S4Screen = None, coordinates : ElementCoordinates = None, input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4Screen(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         input_beam=input_beam)

    def trace_beam(self, **params):
        flag_lost_value = params.get("flag_lost_value", -1)

        footprint = self.get_input_beam().duplicate()

        p, q = self.get_coordinates().get_p_and_q()

        oe = self.get_optical_element()

        if p != 0.0: footprint.retrace(p, resetY=True)

        output_beam = footprint.duplicate()

        apply_crop = True
        negative = oe._i_stop

        # if isinstance(self._beamline_element_syned._optical_element, SyScreen):
        #     apply_crop = False
        # elif isinstance(self._beamline_element_syned._optical_element, SySlit):
        #     apply_crop = True
        #     negative = False
        # elif isinstance(self._beamline_element_syned._optical_element, SyBeamStopper):
        #     apply_crop = True
        #     negative = True
        # elif isinstance(self._beamline_element_syned._optical_element, SyFilter):
        #     apply_crop = True
        #     negative = False
        # elif isinstance(self._beamline_element_syned._optical_element, SyHoledFilter):
        #     apply_crop = True
        #     negative = True

        if apply_crop:
            shape = oe.get_boundary_shape()

            if isinstance(shape, Rectangle):
                x_left, x_right, y_bottom, y_top = shape.get_boundaries()
                output_beam.crop_rectangle(1, x_left, x_right, 3, y_bottom, y_top, negative=negative, flag_lost_value=flag_lost_value)
            elif isinstance(shape, Ellipse):
                a_axis_min, a_axis_max, b_axis_min, b_axis_max = shape.get_boundaries()
                output_beam.crop_ellipse(1, a_axis_min, a_axis_max, 3, b_axis_min, b_axis_max, negative=negative, flag_lost_value=flag_lost_value)
            elif isinstance(shape, MultiplePatch):
                x = output_beam.get_column(1)
                y = output_beam.get_column(3)
                INSIDE = numpy.zeros(x.size, numpy.bool_)
                for i,patch in enumerate(shape.get_patches()):
                    # print(">>> patch: ",patch.info())
                    inside = patch.check_inside_vector(x, y)
                    INSIDE = numpy.logical_or(INSIDE,inside)

                flag = output_beam.get_column(10)
                if negative:
                    flag[numpy.where(INSIDE)] = flag_lost_value
                else:
                    flag[numpy.where(~INSIDE)] = flag_lost_value
                output_beam.rays[:, 9] = flag

            else:
                pass


        # reflectivity calculations
        if oe._i_abs > 0:
            thickness = oe._thick
            # the thickness in Filter syned is ignored. TODO: discuss it it could overwrite
            # thickness = oe._beamline_element_syned._optical_element.get_thickness()
            energy = output_beam.get_column(26)

            if oe._i_abs == 1:  # prerefl preprocessor file

                if oe._file_abs != "":
                    try:
                        pr = PreRefl()
                        pr.read_preprocessor_file(oe._file_abs)
                        print(pr.info())
                    except:
                        raise Exception("Failed to load preprocessor (prerefl) file %s " % oe._file_abs)

                    coeff = pr.get_attenuation_coefficient(energy)
            elif oe._i_abs == 2: # xraylib
                coeff = PreRefl.get_attenuation_coefficient_external_xraylib(photon_energy_ev=energy,
                                                                             material=oe._material,
                                                                             density=oe._density)
            elif oe._i_abs == 3: # dabax
                coeff = PreRefl.get_attenuation_coefficient_external_dabax(photon_energy_ev=energy,
                                                                             material=oe._material,
                                                                             density=oe._density)


            I_over_I0 = numpy.exp(- coeff * thickness * 1e2)
            sqrt_I_over_I0 = numpy.sqrt(I_over_I0)
            output_beam.apply_reflectivities(sqrt_I_over_I0, sqrt_I_over_I0)


        if q != 0.0: output_beam.retrace(q, resetY=True)

        return output_beam, footprint

    def to_python_code(self, **kwargs):
        txt = "\n\n# optical element number XX"
        txt += self.get_optical_element().to_python_code()
        coordinates = self.get_coordinates()
        txt += "\nfrom syned.beamline.element_coordinates import ElementCoordinates"
        txt += "\ncoordinates = ElementCoordinates(p=%g, q=%g, angle_radial=%g, angle_azimuthal=%g, angle_radial_out=%g)" % \
               (coordinates.p(), coordinates.q(), coordinates.angle_radial(), coordinates.angle_azimuthal(), coordinates.angle_radial_out())
        txt += "\nfrom shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement"
        txt += "\nbeamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)"
        txt += "\n\nbeam, footprint = beamline_element.trace_beam()"
        return txt

if __name__ == "__main__":
    o = S4Screen()

    e = S4ScreenElement()

    print(e.to_python_code())
