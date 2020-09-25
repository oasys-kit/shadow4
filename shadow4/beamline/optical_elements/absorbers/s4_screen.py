#
# the screen optical element:
#       deals with screens, slits, beam-stoppers and absorbers (as in shadow3)
#       it is a stand-alone optical element (contrary to shadow4)
#
#
#
import numpy

# from syned.beamline.optical_elements.ideal_elements.screen import Screen as SyScreen
#
#
# from shadow4.syned.absorbers.beam_stopper import BeamStopper as SyBeamStopper   # TODO: syned.beamline.optical_elements.
# from shadow4.syned.absorbers.filter import Filter as SyFilter                   # TODO: syned.beamline.optical_elements.
# from shadow4.syned.absorbers.holed_filter import HoledFilter as SyHoledFilter   # TODO: syned.beamline.optical_elements.
# from shadow4.syned.absorbers.slit import Slit as SySlit                         # TODO: syned.beamline.optical_elements.

from syned.beamline.optical_elements.absorbers.absorber import Absorber

from shadow4.syned.element_coordinates import ElementCoordinates
from shadow4.syned.shape import Rectangle, Ellipse, MultiplePatch # TODO from syned.beamline.shape

from shadow4.physical_models.prerefl.prerefl import PreRefl

from shadow4.beamline.s4_optical_element import S4OpticalElement
from shadow4.beamline.s4_beamline_element import S4BeamlineElement

class S4Screen(Absorber, S4OpticalElement):
    def __init__(self,
                 name="Undefined", boundary_shape=None,  # for syned absorber
                 i_abs=False,  # include absorption
                 i_stop=False, # aperture/stop
                 thick=0.0,    # thickness of the absorber (in SI)
                 file_abs="",  # if i_abs=True, the material file (from prerefl)
                ):
        super().__init__(name=name, boundary_shape=boundary_shape)
        self._i_abs = i_abs
        self._i_stop = i_stop
        self._thick = thick
        self._file_abs = file_abs

class S4ScreenElement(S4BeamlineElement):

    def __init__(self, optical_element=None, coordinates=None):
        super().__init__(optical_element if optical_element is not None else S4Screen(),
                         coordinates if coordinates is not None else ElementCoordinates())


    def trace_beam(self,beam1,flag_lost_value=-1):
        beam0 = beam1.duplicate()

        p,q = self.get_coordinates().get_p_and_q()

        oe = self.get_optical_element()
        coor = self.get_coordinates()

        if p != 0.0:
            beam0.retrace(p, resetY=True)

        beam = beam0.duplicate()

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
                beam.crop_rectangle(1, x_left, x_right, 3, y_bottom, y_top,
                                    negative=negative, flag_lost_value=flag_lost_value)
            elif isinstance(shape, Ellipse):
                a_axis_min, a_axis_max, b_axis_min, b_axis_max = shape.get_boundaries()
                beam.crop_ellipse(1, a_axis_min, a_axis_max, 3, b_axis_min, b_axis_max,
                                  negative=negative, flag_lost_value=flag_lost_value)
            elif isinstance(shape, MultiplePatch):
                x = beam.get_column(1)
                y = beam.get_column(3)
                INSIDE = numpy.zeros(x.size, numpy.bool_)
                for i,patch in enumerate(shape.get_patches()):
                    # print(">>> patch: ",patch.info())
                    inside = patch.check_inside_vector(x, y)
                    INSIDE = numpy.logical_or(INSIDE,inside)

                # print(">>>>",x[0:10],y[0:10],inside[0:10])

                flag = beam.get_column(10)
                if negative:
                    flag[numpy.where(INSIDE)] = flag_lost_value
                else:
                    flag[numpy.where(~INSIDE)] = flag_lost_value
                beam.rays[:, 9] = flag

            else:
                print(">>>>>>>>>>>>  NO CROP !", shape)
                pass


        if oe._i_abs:

            thickness = oe._thick
            # the thickness in Filter syned is ignored. TODO: discuss it it could overwrite
            # thickness = oe._beamline_element_syned._optical_element.get_thickness()

            if oe._file_abs != "":
                try:
                    pr = PreRefl()
                    pr.read_preprocessor_file(oe._file_abs)
                    print(pr.info())
                except:
                    raise Exception("Failed to load preprocessor (prerefl) file %s " % oe._file_abs)

                energy = beam.get_column(26)
                # tmp = pr.get_attenuation_coefficient(energy[0],verbose=1)
                coeff = pr.get_attenuation_coefficient(energy)
                I_over_I0 = numpy.exp(- coeff * thickness * 1e2)
                sqrt_I_over_I0 = numpy.sqrt(I_over_I0)
                print(energy, coeff, I_over_I0)
                beam.apply_reflectivities(sqrt_I_over_I0, sqrt_I_over_I0)


        if q != 0.0:
            beam.retrace(q, resetY=True)

        return beam, beam0



