#
# the screen optical element:
#       deals with screens, slits, beam-stoppers and absorbers (as in shadow3)
#       it is a stand-alone optical element (contrary to shadow4)
#
#
#
import numpy

from syned.beamline.optical_elements.ideal_elements.screen import Screen as SyScreen


from shadow4.syned.absorbers.beam_stopper import BeamStopper as SyBeamStopper   # TODO: syned.beamline.optical_elements.
from shadow4.syned.absorbers.filter import Filter as SyFilter                   # TODO: syned.beamline.optical_elements.
from shadow4.syned.absorbers.holed_filter import HoledFilter as SyHoledFilter   # TODO: syned.beamline.optical_elements.
from shadow4.syned.absorbers.slit import Slit as SySlit                         # TODO: syned.beamline.optical_elements.


from syned.beamline.beamline_element import BeamlineElement
from syned.beamline.element_coordinates import ElementCoordinates
from shadow4.syned.shape import Rectangle, Ellipse, TwoEllipses, MultiplePatch # TODO from syned.beamline.shape

from shadow4.physical_models.prerefl.prerefl import PreRefl

class Screen(object):
    def __init__(self, beamline_element_syned = None,
                 i_abs=False, # include absorption
                 thick=0.0, # thickness of the absorber (in SI)
                 file_abs="", # if i_abs=True, the material file (from prerefl)
                ):

        self._i_abs = i_abs
        self._thick = thick
        self._file_abs = file_abs

        if beamline_element_syned is None:
            self._beamline_element_syned = BeamlineElement(
                SyScreen(name="Undefined"),
                ElementCoordinates(p=0.0, q=0.0, angle_radial=0.0, angle_azimuthal=0.0))
        else:
            ok = False
            for obj in [SyScreen, SySlit, SyBeamStopper, SyFilter, SyHoledFilter]:
                if isinstance(beamline_element_syned._optical_element, obj): ok = True
            if ok:
                self._beamline_element_syned = beamline_element_syned
            else:
                raise Exception("Please initialize shadow4 Screen with syned Screen, Slit, BeamStopper, Filter or HoledFilter")


    def set_positions(self, p, q):
        self._beamline_element_syned.get_coordinates()._p = p
        self._beamline_element_syned.get_coordinates()._q = q
        self._beamline_element_syned.get_coordinates()._angle_radial = 0.0
        self._beamline_element_syned.get_coordinates()._angle_azimuthal = 0.0

    def get_positions(self):
        return self._beamline_element_syned.get_coordinates()._p, \
            self._beamline_element_syned.get_coordinates()._q


    def info(self):
        if self._beamline_element_syned is not None:
            return (self._beamline_element_syned.info())

    def trace_beam(self,beam1,flag_lost_value=-1):
        beam = beam1.duplicate()

        p,q = self.get_positions()

        if p != 0.0:
            beam.retrace(p, resetY=True)

        if isinstance(self._beamline_element_syned._optical_element, SyScreen):
            apply_crop = False
        elif isinstance(self._beamline_element_syned._optical_element, SySlit):
            apply_crop = True
            negative = False
        elif isinstance(self._beamline_element_syned._optical_element, SyBeamStopper):
            apply_crop = True
            negative = True
        elif isinstance(self._beamline_element_syned._optical_element, SyFilter):
            apply_crop = True
            negative = False
        elif isinstance(self._beamline_element_syned._optical_element, SyHoledFilter):
            apply_crop = True
            negative = True

        if apply_crop:
            shape = self._beamline_element_syned._optical_element.get_boundary_shape()
            if isinstance(shape, type(None)):
                pass
            elif isinstance(shape, Rectangle):
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

                # print(">>>>> INSIDE: ",  INSIDE, numpy.where(~INSIDE))

            else:
                raise Exception("Undefined slit shape")


        if self._i_abs:

            thickness = self._thick
            # the thickness in Filter syned is ignored. TODO: discuss it it could overwrite
            # thickness = self._beamline_element_syned._optical_element.get_thickness()

            if self._file_abs != "":
                try:
                    pr = PreRefl()
                    pr.read_preprocessor_file(self._file_abs)
                    print(pr.info())
                except:
                    raise Exception("Failed to load preprocessor (prerefl) file %s " % self._file_abs)

                energy = beam.get_column(26)
                # tmp = pr.get_attenuation_coefficient(energy[0],verbose=1)
                coeff = pr.get_attenuation_coefficient(energy)
                I_over_I0 = numpy.exp(- coeff * thickness * 1e2)
                sqrt_I_over_I0 = numpy.sqrt(I_over_I0)
                print(energy, coeff, I_over_I0)
                beam.apply_reflectivities(sqrt_I_over_I0, sqrt_I_over_I0)


        if q != 0.0:
            beam.retrace(q, resetY=True)

        return beam



