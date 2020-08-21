
import numpy
from collections import OrderedDict
from shadow4.compatibility.beam3 import Beam3

from syned.beamline.optical_elements.ideal_elements.lens import IdealLens as SyIdealLens
from syned.beamline.beamline_element import BeamlineElement
from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.optical_elements.ideal_elements.lens import IdealLens

class LensIdeal(object):
    def __init__(self, beamline_element_syned = None):
        if beamline_element_syned is None:
            self._beamline_element_syned = BeamlineElement(
                SyIdealLens(name="Undefined",focal_x=0.0,focal_y=0.0),
                ElementCoordinates(p=0.0, q=0.0, angle_radial=0.0, angle_azimuthal=0.0))
        else:
            if isinstance(beamline_element_syned._optical_element, SyIdealLens):
                self._beamline_element_syned = beamline_element_syned
            else:
                raise Exception("Please initialize shadow4 IdealLens with syned IdealLens")

    def initialize_from_keywors(self, name=None, focal_x=0.0, focal_z=0.0, p=0.0, q=0.0):
        self._beamline_element_syned._optical_element._focal_x = focal_x
        self._beamline_element_syned._optical_element._focal_y = focal_z
        if name is not None:
            self._name = name
        self.set_positions(p, q)

    def set_positions(self, p, q):
        self._beamline_element_syned.get_coordinates()._p = p
        self._beamline_element_syned.get_coordinates()._q = q
        self._beamline_element_syned.get_coordinates()._angle_radial = 0.0
        self._beamline_element_syned.get_coordinates()._angle_azimuthal = 0.0

    def get_positions(self):
        return self._beamline_element_syned.get_coordinates()._p, \
            self._beamline_element_syned.get_coordinates()._q

    def get_focalX(self):
        return self._beamline_element_syned._optical_element._focal_x

    def get_focalZ(self):
        return self._beamline_element_syned._optical_element._focal_y

    def info(self):
        if self._beamline_element_syned is not None:
            return (self._beamline_element_syned.info())

    def trace_beam(self,beam1):
        beam = beam1.duplicate()

        p,q = self.get_positions()

        if p != 0.0:
            beam.retrace(p,resetY=True)

        # rotate around Z
        if self.get_focalX() != 0.0:
            # whatch out the minus!!
            tan_two_theta = - beam.get_column(1) / self.get_focalX()
            beam.rotate(numpy.arctan(tan_two_theta),axis=3,rad=True)

        # rotate around X
        if self.get_focalZ() != 0.0:
            tan_two_theta = beam.get_column(3) / self.get_focalZ()
            beam.rotate(numpy.arctan(tan_two_theta),axis=1,rad=True)

        if q != 0.0:
            beam.retrace(q,resetY=True)

        return beam

#
# this is not syned compatible!
#
class LensSuperIdeal(object):

    def __init__(self, name="LensIdeal", focal_p_x=0.0, focal_p_z=0.0,
                                        focal_q_x=0.0, focal_q_z=0.0,
                                        p=0.0, q=0.0):
        self._focal_p_x = focal_p_x
        self._focal_p_z = focal_p_z
        self._focal_q_x = focal_q_x
        self._focal_q_z = focal_q_z
        self._name = name
        self._p = p
        self._q = q

    def trace_beam(self,beam1):
        beam = beam1.duplicate()


        if self._p != 0.0:
            beam.retrace(self._p,resetY=True)

        # rotate around Z; watch out the minus!!
        if self._focal_p_x != 0.0:
            tan_theta_p = - beam.get_column(1) / self._focal_p_x
        else:
            tan_theta_p = 0.0

        if self._focal_q_x != 0.0:
            tan_theta_q = - beam.get_column(1) / self._focal_q_x
        else:
            tan_theta_q = 0.0

        two_theta = numpy.arctan(tan_theta_p) + numpy.arctan(tan_theta_q)
        beam.rotate(two_theta,axis=3,rad=True)

        # rotate around X
        if self._focal_p_z != 0.0:
            tan_theta_p = beam.get_column(3) / self._focal_p_z
        else:
            tan_theta_p = 0.0

        if self._focal_q_z != 0.0:
            tan_theta_q = beam.get_column(3) / self._focal_q_z
        else:
            tan_theta_q = 0.0

        two_theta = numpy.arctan(tan_theta_p) + numpy.arctan(tan_theta_q)
        beam.rotate(two_theta,axis=1,rad=True)

        if self._q != 0.0:
            beam.retrace(self._q,resetY=True)

        return beam

