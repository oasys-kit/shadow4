import numpy

from shadow4.syned.element_coordinates import ElementCoordinates
from syned.beamline.optical_elements.ideal_elements.lens import IdealLens
from shadow4.beamline.s4_optical_element import S4OpticalElement
from shadow4.beamline.s4_beamline_element import S4BeamlineElement


class S4IdealLens(IdealLens, S4OpticalElement):
    def __init__(self, name="Undefined",focal_x=0.0,focal_y=0.0):
        super().__init__(name=name,focal_x=focal_x,focal_y=focal_y)


class S4IdealLensElement(S4BeamlineElement):

    def __init__(self, optical_element=None, coordinates=None):
        super().__init__(optical_element if optical_element is not None else S4IdealLens(),
                         coordinates if coordinates is not None else ElementCoordinates())

    def get_focalX(self):
        return self.get_optical_element()._focal_x

    def get_focalZ(self):
        return self.get_optical_element()._focal_y

    def trace_beam(self, beam1, flag_lost_value=-1):
        beam0 = beam1.duplicate()



        p,q = self.get_coordinates().get_p_and_q()
        if p != 0.0:
            beam0.retrace(p,resetY=True)

        beam = beam0.duplicate()

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

        return beam, beam0


class S4SuperIdealLens(IdealLens, S4OpticalElement):
    def __init__(self, name="Undefined",focal_p_x=0.0, focal_p_y=0.0,
                                        focal_q_x=0.0, focal_q_y=0.0,):
        super().__init__(name=name,focal_x=1.0/(1/focal_p_x+1/focal_q_x),focal_y=1.0/(1/focal_p_y+1/focal_q_y))

        self._focal_p_x = focal_p_x
        self._focal_p_y = focal_p_y
        self._focal_q_x = focal_q_x
        self._focal_q_y = focal_q_y


class S4SuperIdealLensElement(S4BeamlineElement):

    def __init__(self, optical_element=None, coordinates=None):
        super().__init__(optical_element if optical_element is not None else S4SuperIdealLens(),
                         coordinates if coordinates is not None else ElementCoordinates())


    def trace_beam(self,beam1):
        beam0 = beam1.duplicate()

        lens = self.get_optical_element()

        p, q = self.get_coordinates().get_p_and_q()

        if p != 0.0:
            beam0.retrace(p,resetY=True)

        beam = beam0.duplicate()

        # rotate around Z; watch out the minus!!
        if lens._focal_p_x != 0.0:
            tan_theta_p = - beam.get_column(1) / lens._focal_p_x
        else:
            tan_theta_p = 0.0

        if lens._focal_q_x != 0.0:
            tan_theta_q = - beam.get_column(1) / lens._focal_q_x
        else:
            tan_theta_q = 0.0

        two_theta = numpy.arctan(tan_theta_p) + numpy.arctan(tan_theta_q)
        beam.rotate(two_theta,axis=3,rad=True)

        # rotate around X
        if lens._focal_p_y != 0.0:
            tan_theta_p = beam.get_column(3) / lens._focal_p_y
        else:
            tan_theta_p = 0.0

        if lens._focal_q_y != 0.0:
            tan_theta_q = beam.get_column(3) / lens._focal_q_y
        else:
            tan_theta_q = 0.0

        two_theta = numpy.arctan(tan_theta_p) + numpy.arctan(tan_theta_q)
        beam.rotate(two_theta,axis=1,rad=True)

        if q != 0.0:
            beam.retrace(q,resetY=True)

        return beam, beam0

