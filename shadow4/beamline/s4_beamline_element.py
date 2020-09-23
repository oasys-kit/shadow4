from shadow4.syned.element_coordinates import ElementCoordinates
from syned.beamline.beamline_element import BeamlineElement


class S4BeamlineElement(BeamlineElement):

    def __init__(self, optical_element=None, coordinates=None):
        super().__init__(optical_element, coordinates)

    def trace_beam(self, beam=None, **params):
        raise NotImplementedError()

    def info(self):
        return self.get_optical_element().info() + "\n" + self.get_coordinates().info()

    def set_positions(self, p=0.0, q=0.0, angle_radial=0.0, angle_radial_out=None, angle_azimuthal=0.0):

        self.get_coordinates()._p = p
        self.get_coordinates()._q = q
        self.get_coordinates()._angle_radial = angle_radial
        self.get_coordinates()._angle_radial_out = angle_radial_out
        self.get_coordinates()._angle_azimuthal = angle_azimuthal

    def get_positions(self):
        coordinates = self.get_coordinates()
        try:
            angle_radial_out = coordinates.angle_radial_out()
        except:
            angle_radial_out = None
        return coordinates.p(), \
            coordinates.q(), \
            coordinates.angle_radial(), \
            angle_radial_out, \
            coordinates.angle_azimuthal()

    def set_p_and_q(self, p=0.0, q=0.0):

        self.get_coordinates()._p = p
        self.get_coordinates()._q = q


    def get_p_and_q(self):
        coordinates = self.get_coordinates()
        return coordinates.p(), \
            coordinates.q()
