"""
Position of a beamline component within a beamline.

"""

from syned.syned_object import SynedObject

class ElementCoordinates(SynedObject):
    def __init__(self, p = 0.0, q = 0.0, angle_radial=0.0, angle_azimuthal=0.0, angle_radial_out=None):
        """

        :param p: distance from previous element.
        :param q: distance to next element.
        :param angle_radial: Radial inclination angle.
        :param angle_azimuthal: Azimuthal inclination angle.
        :return:
        """
        self._p               = p
        self._q               = q
        self._angle_radial    = angle_radial
        self._angle_azimuthal = angle_azimuthal
        self._angle_radial_out = angle_radial_out

        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("p",                "distance from previous continuation plane", "m"    ),
                    ("q",                "distance to next continuation plane",       "m"    ),
                    ("angle_radial",     "incident angle [to normal]",                "rad"  ),
                    ("angle_radial_out", "output angle [to normal]",                  "rad"),
                    ("angle_azimuthal",  "rotation along beam axis",                  "rad"  ),
                ])

    def p(self):
        return self._p

    def q(self):
        return self._q

    def angle_radial(self):
        return self._angle_radial

    def angle_radial_out(self):
        if self._angle_radial_out is None:
            return self.angle_radial()
        else:
            return self._angle_radial_out

    def angle_azimuthal(self):
        return self._angle_azimuthal