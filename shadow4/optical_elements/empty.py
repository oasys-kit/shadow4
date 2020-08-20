
import numpy

#
# this is not syned compatible!
#
class Empty(object):

    def __init__(self, name="Empty", angle_radial=0.0, angle_radial_ref=numpy.pi,
                 angle_azimuthal=0.0, p=0.0, q=0.0):
        self._angle_radial = angle_radial
        self._angle_radial_ref = angle_radial_ref
        self._angle_azimuthal = angle_azimuthal
        self._name = name
        self._p = p
        self._q = q

    def trace_beam(self,beam1,undo_shadow_orientation_angle_rotation=False):
        beam = beam1.duplicate()

        theta_grazing1 = numpy.pi / 2 - self._angle_radial
        theta_grazing2 = numpy.pi / 2 - self._angle_radial_ref
        alpha1 = self._angle_azimuthal

        #
        beam = beam1.duplicate()



        #
        # put beam in mirror reference system
        #
        # print(">>> ",beam.rays[0, 6:9])
        # print(">>>> rotate alpha1 [deg]: ",alpha1*180/numpy.pi)
        beam.rotate(alpha1, axis=2)
        # print(">>> ", beam.rays[0, 6:9])

        # print(">>>> rotate theta_grazing1 [deg]: ", theta_grazing1*180/numpy.pi)
        beam.rotate(theta_grazing1, axis=1)
        # print(">>> ", beam.rays[0, 6:9])

        # print(">>>> translate ",[0.0, -self._p * numpy.cos(theta_grazing1), self._p * numpy.sin(theta_grazing1)])
        # WHY this is not needed?????????
        # beam.translation([0.0, -self._p * numpy.cos(theta_grazing1), self._p * numpy.sin(theta_grazing1)])
        # print(">>> ", beam.rays[0, 6:9])


        #
        # oe does nothing
        #
        # beam.retrace(self._p, resetY=False)


        #
        # from oe reference system to image plane
        #

        beam_out = beam.duplicate()
        # print(">>>> rotate theta_grazing2 deg: ", theta_grazing2*180/numpy.pi)
        beam_out.rotate(theta_grazing2, axis=1)

        if undo_shadow_orientation_angle_rotation:
            beam_out.rotate(-alpha1, axis=2)

        # print(">>>> retrace  ", self._q)
        beam_out.retrace(self._p + self._q, resetY=True)

        return beam_out, beam

