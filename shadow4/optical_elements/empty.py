import numpy
from syned.beamline.optical_elements.ideal_elements.screen import Screen as SyScreen
from syned.beamline.beamline_element import BeamlineElement
from shadow4.syned.element_coordinates import ElementCoordinates # TODO from syned.beamline.element_coordinates

class Empty(object):

    def __init__(self, beamline_element_syned = None):
        if beamline_element_syned is None:
            self._beamline_element_syned = BeamlineElement(
                SyScreen(name="Undefined"),
                ElementCoordinates(p=0.0, q=0.0, angle_radial=0.0, angle_radial_out=numpy.pi, angle_azimuthal=0.0))
        else:
            ok = False
            for obj in [SyScreen]:
                if isinstance(beamline_element_syned._optical_element, obj): ok = True
            if ok:
                self._beamline_element_syned = beamline_element_syned
            else:
                raise Exception("Please initialize shadow4 Empty with syned Screen")

    def set_positions(self, p=0.0, q=0.0, angle_radial=0.0, angle_radial_out=numpy.pi, angle_azimuthal=0.0):
        self._beamline_element_syned.get_coordinates()._p = p
        self._beamline_element_syned.get_coordinates()._q = q
        self._beamline_element_syned.get_coordinates()._angle_radial = angle_radial
        self._beamline_element_syned.get_coordinates()._angle_radial_out = angle_radial_out
        self._beamline_element_syned.get_coordinates()._angle_azimuthal = angle_azimuthal

    def get_positions(self):
        return self._beamline_element_syned.get_coordinates()._p, \
            self._beamline_element_syned.get_coordinates()._q, \
            self._beamline_element_syned.get_coordinates()._angle_radial, \
            self._beamline_element_syned.get_coordinates()._angle_radial_out, \
            self._beamline_element_syned.get_coordinates()._angle_azimuthal

    def trace_beam(self,beam1):
        beam = beam1.duplicate()

        p, q, angle_radial, angle_radial_out, angle_azimuthal = self.get_positions()

        theta_grazing1 = numpy.pi / 2 - angle_radial
        theta_grazing2 = numpy.pi / 2 - angle_radial_out
        alpha1 = angle_azimuthal

        #
        beam = beam1.duplicate()

        #
        # put beam in mirror reference system
        #
        beam.rotate(alpha1, axis=2)
        beam.rotate(theta_grazing1, axis=1)
        beam.translation([0.0, -p * numpy.cos(theta_grazing1), p * numpy.sin(theta_grazing1)])

        #
        # oe does nothing
        #
        pass

        #
        # from oe reference system to image plane
        #
        beam_out = beam.duplicate()
        beam_out.change_to_image_reference_system(theta_grazing2, q)


        return beam_out, beam

