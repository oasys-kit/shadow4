#
# the screen optical element:
#       deals with screens, slits, beam-stoppers and absorbers (as in shadow3)
#       it is a stand-alone optical element (contrary to shadow4)
#
#
#
import numpy
from collections import OrderedDict
from shadow4.compatibility.beam3 import Beam3

from syned.beamline.optical_elements.ideal_elements.screen import Screen as SyScreen
from syned.beamline.optical_elements.absorbers.absorber import Absorber as SyAbsorber
from syned.beamline.optical_elements.absorbers.beam_stopper import BeamStopper as SyBeamStopper
from syned.beamline.optical_elements.absorbers.filter import Filter as SyFilter
from syned.beamline.optical_elements.absorbers.slit import Slit as SySlit

from syned.beamline.beamline_element import BeamlineElement
from syned.beamline.element_coordinates import ElementCoordinates


class Screen(object):
    def __init__(self, beamline_element_syned = None):
        if beamline_element_syned is None:
            self._beamline_element_syned = BeamlineElement(
                SyScreen(name="Undefined"),
                ElementCoordinates(p=0.0, q=0.0, angle_radial=0.0, angle_azimuthal=0.0))
        else:
            ok = False
            for obj in [SyScreen, SyAbsorber, SyBeamStopper, SyFilter, SySlit]:
                if isinstance(beamline_element_syned._optical_element, obj): ok = True
            if ok:
                self._beamline_element_syned = beamline_element_syned
            else:
                raise Exception("Please initialize shadow4 Screen with syned Screen, Absorber, BeamStopper, Filter or Slit")


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

    def trace_beam(self,beam1):
        beam = beam1.duplicate()

        p,q = self.get_positions()

        if p+q != 0.0:
            beam.retrace(p+q,resetY=True)

        #TODO: slit, absorber, etc

        return beam



