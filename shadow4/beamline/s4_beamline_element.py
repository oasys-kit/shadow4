from syned.beamline.beamline_element import BeamlineElement
from syned.beamline.optical_element import OpticalElement
from syned.beamline.element_coordinates import ElementCoordinates
from shadow4.beamline.s4_optical_element import S4OpticalElementDecorator
from shadow4.beam.s4_beam import S4Beam

import copy
class S4BeamlineElement(BeamlineElement):

    def __init__(self,
                 optical_element : OpticalElement = None,
                 coordinates : ElementCoordinates = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element, coordinates)
        self.__input_beam = input_beam

    def get_input_beam(self):
        return self.__input_beam

    def set_input_beam(self, input_beam):
        self.__input_beam = input_beam

    def trace_beam(self, **params):
        raise NotImplementedError()

    def info(self):
        return self.get_optical_element().info() + "\n" + self.get_coordinates().info()

    def to_python_code(self, data=None):
        raise NotImplementedError()

    def duplicate(self):
        return S4BeamlineElement(optical_element=copy.deepcopy(self.get_optical_element()),
                                 coordinates=copy.deepcopy(self.get_coordinates()),
                                 input_beam=self.get_input_beam().duplicate())
