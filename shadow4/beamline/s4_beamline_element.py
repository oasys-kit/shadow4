from syned.beamline.beamline_element import BeamlineElement
from syned.beamline.optical_element import OpticalElement
from syned.beamline.element_coordinates import ElementCoordinates
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

    def duplicate_input_beam(self):
        if self.get_input_beam() is None:
            input_beam = None
        else:
            input_beam = self.get_input_beam().duplicate()
            return input_beam

    def duplicate_optical_element(self):
        optical_element = copy.deepcopy(self.get_optical_element())

    def duplicate_coordinates(self):
        copy.deepcopy(self.get_coordinates())

    def duplicate(self):
        return S4BeamlineElement(optical_element=self.duplicate_coordinates(),
                                coordinates=self.duplicate_coordinates(),
                                input_beam=self.duplicate_input_beam())

