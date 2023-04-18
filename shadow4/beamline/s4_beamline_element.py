from syned.beamline.beamline_element import BeamlineElement
from syned.beamline.optical_element import OpticalElement
from syned.beamline.element_coordinates import ElementCoordinates
from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements

class S4BeamlineElement(BeamlineElement):

    def __init__(self,
                 optical_element : OpticalElement = None,
                 coordinates : ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element, coordinates)
        self.__input_beam = input_beam
        self.__movements = movements

    def get_input_beam(self): return self.__input_beam
    def set_input_beam(self, input_beam): self.__input_beam = input_beam
    def get_movements(self): return self.__movements
    def set_movements(self, movements): self.__movements = movements
    def trace_beam(self, **params): raise NotImplementedError()
    def info(self): return self.get_optical_element().info() + "\n" + self.get_coordinates().info()
    def to_python_code(self, **kwargs): raise NotImplementedError()

    def to_python_code_movements(self):
        movements = self.get_movements()
        if isinstance(movements, S4BeamlineElementMovements):
            txt = "\nfrom shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements"
            txt += "\nmovements = S4BeamlineElementMovements(f_move=%d, offset_x=%g, offset_y=%g, offset_z=%g, rotation_x=%g, rotation_y=%g, rotation_z=%g)" % \
                   (movements.f_move, movements.offset_x, movements.offset_y, movements.offset_z,
                    movements.rotation_x, movements.rotation_y, movements.rotation_z)
        else:
            txt = "\nmovements = None"

        return txt
