from shadow4.syned.beamline_element import BeamlineElement # todo: change when syned upgraded
# from syned.beamline.beamline_element import BeamlineElement


class S4BeamlineElement(BeamlineElement):

    def __init__(self, optical_element=None, coordinates=None):
        super().__init__(optical_element, coordinates)

    def trace_beam(self, beam=None, **params):
        raise NotImplementedError()

    def info(self):
        return self.get_optical_element().info() + "\n" + self.get_coordinates().info()


