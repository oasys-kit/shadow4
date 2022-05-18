from syned.beamline.beamline import Beamline
from shadow4.beamline.s4_beamline_element import S4BeamlineElement


class S4Beamline(Beamline):

    def __init__(self,
                 light_source=None,
                 beamline_elements_list=[]):
        super().__init__(light_source=light_source, beamline_elements_list=beamline_elements_list)


    def to_python_code(self, data=None):
        text_code  =  "# to be implemented "
        return text_code

    def info(self):
        return "Beamline info: to be implemented"

    def run_beamline(self, **params):
        raise NotImplementedError()

    def duplicate(self):
        beamline_elements_list = []
        for beamline_element in self._beamline_elements_list:
            beamline_elements_list.append(beamline_element)

        return S4Beamline(light_source=self._light_source,
                        beamline_elements_list = beamline_elements_list)

    # copied from syned. May be removed when updating syned

    def append_beamline_element(self, beamline_element=S4BeamlineElement()):
        if not isinstance(beamline_element,S4BeamlineElement):
            raise Exception("Input class must be of type: "+S4BeamlineElement.__name__)
        else:
            self._beamline_elements_list.append(beamline_element)




if __name__ == "__main__":
    from shadow4.beamline.optical_elements.mirrors import S4Mirror, S4MirrorElement
    from shadow4.syned.element_coordinates import ElementCoordinates

    m1 = S4Mirror()
    m2 = S4Mirror()

    e1 = S4MirrorElement(m1, ElementCoordinates())
    e2 = S4MirrorElement(m2, ElementCoordinates())

    bl = S4Beamline(beamline_elements_list=[e1,e2])


