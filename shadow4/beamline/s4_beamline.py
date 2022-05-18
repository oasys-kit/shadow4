from syned.beamline.beamline import Beamline
from shadow4.beamline.s4_beamline_element import S4BeamlineElement


class S4Beamline(Beamline):

    def __init__(self,
                 light_source=None,
                 beamline_elements_list=[]):
        super().__init__(light_source=light_source, beamline_elements_list=beamline_elements_list)




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


    def to_python_code(self, data=None):
        script = ''
        try:
            script += self.get_light_source().to_python_code()
        except:
            script +=  "\n\n\n# Error getting python code for S4Beamline S4LightSource "

        for i,element in enumerate(self.get_beamline_elements()):
            try:
                script += element.to_python_code()
            except:
                script += "\n\n\n# Error getting python code for S4Beamline S4BeamlineElement # %d " % (i+1)

        return script



    def run_beamline(self, **params):
        try:
            beam0 = self.get_light_source().get_beam(**params)
            mirr0 = None
        except:
            raise Exception("Error running beamline light source")

        for i,element in enumerate(self.get_beamline_elements()):
            try:
                beam1, mirr1 = element.trace_beam(beam_in=beam0, **params)
            except:
                raise Exception("Error running beamline element # %d" % (i+1) )

            beam0 = beam1
            mirr0 = mirr1

        return beam0, mirr0


if __name__ == "__main__":
    from shadow4.beamline.optical_elements.mirrors.s4_mirror import S4Mirror, S4MirrorElement
    from shadow4.syned.element_coordinates import ElementCoordinates

    m1 = S4Mirror()
    m2 = S4Mirror()

    e1 = S4MirrorElement(m1, ElementCoordinates())
    e2 = S4MirrorElement(m2, ElementCoordinates())

    bl = S4Beamline(beamline_elements_list=[e1,e2])

    print(bl.info())

    print(bl.to_python_code())

    # print(bl.run_beamline())

