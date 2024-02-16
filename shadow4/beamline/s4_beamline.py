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
            beamline_elements_list.append(beamline_element.duplicate())

        return S4Beamline(light_source=self._light_source,
                          beamline_elements_list = beamline_elements_list)

    def append_beamline_element(self, beamline_element: S4BeamlineElement):
        self._beamline_elements_list.append(beamline_element)

    def to_python_code(self, **kwargs):
        script = "from shadow4.beamline.s4_beamline import S4Beamline"
        script += "\n\nbeamline = S4Beamline()\n"
        try:
            script += self.get_light_source().to_python_code()
            script += "\n\nbeamline.set_light_source(light_source)"
        except:
            script +=  "\n\n\n# Error getting python code for S4Beamline S4LightSource "

        for i,element in enumerate(self.get_beamline_elements()):
            try:
                script += element.to_python_code()
                script += "\n\nbeamline.append_beamline_element(beamline_element)"
            except:
                script += "\n\n\n# Error getting python code for S4Beamline S4BeamlineElement # %d  :" % (i+1)
                script += "\n#       %s " % (str(element))

        return script

    def run_beamline(self, **params):
        try:
            output_beam = self.get_light_source().get_beam(**params)
            output_mirr = None
        except:
            raise Exception("Error running beamline light source")

        for i, element in enumerate(self.get_beamline_elements()):
            try:
                element.set_input_beam(output_beam)
                output_beam, output_mirr = element.trace_beam(**params)
            except:
                raise Exception("Error running beamline element # %d" % (i+1) )

        return output_beam, output_mirr


if __name__ == "__main__":
    from shadow4.beamline.optical_elements.mirrors.s4_mirror import S4Mirror, S4MirrorElement
    from syned.beamline.element_coordinates import ElementCoordinates

    m1 = S4Mirror()
    m2 = S4Mirror()

    e1 = S4MirrorElement(m1, ElementCoordinates())
    e2 = S4MirrorElement(m2, ElementCoordinates())

    bl = S4Beamline(beamline_elements_list=[e1,e2])

    print(bl.info())

    print(bl.to_python_code())

    # print(bl.run_beamline())

