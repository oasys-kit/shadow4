from syned.beamline.beamline import Beamline as SynedBeamline
from syned.storage_ring.light_source import LightSource as SynedLightSource




class Beamline(SynedBeamline):

    def __init__(self,
                 light_source=None,
                 beamline_elements_list=[]):
        if light_source is None: light_source=SynedLightSource()
        super().__init__(light_source=light_source, beamline_elements_list=beamline_elements_list)


    def to_python_code(self, data=None):
        text_code  =  "# to be implemented "
        return text_code

    def info(self):
        return "Beamline info: to be implemented"

    def run_beamline(self):
        raise NotImplementedError()



if __name__ == "__main__":
    from shadow4.optical_elements.mirror import Mirror
    from syned.beamline.beamline_element import BeamlineElement
    from syned.beamline.element_coordinates import ElementCoordinates

    m1 = Mirror()
    m2 = Mirror()

    e1 = BeamlineElement(m1, ElementCoordinates())
    e2 = BeamlineElement(m2, ElementCoordinates())


    bl = Beamline(beamline_elements_list=[e1,e2])


