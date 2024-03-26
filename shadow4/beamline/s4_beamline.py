"""
Defines the S4 beamline.

As a reminder, and following the SYNED philosophy, the S4 beamline is a container of the S4 light source and a
list of S4 beamline elements.

"""
from syned.beamline.beamline import Beamline
from shadow4.beamline.s4_beamline_element import S4BeamlineElement


class S4Beamline(Beamline):
    """
    Constructor.

    Parameters
    ----------
    light_source : instance of LightSource
        The light source
    beamline_elements_list : list
        The beamline elements (each one an instance of S4BeamlineElement).
    """
    def __init__(self,
                 light_source=None,
                 beamline_elements_list=None):
        super().__init__(light_source=light_source, beamline_elements_list=beamline_elements_list)

    def duplicate(self):
        """
        Returns a copy of the S4 beamline instance.

        Returns
        -------
        S4BeamlineElement  instance
            A copy of the object instance.

        """
        if self.get_beamline_elements_number() == 0:
            beamline_elements_list = None
        else:
            beamline_elements_list = []
            for beamline_element in self._beamline_elements_list:
                beamline_elements_list.append(beamline_element.duplicate())

        return S4Beamline(light_source=self._light_source,
                          beamline_elements_list = beamline_elements_list)

    def append_beamline_element(self, beamline_element: S4BeamlineElement):
        """
        Appends a S4 beamline element.

        Parameters
        ----------
        beamline_element : instance of S4BeamlineElement.
            The beamline element to append.
        """
        self._beamline_elements_list.append(beamline_element)

    def to_python_code(self, **kwargs):
        """
        Returns the python code to create the beamline.

        Parameters
        ----------
        **kwargs
            Passed arguments

        Returns
        -------
        str
            The python code.
        """
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
        """
        Runs (performs the ray tracing) of the full beamline.

        Parameters
        ----------
        **params
            Passed params.

        Returns
        -------
        tuple
            (output_beam, output_mirr) the traced beam and footprint (after the last beamline element).
        """
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
    from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical
    from shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirror, S4PlaneMirrorElement
    from syned.beamline.element_coordinates import ElementCoordinates

    light_source = SourceGeometrical(name='SourceGeometrical', nrays=10000, seed=5676561)

    m1 = S4PlaneMirror()
    m2 = S4PlaneMirror()

    e1 = S4PlaneMirrorElement(m1, ElementCoordinates())
    e2 = S4PlaneMirrorElement(m2, ElementCoordinates())

    bl = S4Beamline(light_source=light_source, beamline_elements_list=[e1,e2])

    print(bl.info())

    print(bl.to_python_code())

    output_beam, output_mirr = bl.run_beamline()

    # test plot
    from srxraylib.plot.gol import plot_scatter
    rays = output_beam.get_rays()
    plot_scatter(1e6 * rays[:, 3], 1e6 * rays[:, 5], title='(Xp,Zp) in microns')

