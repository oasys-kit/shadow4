import numpy
from syned.beamline.optical_elements.ideal_elements.screen import Screen
from syned.beamline.element_coordinates import ElementCoordinates

from shadow4.beamline.s4_optical_element_decorators import S4OpticalElementDecorator
from shadow4.beamline.s4_beamline_element import S4BeamlineElement
from shadow4.beam.s4_beam import S4Beam

class S4Empty(Screen, S4OpticalElementDecorator):
    def __init__(self, name="Undefined"):
        super().__init__(name=name)

    def get_info(self):
        """
        Returns the specific information of the S4 empty optical element.

        Returns
        -------
        str
        """
        txt = "\n\n"
        txt += "EMPTY\n"

        return txt

    def to_python_code(self, **kwargs):
        """
        Creates the python code for defining the element.

        Parameters
        ----------
        **kwargs

        Returns
        -------
        str
            Python code.
        """
        txt_pre = """

from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4Empty
optical_element = S4Empty(name='{name:s}')
"""
        txt = txt_pre.format(**{'name': self.get_name()})
        return txt

class S4EmptyElement(S4BeamlineElement):
    """
    Constructor.

    Parameters
    ----------
    optical_element : instance of OpticalElement, optional
        The syned optical element.
    coordinates : instance of ElementCoordinates, optional
        The syned element coordinates.
    input_beam : instance of S4Beam, optional
        The S4 incident beam.

    Returns
    -------
    instance of S4EmptyElement.

    """
    def __init__(self,
                 optical_element : S4Empty = None,
                 coordinates : ElementCoordinates = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4Empty(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         input_beam=input_beam)

    def trace_beam(self, **params):
        """
        Runs (ray tracing) the input beam through the element.

        Parameters
        ----------
        **params

        Returns
        -------
        tuple
            (output_beam, footprint) instances of S4Beam.
        """
        p, q, angle_radial, angle_radial_out, angle_azimuthal = self.get_coordinates().get_positions()

        theta_grazing1 = numpy.pi / 2 - angle_radial
        theta_grazing2 = numpy.pi / 2 - angle_radial_out
        alpha1 = angle_azimuthal

        #
        input_beam = self.get_input_beam().duplicate()

        #
        # put beam in mirror reference system
        #
        input_beam.rotate(alpha1, axis=2)
        input_beam.rotate(theta_grazing1, axis=1)
        input_beam.translation([0.0, -p * numpy.cos(theta_grazing1), p * numpy.sin(theta_grazing1)])

        #
        # oe does nothing
        #
        pass

        #
        # from oe reference system to image plane
        #
        output_beam = input_beam.duplicate()
        output_beam.change_to_image_reference_system(theta_grazing2, q)


        return output_beam, input_beam

    def to_python_code(self, **kwargs):
        """
        Creates the python code for defining the element.

        Parameters
        ----------
        **kwargs

        Returns
        -------
        str
            Python code.
        """
        txt = "\n\n# optical element number XX"
        txt += self.get_optical_element().to_python_code()
        txt += self.to_python_code_coordinates()
        txt += "\nfrom shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4EmptyElement"
        txt += "\nbeamline_element = S4EmptyElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)"
        txt += "\n\nbeam, mirr = beamline_element.trace_beam()"
        return txt

if __name__ == "__main__":
    # a = S4Empty()
    # print(a.to_python_code())
    b = S4EmptyElement()
    print(b.to_python_code())