"""
Defines the S4 beamline element.

The S4 beamline element is a container with:
- The (SYNED) optical element (OpticalElement)
- The (SYNED) element coordinates (ElementCoordinates)
- The S4 element movements (S4BeamlineElementMovements)
- The S4 beam to be received by te element (S4Beam)

"""
from syned.beamline.beamline_element import BeamlineElement
from syned.beamline.optical_element import OpticalElement
from syned.beamline.element_coordinates import ElementCoordinates
from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements

class S4BeamlineElement(BeamlineElement):
    """
    Constructor.

    Parameters
    ----------
    optical_element : instance of OpticalElement, optional
        The syned optical element.
    coordinates : instance of ElementCoordinates, optional
        The syned element coordinates.
    movements : instance of S4BeamlineElementMovements, optional
        The S4 element movements.
    input_beam : instance of S4Beam, optional
        The S4 incident beam.
    """
    def __init__(self,
                 optical_element : OpticalElement = None,
                 coordinates : ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element, coordinates)
        self.__input_beam = input_beam
        self.__movements = movements

    def get_input_beam(self):
        """
        Returns the input beam.

        Returns
        -------
        instance of S4Beam
        """
        return self.__input_beam

    def set_input_beam(self, input_beam):
        """
        Sets the input beam.
        Parameters
        ----------
        input_beam : instance of S4Beam
            The input beam.
        """
        self.__input_beam = input_beam

    def get_movements(self):
        """
        Returns the element movements.

        Returns
        -------
        instance of S4BeamlineElementMovements
        """
        return self.__movements

    def set_movements(self, movements):
        """
        Sets the element movements.

        Parameters
        ----------
        movements : instance of S4BeamlineElementMovements
            The descriptod of the movements.
        """
        self.__movements = movements

    def trace_beam(self, **params):
        """
        Performs the ray tracing (trace). To be implemented in the derived classes.

        Parameters
        ----------
        **params

        Raises
        ------
        NotImplementedError
        """
        raise NotImplementedError()

    def info(self):
        """
        Gets the information text (syned doc).

        Returns
        -------
        str
            A info text for the S4 beamline element.
        """
        return self.get_optical_element().info() + "\n" + self.get_coordinates().info()

    def to_python_code(self, **kwargs):
        """
        Returns the python code to describe and run th eelement. To be implemented in the derived classes.

        Parameters
        ----------
        **kwargs
            The arguments to pass.

        Raises
        ------
        NotImplementedError
        """
        raise NotImplementedError()

    def to_python_code_movements(self):
        """
        Returns a code block to implement the element movements (to be used by the upper classes).

        Returns
        -------
        str
            A code block.
        """
        movements = self.get_movements()
        if isinstance(movements, S4BeamlineElementMovements):
            txt = "\nfrom shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements"
            txt += "\nmovements = S4BeamlineElementMovements(f_move=%d, offset_x=%g, offset_y=%g, offset_z=%g, rotation_x=%g, rotation_y=%g, rotation_z=%g)" % \
                   (movements.f_move, movements.offset_x, movements.offset_y, movements.offset_z,
                    movements.rotation_x, movements.rotation_y, movements.rotation_z)
        else:
            txt = "\nmovements = None"

        return txt

    def to_python_code_coordinates(self):
        """
        Returns a code block to implement the element coordinates (to be used by the upper classes).

        Returns
        -------
        str
            A code block.
        """
        coordinates = self.get_coordinates()
        txt = "\nfrom syned.beamline.element_coordinates import ElementCoordinates"
        txt += "\ncoordinates = ElementCoordinates(p=%.10g, q=%.10g, angle_radial=%.10g, angle_azimuthal=%.10g, angle_radial_out=%.10g)" % \
               (coordinates.p(), coordinates.q(), coordinates.angle_radial(), coordinates.angle_azimuthal(), coordinates.angle_radial_out())
        return txt
