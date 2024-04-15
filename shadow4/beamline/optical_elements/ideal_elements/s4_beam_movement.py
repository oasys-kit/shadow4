from syned.beamline.optical_element import OpticalElement
from syned.beamline.element_coordinates import ElementCoordinates

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.s4_optical_element_decorators import S4OpticalElementDecorator
from shadow4.beamline.s4_beamline_element import S4BeamlineElement
from shadow4.beam.s4_beam import S4Beam
from shadow4.tools.logger import is_verbose, is_debug

class S4BeamMovement(OpticalElement, S4OpticalElementDecorator):
    def __init__(self, name="Undefined",
                 apply_flag=0, # 0=movement off, 1=movement on
                 translation_x=0.0,
                 translation_y=0.0,
                 translation_z=0.0,
                 rotation_x=0.0,
                 rotation_y=0.0,
                 rotation_z=0.0,
                 ):
        super().__init__(name=name)
        self.apply_flag    = apply_flag
        self.translation_x = translation_x
        self.translation_y = translation_y
        self.translation_z = translation_z
        self.rotation_x    = rotation_x
        self.rotation_y    = rotation_y
        self.rotation_z    = rotation_z

        self.__inputs = {
            "name": name,
            "apply_flag": apply_flag,
            "translation_x": translation_x,
            "translation_y": translation_y,
            "translation_z": translation_z,
            "rotation_x": rotation_x,
            "rotation_y": rotation_y,
            "rotation_z": rotation_z,
        }

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

from shadow4.beamline.optical_elements.ideal_elements.s4_beam_movement import S4BeamMovement
optical_element = S4BeamMovement(name='{name:s}', apply_flag={apply_flag:d},
    translation_x={translation_x:g}, translation_y={translation_y:g}, translation_z={translation_z:g},
    rotation_x={rotation_x:g}, rotation_y={rotation_y:g}, rotation_z={rotation_z:g})
"""
        txt = txt_pre.format(**self.__inputs)
        return txt

class S4BeamMovementElement(S4BeamlineElement):
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
        instance of S4BeamMovementElement.

    """
    def __init__(self,
                 optical_element : S4BeamMovement = None,
                 coordinates : ElementCoordinates = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4BeamMovement(),
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
        verbose = 0
        #
        input_beam = self.get_input_beam().duplicate()

        #
        # apply movement to beam
        #
        oe = self.get_optical_element()

        if oe.apply_flag:
            if is_verbose():
                print("Beam translated [%f,%f,%f] um" % (1e6 * oe.translation_x,
                                                         1e6 * oe.translation_y,
                                                         1e6 * oe.translation_z, ))
            input_beam.translation([oe.translation_x, oe.translation_y, oe.translation_z])
            if oe.rotation_x != 0.0:
                if is_verbose(): print("Beam rotated %f urad along X" % (1e6 * self.rotation_x))
                input_beam.rotate(oe.rotation_x, axis=1, rad=True)
            if oe.rotation_y != 0.0:
                if is_verbose(): print("Beam rotated %f urad along Y" % (1e6 * oe.rotation_y))
                input_beam.rotate(oe.rotation_y, axis=2, rad=True)
            if oe.rotation_z != 0.0:
                if is_verbose(): print("Beam rotated %f urad along Z" % (1e6 * oe.rotation_z))
                input_beam.rotate(oe.rotation_z, axis=3, rad=True)

        #
        # from oe reference system to image plane  (nothing to do...)
        #
        output_beam = input_beam.duplicate()
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
        txt += "\nfrom shadow4.beamline.optical_elements.ideal_elements.s4_beam_movement import S4BeamMovementElement"
        txt += "\nbeamline_element = S4BeamMovementElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)"
        txt += "\n\nbeam, mirr = beamline_element.trace_beam()"
        return txt

if __name__ == "__main__":
    # a = S4BeamMovement()
    # print(a.to_python_code())

    b = S4BeamMovementElement()
    print(b.to_python_code())