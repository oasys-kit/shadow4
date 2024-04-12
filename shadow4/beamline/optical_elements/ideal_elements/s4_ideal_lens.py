import numpy

from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.optical_elements.ideal_elements.ideal_lens import IdealLens

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.s4_optical_element_decorators import S4OpticalElementDecorator
from shadow4.beamline.s4_beamline_element import S4BeamlineElement


class S4IdealLens(IdealLens, S4OpticalElementDecorator):
    def __init__(self, name="Undefined",focal_x=0.0,focal_y=0.0):
        super().__init__(name=name,focal_x=focal_x,focal_y=focal_y)

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

from shadow4.beamline.optical_elements.ideal_elements.s4_ideal_lens import S4IdealLens
optical_element = S4IdealLens(name='{name:s}', focal_x={focal_x:g},focal_y={focal_x:g})
"""
        txt = txt_pre.format(**{'name':self.get_name(), 'focal_x': self._focal_x, 'focal_y': self._focal_y})
        return txt

class S4IdealLensElement(S4BeamlineElement):
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
    instance of S4IdealLensElement.
    """
    def __init__(self,
                 optical_element : S4IdealLens = None,
                 coordinates : ElementCoordinates = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4IdealLens(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         input_beam=input_beam)

    def get_focalX(self):
        return self.get_optical_element()._focal_x

    def get_focalZ(self):
        return self.get_optical_element()._focal_y

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
        footprint = self.get_input_beam().duplicate()

        p, q = self.get_coordinates().get_p_and_q()
        if p != 0.0: footprint.retrace(p,resetY=True)

        output_beam = footprint.duplicate()

        # rotate around Z
        if self.get_focalX() != 0.0:
            # whatch out the minus!!
            tan_two_theta = - output_beam.get_column(1) / self.get_focalX()
            output_beam.rotate(numpy.arctan(tan_two_theta),axis=3,rad=True)

        # rotate around X
        if self.get_focalZ() != 0.0:
            tan_two_theta = output_beam.get_column(3) / self.get_focalZ()
            output_beam.rotate(numpy.arctan(tan_two_theta),axis=1,rad=True)

        if q != 0.0:
            output_beam.retrace(q,resetY=True)

        return output_beam, footprint

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
        txt += "\nfrom shadow4.beamline.optical_elements.ideal_elements.s4_ideal_lens import S4IdealLensElement"
        txt += "\nbeamline_element = S4IdealLensElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)"
        txt += "\n\nbeam, mirr = beamline_element.trace_beam()"
        return txt

class S4SuperIdealLens(IdealLens, S4OpticalElementDecorator):
    def __init__(self, name="Undefined",focal_p_x=1.0, focal_p_y=1.0,
                                        focal_q_x=1.0, focal_q_y=1.0,):
        super().__init__(name=name,focal_x=1.0/(1/focal_p_x+1/focal_q_x),focal_y=1.0/(1/focal_p_y+1/focal_q_y))

        self._focal_p_x = focal_p_x
        self._focal_p_y = focal_p_y
        self._focal_q_x = focal_q_x
        self._focal_q_y = focal_q_y

        self.__inputs = {
            "name": name,
            "focal_p_x": focal_p_x,
            "focal_q_x": focal_q_x,
            "focal_p_y": focal_p_y,
            "focal_q_y": focal_q_y,
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

from shadow4.beamline.optical_elements.ideal_elements.s4_ideal_lens import S4SuperIdealLens
optical_element = S4SuperIdealLens(name='{name:s}', 
    focal_p_x={focal_p_x:g}, focal_q_x={focal_q_x:g},
    focal_p_y={focal_p_y:g}, focal_q_y={focal_q_y:g})
"""
        txt = txt_pre.format(**self.__inputs)
        return txt

class S4SuperIdealLensElement(S4BeamlineElement):
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
    instance of S4SuperIdealLensElement.
    """
    def __init__(self,
                 optical_element : S4SuperIdealLens = None,
                 coordinates : ElementCoordinates = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4SuperIdealLens(),
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
        footprint = self.get_input_beam().duplicate()

        p, q = self.get_coordinates().get_p_and_q()
        if p != 0.0: footprint.retrace(p,resetY=True)

        lens = self.get_optical_element()

        output_beam = footprint.duplicate()

        # rotate around Z; watch out the minus!!
        if lens._focal_p_x != 0.0:
            tan_theta_p = - output_beam.get_column(1) / lens._focal_p_x
        else:
            tan_theta_p = 0.0

        if lens._focal_q_x != 0.0:
            tan_theta_q = - output_beam.get_column(1) / lens._focal_q_x
        else:
            tan_theta_q = 0.0

        two_theta = numpy.arctan(tan_theta_p) + numpy.arctan(tan_theta_q)
        output_beam.rotate(two_theta,axis=3,rad=True)

        # rotate around X
        if lens._focal_p_y != 0.0:
            tan_theta_p = output_beam.get_column(3) / lens._focal_p_y
        else:
            tan_theta_p = 0.0

        if lens._focal_q_y != 0.0:
            tan_theta_q = output_beam.get_column(3) / lens._focal_q_y
        else:
            tan_theta_q = 0.0

        two_theta = numpy.arctan(tan_theta_p) + numpy.arctan(tan_theta_q)
        output_beam.rotate(two_theta,axis=1,rad=True)

        if q != 0.0:
            output_beam.retrace(q,resetY=True)

        return output_beam, footprint

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
        txt += "\nfrom shadow4.beamline.optical_elements.ideal_elements.s4_ideal_lens import S4IdealLensElement"
        txt += "\nbeamline_element = S4IdealLensElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)"
        txt += "\n\nbeam, mirr = beamline_element.trace_beam()"
        return txt

if __name__ == "__main__":
    # a = S4IdealLens()
    # print(a.to_python_code())
    # a = S4SuperIdealLens()
    # print(a.to_python_code())
    b = S4IdealLensElement()
    print(b.to_python_code())
    # b = S4SuperIdealLensElement()
    # print(b.to_python_code())