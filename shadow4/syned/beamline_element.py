"""
Base class for all beamline components.
Enforce to name every every component.
Every beamline component can store settings
"""

from syned.syned_object import SynedObject
from syned.beamline.optical_element import OpticalElement
# from syned.beamline.element_coordinates import ElementCoordinates
from shadow4.syned.element_coordinates import ElementCoordinates # TODO: update syned

class BeamlineElement(SynedObject):
    def __init__(self, optical_element=OpticalElement(), coordinates=ElementCoordinates()):
        self._optical_element = optical_element
        self._coordinates = coordinates
        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("optical_element",       "Optical Element",      ""),
                    ("coordinates",           "Element coordinates",  ""),
            ] )

    def get_optical_element(self):
        return self._optical_element

    def get_coordinates(self):
        return self._coordinates

    def set_optical_element(self, value):
        if isinstance(value, OpticalElement):
            self._optical_element = value
        else:
            raise Exception("entry is not an instance of OpticalElement")

    def set_coordinates(self, value):
        if isinstance(value, ElementCoordinates):
            self._coordinates = value
        else:
            raise Exception("entry is not an instance of ElementCoordinates")