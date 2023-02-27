
from syned.beamline.optical_element import OpticalElement
from syned.beamline.shape import BoundaryShape

class Absorber(OpticalElement):
    def __init__(self, name="Undefined", boundary_shape=None):
        if boundary_shape is None: boundary_shape = BoundaryShape()
        OpticalElement.__init__(self, name=name, boundary_shape=boundary_shape)
