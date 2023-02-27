
from syned.beamline.shape import BoundaryShape
from syned.beamline.shape import Rectangle, Ellipse

from syned.beamline.optical_elements.absorbers.absorber import Absorber

class BeamStopper(Absorber):
    def __init__(self, name="Undefined", boundary_shape=None):
        if boundary_shape is None: boundary_shape = BoundaryShape()
        Absorber.__init__(self, name=name, boundary_shape=boundary_shape)

    def set_rectangle(self,width=3e-3,height=4e-3):
        self._boundary_shape=Rectangle(-0.5*width,0.5*width,-0.5*height,0.5*height)

    def set_circle(self,radius=3e-3):
        self._boundary_shape=Ellipse(-0.5*radius,0.5*radius,-0.5*radius,0.5*radius)
