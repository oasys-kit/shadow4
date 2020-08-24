
from syned.beamline.shape import BoundaryShape
from syned.beamline.shape import Rectangle, Ellipse, Circle

from syned.beamline.optical_elements.absorbers.absorber import Absorber

class Slit(Absorber):
    def __init__(self, name="Undefined", boundary_shape=None):
        if boundary_shape is None:
            boundary_shape = BoundaryShape()
        Absorber.__init__(self, name=name, boundary_shape=boundary_shape)

    def set_rectangle(self,width=3e-3,height=4e-3,center_x=0.0,center_y=0.0):
        self._boundary_shape=Rectangle(-0.5*width+center_x,0.5*width+center_x,-0.5*height+center_y,0.5*height+center_y)

    def set_circle(self,radius=3e-3,center_x=0.0,center_y=0.0):
        self._boundary_shape=Circle(radius,center_x,center_y)

    def set_ellipse(self,width=3e-3,height=4e-3,center_x=0.0,center_y=0.0):
        self._boundary_shape=Ellipse(-0.5*width+center_x,0.5*width+center_x,-0.5*height+center_y,0.5*height+center_y)