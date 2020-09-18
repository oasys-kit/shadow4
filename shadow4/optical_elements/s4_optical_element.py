

class S4OpticalElement(object):

    def __init__(self):
        pass

    def info(self):
        raise NotImplementedError()

    def to_python_code(self, data=None):
        raise NotImplementedError()


class SurfaceCalculation:
    INTERNAL = 0
    EXTERNAL = 1

class S4CurvedOpticalElement(object):

    def __init__(self,
                 surface_calculation=SurfaceCalculation.INTERNAL,
                 is_cylinder=False,
                 ):
        self._surface_calculation = surface_calculation
        self._is_cylinder = is_cylinder

class S4PlaneOpticalElement(object):
    pass

class S4EllipsoidOpticalElement(S4CurvedOpticalElement):
    pass

class S4HyperboloidOpticalElement(S4CurvedOpticalElement):
    pass

class S4SphereOpticalElement(S4CurvedOpticalElement):
    pass

class S4ToroidalOpticalElement(S4CurvedOpticalElement):
    pass

class S4ConicOpticalElement(S4CurvedOpticalElement):
    pass
