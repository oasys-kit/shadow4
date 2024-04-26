# TODO: move this to syned

from syned.beamline.shape import SurfaceShape
from syned.beamline.optical_element_with_surface_shape import OpticalElementsWithSurfaceShape

class Multilayer(OpticalElementsWithSurfaceShape):
    """
    Constructor.

    Parameters
    ----------
    name : str
        The name of the optical element.
    surface_shape : instance of SurfaceShape, optional
        The geometry of the crystal surface. if None, it is initialized to SurfaceShape().
    boundary_shape : instance of BoundaryShape, optional
        The geometry of the slit aperture. if None, it is initialized to BoundaryShape().
    structure : str, optional
        The multilayer structure e.g. [B,W]x50+Si.
    period : float, optional
        The period of the repeated bilayer in A.
    Gamma : float, optional
        The gamma factor.
    """
    def __init__(self,
                 name="Undefined",
                 surface_shape=SurfaceShape(),
                 boundary_shape=None,
                 structure='[B/W]x50+Si',
                 period=25.0,
                 Gamma=0.5,
                 ):

        super().__init__(name, surface_shape, boundary_shape)
        self._structure = structure
        self._period = period
        self._Gamma = Gamma
        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("name",                "Name" ,                                          "" ),
                    ("surface_shape",       "Surface shape",                                  "" ),
                    ("boundary_shape",      "Boundary shape",                                 "" ),
                    ("structure",           "structure ([Odd,Even]xN+Substrate)",             "" ),
                    ("period",              "period of the repeated structure",               "A"),
                    ("Gamma",               "Gamma factor [thickness ratio Even)/(Odd+Even)", ""),
            ] )

