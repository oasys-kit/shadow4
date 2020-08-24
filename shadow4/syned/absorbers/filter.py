
from syned.beamline.optical_elements.absorbers.absorber import Absorber

class Filter(Absorber):
    def __init__(self,
                 name="Undefined",
                 material="Be",
                 thickness=1e-3,
                 boundary_shape=None):
        Absorber.__init__(self, name=name, boundary_shape=boundary_shape)
        self._material = material
        self._thickness = thickness


        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("material"      , "Material (symbol, formula or name)",    "" ),
                    ("thickness"     , "Thickness ",                            "m" ),
            ] )
    def get_material(self):
        return self._material

    def get_thickness(self):
        return self._thickness
