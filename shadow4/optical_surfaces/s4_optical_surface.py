
# abstract class defining the interfaces of the optical surfaces to be impremented in the clildren classes


class S4OpticalSurface(object):

    def info(self):
        raise NotImplementedError

    def duplicate(self):
        raise NotImplementedError

    def surface_height(self, x, y, **kwargs):
        raise NotImplementedError

    def get_normal(self, x, **kwargs):
        raise NotImplementedError

    def calculate_intercept(self, XIN, VIN, **kwargs):
        raise NotImplementedError

    def apply_specular_reflection_on_beam(self, beam, **kwargs):
        raise NotImplementedError

    def apply_refraction_on_beam(self, beam, **kwargs):
        raise NotImplementedError

    def apply_crystal_diffraction_bragg_symmetric_on_beam(self, beam, **kwargs):
        raise NotImplementedError



if __name__ == "__main__":
    a = S4OpticalSurface()
    a.info()