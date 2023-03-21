from syned.storage_ring.empty_light_source import EmptyLightSource

#

# this is an abstract class to support non synchrotron light sources (e.g. geometrical source)

class S4LightSourceBase(EmptyLightSource):

    def __init__(self, name="Undefined", nrays=5000, seed=1234567):
        super().__init__(name=name)
        self.__nrays = nrays
        self.__seed = seed

    def set_nrays(self, nrays):
        self.__nrays = nrays

    def get_nrays(self):
        return self.__nrays

    def set_seed(self, seed):
        self.__seed = seed

    def get_seed(self):
        return self.__seed

    def to_python_code(self, **kwargs):
        raise NotImplementedError()

    def get_beam(self, **params):
        raise NotImplementedError()

    def calculate_spectrum(self, **params):
        raise NotImplementedError()




if __name__ == "__main__":
    a = S4LightSourceBase()
    print(a.info())
    # print(a.to_python_code())

    from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical
    b = SourceGeometrical(seed=0)

