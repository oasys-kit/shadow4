from syned.storage_ring.light_source import LightSource

class S4LightSource(LightSource):

    def __init__(self,
                 name="Undefined",
                 electron_beam=None,
                 magnetic_structure=None,
                 nrays=5000,
                 seed=12345,
                 ):
        super().__init__(name=name, electron_beam=electron_beam, magnetic_structure=magnetic_structure)

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
        script = ''
        try:
            script += self.get_electron_beam().to_python_code()
        except:
            script += "\n\n#Error retrieving electron_beam code"

        try:
            script += self.get_magnetic_structure().to_python_code()
        except:
            script += "\n\n#Error retrieving magnetic structure code"


        script += "\n\n\nfrom shadow4.sources.s4_light_source import S4LightSource"
        script += "\nlight_source = S4LightSource(name='%s', electron_beam=electron_beam, magnetic_structure=source, nrays=%s, seed=%s)" % \
                                (self.get_name(), self.get_nrays(), self.get_seed())
        return script

    def get_beam(self, **params):
        raise NotImplementedError()

    def calculate_spectrum(self, **params):
        raise NotImplementedError()

if __name__ == "__main__":
    a = S4LightSource()
    print(a.info())
    print(a.to_python_code())
