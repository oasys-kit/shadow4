"""
Defines the a BaseLightSource to support non-synchrotron sources.
"""
from syned.storage_ring.empty_light_source import EmptyLightSource

class S4LightSourceBase(EmptyLightSource):
    """
    Abstract class to support non-synchrotron light sources (e.g. geometrical source)

    Parameters
    ----------
    name : str
        A name.
    nrays : int
        The number of rays.
    seed : int
        The Monte Carlo seed.

    """
    def __init__(self, name="Undefined", nrays=5000, seed=1234567):
        super().__init__(name=name)
        self.__nrays = nrays
        self.__seed = seed

    def set_nrays(self, nrays):
        """
        Defines the number of rays to be used.

        Parameters
        ----------
        nrays : int
            The number of rays.

        """
        self.__nrays = nrays

    def get_nrays(self):
        """
        Returns the stored number of rays.

        Returns
        -------
        int
            The number of rays.

        """
        return self.__nrays

    def set_seed(self, seed):
        """
        Defines the Monte Carlo seed.

        Parameters
        ----------
        seed : int
            The seed.

        """
        self.__seed = seed

    def get_seed(self):
        """
        Returns the stored seed.

        Returns
        -------
        int
            The seed.
        """
        return self.__seed

    def to_python_code(self, **kwargs):
        """
        To be implemented in a derived class.

        Raises
        ------
        NotImplementedError()

        """
        raise NotImplementedError()

    def get_beam(self, **params):
        raise NotImplementedError()

    def calculate_spectrum(self, **params):
        """
        To be implemented in a derived class.

        Raises
        ------
        NotImplementedError()

        """
        raise NotImplementedError()

if __name__ == "__main__":
    a = S4LightSourceBase()
    print(a.info())

    from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical
    b = SourceGeometrical(seed=0)

