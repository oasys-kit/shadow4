"""
Defines the shadow4 synchrotron LightSource (the base class of beanding magnet, wiggler and undulator lightsources).

It contains an electron beam and a magnetic structure and some parameters like nrays and seed.
"""
from syned.storage_ring.light_source import LightSource

class S4LightSource(LightSource):
    """
    Constructor

    Parameters
    ----------
    name : str
        A name.
    electron_beam : instance of S4ElectronBeam
        The electron beam.
    magnetic_structure : instance of a class derived from SYNED MagneticStructre.
        The magnetic structure (bending magnet, wiggler or undulator).
    nrays : int
        The number of rays.
    seed : int
        The seed.
    """
    def __init__(self,
                 name="Undefined",
                 electron_beam=None,
                 magnetic_structure=None,
                 nrays=5000,
                 seed=12345,
                 ):
        super().__init__(name=name, electron_beam=electron_beam, magnetic_structure=magnetic_structure)

        self._nrays = nrays
        self._seed = seed

        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._add_support_text([
            ("name","Name",""),
            ("nrays","Number of rays to be generated",""),
            ("seed","Seed for the Monte Carlo generator",""),
            ] )

    def set_nrays(self, nrays):
        """
        Defines the number of rays to be used.

        Parameters
        ----------
        nrays : int
            The number of rays.

        """
        self._nrays = nrays

    def get_nrays(self):
        """
        Returns the stored number of rays.

        Returns
        -------
        int
            The number of rays.

        """
        return self._nrays

    def set_seed(self, seed):
        """
        Defines the Monte Carlo seed.

        Parameters
        ----------
        seed : int
            The seed.

        """
        self._seed = seed

    def get_seed(self):
        """
        Returns the stored seed.

        Returns
        -------
        int
            The seed.
        """
        return self._seed

    def to_python_code(self, **kwargs):
        """
        Returns the python code to create the light source. To be fully defined in the derived classes.

        Parameters
        ----------
        **kwargs
            Passed arguments

        Returns
        -------
        str
            The python code.
        """
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
        """
        To be implemented in a derived class.

        Raises
        ------
        NotImplementedError()

        """
        raise NotImplementedError()

    def calculate_spectrum(self, **params):
        """
        To be implemented in a derived class.

        Raises
        ------
        NotImplementedError()

        """
        raise NotImplementedError()

    def get_info(self):
        """
        Returns specific information.

        Returns
        -------
        str
        """
        return "\n\nSpecific information not available (get_info() method not overloaded).\n\n"

if __name__ == "__main__":
    a = S4LightSource()
    print(a.info())
    print(a.to_python_code())
    print(a.get_info())
    print(a.info())