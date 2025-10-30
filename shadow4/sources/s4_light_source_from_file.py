"""
Defines the a LightSource with a beam defined in a HDF5 file.
"""
from syned.storage_ring.empty_light_source import EmptyLightSource
from shadow4.beam.s4_beam import S4Beam

class S4LightSourceFromFile(EmptyLightSource):
    """
    Class to create a light source from a S4Beam in an H5 file.

    Parameters
    ----------
    name : str, optional
        A name.
    file_name : str, optional
        The name of the H5 file.
    simulation_name : str, optional
        A name or key to define the simulation within the H5 file.
    beam_name : str, optional
        A name or key to define the name of the beam with the simulation.

    """
    def __init__(self, name="Undefined", file_name="", simulation_name='run001', beam_name='begin'):
        super().__init__(name=name)
        self._file_name = file_name
        self._simulation_name = simulation_name
        self._beam_name = beam_name

        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
            ("name","Name",""),
            ("file_name","HDF5 file name",""),
            ("simulation_name","name of simulation",""),
            ("beam_name", "name of the beam in a simulation", ""),
            ] )

        ierr = self._load()
        if ierr == 1: print("Error loading data in: %s::/%s/%s/" % (self._file_name, self._simulation_name, self._beam_name))


    def _load(self):
        try:
            self._beam = S4Beam.load_h5(self._file_name, simulation_name=self._simulation_name, beam_name=self._beam_name)
            return 0
        except:
            self._beam = None
            return 1

    def get_beam(self, copy=0):
        """
        Retirns the S4 beam.

        Parameters
        ----------
        copy : int
            Returns the beam stored in the class (0) or a copy of it (1).

        Returns
        -------
            S4Beam instance
            The S4 beam.

        """
        if copy:
            return self._beam.duplicate()
        else:
            return self._beam

    def to_python_code(self, **kwargs):
        """
        Returns the python code for calculating the geometrical source.

        Returns
        -------
        str
            The python code.
        """
        txt = ""
        txt += "\n#\n#\n#"
        txt += "\nfrom shadow4.sources.s4_light_source_from_file import S4LightSourceFromFile"
        txt += "\nlight_source = S4LightSourceFromFile(name='%s', file_name='%s', simulation_name='%s', beam_name='%s')" % \
               (self.get_name(), self._file_name, self._simulation_name, self._beam_name)
        txt += "\nbeam = light_source.get_beam()"
        return txt

    def get_info(self):
        """
        Returns specific information.

        Returns
        -------
        str
        """
        if self._beam is None:
            return "\n\nEmpty beam.\n\n"
        else:
            return "\n\nBeam from %s::/%s/%s/ with %d rays" % \
                   (self._file_name, self._simulation_name, self._beam_name, self.get_beam().N)

if __name__ == "__main__":
    a = S4LightSourceFromFile(file_name="/nobackup/gurb1/srio/Oasys/tmp4.h5")
    print(a.get_beam(copy=0), a.get_beam(copy=1))
    print(a.info())
    print('beam: ', a.get_beam())
    print('script: ', a.to_python_code())
