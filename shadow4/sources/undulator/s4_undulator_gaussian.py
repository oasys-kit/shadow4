"""
S4 Undulator magnetic structure.
"""

from syned.storage_ring.magnetic_structures.undulator import Undulator


class S4UndulatorGaussian(Undulator):
    """
    Defines a shadow4 undulator magnetic structure.

    Parameters
    ----------
    K_vertical : float, optional
        the K value.
    period_length : float, optional
        for magnetic_field_periodic=1, the period in m.
    number_of_periods : float, optional
        for magnetic_field_periodic=1, the number of periods.
    emin : float, optional
        minimum photon energy in eV.
    emax : float, optional
        maximum photon energy in eV.
    ng_e : int, optional
        Number of points in energy. This is used for the calculation of the spectrum, and also for
        sampling rays if flag_interpolation=1.
    flag_emittance : int, optional
        Flag: 0=Zero emittance (filament beam), 1=Use emittance.
    flag_energy_spread : int, optional
        Flag: 0=do not include energy spread effects, 1=do include energy spread effects
        (the value will be available from the electron beam).
    harmonic_number : int, optional
        The harmonic in use (used if flag_energy_spread=1).
    flag_autoset_flux_central_cone : int, optional
        Autocalculate the flux in central cone (to be used for calculating the spectrum).
    flux_central_cone : float, optional
        The value of the flux in the central cone in photons/s/0.1%bw (overwritten if flag_autoset_flux_central_cone=1).
    """
    def __init__(self,
                 period_length=0.01,                 # syned Undulator parameter
                 number_of_periods=100,              # syned Undulator parameter
                 # K_vertical=None,  # syned Undulator parameter NOT USED
                 photon_energy=10000.0,                       # Photon energy scan from energy (in eV)
                 delta_e=0.0,                        # Photon energy width (in eV)
                 ng_e=11,                            # Photon energy scan number of points (for lightsource.calulate_spectrum)
                 flag_emittance=0,  # when sampling rays: Use emittance (0=No, 1=Yes)
                 flag_energy_spread=0,
                 harmonic_number=1,
                 flag_autoset_flux_central_cone=0, # 0=No, 1=Yes
                 flux_central_cone=1e10,
                 ):
        super().__init__(K_vertical=None,
                         K_horizontal=None,
                         period_length=period_length,
                         number_of_periods=number_of_periods,
                         )

        # Photon energy scan
        self._photon_energy = photon_energy
        self._delta_e = delta_e
        # Photon energy scan number of points
        self._NG_E = ng_e

        self._flag_emittance =  flag_emittance # Yes  # Use emittance (0=No, 1=Yes)
        self._flag_energy_spread =  flag_energy_spread # Yes  # Use energy spread (0=No, 1=Yes)
        self._harmonic_number = harmonic_number
        self._flag_autoset_flux_central_cone = flag_autoset_flux_central_cone
        self._flux_central_cone = flux_central_cone


        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._add_support_text([
            ("photon_energy",                    "photon energy in eV",                    "eV" ),
            ("delta_e",                          "photon energy width",                    "eV"),
            ("NG_E",                             "number of energy points",                ""),
            ("flag_emittance",                   "Use emittance (0=No, 1=Yes)",            ""),
            ("flag_energy_spread",               "Use energy spread (0=No,1=Yes)",         ""),
            ("harmonic_number",                  "Number of harmonic",                     ""),
            ("flux_central_cone",                "Value to set flux in central cone",      ""),
            ("flag_autoset_flux_central_cone",   "Automatic set flux in central cone",     ""),

        ] )

    def __script_dict(self):
        return {
            "period_length"                    : self.period_length(),
            "number_of_periods"                : self.number_of_periods(),
            "photon_energy"                    : self._photon_energy,
            "delta_e"                          : self._delta_e,
            "ng_e"                             : self._NG_E            ,
            "flag_emittance"                   : self._flag_emittance  ,
            "flag_energy_spread"               : self._flag_energy_spread,
            "harmonic_number"                  : self._harmonic_number,
            "flux_central_cone"                : self._flux_central_cone,
            "flag_autoset_flux_central_cone"   : self._flag_autoset_flux_central_cone,
        }


    def get_info(self):
        """
        Returns the specific information of the undulator magnetic structure in Gaussian approximation.

        Returns
        -------
        str
        """
        txt = "Undulator source in GAUSSIAN approximation\n"
        txt += "\n"

        txt += "Input Undulator parameters: \n"
        txt += "        period: %f m\n"%self.period_length()
        txt += "        number of periods: %d\n"%self.number_of_periods()
        if self.K_vertical() is not None:
            txt += "        K-value: %f\n"%self.K_vertical()

        txt += "\nEnergy Grid: \n"
        if self.is_monochromatic():
            txt += "        (monochromatic) photon energy %f eV\n" % (self._photon_energy)
        else:
            txt += "        (polychromatic) photon energy from %10.3f eV to %10.3f eV\n" % (\
                self._photon_energy - 0.5 * self._delta_e, self._photon_energy - 0.5 * self._delta_e)


        txt += "\n"
        if self._flag_emittance:
            txt += "Using electron beam emittance\n"
        else:
            txt += "Using ZERO electron beam emittance\n"

        txt += "\n"
        txt += "When sampling rays: use e- energy spread (0=No, 1=Yes): %d \n" % self._flag_energy_spread
        if self._flag_energy_spread:
            txt += "                  harmonic_number = %d \m" % self._harmonic_number

        txt += "\n"
        txt += "Flux peak: "
        if self._flag_autoset_flux_central_cone:
            txt += " calculated \n"
        else:
            txt += " user-defined (%g photons/s/0.1percent bw) \n" % self._flux_central_cone

        txt += "\n"
        txt += "Number of points in energy (for flux calculation): %d \n" % self._NG_E
        return txt

    def get_flag_emittance(self):
        """
        Returns the flag for considering electron beam emittance.

        Returns
        -------
        int
            0=No, 1=Yes.
        """
        return self._flag_emittance

    def get_flag_energy_spread(self):
        """
        Returns the flag for considering electron energy spread.

        Returns
        -------
        int
            0=No, 1=Yes.
        """
        return self._flag_energy_spread

    def get_number_of_energy_points(self):
        """
        Returns the number of photon energy points.

        Returns
        -------
        int
        """

        return self._NG_E

    def set_energy_monochromatic(self, photon_energy):
        """
        Sets a single energy line for the source (monochromatic).

        Parameters
        ----------
        emin : float
            photon energy in eV.
        """
        self._photon_energy = photon_energy
        self._delta_e = 0.0


    def set_energy_box(self, photon_energy, delta_e=0.0):
        """
        Sets a box for photon energy distribution for the source.

        Parameters
        ----------
        photon_energy : float
            photon energy in eV.
        delta_e : float, optional
            interval for photon energy in eV.
        """
        self._photon_energy = photon_energy
        self._delta_e = delta_e

    def get_energy_box(self):
        """
        Returns the limits of photon energy distribution for the source and the number of points.

        Returns
        -------
        tuple
            (photon_energy, delta_e, self._NG_E)
        """
        return self._photon_energy, self._delta_e, self._NG_E

    def is_monochromatic(self):
        """
        Returns a flag indicating if the source is monochromatic (True) or polychromatic (False).

        Returns
        -------
        boolean
        """
        if self._delta_e == 0.0: return True
        return False

    def to_python_code(self):
        """
        Returns the python code for defining the wiggler magnetic structure.

        Returns
        -------
        str
            The python code.
        """
        script_template = """

# magnetic structure
from shadow4.sources.undulator.s4_undulator_gaussian import S4UndulatorGaussian
source = S4UndulatorGaussian(
    period_length     = {period_length},     # syned Undulator parameter (length in m)
    number_of_periods = {number_of_periods}, # syned Undulator parameter
    photon_energy     = {photon_energy}, # Photon energy (in eV)
    delta_e           = {delta_e}, # Photon energy width (in eV)
    ng_e              = {ng_e}, # Photon energy scan number of points
    flag_emittance    = {flag_emittance}, # when sampling rays: Use emittance (0=No, 1=Yes)
    flag_energy_spread = {flag_energy_spread}, # when sampling rays: Use e- energy spread (0=No, 1=Yes)
    harmonic_number    = {harmonic_number}, # harmonic number
    flag_autoset_flux_central_cone  = {flag_autoset_flux_central_cone}, # value to set the flux peak
    flux_central_cone  = {flux_central_cone}, # value to set the flux peak
    )"""


        script = script_template.format_map(self.__script_dict())

        return script

if __name__ == "__main__":
    u = S4UndulatorGaussian()
    # print(u.info())
    print(u.get_info())
    # print(u.to_python_code())