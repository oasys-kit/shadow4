"""
S4 Undulator magnetic structure.
"""

from syned.storage_ring.magnetic_structures.undulator import Undulator



class S4Undulator(Undulator):
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
    maxangle
    ng_t : int, optional
        the number of points for angle theta.
    ng_p : int, optional
        the number of points for angle phi.
    ng_j : int, optional
        Number of points per period in calculating the electron trajectory.
    code_undul_phot : str, optional
        code to be used for calculating the emission vs angle: internal, pysru, SRW
    flag_emittance : int, optional
        Flag: 0=Zero emittance (filament beam), 1=Use emittance.
    flag_size : int, optional
        Method used for sampling the source size:
        0=point, 1=Gaussian, 2=FT(Divergences).
    use_gaussian_approximation : int, optional
        0=No, 1=Yes.
    """
    def __init__(self,
                 K_vertical=1.0,               # syned Undulator parameter
                 period_length=0.01,           # syned Undulator parameter
                 number_of_periods=100,        # syned Undulator parameter
                 emin=10000.0,                 # Photon energy scan from energy (in eV)
                 emax=11000.0,                 # Photon energy scan to energy (in eV)
                 ng_e=11,                      # Photon energy scan number of points
                 maxangle=50e-6,               # Maximum radiation semiaperture in RADIANS
                 ng_t=31,                      # Number of points in angle theta
                 ng_p=21,                      # Number of points in angle phi
                 ng_j=20,                      # Number of points in electron trajectory (per period) for internal calculation only
                 code_undul_phot="internal",   # internal, pysru, srw # todo: remove
                 flag_emittance=0,             # when sampling rays: Use emittance (0=No, 1=Yes)
                 flag_size=0,                  # when sampling rays: 0=point,1=Gaussian,2=FT(Divergences)
                 use_gaussian_approximation=0, # use Gaussian approximation for generating simplified beam
                 distance=100,
                 srw_range=0.05,
                 srw_resolution=50.0,
                 srw_semianalytical=0,
                 magnification=0.01,
                 pysru_source=0, # for pysru/wofry backpropagation: source interpolated from polar (0) or recalculated (1)
                 ):
        super().__init__(K_vertical=K_vertical,
                 K_horizontal = 0.0,
                 period_length = period_length,
                 number_of_periods = number_of_periods)

        # Photon energy scan
        self._EMIN            = emin   # Photon energy scan from energy (in eV)
        self._EMAX            = emax   # Photon energy scan to energy (in eV)
        # Photon energy scan number of points
        if emin == emax: self._NG_E = 1
        else:            self._NG_E            = ng_e

        # Geometry
        self._MAXANGLE        = maxangle   # Maximum radiation semiaperture in RADIANS
        self._NG_T            = ng_t       # Number of points in angle theta
        self._NG_P            = ng_p       # Number of points in angle phi
        self._NG_J            = ng_j       # Number of points in electron trajectory (per period)


        self.code_undul_phot = code_undul_phot

        self._FLAG_EMITTANCE  =  flag_emittance # Yes  # Use emittance (0=No, 1=Yes)
        self._FLAG_SIZE  =  flag_size # 0=point,1=Gaussian,2=backpropagate Divergences
        self._use_gaussian_approximation = use_gaussian_approximation

        # specific parameters
        self._distance = distance
        self._srw_range = srw_range
        self._srw_resolution = srw_resolution
        self._srw_semianalytical = srw_semianalytical
        self._magnification = magnification
        self._pysru_source = pysru_source


        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._add_support_text([
                    ("EMIN", "minimum photon energy", "eV" ),
                    ("EMAX", "maximum photon energy", "eV"),
                    ("NG_E", "number of energy points", ""),
                    ("MAXANGLE", "Maximum radiation semiaperture", "rad"),
                    ("NG_T", "Number of points in angle theta", ""),
                    ("NG_P", "Number of points in angle phi", ""),
                    ("NG_J", "number of points of the electron trajectory", ""),
                    ("FLAG_EMITTANCE", "Use emittance (0=No, 1=Yes)", ""),
                    ("FLAG_SIZE", "Size philament beam 0=point,1=Gaussian,2=backpropagate Divergences", ""),
                    ("use_gaussian_approximation", "0=No, 1=Yes", ""),
            ("distance", "Distance to far field plane", "m"),
            ("srw_range", "for SRW (range magnification)", ""),
            ("srw_resolution", "for SRW (resolution factor)", ""),
            ("srw_semianalytical", "for SRW (semianalytical treatment of phase)", ""),
            ("magnification", "for internal/wofry magnification in propagation", ""),
            ("pysru_source", "for pysru/wofry backpropagation: source interpolated from polar (0) or recalculated (1)", ""),
            ] )


    def use_gaussian_approximation(self):
        """
        Returns a flag to indicate if the computations use the Gaussian approximation.

        Returns
        -------
        int
            0=No, 1=Yes.
        """
        return self._use_gaussian_approximation

    def get_info(self):
        """
        Returns the specific information of the wiggler magnetic structure.

        Returns
        -------
        str
        """
        txt = ""

        txt += "Input Undulator parameters: \n"
        txt += "        period: %f m\n"%self.period_length()
        txt += "        number of periods: %d\n"%self.number_of_periods()
        txt += "        K-value: %f\n"%self.K_vertical()

        # txt += "-----------------------------------------------------\n"

        txt += "\nEnergy Grid: \n"
        if self._NG_E == 1:
            txt += "        photon energy %f eV\n" % (self._EMIN)
        else:
            txt += "        photon energy from %10.3f eV to %10.3f eV\n" % (self._EMIN, self._EMAX)


        if self.use_gaussian_approximation():
            txt += "\nUsing Gaussian Approximation!!!\n"
        else:
            txt += "\nOther Grids: \n"
            txt += "        number of points for the trajectory: %d\n"%(self._NG_J)
            txt += "        number of energy points: %d\n"%(self._NG_E)
            txt += "        maximum elevation angle: %f urad\n"%(1e6*self._MAXANGLE)
            txt += "        number of angular elevation points: %d\n"%(self._NG_T)
            txt += "        number of angular azimuthal points: %d\n"%(self._NG_P)
            # txt += "        number of rays: %d\n"%(self.NRAYS)
            # txt += "        random seed: %d\n"%(self.SEED)
            txt += "-----------------------------------------------------\n"

            txt += "calculation code: %s\n"%self.code_undul_phot
            # if self._result_radiation is None:
            #     txt += "radiation: NOT YET CALCULATED\n"
            # else:
            #     txt += "radiation: CALCULATED\n"
            txt += "Sampling: \n"
            if self._FLAG_SIZE == 0:
                flag = "point"
            elif self._FLAG_SIZE == 1:
                flag = "Gaussian"
            elif self._FLAG_SIZE == 2:
                flag = "Far field backpropagated"

            txt += "        Photon source size sampling flag: %d (%s)\n"%(self._FLAG_SIZE,flag)
            # if self._FLAG_SIZE == 1:
            #     if self._result_photon_size_sigma is not None:
            #         txt += "        Photon source size sigma (Gaussian): %6.3f um \n"%(1e6*self._result_photon_size_sigma)

            # txt += "-----------------------------------------------------\n"
        return txt

    def get_flag_emittance(self):
        """
        Returns the flag for considering electron beam emittance.

        Returns
        -------
        int
            0=No, 1=Yes.
        """
        return self._FLAG_EMITTANCE

    def set_energy_monochromatic(self, emin):
        """
        Sets a single energy line for the source (monochromatic).

        Parameters
        ----------
        emin : float
            photon energy in eV.
        """
        self._EMIN = emin
        self._EMAX = emin
        self._NG_E = 1


    def set_energy_box(self, emin, emax, npoints=None):
        """
        Sets a box for photon energy distribution for the source.

        Parameters
        ----------
        emin : float
            minimum photon energy in eV.
        emax : float
            maximum photon energy in eV.
        npoints : int or None, optional
            Number of points in energy. If None, it does not modify the number of points already stored.
        """
        self._EMIN = emin
        self._EMAX = emax
        if npoints != None:
            self._NG_E = npoints

    def set_maxangle(self, maxangle):
        """
        Defines the maximum angle (in rads) for radiation calculations.

        Parameters
        ----------
        maxangle : float
            The angle in rads.
        """
        self._MAXANGLE = maxangle

    def get_energy_box(self):
        """
        Returns the limits of photon energy distribution for the source and the number of points.

        Returns
        -------
        tuple
            (emin, emax, npoints)
        """
        return self._EMIN, self._EMAX, self._NG_E

    def is_monochromatic(self):
        """
        Returns a flag indicating if the source is monochromatic (True) or polychromatic (False).

        Returns
        -------
        boolean
        """
        if self._NG_E == 1: return True
        if self._EMAX == self._EMIN: return True
        return False

    def to_python_code(self):
        """
        returns the python code for defining the wiggler magnetic structure.

        Returns
        -------
        str
            The python code.
        """
        script_template = """

# magnetic structure
from shadow4.sources.undulator.s4_undulator import S4Undulator
source = S4Undulator(
    K_vertical        = {K_vertical}, # syned Undulator parameter
    period_length     = {period_length}, # syned Undulator parameter
    number_of_periods = {number_of_periods}, # syned Undulator parameter
    emin              = {emin}, # Photon energy scan from energy (in eV)
    emax              = {emax}, # Photon energy scan to energy (in eV)
    ng_e              = {ng_e}, # Photon energy scan number of points
    maxangle          = {maxangle}, # Maximum radiation semiaperture in RADIANS
    ng_t              = {ng_t}, # Number of points in angle theta
    ng_p              = {ng_p}, # Number of points in angle phi
    ng_j              = {ng_j}, # Number of points in electron trajectory (per period) for internal calculation only
    code_undul_phot   = '{code_undul_phot}', # internal, pysru, srw
    flag_emittance    = {flag_emittance}, # when sampling rays: Use emittance (0=No, 1=Yes)
    flag_size         = {flag_size}, # when sampling rays: 0=point,1=Gaussian,2=FT(Divergences)
    use_gaussian_approximation = {use_gaussian_approximation}, # use Gaussian approximation for generating simplified beam
    distance          = {distance}, # distance to far field plane
    srw_range         = {srw_range}, # for SRW backpropagation, the range factor
    srw_resolution    = {srw_resolution}, # for SRW backpropagation, the resolution factor
    srw_semianalytical= {srw_semianalytical}, # for SRW backpropagation, use semianalytical treatement of phase
    magnification     = {magnification}, # for internal/wofry backpropagation, the magnification factor
    pysru_source      = {pysru_source}, # for pysru/wofry backpropagation: source interpolated from polar (0) or recalculated (1)
    )"""


        script_dict = {
            "K_vertical"                 : self.K_vertical(),
            "period_length"              : self.period_length(),
            "number_of_periods"          : self.number_of_periods(),
            "emin"                       : self._EMIN            ,
            "emax"                       : self._EMAX            ,
            "ng_e"                       : self._NG_E            ,
            "maxangle"                   : self._MAXANGLE,
            "ng_t"                       : self._NG_T,
            "ng_p"                       : self._NG_P,
            "ng_j"                       : self._NG_J,
            "code_undul_phot"            : self.code_undul_phot,
            "flag_emittance"             : self._FLAG_EMITTANCE  ,
            "flag_size"                  : self._FLAG_SIZE,
            "use_gaussian_approximation" : self._use_gaussian_approximation,
            "distance"                   : self._distance,
            "srw_range"                  : self._srw_range,
            "srw_resolution"             : self._srw_resolution,
            "srw_semianalytical"         : self._srw_semianalytical,
            "magnification"              : self._magnification,
            "pysru_source"               : self._pysru_source,
        }

        script = script_template.format_map(script_dict)

        return script

if __name__ == "__main__":
    u = S4Undulator(use_gaussian_approximation=1)
    print(u.info())
    print(u.get_info())
    print(u.to_python_code())


