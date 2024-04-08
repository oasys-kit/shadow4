"""
S4 Undulator magnetic structure.
"""

from syned.storage_ring.magnetic_structures.undulator import Undulator
import numpy
import scipy.constants as codata

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
    maxangle: float, optional
        Maximum radiation semiaperture in RADIANS.
    ng_t : int, optional
        the number of points for angle theta.
    ng_p : int, optional
        the number of points for angle phi.
    ng_j : int, optional
        Number of points PER PERIOD in calculating the electron trajectory.
    code_undul_phot : str, optional
        code to be used for calculating the emission vs angle: internal, pysru, SRW
    flag_emittance : int, optional
        Flag: 0=Zero emittance (filament beam), 1=Use emittance.
    flag_size : int, optional
        Method used for sampling the source size:
        0=point, 1=Gaussian, 2=Backpropagated far field.
    use_gaussian_approximation : int, optional
        0=No, 1=Yes.
    srw_range: float, optional
        for code_undul_phot="srw" the range factor.
    srw_resolution: float, optional
        for code_undul_phot="srw" the resolution factor.
    srw_semianalytical: int, optional
        for code_undul_phot="srw" flag to apply the semianalytical phase treatement (0=No, 1=Yes).
    magnification: float, optional
        for code_undul_phot in ["internal", "pysru"]: the magnification factor.
    flag_backprop_recalculate_source: int, optional
        for code_undul_phot in ["internal", "pysru"] and flag_size=2: apply Gaussian weight to backprop amplitudes (0=No, 1=Yes).
    flag_backprop_weight: int, optional
        for code_undul_phot in ["internal", "pysru"] and flag_size=2: apply Gaussian weight to backprop amplitudes.
    weight_ratio: float, optional
        for flag_backprop_weight=1: the Gaussian sigma in units of r.max().
    flag_energy_spread : int, optional
        Flag: 0=do not include energy spread effects, 1=do include energy spread effects
        (the value will be available from the electron beam).
    """
    def __init__(self,
                 K_vertical=1.0,                     # syned Undulator parameter
                 period_length=0.01,                 # syned Undulator parameter
                 number_of_periods=100,              # syned Undulator parameter
                 emin=10000.0,                       # Photon energy scan from energy (in eV)
                 emax=11000.0,                       # Photon energy scan to energy (in eV)
                 ng_e=11,                            # Photon energy scan number of points
                 maxangle=50e-6,                     # Maximum radiation semiaperture in RADIANS
                 ng_t=31,                            # Number of points in angle theta
                 ng_p=21,                            # Number of points in angle phi
                 ng_j=20,                            # Number of points in electron trajectory (per period)
                 code_undul_phot="internal",         # code used for far field and backpropagation: internal, pysru, srw
                 flag_emittance=0,                   # when sampling rays: Use emittance (0=No, 1=Yes)
                 flag_size=0,                        # when sampling rays: 0=point,1=Gaussian, 2=Backpropagation (Divergences)
                 distance=100,                       # distance to the image plane (for far field calculation)
                 srw_range=0.05,                     # for code_undul_phot="srw" the range factor.
                 srw_resolution=50.0,                # for code_undul_phot="srw" the resolution factor.
                 srw_semianalytical=0,               # for code_undul_phot="srw" flag to apply the semianalytical phase treatement
                 magnification=0.01,                 # for code_undul_phot in ["internal", "pysru"]: the magnification factor
                 flag_backprop_recalculate_source=0, # for code_undul_phot in ["internal", "pysru"] and flag_size=2: source reused (0) or recalculated (1)
                 flag_backprop_weight=0,             # for code_undul_phot in ["internal", "pysru"] and flag_size=2: apply Gaussian weight to backprop amplitudes
                 weight_ratio=0.5,                   # for flag_backprop_weight=1: the Gaussian sigma in units of r.max()
                 flag_energy_spread=0,               # include (1) or not (0) energy spread effects. Only valid for monochromatoc simulations.
                 ):
        super().__init__(K_vertical=K_vertical,
                 K_horizontal = 0.0,
                 period_length = period_length,
                 number_of_periods = number_of_periods)

        # Photon energy scan
        self._emin = emin   # Photon energy scan from energy (in eV)
        self._emax = emax   # Photon energy scan to energy (in eV)
        # Photon energy scan number of points
        if emin == emax: self._ng_e = 1
        else:            self._ng_e = ng_e

        # Geometry
        self._maxangle = maxangle   # Maximum radiation semiaperture in RADIANS
        self._ng_t = ng_t       # Number of points in angle theta
        self._ng_p = ng_p       # Number of points in angle phi
        self._ng_j = ng_j       # Number of points in electron trajectory (per period)


        self.code_undul_phot = code_undul_phot

        self._flag_emittance             =  flag_emittance # Yes  # Use emittance (0=No, 1=Yes)
        self._flag_size                  =  flag_size # 0=point,1=Gaussian,2=backpropagate Divergences

        # specific parameters
        self._distance                         = distance
        self._srw_range                        = srw_range
        self._srw_resolution                   = srw_resolution
        self._srw_semianalytical               = srw_semianalytical
        self._magnification                    = magnification

        self._flag_backprop_recalculate_source = flag_backprop_recalculate_source
        self._flag_backprop_weight             = flag_backprop_weight
        self._weight_ratio                     = weight_ratio

        self._flag_energy_spread               = flag_energy_spread


        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._add_support_text([
            ("emin",                             "minimum photon energy",                                                                        "eV" ),
            ("emax",                             "maximum photon energy",                                                                        "eV"),
            ("ng_e",                             "number of energy points",                                                                      ""),
            ("maxangle",                         "Maximum radiation semiaperture",                                                               "rad"),
            ("ng_t",                             "Number of points in angle theta",                                                              ""),
            ("ng_p",                             "Number of points in angle phi",                                                                ""),
            ("ng_j",                             "number of points of the electron trajectory",                                                  ""),
            ("flag_emittance",                   "Use emittance (0=No, 1=Yes)",                                                                  ""),
            ("flag_size",                        "Size philament beam 0=point,1=Gaussian,2=backpropagate Divergences",                           ""),
            ("distance",                         "Distance to far field plane",                                                                  "m"),
            ("srw_range",                        "for SRW (range magnification)",                                                                ""),
            ("srw_resolution",                   "for SRW (resolution factor)",                                                                  ""),
            ("srw_semianalytical",               "for SRW (semianalytical treatment of phase)",                                                  ""),
            ("magnification",                    "for internal/wofry magnification in propagation",                                              ""),
            ("flag_backprop_recalculate_source", "for internal or pysru/wofry backpropagation: source reused (0) or recalculated (1)",           ""),
            ("flag_backprop_weight",             "for internal or pysru/wofry backpropagation: apply Gaussian weight to amplitudes 0=No, 1=Yes", ""),
            ("weight_ratio",                     "for flag_backprop_recalculate_source=1, the ratio value sigma/window halfwidth",               ""),
            ("flag_energy_spread",               "for monochromatod sources, apply (1) or not (0) electron energy spread correction", ""),
        ] )

    def get_info(self, debug=False):
        """
        Returns the specific information for the undulator magnetic structure.

        Returns
        -------
        str
        """
        undulator = self

        txt = ""
        txt = "Undulator source\n"
        txt += "\n"
        txt += "Input Undulator parameters: \n"
        txt += "        period: %f m\n"%undulator.period_length()
        txt += "        number of periods: %d\n"%undulator.number_of_periods()
        txt += "        K-value: %f\n"%undulator.K_vertical()

        txt += "Undulator length: %f m\n"%(undulator.period_length()*undulator.number_of_periods())
        K_to_B = (2.0 * numpy.pi / undulator.period_length()) * codata.m_e * codata.c / codata.e

        txt += "Undulator peak magnetic field: %f T\n"%(K_to_B*undulator.K_vertical())

        txt += "\nGrids: \n"
        if undulator.is_monochromatic():
            txt += "        photon energy %f eV\n"%(undulator._emin)
        else:
            txt += "        photon energy from %10.3f eV to %10.3f eV\n"%(undulator._emin,undulator._emax)
        txt += "        number of points for the trajectory: %d\n"%(undulator._ng_j)
        txt += "        number of energy points: %d\n"%(undulator._ng_e)
        txt += "        maximum elevation angle: %f urad\n"%(1e6*undulator._maxangle)
        txt += "        number of angular elevation points: %d\n"%(undulator._ng_t)
        txt += "        number of angular azimuthal points: %d\n"%(undulator._ng_p)

        txt += "\n"
        txt += "calculation code: %s\n"%undulator.code_undul_phot

        txt += "\n"
        txt += "Sampling: \n"
        if undulator._flag_size == 0:
            flag = "point"
        elif undulator._flag_size == 1:
            flag = "Gaussian"
        elif undulator._flag_size == 2:
            flag = "Far field backpropagated"

        txt += "        Photon source size sampling flag: %d (%s)\n"%(undulator._flag_size, flag)
        return txt

    def __script_dict(self):
        return {
            "K_vertical"                       : self.K_vertical(),
            "period_length"                    : self.period_length(),
            "number_of_periods"                : self.number_of_periods(),
            "emin"                             : self._emin,
            "emax"                             : self._emax,
            "ng_e"                             : self._ng_e,
            "maxangle"                         : self._maxangle,
            "ng_t"                             : self._ng_t,
            "ng_p"                             : self._ng_p,
            "ng_j"                             : self._ng_j,
            "code_undul_phot"                  : self.code_undul_phot,
            "flag_emittance"                   : self.get_flag_emittance()  ,
            "flag_size"                        : self._flag_size,
            "distance"                         : self._distance,
            "srw_range"                        : self._srw_range,
            "srw_resolution"                   : self._srw_resolution,
            "srw_semianalytical"               : self._srw_semianalytical,
            "magnification"                    : self._magnification,
            "flag_backprop_recalculate_source" : self._flag_backprop_recalculate_source,
            "flag_backprop_weight"             : self._flag_backprop_weight,
            "weight_ratio"                     : self._weight_ratio,
            "flag_energy_spread"               : self._flag_energy_spread,
        }

    def get_flag_emittance(self):
        """
        Returns the flag for considering electron beam emittance.

        Returns
        -------
        int
            0=No, 1=Yes.
        """
        return self._flag_emittance

    def get_number_of_energy_points(self):
        """
        Returns the number of photon energy points.

        Returns
        -------
        int
        """
        if self.is_monochromatic():
            return 1
        else:
            return self._ng_e

    def set_energy_monochromatic(self, emin):
        """
        Sets a single energy line for the source (monochromatic).

        Parameters
        ----------
        emin : float
            photon energy in eV.
        """
        self._emin = emin
        self._emax = emin
        self._ng_e = 1


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
        self._emin = emin
        self._emax = emax
        if npoints != None:
            self._ng_e = npoints

    def set_maxangle(self, maxangle):
        """
        Defines the maximum angle (in rads) for radiation calculations.

        Parameters
        ----------
        maxangle : float
            The angle in rads.
        """
        self._maxangle = maxangle

    def get_energy_box(self):
        """
        Returns the limits of photon energy distribution for the source and the number of points.

        Returns
        -------
        tuple
            (emin, emax, npoints)
        """
        if self.is_monochromatic():
            return self._emin, self._emin, 1
        else:
            return self._emin, self._emax, self._ng_e

    def is_monochromatic(self):
        """
        Returns a flag indicating if the source is monochromatic (True) or polychromatic (False).

        Returns
        -------
        boolean
        """
        if self._ng_e == 1: return True
        if self._emax == self._emin: return True
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
    distance          = {distance}, # distance to far field plane
    srw_range         = {srw_range}, # for SRW backpropagation, the range factor
    srw_resolution    = {srw_resolution}, # for SRW backpropagation, the resolution factor
    srw_semianalytical= {srw_semianalytical}, # for SRW backpropagation, use semianalytical treatement of phase
    magnification     = {magnification}, # for internal/wofry backpropagation, the magnification factor
    flag_backprop_recalculate_source      = {flag_backprop_recalculate_source}, # for internal or pysru/wofry backpropagation: source reused (0) or recalculated (1)
    flag_backprop_weight = {flag_backprop_weight}, # for internal or pysru/wofry backpropagation: apply Gaussian weight to amplitudes
    weight_ratio         = {weight_ratio}, # for flag_backprop_recalculate_source=1, the ratio value sigma/window halfwidth
    flag_energy_spread   = {flag_energy_spread}, # for monochromatod sources, apply (1) or not (0) electron energy spread correction
    )"""

        script = script_template.format_map(self.__script_dict())

        return script

if __name__ == "__main__":
    u = S4Undulator()
    # print(u.info())
    print(u.get_info())
    # print(u.to_python_code())


