"""
Bending magnet magnetic structure.
"""
import numpy
from syned.storage_ring.magnetic_structures.bending_magnet import BendingMagnet

class S4BendingMagnet(BendingMagnet):
    """
    Defines a shadow4 bending magnet magnetic structure.

    Parameters
    ----------
    radius : float, optional
        Physical Radius/curvature of the magnet in m.
    magnetic_field : float, optional
         Magnetic field strength in T.
    length : float, optional
        physical length of the bending magnet (along the arc) in m.
    emin : float, optional
        minimum photon energy in eV.
    emax : float, optional
        maximum photon energy in eV.
    ng_e : int, optional
        Number of points in energy.
    flag_emittance : int, optional
        Flag: 0=Zero emmitance (filament beam), 1=Use emittance.
    """
    def __init__(self,
                 radius=1.0,         # syned BM
                 magnetic_field=1.0, # syned BM
                 length=1.0,         # syned BM
                 emin=10000.0,       # Photon energy scan from energy (in eV)
                 emax=11000.0,       # Photon energy scan to energy (in eV)
                 ng_e=11,            # Photon energy scan number of points
                 flag_emittance=0,   # when sampling rays: Use emittance (0=No, 1=Yes)
                 ):
        super().__init__(radius=radius, magnetic_field=magnetic_field, length=length)

        # Photon energy scan
        # todo: rename to low case?
        self._EMIN            = emin   # Photon energy scan from energy (in eV)
        self._EMAX            = emax   # Photon energy scan to energy (in eV)
        self._NG_E            = ng_e   # Photon energy scan number of points
        self._FLAG_EMITTANCE  =  flag_emittance # Yes  # Use emittance (0=No, 1=Yes) #todo kw in calculate rays

        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._add_support_text([
                    ("EMIN", "minimum photon energy", "eV" ),
                    ("EMAX", "maximum photon energy", "eV"),
                    ("NG_E", "number of energy points", ""),
                    ("FLAG_EMITTANCE", "Use emittance (0=No, 1=Yes)", "" ),
            ] )

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
        npoints : int, optional
            Number of points in energy.
        """
        self._EMIN = emin
        self._EMAX = emax
        if npoints != None: self._NG_E = npoints

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

    @classmethod
    def initialize_from_magnetic_field_divergence_and_electron_energy(cls,
                magnetic_field=1.0,
                divergence=1e-3,
                electron_energy_in_GeV=1.0,
                emin=10000.0,
                emax=11000.0,
                ng_e=11,  # Photon energy scan number of points
                flag_emittance=0,
                ):
        """
        Constructor from magnetic field divergence and electron energy.

        Parameters
        ----------
        magnetic_field : float, optional
             Magnetic field strength in T.
        divergence : float, optional
            The accepted divergence in rad.
        electron_energy_in_GeV : float, optional
            The electron beam energy in GeV.
        emin : float, optional
            minimum photon energy in eV.
        emax : float, optional
            maximum photon energy in eV.
        ng_e : int, optional
            Number of points in energy.
        flag_emittance : int, optional
            Flag: 0=Zero emmitance (filament beam), 1=Use emittance.

        Returns
        -------
        instance of S4BendingMagnet
        """
        magnetic_radius = cls.calculate_magnetic_radius(magnetic_field, electron_energy_in_GeV)
        return S4BendingMagnet(radius=magnetic_radius,
                               magnetic_field=magnetic_field,
                               emin=emin,
                               emax=emax,
                               ng_e=ng_e,
                               flag_emittance=flag_emittance)

    @classmethod
    def initialize_from_magnetic_radius_divergence_and_electron_energy(cls,
                magnetic_radius=10.0,
                divergence=1e-3,
                electron_energy_in_GeV=1.0,
                emin=10000.0,
                emax=11000.0,
                ng_e=11,  # Photon energy scan number of points
                flag_emittance=0,
                ):
        """
        Constructor from  magnetic radius, divergence and electron energy.

        Parameters
        ----------
        magnetic_radius : float, optional
            The magnetic radius in m.
        divergence : float, optional
            The accepted divergence in rad.
        electron_energy_in_GeV : float, optional
            The electron beam energy in GeV.
        emin : float, optional
            minimum photon energy in eV.
        emax : float, optional
            maximum photon energy in eV.
        ng_e : int, optional
            Number of points in energy.
        flag_emittance : int, optional
            Flag: 0=Zero emmitance (filament beam), 1=Use emittance.

        Returns
        -------
        instance of S4BendingMagnet
        """
        magnetic_field = cls.calculate_magnetic_field(magnetic_radius, electron_energy_in_GeV)
        return S4BendingMagnet(radius=magnetic_radius,
                               magnetic_field=magnetic_field,
                               emin=emin,
                               emax=emax,
                               ng_e=ng_e,
                               flag_emittance=flag_emittance)

    def to_python_code(self):
        """
        returns the python code for defining the bending magnet magnetic structure.

        Returns
        -------
        str
            The python code.
        """
        script_template = """

#magnetic structure
from shadow4.sources.bending_magnet.s4_bending_magnet import S4BendingMagnet
source = S4BendingMagnet(
                 radius={radius}, # from syned BM, can be obtained as S4BendingMagnet.calculate_magnetic_radius({magnetic_field}, electron_beam.energy())
                 magnetic_field={magnetic_field}, # from syned BM
                 length={length}, # from syned BM = abs(BM divergence * magnetic_radius)
                 emin={emin},     # Photon energy scan from energy (in eV)
                 emax={emax},     # Photon energy scan to energy (in eV)
                 ng_e={ng_e},     # Photon energy scan number of points
                 flag_emittance={flag_emittance}, # when sampling rays: Use emittance (0=No, 1=Yes)
                 )
"""

        script_dict = {
            "radius": self.radius(),
            "magnetic_field": self.magnetic_field(),
            "length": self.length(),
            "emin": self._EMIN,
            "emax": self._EMAX,
            "ng_e": self._NG_E,
            "flag_emittance": self._FLAG_EMITTANCE,
        }

        script = script_template.format_map(script_dict)

        return script

if __name__ == "__main__":

    emin = 5000.0                # Photon energy scan from energy (in eV)
    emax = 100000.0              # Photon energy scan to energy (in eV)
    ng_e = 51                    # Photon energy scan number of points
    flag_emittance = 1           # when sampling rays: Use emittance (0=No, 1=Yes)



    bm = S4BendingMagnet(
                 radius=25.1772,
                 magnetic_field=0.8,
                 length=25.1772 * 0.001,
                 emin=emin,               # Photon energy scan from energy (in eV)
                 emax=emax,               # Photon energy scan to energy (in eV)
                 ng_e=ng_e,                    # Photon energy scan number of points
                 flag_emittance=flag_emittance,           # when sampling rays: Use emittance (0=No, 1=Yes)
                )

    print(bm.info())

    print(bm.to_python_code())

