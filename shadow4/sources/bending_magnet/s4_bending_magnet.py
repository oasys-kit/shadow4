import numpy


from syned.storage_ring.magnetic_structures.bending_magnet import BendingMagnet

class S4BendingMagnet(BendingMagnet):
    def __init__(self,
                 radius=1.0, magnetic_field=1.0, length=1.0, # syned BM
                 emin=10000.0,               # Photon energy scan from energy (in eV)
                 emax=11000.0,               # Photon energy scan to energy (in eV)
                 ng_e=11,                    # Photon energy scan number of points
                 ng_j=20,                    # Number of points in electron trajectory (per period) for internal calculation only
                 flag_emittance=0,           # when sampling rays: Use emittance (0=No, 1=Yes)
                 ):

        super().__init__(radius=radius, magnetic_field=magnetic_field, length=length)

        # Photon energy scan
        self._EMIN            = emin   # Photon energy scan from energy (in eV)
        self._EMAX            = emax   # Photon energy scan to energy (in eV)
        self._NG_E            = ng_e   # Photon energy scan number of points

        self._NG_J            = ng_j       # Number of points in electron trajectory (per period)

        self._FLAG_EMITTANCE  =  flag_emittance # Yes  # Use emittance (0=No, 1=Yes) #todo kw in calculate rays

    def info(self, debug=False):

        txt = ""

        txt += "-----------------------------------------------------\n"
        txt += "Grids: \n"
        if self._NG_E == 1:
            txt += "        photon energy %f eV\n"%(self._EMIN)
        else:
            txt += "        photon energy from %10.3f eV to %10.3f eV\n"%(self._EMIN,self._EMAX)
        txt += "        number of energy points: %d\n"%(self._NG_E)
        txt += "        number of points for the trajectory: %d\n"%(self._NG_J)


        txt += "-----------------------------------------------------\n"


        return super(S4BendingMagnet, self).info() + "\n" + txt

        return txt


    def set_energy_monochromatic(self,emin):
        """
        Sets a single energy line for the source (monochromatic)
        :param emin: the energy in eV
        :return:
        """
        self._EMIN = emin
        self._EMAX = emin
        self._NG_E = 1


    def set_energy_box(self,emin,emax,npoints=None):
        """
        Sets a box for photon energy distribution for the source
        :param emin:  Photon energy scan from energy (in eV)
        :param emax:  Photon energy scan to energy (in eV)
        :param npoints:  Photon energy scan number of points (optinal, if not set no changes)
        :return:
        """

        self._EMIN = emin
        self._EMAX = emax
        if npoints != None:
            self._NG_E = npoints

    def get_energy_box(self):
        """
        Gets the limits of photon energy distribution for the source
        :return: emin,emax,number_of_points
        """
        return self._EMIN,self._EMAX,self._NG_E

    def is_monochromatic(self):
        if self._NG_E == 1:
            return True
        if self._EMAX == self._EMIN:
            return True
        return False

    @classmethod
    def initialize_from_magnetic_field_divergence_and_electron_energy(cls,
                magnetic_field=1.0,divergence=1e-3,electron_energy_in_GeV=1.0,
                emin=10000.0,
                emax=11000.0,
                ng_e=11,  # Photon energy scan number of points
                ng_j=20,
                flag_emittance=0,
                ):
        """
        Constructor from  magnetic field divergence and electron energy
        :param magnetic_field: in T
        :param divergence: in rad
        :param electron_energy_in_GeV: in GeV
        :return:
        """
        magnetic_radius = cls.calculate_magnetic_radius(magnetic_field, electron_energy_in_GeV)
        return S4BendingMagnet(magnetic_radius,magnetic_field,numpy.abs(divergence * magnetic_radius),
                               emin=emin, emax=emax, ng_e=ng_e, ng_j=ng_j, flag_emittance=flag_emittance)

    @classmethod
    def initialize_from_magnetic_radius_divergence_and_electron_energy(cls,
                magnetic_radius=10.0,divergence=1e-3,electron_energy_in_GeV=1.0,
                emin=10000.0,
                emax=11000.0,
                ng_e=11,  # Photon energy scan number of points
                ng_j=20,
                flag_emittance=0,
                ):

        """
        Constructor from  magnetic radius, divergence and electron energy
        :param magnetic_radius: in m
        :param divergence: in rad
        :param electron_energy_in_GeV: in GeV
        :return:
        """
        magnetic_field = cls.calculate_magnetic_field(magnetic_radius, electron_energy_in_GeV)
        return S4BendingMagnet(magnetic_radius,magnetic_field,numpy.abs(divergence * magnetic_radius),
                               emin=emin, emax=emax, ng_e=ng_e, ng_j=ng_j, flag_emittance=flag_emittance)

    def to_python_code(self):
        script_template = """

#magnetic structure
from shadow4.sources.bending_magnet.s4_bending_magnet import S4BendingMagnet
source = S4BendingMagnet(
                 radius={radius}, # from syned BM, can be obtained as S4BendingMagnet.calculate_magnetic_radius({magnetic_field}, electron_beam.energy())
                 magnetic_field={magnetic_field}, # from syned BM
                 length={length}, # from syned BM = abs(BM divergence * magnetic_field)
                 emin={emin},     # Photon energy scan from energy (in eV)
                 emax={emax},     # Photon energy scan to energy (in eV)
                 ng_e={ng_e},     # Photon energy scan number of points
                 ng_j={ng_j},     # Number of points in electron trajectory (per period) for internal calculation only
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
            "ng_j": self._NG_J,
            "flag_emittance": self._FLAG_EMITTANCE,
        }

        script = script_template.format_map(script_dict)

        return script

if __name__ == "__main__":

    emin = 5000.0                # Photon energy scan from energy (in eV)
    emax = 100000.0              # Photon energy scan to energy (in eV)
    ng_e = 51                    # Photon energy scan number of points
    ng_j = 20                    # Number of points in electron trajectory (per period) for internal calculation only
    flag_emittance = 1           # when sampling rays: Use emittance (0=No, 1=Yes)



    bm = S4BendingMagnet(
                 radius=25.1772, magnetic_field=0.8, length=25.1772 * 0.001,
                 emin=emin,               # Photon energy scan from energy (in eV)
                 emax=emax,               # Photon energy scan to energy (in eV)
                 ng_e=ng_e,                    # Photon energy scan number of points
                 ng_j=ng_j,                    # Number of points in electron trajectory (per period) for internal calculation only
                 flag_emittance=flag_emittance,           # when sampling rays: Use emittance (0=No, 1=Yes)
                )

    print(bm.info())

    print(bm.to_python_code())

