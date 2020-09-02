"""
TODO: remove this file
"""
#
# This is just a trial for using in shadow4 only syned objects as inputs
# The parameters that are not inside the standard syned objects are stored
# in ad-hoc created syned-heritated objects (like InputSourceBM)
#
# It is not used for the moment
#

from syned.storage_ring.magnetic_structures.bending_magnet import BendingMagnet
from syned.storage_ring.light_source import LightSource
from syned.storage_ring.electron_beam import ElectronBeam
from syned.syned_object import SynedObject

class InputSourceBM(SynedObject):
    def __init__(self,
        emin=5000.0,  # Photon energy scan from energy (in eV)
        emax = 100000.0,  # Photon energy scan to energy (in eV)
        ng_e = 51,  # Photon energy scan number of points
        ng_j = 20,  # Number of points in electron trajectory (per period) for internal calculation only
        flag_emittance = 1,  # when sampling rays: Use emittance (0=No, 1=Yes)
        ):

        self._emin = emin
        self._emax = emax
        self._ng_e = ng_e
        self._ng_j = ng_j
        self._flag_emittance = flag_emittance

        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("emin", "Minimum photon energy" , "eV"    ),
                    ("emax", "Maximum photon energy", "eV"),
                    ("ng_e", "Number of points in energy", ""),
                    ("ng_j", "Number of points in trajectory", ""),
                    ("flag_emittance", "Use or not emittance", ""),

            ] )

class LightSourceBM(SynedObject):
    def __init__(self,
                 syned_electron_beam=None,
                 syned_bending_magnet=None,
                 input_source_bm=None
                 ):
        self._syned_electron_beam = syned_electron_beam
        self._syned_bending_magnet = syned_bending_magnet
        self._input_source_bm = input_source_bm
        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("syned_electron_beam",     "Electron Beam",        ""),
                    ("syned_bending_magnet",    "Bending Magnet",       ""),
                    ("input_source_bm",         "Specific BM",          ""),
            ] )


if __name__ == "__main__":

    from srxraylib.plot.gol import plot_scatter

    syned_electron_beam = ElectronBeam(energy_in_GeV=6.04,current=0.2,
                                       moment_xx=(0.0078e-2)**2,
                                       moment_xpxp=(3.8e-07/0.0078)**2,
                                       moment_yy=(0.0036*1e-2)**2,
                                       moment_ypyp=(3.8e-09/0.0036)**2,
                                       )

    syned_bending_magnet = BendingMagnet(radius=25.1772,magnetic_field=0.8,length=25.1772*0.001)

    emin = 5000.0                # Photon energy scan from energy (in eV)
    emax = 100000.0              # Photon energy scan to energy (in eV)
    ng_e = 51                    # Photon energy scan number of points
    ng_j = 20                    # Number of points in electron trajectory (per period) for internal calculation only
    flag_emittance = 1           # when sampling rays: Use emittance (0=No, 1=Yes)



    bm_input = InputSourceBM(
                 emin=emin,               # Photon energy scan from energy (in eV)
                 emax=emax,               # Photon energy scan to energy (in eV)
                 ng_e=ng_e,                    # Photon energy scan number of points
                 ng_j=ng_j,                    # Number of points in electron trajectory (per period) for internal calculation only
                 flag_emittance=flag_emittance,           # when sampling rays: Use emittance (0=No, 1=Yes)
                )

    ls = LightSource(syned_electron_beam,syned_bending_magnet)
    ls1 = LightSourceBM(syned_electron_beam,syned_bending_magnet,bm_input)
    # print(bm_input.info())
    print(ls1.info())

    ls1.to_json("tmp.json")

