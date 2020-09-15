import numpy

from syned.storage_ring.light_source import LightSource as SynedLightSource
from syned.storage_ring.electron_beam import ElectronBeam as SynedElectronBeam


from shadow4.sources.bending_magnet.bending_magnet import BendingMagnet
from shadow4.beam.beam import Beam
from shadow4.sources.light_source import LightSource

class BendingMagnetLightSource(SynedLightSource, LightSource):

    def __init__(self, name="Undefined", electron_beam=None, bending_magnet_magnetic_structure=None):

        if electron_beam is None:
            electron_beam = SynedElectronBeam()

        if bending_magnet_magnetic_structure is None:
            bending_magnet_magnetic_structure = BendingMagnet()

        super().__init__(name,
                         electron_beam=electron_beam,
                         magnetic_structure=bending_magnet_magnetic_structure)

    def get_beam(self, F_COHER=0, NRAYS=5000, SEED=123456,
                       EPSI_DX=0.0, EPSI_DZ=0.0,
                       psi_interval_in_units_one_over_gamma=None,
                       psi_interval_number_of_points=1001,
                       verbose=False):

        rays = self.get_magnetic_structure().calculate_rays(
                       self.get_electron_beam(),
                       F_COHER=F_COHER, NRAYS=NRAYS, SEED=SEED,
                       EPSI_DX=EPSI_DX, EPSI_DZ=EPSI_DZ,
                       psi_interval_in_units_one_over_gamma=psi_interval_in_units_one_over_gamma,
                       psi_interval_number_of_points=psi_interval_number_of_points,
                       verbose=verbose)

        return Beam.initialize_from_array(rays)

    def info(self):
        syned_electron_beam = self.get_electron_beam()
        magnetic_structure = self.get_magnetic_structure()

        txt = ""
        txt += "-----------------------------------------------------\n"

        txt += "Input Electron parameters: \n"
        txt += "        Electron energy: %f geV\n"%syned_electron_beam.energy()
        txt += "        Electron current: %f A\n"%syned_electron_beam.current()
        if magnetic_structure._FLAG_EMITTANCE:
            sigmas = syned_electron_beam.get_sigmas_all()
            txt += "        Electron sigmaX: %g [um]\n"%(1e6*sigmas[0])
            txt += "        Electron sigmaZ: %g [um]\n"%(1e6*sigmas[2])
            txt += "        Electron sigmaX': %f urad\n"%(1e6*sigmas[1])
            txt += "        Electron sigmaZ': %f urad\n"%(1e6*sigmas[3])
        txt += "Input Bending Magnet parameters: \n"
        txt += "        radius: %f m\n"%magnetic_structure._radius
        txt += "        magnetic field: %f T\n"%magnetic_structure._magnetic_field
        txt += "        length: %f m\n"%magnetic_structure._length


        txt += "-----------------------------------------------------\n"

        txt += "Lorentz factor (gamma): %f\n"%syned_electron_beam.gamma()

        txt2 = magnetic_structure.info()
        return (txt + "\n\n" + txt2)

    def to_python_code(self):
        return "# to be implemented..."
