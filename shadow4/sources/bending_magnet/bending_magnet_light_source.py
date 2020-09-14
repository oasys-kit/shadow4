from wofrysrw.storage_ring.srw_light_source import SRWLightSource
from wofrysrw.storage_ring.magnetic_structures.srw_bending_magnet import SRWBendingMagnet
from wofrysrw.storage_ring.srw_electron_beam import SRWElectronBeam

from syned.storage_ring.light_source import LightSource as SynedLightSource
from syned.storage_ring.electron_beam import ElectronBeam as SynedElectronBeam
from syned.storage_ring.magnetic_structures.bending_magnet import BendingMagnet as SynedBendingMagnet

from shadow4.sources.bending_magnet.bending_magnet import BendingMagnet


class BendingMagnetLightSource(SynedLightSource):

    def __init__(self,
                 name="Undefined",
                 electron_beam=SynedElectronBeam(),
                 bending_magnet_magnetic_structure=BendingMagnet()):


        super().__init__(name,
                         electron_beam=electron_beam,
                         magnetic_structure=bending_magnet_magnetic_structure)

