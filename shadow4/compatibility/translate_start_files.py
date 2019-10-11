from syned.storage_ring.electron_beam import ElectronBeam
from syned.storage_ring.magnetic_structures.bending_magnet import BendingMagnet
import Shadow


from shadow4.sources.bending_magnet.source_bending_magnet import SourceBendingMagnet
from shadow4.compatibility.source import Source as Source4

def is_bending_magnet(oe0):
    if oe0.F_WIGGLER == 0 and (oe0.FDISTR == 4 or oe0.FDISTR == 6):
        return True
    else:
        return False

def start00_to_bending_magnet(input,to_meters=1.0):

    print(">>>>>>>",input,type(input))
    if isinstance(input,type("")):
        # # old
        # oe0 = Shadow.Source()
        # oe0.load(input)

        #new

        oe0 = Source4()
        oe0.load(input)
    else:
        oe0 = input

    if not is_bending_magnet(oe0):
        raise Exception("Not a bending magnet source")


    syned_electron_beam = ElectronBeam(energy_in_GeV=oe0.BENER,current=0.4,
                                       moment_xx   = (oe0.SIGMAX * to_meters)**2,
                                       moment_xpxp = (oe0.EPSI_X / oe0.SIGMAX)**2,
                                       moment_yy   = (oe0.SIGMAZ * to_meters)**2,
                                       moment_ypyp = (oe0.EPSI_Z / oe0.SIGMAZ)**2,
                                       )

    syned_bending_magnet = BendingMagnet.initialize_from_magnetic_radius_divergence_and_electron_energy(
        magnetic_radius=oe0.R_MAGNET,divergence=oe0.HDIV1+oe0.HDIV2,electron_energy_in_GeV=oe0.BENER)


    bm = SourceBendingMagnet(syned_electron_beam=syned_electron_beam,
                 syned_bending_magnet=syned_bending_magnet,
                 emin=oe0.PH1,               # Photon energy scan from energy (in eV)
                 emax=oe0.PH2,               # Photon energy scan to energy (in eV)
                 ng_e=200,                   # Photon energy scan number of points
                 ng_j=100,                   # Number of points in electron trajectory (per period) for internal calculation only
                 flag_emittance=True,        # when sampling rays: Use emittance (0=No, 1=Yes)
                )

    return bm





