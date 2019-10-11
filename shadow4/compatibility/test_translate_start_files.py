

from syned.storage_ring.magnetic_structures.bending_magnet import BendingMagnet
import Shadow
from shadow4.compatibility.translate_start_files import start00_to_bending_magnet

from shadow4.sources.bending_magnet.source_bending_magnet import SourceBendingMagnet
from shadow4.compatibility.source import Source as Source4

def example_bm(iwrite=0):

    oe0 = Shadow.Source()


    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.BENER = 1.9
    oe0.EPSI_X = 1.989e-09
    oe0.EPSI_Z = 3.007e-11
    oe0.FDISTR = 6
    oe0.FSOURCE_DEPTH = 4
    oe0.F_COLOR = 3
    oe0.F_PHOT = 0
    oe0.HDIV1 = 0.035
    oe0.HDIV2 = 0.035
    oe0.ISTAR1 = 5676561
    oe0.NCOL = 0
    oe0.NPOINT = 10000
    oe0.N_COLOR = 0
    oe0.PH1 = 0.4
    oe0.PH2 = 0.401
    oe0.POL_DEG = 0.0
    oe0.R_ALADDIN = -5.0
    oe0.R_MAGNET = -5.0
    oe0.SIGDIX = 0.0
    oe0.SIGDIZ = 0.0
    oe0.SIGMAX = 3.9e-05
    oe0.SIGMAY = 0.0
    oe0.SIGMAZ = 3.1e-05
    oe0.VDIV1 = 0.05
    oe0.VDIV2 = 0.05
    oe0.WXSOU = 0.0
    oe0.WYSOU = 0.0
    oe0.WZSOU = 0.0


    if iwrite:
        oe0.write("start.00")

    return oe0

def test_start00_to_bending_magnet(iwrite=0,run=False):


    oe0 = example_bm(iwrite=iwrite)

    if iwrite:
        bm = start00_to_bending_magnet("start.00")
    else:
        bm = start00_to_bending_magnet(oe0)

    print(bm.info())

    if run:
        print(bm.info())

        rays = bm.calculate_rays(F_COHER=oe0.F_COHER,NRAYS=oe0.NPOINT,SEED=oe0.ISTAR1,
                                 EPSI_DX=oe0.EPSI_DX,EPSI_DZ=oe0.EPSI_DZ,verbose=False)

        from srxraylib.plot.gol import plot_scatter,set_qt
        set_qt()

        plot_scatter(rays[:,0]*1e6,rays[:,2]*1e6,xtitle="X um",ytitle="Z um")
        plot_scatter(rays[:,1],rays[:,0]*1e6,xtitle="Y m",ytitle="X um")
        plot_scatter(rays[:,1],rays[:,2]*1e6,xtitle="Y m",ytitle="Z um")
        plot_scatter(rays[:,3]*1e6,rays[:,5]*1e6,xtitle="X' urad",ytitle="Z' urad")




if __name__ == "__main__":


    test_start00_to_bending_magnet(iwrite=0,run=False)
    test_start00_to_bending_magnet(iwrite=1,run=True)