import Shadow
import numpy

from minishadow.bending_magnet.source_bending_magnet import SourceBendingMagnet

def run_bm_shadow3(iwrite=0):
    #
    # Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
    #


    # write (1) or not (0) SHADOW files start.xx end.xx star.xx


    #
    # initialize shadow3 source (oe0) and beam
    #
    beam = Shadow.Beam()
    oe0 = Shadow.Source()

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.BENER = 6.04
    oe0.EPSI_X = 3.8e-07
    oe0.EPSI_Z = 3.8e-09
    oe0.FDISTR = 4
    oe0.FSOURCE_DEPTH = 4
    oe0.F_COLOR = 3
    oe0.F_PHOT = 0
    oe0.HDIV1 = 0.0005
    oe0.HDIV2 = 0.0005
    oe0.NCOL = 0
    oe0.N_COLOR = 0
    oe0.PH1 = 5000.0
    oe0.PH2 = 100000.0
    oe0.POL_DEG = 0.0
    oe0.R_ALADDIN = 2517.72
    oe0.R_MAGNET = 25.1772
    oe0.SIGDIX = 0.0
    oe0.SIGDIZ = 0.0
    oe0.SIGMAX = 0.0078
    oe0.SIGMAY = 0.0
    oe0.SIGMAZ = 0.0036
    oe0.VDIV1 = 1.0
    oe0.VDIV2 = 1.0
    oe0.WXSOU = 0.0
    oe0.WYSOU = 0.0
    oe0.WZSOU = 0.0



    #Run SHADOW to create the source

    if iwrite:
        oe0.write("start.00")

    beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")


    Shadow.ShadowTools.plotxy(beam,1,3,nbins=101,nolost=1,title="Real space")
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")

    return beam

if __name__ == "__main__":

    shadow3_beam = run_bm_shadow3()

    bm = SourceBendingMagnet()