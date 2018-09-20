import numpy
from collections import OrderedDict

from minishadow.source_geometrical.gaussian import SourceGaussian
from minishadow.beam.beam import Beam
from minishadow.optical_elements.lens_ideal import LensIdeal, LensSuperIdeal

try:
    import Shadow
    HAVE_SHADOW3 = True
except:
    HAVE_SHADOW3 = False

f2dot35 = 2*numpy.sqrt(2*numpy.log(2))

def get_sigmas_radiation(photon_energy,undulator_length):
    import scipy.constants as codata
    lambdan = 1e-10 * codata.h*codata.c/codata.e*1e10 / photon_energy # in m
    print("wavelength in m",lambdan)
    return 1e6*2.740/4/numpy.pi*numpy.sqrt(lambdan*undulator_length),1e6*0.69*numpy.sqrt(lambdan/undulator_length)

def id16ni_source(do_plot=False):

    ESRF =  {'sigmaX':387.8,"sigmaZ":3.5,"sigmaX'":10.3,"sigmaZ'":1.2}
    EBS =  {'sigmaX':27.2, "sigmaZ":3.4,"sigmaX'":5.2,"sigmaZ'":1.4}
    m = ESRF

    photon_energy = 17000.0
    undulator_length = 1.4


    sr,srp = get_sigmas_radiation(photon_energy,undulator_length)

    print("radiation sigmas: ",sr,srp)

    #


    sx,sz,sxp,szp = m['sigmaX'],m['sigmaZ'],m["sigmaX'"],m["sigmaZ'"]

    Sx  = 1e-6 * numpy.sqrt( sx**2 + sr**2)
    Sz  = 1e-6 * numpy.sqrt( sz**2 + sr**2)
    Sxp = 1e-6 * numpy.sqrt( sxp**2 + srp**2)
    Szp = 1e-6 * numpy.sqrt( szp**2 + srp**2)

    print("Gaussian source dimensions (rms) x z x' z'",1e6*Sx,1e6*Sz,1e6*Sxp,1e6*Szp)
    print("Gaussian source dimensions (FWHM) x z x' z'",f2dot35*1e6*Sx,f2dot35*1e6*Sz,f2dot35*1e6*Sxp,f2dot35*1e6*Szp)

    src = SourceGaussian.initialize_from_keywords(number_of_rays=500000,
                                                  sigmaX=Sx,sigmaZ=Sz,
                                                  sigmaXprime=Sxp,sigmaZprime=Szp)

    beam = Beam()

    beam.genSource(src)

    # point source
    # beam.set_column(1,0.0)
    # beam.set_column(3,0.0)
    # beam.set_column(4,0.0)
    # beam.set_column(5,1.0)
    # beam.set_column(6,0.0)

    print(beam.info())
    SX, SZ = (1e6*beam.get_standard_deviation(1),1e6*beam.get_standard_deviation(3))

    if do_plot:
        if HAVE_SHADOW3:
            Shadow.ShadowTools.plotxy(beam.get_shadow3_beam(),1,3,nbins=100)
        else:
            print("Error: cannot plot right now (need import shadow3)")

    print("----------------- Source---------------------")
    print("Source dimensions (rms): %f %f um"%(SX,SZ))
    print("Source dimensions (2.35*rms): %f %f um"%(f2dot35*SX,f2dot35*SZ))

    return beam

def id16ni_multilayer_as_ideal_lens(beam,do_plot=False):

    SX, SZ = (1e6*beam.get_standard_deviation(1),1e6*beam.get_standard_deviation(3))

    # multilayer
    p = 28.3
    q = 11.70
    F = 1/(1/p+1/q)
    lens1 = LensIdeal("ML",focal_x=F,focal_z=0,p=p,q=q)
    beam.traceOE(lens1,1,overwrite=True)


    FX, FZ = (1e6*beam.get_standard_deviation(1),1e6*beam.get_standard_deviation(3))
    print("----------------- Secondary source---------------------")
    print("Focal dimensions (rms): %f %f um"%(FX,FZ))
    print("Focal dimensions (2.35*rms): %f %f um"%(f2dot35*FX,f2dot35*FZ))
    print("Demagnification: H:%g V:%g (theoretical: %g) "%(SX/FX,SZ/FZ,p/q))

    if do_plot:
        Shadow.ShadowTools.plotxy(beam.get_shadow3_beam(),1,3,nbins=300,title="Secondary source")


    return beam

def id16ni_kb_as_ideal_lenses(beam,do_plot=False):

    SX, SZ = (1e6*beam.get_standard_deviation(1),1e6*beam.get_standard_deviation(3))

    # first KB mirror
    p = 144.90
    q = 0.025
    F = 1.0/(1/184.90+1/0.10)
    lens2 = LensIdeal("KBV",focal_x=0,focal_z=F,p=p,q=q)
    # lens2 = OpticalElemenLensSuperIdeal("KBV",focal_p_x=0,focal_q_x=0,focal_p_z=184.90,focal_q_z=0.10,p=p,q=q)
    beam.traceOE(lens2,1,overwrite=True)

    # second KB mirror
    p = 0.025
    q = 0.05
    F = 1.0/(1/144.95+1/0.05)
    lens3 = LensIdeal("KBH",focal_x=F,focal_z=0,p=p,q=q)
    # lens3 = OpticalElemenLensSuperIdeal("KBH",focal_p_x=144.95,focal_q_x=0.05,focal_p_z=0,focal_q_z=0,p=p,q=q)
    beam.traceOE(lens3,1,overwrite=True)

    #

    tkt = beam.histo2(1,3,nbins=300,xrange=[-0.0000005,0.0000005],yrange=[-0.0000005,0.0000005])
    if do_plot:
        if HAVE_SHADOW3:
            Shadow.ShadowTools.plotxy(beam.get_shadow3_beam(),1,3,nbins=300,xrange=[-0.0000005,0.0000005],yrange=[-0.0000005,0.0000005],title="FOCAL PLANE")
        else:
            print("Error: cannot plot right now (need import shadow3)")

    print(tkt['fwhm_h'],tkt['fwhm_v'])

    FX, FZ = (1e6*beam.get_standard_deviation(1),1e6*beam.get_standard_deviation(3))

    demagX = [2.42,2899]
    demagZ = 1849

    print("----------------- Focal position ---------------------")
    print("Source dimensions (rms): %f %f um"%(SX,SZ))
    print("Source dimensions (2.35*rms): %f %f um"%(f2dot35*SX,f2dot35*SZ))
    print("Focal dimensions (rms): %f %f um"%(FX,FZ))
    print("Focal dimensions (2.35*rms): %f %f um"%(f2dot35*FX,f2dot35*FZ))
    print("Focal dimensions (FWHM HISTO): %f %f nm"%(1e9*tkt['fwhm_h'],1e9*tkt['fwhm_v']))
    print("Demagnification (StDev) H:%g V:%g (theoretical: %g,%g) "%(SX/FX,SZ/FZ,demagX[0]*demagX[1],demagZ))
    print("Demagnification:(HISTO) H:%g V:%g (theoretical: %g,%g) "%(f2dot35*SX/(f2dot35*1e6*tkt['fwhm_h']),SZ/(1e6*tkt['fwhm_v']),demagX[0]*demagX[1],demagZ))

if __name__ == "__main__":
    beam = id16ni_source(do_plot=False)
    beam = id16ni_multilayer_as_ideal_lens(beam,do_plot=False)
    beam = id16ni_kb_as_ideal_lenses(beam,do_plot=True)
