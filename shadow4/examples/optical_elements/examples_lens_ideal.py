
import numpy
from shadow4.beam.beam import Beam
from shadow4.sources.source_geometrical.source_gaussian import SourceGaussian
from shadow4.optical_elements.s4_lens_ideal import S4LensIdeal, S4LensIdealElement
from shadow4.optical_elements.s4_lens_ideal import S4LensSuperIdeal, S4LensSuperIdealElement


from shadow4.compatibility.beam3 import Beam3
from Shadow.ShadowTools import plotxy

# from syned.beamline.optical_elements.ideal_elements.lens import IdealLens as SyIdelLens
from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.beamline_element import BeamlineElement


def get_sigmas_radiation(photon_energy,undulator_length):
    import scipy.constants as codata
    lambdan = 1e-10 * codata.h*codata.c/codata.e*1e10 / photon_energy # in m
    print("wavelength in m",lambdan)
    return 1e6*2.740/4/numpy.pi*numpy.sqrt(lambdan*undulator_length),1e6*0.69*numpy.sqrt(lambdan/undulator_length)


def test_with_collimated_beam(do_plot=True,interface='new'):

    #
    # collimated source
    #
    src = SourceGaussian.initialize_collimated_source(number_of_rays=10000,sigmaX=1e-6,sigmaZ=1e-6)



    if interface == 'new':
        beam = src.get_beam()
    elif interface == "old":
        beam = Beam()
        beam.genSource(src)
    print(beam.info())
    SX, SZ = (1e6*beam.get_standard_deviation(1),1e6*beam.get_standard_deviation(3))

    if do_plot:
        beam3 = Beam3.initialize_from_shadow4_beam(beam)
        plotxy(beam3,1,3,nbins=100,title="SOURCE")

    #
    # lens definition
    #



    lens1e = S4LensIdealElement(optical_element=S4LensIdeal(name="Undefined", focal_x=10.0, focal_y=10.0),
                                coordinates=ElementCoordinates(p=100.0, q=10.0))

    # lens1 = LensIdeal(beamline_element_syned=beamline_element_syned)

    print(lens1e.info())

    #
    # trace
    #
    if interface == 'new':
        beam2 = lens1e.trace_beam(beam)
    elif interface == 'old':
        beam2 = beam.traceOE(lens1, 1, overwrite=False)

    #
    if do_plot:
        beam3 = Beam3.initialize_from_shadow4_beam(beam2)
        plotxy(beam3,1,3,nbins=100,title="FOCAL PLANE")

    FX, FZ = (1e6*beam2.get_standard_deviation(1),1e6*beam2.get_standard_deviation(3))
    print("Source dimensions: %f %f um"%(SX,SZ))
    print("Focal dimensions: %f %f um"%(FX,FZ))
    print("Demagnification: %g %g"%(SX/FX,SX/FZ))


def test_with_divergent_beam_super_ideal(do_plot=True):

    src = SourceGaussian.initialize_from_keywords(number_of_rays=100000,
                                                  sigmaX=     20e-6/2.35,sigmaZ=     10e-6/2.35,
                                                  sigmaXprime=50e-6/2.35,sigmaZprime=10e-6/2.35)

    beam = src.get_beam()

    # point source
    # beam.set_column(1,0.0)
    # beam.set_column(3,0.0)

    print(beam.info())
    SX, SZ = (1e6*beam.get_standard_deviation(1),1e6*beam.get_standard_deviation(3))

    if do_plot:
        beam3 = Beam3.initialize_from_shadow4_beam(beam)
        plotxy(beam3,4,6,nbins=100,title="SOURCE DIVERGENCES")

    beam_tmp = beam.duplicate()
    beam_tmp.retrace(100.0)
    if do_plot:
        beam3 = Beam3.initialize_from_shadow4_beam(beam_tmp)
        plotxy(beam3,1,3,nbins=100,title="SOURCE AFTER 10m")

    p = 10.0
    q = 10.0

    lens1e = S4LensSuperIdealElement(
                        optical_element=S4LensSuperIdeal(name="test1", focal_p_x=p,focal_p_y=p,focal_q_x=q,focal_q_y=q),
                        coordinates=ElementCoordinates(p=p, q=q))


    beam2 = lens1e.trace_beam(beam)

    X = beam2.get_column(1)
    Y = beam2.get_column(3)
    print("X: ",X.min(),X.max(),X.std())
    print("Y: ",Y.min(),Y.max(),Y.std())

    if do_plot:
        beam3 = Beam3.initialize_from_shadow4_beam(beam2)
        plotxy(beam3,1,3,nbins=100,title="FOCAL PLANE",xrange=[-5e-9,5e-9],yrange=[-5e-9,5e-9])

    FX, FZ = (1e6*beam2.get_standard_deviation(1),1e6*beam2.get_standard_deviation(3))
    print("Source dimensions (rms): %f %f um"%(SX,SZ))
    print("Focal dimensions (rms): %f %f um"%(FX,FZ))
    print("Demagnification: H:%g V:%g (theoretical: %g) "%(SX/FX,SZ/FZ,p/q))


def test_with_divergent_beam(do_plot=True):



    src = SourceGaussian.initialize_from_keywords(number_of_rays=100000,
                                                  sigmaX=     20e-6/2.35,sigmaZ=     10e-6/2.35,
                                                  sigmaXprime=50e-6/2.35,sigmaZprime=10e-6/2.35)

    beam = src.get_beam()

    # point source
    # beam.set_column(1,0.0)
    # beam.set_column(3,0.0)

    print(beam.info())
    SX, SZ = (1e6*beam.get_standard_deviation(1),1e6*beam.get_standard_deviation(3))

    if do_plot:
        beam3 = Beam3.initialize_from_shadow4_beam(beam)
        plotxy(beam3,4,6,nbins=100,title="SOURCE DIVERGENCES")

    beam_tmp = beam.duplicate()
    beam_tmp.retrace(100.0)
    if do_plot:
        beam3 = Beam3.initialize_from_shadow4_beam(beam_tmp)
        plotxy(beam3,1,3,nbins=100,title="SOURCE AFTER 10m")


    #
    # lens definition
    #

    p = 10.0
    q = 10.0
    F = 1.0 / (1/p + 1/q)

    lens1e = S4LensIdealElement(optical_element=S4LensIdeal(name="Undefined", focal_x=F, focal_y=F),
                                coordinates=ElementCoordinates(p=p, q=q))

    print(lens1e.info())
    beam2 = lens1e.trace_beam(beam)


    X = beam2.get_column(1)
    Y = beam2.get_column(3)
    print("X: ",X.min(),X.max(),X.std())
    print("Y: ",Y.min(),Y.max(),Y.std())


    if do_plot:
        beam3 = Beam3.initialize_from_shadow4_beam(beam2)
        plotxy(beam3,1,3,nbins=100,title="FOCAL PLANE",xrange=[-5e-9,5e-9],yrange=[-5e-9,5e-9])

    FX, FZ = (1e6*beam2.get_standard_deviation(1),1e6*beam2.get_standard_deviation(3))
    print("Source dimensions (rms): %f %f um"%(SX,SZ))
    print("Focal dimensions (rms): %f %f um"%(FX,FZ))
    print("Demagnification: H:%g V:%g (theoretical: %g) "%(SX/FX,SZ/FZ,p/q))


def test_id16ni(do_plot=True):

    #
    # source
    #
    ESRF =  {'sigmaX':387.8,"sigmaZ":3.5,"sigmaX'":10.3,"sigmaZ'":1.2}
    EBS =  {'sigmaX':27.2, "sigmaZ":3.4,"sigmaX'":5.2,"sigmaZ'":1.4}
    m = EBS

    photon_energy = 17000.0
    undulator_length = 1.4

    sr,srp = get_sigmas_radiation(photon_energy,undulator_length)
    print("radiation sigmas: ",sr,srp)

    demagX = [2.42,2899]
    demagZ = 1849

    #
    f2dot35 = 2*numpy.sqrt(2*numpy.log(2))
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

    # src = SourceGaussian.initialize_from_keywords(number_of_rays=10000,
    #                                               sigmaX=     1000e-6/2.35,sigmaZ=   10e-6/2.35,
    #                                               sigmaXprime=30e-6/2.35,sigmaZprime=15e-6/2.35)
    beam = src.get_beam()
    print(beam.info())
    SX, SZ = (1e6*beam.get_standard_deviation(1),1e6*beam.get_standard_deviation(3))

    if do_plot:
        beam3 = Beam3.initialize_from_shadow4_beam(beam)
        plotxy(beam3,1,3,nbins=100,title="SOURCE")


    #
    # multilayer
    #
    p = 28.3
    q = 11.70
    F = 1/(1/p+1/q)


    lens1e = S4LensIdealElement(optical_element=S4LensIdeal(name="ML", focal_x=F, focal_y=F),
                                coordinates=ElementCoordinates(p=p, q=q))

    beam = lens1e.trace_beam(beam)




    FX, FZ = (1e6*beam.get_standard_deviation(1),1e6*beam.get_standard_deviation(3))
    print("----------------- Secondary source---------------------")
    print("Source dimensions (rms): %f %f um"%(SX,SZ))
    print("Focal dimensions (rms): %f %f um"%(FX,FZ))
    print("Focal dimensions (2.35*rms): %f %f um"%(f2dot35*FX,f2dot35*FZ))
    print("Demagnification: H:%g V:%g (theoretical: %g) "%(SX/FX,SZ/FZ,p/q))


    # first KB mirror
    p = 144.90
    q = 0.025
    F = 1.0/(1/184.90+1/0.10)

    lens2e = S4LensIdealElement(optical_element=S4LensIdeal(name="KBV", focal_x=0, focal_y=F),
                                coordinates=ElementCoordinates(p=p, q=q))


    print(lens2e.info())
    beam = lens2e.trace_beam(beam)



    # second KB mirror
    p = 0.025
    q = 0.05
    F = 1.0/(1/144.95+1/0.05)

    lens3e = S4LensIdealElement(optical_element=S4LensIdeal(name="KBH", focal_x=F, focal_y=0),
                                coordinates=ElementCoordinates(p=p, q=q))

    print(lens3e.info())
    beam = lens3e.trace_beam(beam)


    #
    #
    #
    beam3 = Beam3.initialize_from_shadow4_beam(beam)
    if do_plot:
        tkt = plotxy(beam3,1,3,nbins=300,xrange=[-0.0000005,0.0000005],yrange=[-0.0000005,0.0000005],title="FOCAL PLANE")

        print(tkt['fwhm_h'],tkt['fwhm_v'])

    FX, FZ = (1e6*beam.get_standard_deviation(1),1e6*beam.get_standard_deviation(3))
    print("----------------- Focal position ---------------------")
    print("Source dimensions (rms): %f %f um"%(SX,SZ))
    print("Source dimensions (2.35*rms): %f %f um"%(f2dot35*SX,f2dot35*SZ))
    print("Focal dimensions (rms): %f %f um"%(FX,FZ))
    print("Focal dimensions (2.35*rms): %f %f um"%(f2dot35*FX,f2dot35*FZ))
    if do_plot: print("Focal dimensions (FWHM HISTO): %f %f nm"%(1e9*tkt['fwhm_h'],1e9*tkt['fwhm_v']))
    print("Demagnification (StDev) H:%g V:%g (theoretical: %g,%g) "%(SX/FX,SZ/FZ,demagX[0]*demagX[1],demagZ))
    if do_plot: print("Demagnification:(HISTO) H:%g V:%g (theoretical: %g,%g) "%(f2dot35*SX/(f2dot35*1e6*tkt['fwhm_h']),SZ/(1e6*tkt['fwhm_v']),demagX[0]*demagX[1],demagZ))

if __name__ == "__main__":
    test_with_collimated_beam(do_plot=False,interface='new')
    test_with_divergent_beam(do_plot=False)
    test_id16ni(do_plot=False)
    test_with_divergent_beam_super_ideal(do_plot=False)