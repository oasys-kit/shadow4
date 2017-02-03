
import numpy
from collections import OrderedDict

class OpticalElemenLensIdeal(object):

    def __init__(self, name="LensIdeal", focal_x=0.0, focal_z=0.0, p=0.0, q=0.0):
        self._focal_x = focal_x
        self._focal_z = focal_z
        self._name = name
        self._p = p
        self._q = q

    def get_focalX(self):
        return self._focal_x

    def get_focalZ(self):
        return self._focal_z

    def to_dictionary(self):
        #returns a dictionary with the variable names as keys, and a tuple with value, unit and doc string
        mytuple = [ ("focal_x"   ,( self._focal_x ,"m",  "Ideal lens focal length (horizontal)" ) ),
                    ("focal_y"   ,( self._focal_y ,"m",  "Ideal lens focal length (vertical)"  ) )]
        return(OrderedDict(mytuple))

    def trace_beam(self,beam1):
        beam = beam1.duplicate()


        if self._p != 0.0:
            beam.retrace(self._p,resetY=True)

        # rotate around Z
        if self._focal_x != 0.0:
            # whatch out the minus!!
            tan_two_theta = - beam.get_column(1) / self._focal_x
            beam.rotate(numpy.arctan(tan_two_theta),axis=3,rad=True)

        # rotate around X
        if self._focal_z != 0.0:
            tan_two_theta = beam.get_column(3) / self._focal_z
            beam.rotate(numpy.arctan(tan_two_theta),axis=1,rad=True)

        if self._q != 0.0:
            beam.retrace(self._q,resetY=True)

        return beam


def test_with_collimated_beam():


    from SourceGaussian import SourceGaussian
    from Beam import Beam
    from Shadow.ShadowTools import plotxy

    src = SourceGaussian.initialize_collimated_source(number_of_rays=10000,sigmaX=1e-6,sigmaZ=1e-6)

    beam = Beam()

    beam.genSource(src)

    print(beam.info())
    SX, SZ = (1e6*beam.get_standard_deviation(1),1e6*beam.get_standard_deviation(3))

    plotxy(beam.get_shadow3_beam(),1,3,nbins=100,title="SOURCE")

    lens1 = OpticalElemenLensIdeal("test",focal_x=10.0,focal_z=10.0,p=100.0,q=10.0)

    method = 2 # 0:direct, 1:interface with overwrite, 2: no overwrite
    if method == 0:
        beam2 = lens1.trace_beam(beam)
    elif method == 1:
        beam.traceOE(lens1,1,overwrite=True)
        beam2 = beam
    elif method == 2:
        beam2 = beam.traceOE(lens1,1,overwrite=True)
    else:
        raise Exception("Undefined method")

    #

    plotxy(beam2.get_shadow3_beam(),1,3,nbins=100,title="FOCAL PLANE")
    FX, FZ = (1e6*beam2.get_standard_deviation(1),1e6*beam2.get_standard_deviation(3))
    print("Source dimensions: %f %f um"%(SX,SZ))
    print("Focal dimensions: %f %f um"%(FX,FZ))
    print("Demagnification: %g %g"%(SX/FX,SX/FZ))


def test_with_divergent_beam(p=30.0,q=5.0):

    from SourceGaussian import SourceGaussian
    from Beam import Beam
    from Shadow.ShadowTools import plotxy

    src = SourceGaussian.initialize_from_keywords(number_of_rays=10000,sigmaX=20e-6,sigmaZ=10e-6,
                                                  sigmaXprime=50e-6,sigmaZprime=10e-6)

    beam = Beam()

    beam.genSource(src)

    print(beam.info())
    SX, SZ = (1e6*beam.get_standard_deviation(1),1e6*beam.get_standard_deviation(3))

    # plotxy(beam.get_shadow3_beam(),1,3,nbins=100,title="SOURCE")
    plotxy(beam.get_shadow3_beam(),1,4,nbins=100,title="SOURCE PHASE SPACE")

    F = 1.0 / (1/p + 1/q)
    lens1 = OpticalElemenLensIdeal("test",focal_x=F,focal_z=F,p=p,q=q)

    method = 2 # 0:direct, 1:interface with overwrite, 2: no overwrite
    if method == 0:
        beam2 = lens1.trace_beam(beam)
    elif method == 1:
        beam.traceOE(lens1,1,overwrite=True)
        beam2 = beam
    elif method == 2:
        beam2 = beam.traceOE(lens1,1,overwrite=True)
    else:
        raise Exception("Undefined method")

    #

    # plotxy(beam2.get_shadow3_beam(),1,3,nbins=100,title="FOCAL PLANE")
    plotxy(beam2.get_shadow3_beam(),1,4,nbins=100,title="FOCAL PLANE PHASE SPACE")

    FX, FZ = (1e6*beam2.get_standard_deviation(1),1e6*beam2.get_standard_deviation(3))
    print("Source dimensions (rms): %f %f um"%(SX,SZ))
    print("Focal dimensions (rms): %f %f um"%(FX,FZ))
    print("Demagnification: H:%g V:%g (theoretical: %g) "%(SX/FX,SZ/FZ,p/q))



def get_sigmas_radiation(photon_energy,undulator_length):
    import scipy.constants as codata
    lambdan = 1e-10 * codata.h*codata.c/codata.e*1e10 / photon_energy # in m
    print("wavelength in m",lambdan)
    return 1e6*2.740/4/numpy.pi*numpy.sqrt(lambdan*undulator_length),1e6*0.69*numpy.sqrt(lambdan/undulator_length)

def test_id16ni():

    from SourceGaussian import SourceGaussian
    from Beam import Beam
    from Shadow.ShadowTools import plotxy


    ESRF =  {'sigmaX':387.8,"sigmaZ":3.5,"sigmaX'":10.3,"sigmaZ'":1.2}
    EBS =  {'sigmaX':27.2, "sigmaZ":3.4,"sigmaX'":5.2,"sigmaZ'":1.4}
    m = ESRF

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

    src = SourceGaussian.initialize_from_keywords(number_of_rays=10000,sigmaX=Sx,sigmaZ=Sz,
                                                  sigmaXprime=Sxp,sigmaZprime=Szp)

    beam = Beam()

    beam.genSource(src)

    print(beam.info())
    SX, SZ = (1e6*beam.get_standard_deviation(1),1e6*beam.get_standard_deviation(3))

    plotxy(beam.get_shadow3_beam(),1,4,nbins=100,title="SOURCE H phase space")
    # plotxy(beam.get_shadow3_beam(),1,3,nbins=100)

    # multilayer
    p = 28.3
    q = 11.70
    F = 1/(1/p+1/q)
    lens1 = OpticalElemenLensIdeal("ML",focal_x=F,focal_z=0,p=p,q=q)
    beam.traceOE(lens1,1,overwrite=True)

    # TODO: problem with divergences
    plotxy(beam.get_shadow3_beam(),1,4,nbins=100,title="H phase space")


    FX, FZ = (1e6*beam.get_standard_deviation(1),1e6*beam.get_standard_deviation(3))
    print("----------------- Secondary source---------------------")
    print("Source dimensions (rms): %f %f um"%(SX,SZ))
    print("Focal dimensions (rms): %f %f um"%(FX,FZ))
    print("Focal dimensions (FWHM): %f %f um"%(f2dot35*FX,f2dot35*FZ))
    print("Demagnification: H:%g V:%g (theoretical: %g) "%(SX/FX,SZ/FZ,p/q))


    # first KB mirror
    p = 144.90
    q = 0.025
    F = 1.0/(1/184.90+1/0.10)
    lens2 = OpticalElemenLensIdeal("KBV",focal_x=0,focal_z=F,p=p,q=q)
    beam.traceOE(lens2,1,overwrite=True)

    # second KB mirror
    p = 0.025
    q = 0.05
    F = 1.0/(1/144.95+1/0.05)
    lens3 = OpticalElemenLensIdeal("KBH",focal_x=F,focal_z=0,p=p,q=q)
    beam.traceOE(lens3,1,overwrite=True)

    #

    # plotxy(beam2.get_shadow3_beam(),1,3,nbins=100)
    FX, FZ = (1e6*beam.get_standard_deviation(1),1e6*beam.get_standard_deviation(3))
    print("----------------- Focal position ---------------------")
    print("Source dimensions (rms): %f %f um"%(SX,SZ))
    print("Focal dimensions (rms): %f %f um"%(FX,FZ))
    print("Focal dimensions (FWHM): %f %f um"%(f2dot35*FX,f2dot35*FZ))
    print("Demagnification: H:%g V:%g (theoretical: %g) "%(SX/FX,SZ/FZ,p/q))


if __name__ == "__main__":
    # test_with_collimated_beam()
    # todo check with two elements
    test_with_divergent_beam()
    # test_id16ni()