
import numpy
from collections import OrderedDict

class LensIdeal(object):

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

class LensSuperIdeal(object):

    def __init__(self, name="LensIdeal", focal_p_x=0.0, focal_p_z=0.0,
                                        focal_q_x=0.0, focal_q_z=0.0,
                                        p=0.0, q=0.0):
        self._focal_p_x = focal_p_x
        self._focal_p_z = focal_p_z
        self._focal_q_x = focal_q_x
        self._focal_q_z = focal_q_z
        self._name = name
        self._p = p
        self._q = q

    # def get_focalX(self):
    #     return self._focal_x
    #
    # def get_focalZ(self):
    #     return self._focal_z
    #
    # def to_dictionary(self):
    #     #returns a dictionary with the variable names as keys, and a tuple with value, unit and doc string
    #     mytuple = [ ("focal_x"   ,( self._focal_x ,"m",  "Ideal lens focal length (horizontal)" ) ),
    #                 ("focal_y"   ,( self._focal_y ,"m",  "Ideal lens focal length (vertical)"  ) )]
    #     return(OrderedDict(mytuple))

    def trace_beam(self,beam1):
        beam = beam1.duplicate()


        if self._p != 0.0:
            beam.retrace(self._p,resetY=True)

        # rotate around Z; watch out the minus!!
        if self._focal_p_x != 0.0:
            tan_theta_p = - beam.get_column(1) / self._focal_p_x
        else:
            tan_theta_p = 0.0

        if self._focal_q_x != 0.0:
            tan_theta_q = - beam.get_column(1) / self._focal_q_x
        else:
            tan_theta_q = 0.0

        two_theta = numpy.arctan(tan_theta_p) + numpy.arctan(tan_theta_q)
        beam.rotate(two_theta,axis=3,rad=True)

        # rotate around X
        if self._focal_p_z != 0.0:
            tan_theta_p = beam.get_column(3) / self._focal_p_z
        else:
            tan_theta_p = 0.0

        if self._focal_q_z != 0.0:
            tan_theta_q = beam.get_column(3) / self._focal_q_z
        else:
            tan_theta_q = 0.0

        two_theta = numpy.arctan(tan_theta_p) + numpy.arctan(tan_theta_q)
        beam.rotate(two_theta,axis=1,rad=True)

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

    lens1 = LensIdeal("test",focal_x=10.0,focal_z=10.0,p=100.0,q=10.0)

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


def test_with_divergent_beam():

    from SourceGaussian import SourceGaussian
    from SourceGridCartesian import SourceGridCartesian
    from Beam import Beam
    from Shadow.ShadowTools import plotxy

    src = SourceGaussian.initialize_from_keywords(number_of_rays=100000,
                                                  sigmaX=     20e-6/2.35,sigmaZ=     10e-6/2.35,
                                                  sigmaXprime=50e-6/2.35,sigmaZprime=10e-6/2.35)

    # src = SourceGridCartesian.initialize_point_source(
    #             direction_space        = [2*numpy.arccos(numpy.pi/4),2*numpy.arccos(numpy.pi/4)],
    #             direction_space_points = [20,  20],
    #             direction_space_center = [0.0, 0.0] )

    beam = Beam()

    beam.genSource(src)

    # point source
    beam.set_column(1,0.0)
    beam.set_column(3,0.0)

    print(beam.info())
    SX, SZ = (1e6*beam.get_standard_deviation(1),1e6*beam.get_standard_deviation(3))

    plotxy(beam.get_shadow3_beam(),4,6,nbins=100,title="SOURCE DIVERGENCES")

    beam_tmp = beam.duplicate()
    beam_tmp.retrace(100.0)
    plotxy(beam_tmp.get_shadow3_beam(),1,3,nbins=100,title="SOURCE AFTER 10m")
    # plotxy(beam.get_shadow3_beam(),1,4,nbins=100,title="SOURCE PHASE SPACE")

    # p = 30
    # q = 5.0
    p = 10.0
    q = 10.0
    F = 1.0 / (1/p + 1/q)
    # lens1 = LensIdeal("test1",focal_x=F,focal_z=F,p=p,q=q)
    lens1 = LensSuperIdeal("test1",focal_p_x=p,focal_p_z=p,focal_q_x=q,focal_q_z=q,p=p,q=q)

    # p = 5.0
    # q = 30
    # F = 1.0 / (1/p + 1/q)
    # lens2 = LensIdeal("test2",focal_x=F,focal_z=F,p=p,q=q)

    method = 2 # 0:direct, 1:interface with overwrite, 2: no overwrite
    if method == 0:
        beam2 = lens1.trace_beam(beam)
    elif method == 1:
        beam.traceOE(lens1,1,overwrite=True)
        beam2 = beam
    elif method == 2:
        beam2 = beam.traceOE(lens1,1,overwrite=True)
        # beam2 = beam.traceOE(lens2,2,overwrite=True)
    else:
        raise Exception("Undefined method")

    #

    X = beam2.get_column(1)
    Y = beam2.get_column(3)
    print("X: ",X.min(),X.max(),X.std())
    print("Y: ",Y.min(),Y.max(),Y.std())
    # from srxraylib.plot.gol import plot_scatter
    # plot_scatter(X,Y)


    plotxy(beam2.get_shadow3_beam(),1,3,nbins=100,title="FOCAL PLANE",xrange=[-5e-9,5e-9],yrange=[-5e-9,5e-9])
    # plotxy(beam2.get_shadow3_beam(),1,4,nbins=100,title="FOCAL PLANE PHASE SPACE")

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

    from minishadow.source_geometrical.gaussian import SourceGaussian
    from minishadow.beam.beam import Beam
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

    src = SourceGaussian.initialize_from_keywords(number_of_rays=500000,
                                                  sigmaX=Sx,sigmaZ=Sz,
                                                  sigmaXprime=Sxp,sigmaZprime=Szp)

    # src = SourceGaussian.initialize_from_keywords(number_of_rays=10000,
    #                                               sigmaX=     1000e-6/2.35,sigmaZ=   10e-6/2.35,
    #                                               sigmaXprime=30e-6/2.35,sigmaZprime=15e-6/2.35)
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

    # plotxy(beam.get_shadow3_beam(),1,4,nbins=100,title="SOURCE H phase space")
    plotxy(beam.get_shadow3_beam(),1,3,nbins=100)

    # multilayer
    p = 28.3
    q = 11.70
    F = 1/(1/p+1/q)
    lens1 = LensIdeal("ML",focal_x=F,focal_z=0,p=p,q=q)
    # lens1 = LensSuperIdeal("ML",focal_p_x=p,focal_q_x=q,focal_p_z=0,focal_q_z=0,p=p,q=q)
    beam.traceOE(lens1,1,overwrite=True)

    # plotxy(beam.get_shadow3_beam(),1,4,nbins=100,title="H phase space")
    # plotxy(beam.get_shadow3_beam(),3,6,nbins=100,title="V phase space")


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
    lens2 = LensIdeal("KBV",focal_x=0,focal_z=F,p=p,q=q)
    # lens2 = LensSuperIdeal("KBV",focal_p_x=0,focal_q_x=0,focal_p_z=184.90,focal_q_z=0.10,p=p,q=q)
    beam.traceOE(lens2,1,overwrite=True)

    # second KB mirror
    p = 0.025
    q = 0.05
    F = 1.0/(1/144.95+1/0.05)
    lens3 = LensIdeal("KBH",focal_x=F,focal_z=0,p=p,q=q)
    # lens3 = LensSuperIdeal("KBH",focal_p_x=144.95,focal_q_x=0.05,focal_p_z=0,focal_q_z=0,p=p,q=q)
    beam.traceOE(lens3,1,overwrite=True)


    #

    tkt = plotxy(beam.get_shadow3_beam(),1,3,nbins=300,xrange=[-0.0000005,0.0000005],yrange=[-0.0000005,0.0000005],title="FOCAL PLANE")

    print(tkt['fwhm_h'],tkt['fwhm_v'])

    FX, FZ = (1e6*beam.get_standard_deviation(1),1e6*beam.get_standard_deviation(3))
    print("----------------- Focal position ---------------------")
    print("Source dimensions (rms): %f %f um"%(SX,SZ))
    print("Source dimensions (2.35*rms): %f %f um"%(f2dot35*SX,f2dot35*SZ))
    print("Focal dimensions (rms): %f %f um"%(FX,FZ))
    print("Focal dimensions (2.35*rms): %f %f um"%(f2dot35*FX,f2dot35*FZ))
    print("Focal dimensions (FWHM HISTO): %f %f nm"%(1e9*tkt['fwhm_h'],1e9*tkt['fwhm_v']))
    print("Demagnification (StDev) H:%g V:%g (theoretical: %g,%g) "%(SX/FX,SZ/FZ,demagX[0]*demagX[1],demagZ))
    print("Demagnification:(HISTO) H:%g V:%g (theoretical: %g,%g) "%(f2dot35*SX/(f2dot35*1e6*tkt['fwhm_h']),SZ/(1e6*tkt['fwhm_v']),demagX[0]*demagX[1],demagZ))

if __name__ == "__main__":
    # test_with_collimated_beam()
    # todo check with two elements
    # test_with_divergent_beam()
    test_id16ni()