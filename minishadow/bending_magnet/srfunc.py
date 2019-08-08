r"""

srfunc: calculates synchrotron radiation emission (radiation and angle 
        distributions) 

        functions: 

             basic:

             fintk53: integral of the Bessel function K5/3(x)
             sync_g1: energy spectrum integrated over the full vertical angle
             sync_f:  angular dependence of synchrotron radiation
             sync_hi: function Hi(x) = x^i * BeselK(x/2,2/3) 

             bending magnet radiation:

             sync_ang: angular distributions
             sync_ene: energy distributions


             wiggler radiation:

             wiggler_trajectory: computes the electron trajectroy in a
                                 magnetic field (sinusoidal or from a file
                                 with B[T] or its harmonic decomposition).
             wiggler_spectrum: computes the wiggler spectrum (full emission)
                               given the trajectory.

             wiggler_nphoton: computes the number of photons emtted versus
                              bending radius per (1mA 1mrad (horizontal) 0.1% bandwidth)
             wiggler_harmonics: computes the harmonic decomposition of the 
                              magnetic field B(s)

             wiggler_cdf: computes the cumulative distribution function from radiation

TODO:
- elliptical wigglers
- vectorize some loops
- remove scipy dependency

"""


__author__ = "Manuel Sanchez del Rio"
__contact__ = "srio@esrf.eu"
__copyright = "ESRF, 2002-2014"


import numpy
import scipy.special
import scipy.constants.codata
import os
import sys

#
#----------------------------  GLOBAL NAMES ------------------------------------
#

#Physical constants (global, by now)
try:
    import scipy.constants as codata
    codata_c = codata.c
    codata_mee = 1e-6 * codata.m_e * codata.c**2 / codata.e
    codata_me = codata.m_e
    codata_h = codata.h
    codata_ec = codata.e

except ImportError:
    #TODO update these constants
    print("Failed to import scipy. Finding alternative ways.")
    codata_c = numpy.array(299792458.0)
    codata_me = numpy.array(9.10938291e-31)
    codata_mee = numpy.array(0.510999)
    codata_h = numpy.array(6.62606957e-34)
    codata_ec = numpy.array(1.602176565e-19)

m2ev = codata_c*codata_h/codata_ec      # lambda(m)  = m2eV / energy(eV)


def fintk53(xd):
    r"""
     Calculates the integral from x to infinity of the Bessel function K5/3(x) 

      NAME: 
            fintk53 
      
      PURPOSE: 
        Calculates the function consisting of the integral, from x to infinity,  
         of the Bessel function K5/3(x).  
         g_one=fintk53*x is the universal curve from which the energy spectrum  
         of the synchrotron bending magnet is calculated. 
      
      CATEGORY: 
            Mathematics. 
          
      CALLING SEQUENCE: 
            Result = fintk53(x) 
      
      INPUTS: 
            x: the argument of the function. All calculations are done in doble 
                precision. 
      
      OUTPUTS: 
            returns the value  of the fintk53 function 
      
      PROCEDURE: 
            Translated from a Fortran program, original from Umstatter.  
          C 
          C Routines taken from  
          C 
          C http://www.slac.stanford.edu/grp/arb/tn/arbvol2/ARDB162.pdf 
          C The reference  1981 CERN/PS/SM/81-13  is cited in 
        C 'Synchrotron Radiation Spectra' by G.Brown and W. 
          C Lavender pp.37-61 Handbook on Synchrotron Radiation, 
        C vol. 3 edited by G.S. Brown and D.E. Moncton (Elsevier 
          C Science Publishers B.V. 1991 North-Holland, Amsterdam)  
          C  
      
            I have performed a comparison of the result with Mathematica 
            with very good agreement (note that Mathematica values diverge 
            for x> 20. I do not know why): 
            (Mathematica evaluation N[g_one[0.001],20]) 
            (x = 5.0 & print,fintk53(x),x*fintk53(x),format='(2G32.18)') 
      
          x      mathematica         idl (x*fintk53(x))       python x*fintk3(x)
        0.001 0.2131390650914501     0.213139096577768417     2.13139066e-01
        0.01  0.4449725041142102     0.444972550630643671     4.44972505e-01
        0.1   0.818185534872854      0.818185588215680770     8.18185536e-01
        1.0   0.6514228153553639697  0.651422821506926542     6.51422815e-01
        5.0   0.021248129774982      0.0212481300729910755    2.12481298e-02
        10.0  0.00019223826428       0.000192238266987909711  1.92238264e-04
        20.0  5.960464477539063E-7   1.19686346217633044E-08  1.19686345e-08
        50.0  6.881280000000002E7    1.73478522828932108E-21  1.73478520e-21
        100.0 4.642275147320176E29   4.69759373162073832E-43  4.69759367e-43
        1000.0 -1.7E424              Floating underflow (<620 OK) 0.00000000e+00
      
      
      MODIFICATION HISTORY: 
            Written by:     M. Sanchez del Rio, srio@esrf.fr, 2002-04-22 
            20120208 srio@esrf.eu: python version
      
    """ 
    #
    #C
    #C Computes Integral (from x to infinity) {K5/3(y) dy}
    #C
    xd = numpy.array(xd)
    oldshape=xd.shape
    xd.shape=-1
    a=1.0
    b=1.0
    p=1.0 
    q=1.0 
    x=1.0 
    y=1.0 
    z=1.0 
    
    xi = numpy.where(xd >= 5.0)
    xi = numpy.array(xi)
    count1 = xi.size
    
    fintk53=xd*0.0
    
    if (count1 > 0): 
        x = xd[xi]
        z=20./x-2.
        a= 0.0000000001
        b=z*a - 0.0000000004
        a=z*b-a + 0.0000000020
        b=z*a-b - 0.0000000110
        a=z*b-a + 0.0000000642
        b=z*a-b - 0.0000004076
        a=z*b-a + 0.0000028754
        b=z*a-b - 0.0000232125
        a=z*b-a + 0.0002250532
        b=z*a-b - 0.0028763680
        a=z*b-a + 0.0623959136
        p=0.5*z*a-b + 1.0655239080
        p=p* numpy.power(numpy.pi/2/x,0.5)/numpy.exp(x)
        fintk53[xi]=p
    
    xi = numpy.where(xd < 5.0)
    xi = numpy.array(xi)
    count2 = xi.size
    
    if ((count1+count2) != xd.size):
        print('Error: (count1+count2) NE N_Elements(xd)')
        print(count1)
        print(count2)
        raise ValueError("Error: (count1+count2) != size(xd)=%1" % xd.size)
    
    if (count2 > 0):
        x = xd[xi]
        z=numpy.power(x,2)/16.-2.
        a= 0.0000000001
        b=z*a + 0.0000000023
        a=z*b-a + 0.0000000813
        b=z*a-b + 0.0000024575
        a=z*b-a + 0.0000618126
        b=z*a-b + 0.0012706638
        a=z*b-a + 0.0209121680
        b=z*a-b + 0.2688034606
        a=z*b-a + 2.6190218379
        b=z*a-b + 18.6525089687
        a=z*b-a + 92.9523266592
        b=z*a-b + 308.1591941313
        a=z*b-a + 644.8697965824
        p=0.5*z*a-b + 414.5654364883
        a= 0.0000000012
        b=z*a + 0.0000000391
        a=z*b-a + 0.0000011060
        b=z*a-b + 0.0000258145
        a=z*b-a + 0.0004876869
        b=z*a-b + 0.0072845620
        a=z*b-a + 0.0835793546
        b=z*a-b + 0.7103136120
        a=z*b-a + 4.2678026127
        b=z*a-b + 17.0554078580
        a=z*b-a + 41.8390348678
        q=0.5*z*a-b + 28.4178737436
        y=numpy.power(x,0.666666667)
        p=(p/y-q*y-1.)*1.8137993642
        fintk53[xi]=p

        fintk53.shape=oldshape
    return fintk53


def sync_g1(x,polarization=0):
    r"""
    calculates the synchrotron radiation g1 function
    
      NAME:
            sync_g1
     
      PURPOSE:
            Calculates the functions used for calculating synchrotron
         radiation energy spectrum integrated over the full vertical
         angle.
     
      CATEGORY:
            Mathematics.
     
      CALLING SEQUENCE:
            Result = sync_g1(x)
     
      INPUTS:
            x:      the argument of the function. It is converted to double
             precision for calculations. 
     
      KEYWORD PARAMETERS:
         POLARIZATION: 0 Total 
                   1 Parallel       (l2=1, l3=0, in Sokolov&Ternov notation)
                   2 Perpendicular  (l2=0, l3=1)
      OUTPUTS:
            returns the value  of the sync_g1 function
     
      PROCEDURE:
            The number of emitted photons versus energy is:
         N(E) = 2.4605e13 I[A] Ee[Gev] Theta[mrad] Sync_G1(E/Ec]
            Where: 
             I is the storage ring intensity in A
             Ee is the energy of the electrons in the storage ring 
             E is the photon energy
             Ec is the critical energy
             The value Sync_G1 returned by this function is:
                 sync_g1(x) (total polarization):
                 x* Integrate[BeselK[x,5/3],{x,y,Infinity}]
                 sync_g1(x,Pol=1) (parallel polarization):
                 (1/2)* [x* Integrate[BesselK[x,5/3],{x,y,Infinity}] + 
                 x*BesselK(x,2/3)]
                 sync_g1(x,Pol=2) (perpendicular polarization):
                 (1/2)* [x* Integrate[BesselK[x,5/3],{x,y,Infinity}] -
                 x*BesselK(x,2/3)]
     
         For calculating the Integrate[BeselK[x,5/3],{x,y,Infinity}]
                 function, the function fintk53 is used. 
     
         Reference: A A Sokolov and I M Ternov, Synchrotron Radiation, 
                 Akademik-Verlag, Berlin, 1968, Formula 5.19, 
                 pag 32.
     
      MODIFICATION HISTORY:
            Written by:     M. Sanchez del Rio, srio@esrf.fr, 2002-05-24
            20120208 srio@esrf.eu: python version
      
    """
    #y = fintk53(x)*x
    x = numpy.array(x)
    y = fintk53(x)
    y = y*x
    if polarization == 0:
        return y
    
    if polarization == 1:
        #return 0.5*(y+(x*BeselK(x,2.0/3.0)))
        return 0.5*(y+(x*scipy.special.kv(2.0/3.0,x)))
    
    if polarization == 2:
        #return 0.5*(y-(x*BeselK(x,2.0/3.0)))
        return 0.5*(y-(x*scipy.special.kv(2.0/3.0,x)))
    
    raise ValueError("invalid polarization=: %s" % polarization)


def sync_f(rAngle,rEnergy=None,polarization=0,gauss=0,l2=1,l3=0 ):
    r""" angular dependency of synchrotron radiation emission

      NAME:
            sync_f
     
      PURPOSE:
            Calculates the function used for calculating the angular 
         dependence of synchrotron radiation. 
     
      CATEGORY:
            Mathematics.
     
      CALLING SEQUENCE:
            Result = sync_f(rAngle [,rEnergy] )
     
      INPUTS:
            rAngle:  the reduced angle, i.e., angle[rads]*Gamma. It can be a
             scalar or a vector.
      OPTIONAL INPUTS:
            rEnergy:  a value for the reduced photon energy, i.e., 
             energy/critical_energy. It can be an scalar or a verctor. 
             If this input is present, the calculation is done for this 
             energy. Otherwise, the calculation results is the integration 
             over all photon energies.
     
      KEYWORD PARAMETERS:
         POLARIZATION: 0 Total 
                   1 Parallel       (l2=1, l3=0, in Sokolov&Ternov notation)
                   2 Perpendicular  (l2=0, l3=1)
                   3 Any            (define l2 and l3)
     
         l2: The polarization value of L2
         l3: The polarization value of L3
             Note: If using L2 and L3, both L2 and L3 must be defined.
                   In this case, the Pol keyword is ignored.
     
         GAUSS: When this keyword is set, the "Gaussian" approximaxion 
             instead of the full calculation is used. 
             Only valid for integrated flux aver all photon energies.
      OUTPUTS:
            returns the value  of the sync_f function
             It is a scalar if both inputs are scalar. If one input
             is an array, the result is an array of the same dimension. 
             If both inputs are arrays, the resulting array has dimension
             NxM, N=Dim(rAngle) and M=Dim(rEnergy)
     
      PROCEDURE:
            The number of emitted photons versus vertical angle Psi is
         proportional to sync_f, which value is given by the formulas
         in the references.
     
         For angular distribution integrated over full photon energies (rEnergy 
         optional input not present) we use the Formula 9, pag 4 in Green. 
         For its gaussian approximation (in this case the polarization keyword 
         has no effect) we use for 87 in pag 32 in Green.
     
         For angular distribution at a given photon energy (rEnergy 
         optional input not present) we use the Formula 11, pag 6 in Green. 
     
     
         References: 
             G K Green, "Spectra and optics of synchrotron radiation" 
                 BNL 50522 report (1976)
             A A Sokolov and I M Ternov, Synchrotron Radiation, 
                 Akademik-Verlag, Berlin, 1968
     
      OUTPUTS:
            returns the value  of the sync_hi function
     
      PROCEDURE:
            Uses IDL's BeselK() function
     
      MODIFICATION HISTORY:
            Written by:     M. Sanchez del Rio, srio@esrf.fr, 2002-05-23
         2002-07-12 srio@esrf.fr adds circular polarization term for 
             wavelength integrated spectrum (S&T formula 5.25)
         2012-02-08 srio@esrf.eu: python version
      
    """
    # auto-call for total polarization
    if polarization == 0:
        return sync_f(rAngle,rEnergy,polarization=1)+ \
               sync_f(rAngle,rEnergy,polarization=2)

    rAngle=numpy.array(rAngle)
    rAngle.shape=-1
    
    if polarization == 1:
        l2=1.0
        l3=0.0

    if polarization == 2:
        l2=0.0
        l3=1.0
    
    #;
    #; angle distribution integrated over full energies
    #;
    if (rEnergy is None):
        if gauss == 1:
            #; Formula 87 in Pag 32 in Green 1975
            efe = 0.4375*numpy.exp(-0.5* numpy.power(rAngle/0.608,2) )
            return efe
    
        #if polarization == 0:
        #    return sync_f(rAngle,polarization=1)+sync_f(rAngle,polarization=2)
        #
        #; For 9 in Pag 4 in Green 1975
        #; The two summands correspond to the par and per polarizations, as 
        #; shown in Sokolov&Ternov formulas (3.31) and 5.26)
        #; However, for circular polarization a third term (S&T 5.25) 
        #; must also be used
        efe = (7.0/16.0)*l2*l2+ \
        (5.0/16.0)*(rAngle*rAngle/(1.0+rAngle*rAngle))*l3*l3 + \
        (64.0/16.0/numpy.pi/numpy.sqrt(3.0))* \
        (rAngle/numpy.power(1+rAngle*rAngle,0.5))*l2*l3
        efe = efe * ( numpy.power(1.0+rAngle*rAngle,-5.0/2.0) )
        return efe

    #;
    #; angle distribution for given energy/ies
    #;
    rEnergy=numpy.array(rEnergy)
    rEnergy.shape=-1
    #
    #; For 11 in Pag 6 in Green 1975
    #
    ji = numpy.sqrt( numpy.power(1.0+numpy.power(rAngle,2),3) )
    ji = numpy.outer(ji,rEnergy/2.0)
    rAngle2 = numpy.outer(rAngle,(rEnergy*0.0+1.0))
    efe = l2*scipy.special.kv(2.0/3.0,ji)+ \
            l3* rAngle2*scipy.special.kv(1.0/3.0,ji)/ \
            numpy.sqrt(1.0+numpy.power(rAngle2,2))
    efe = efe* (1.0+numpy.power(rAngle2,2))
    efe = efe*efe

    return efe


def sync_hi(x,i=2,polarization=0): 
    r""" calculates the function Hi(x) used for Synchrotron radiation

      NAME:
            sync_hi
     
      PURPOSE:
            Calculates the function Hi(x) used for Synchrotron radiation 
         Hi(x) = x^i * BesselK(x/2,2/3) (for total polarization)
     
     
      CATEGORY:
            Mathematics.
     
      CALLING SEQUENCE:
            Result = sync_hi(x [,i] )
     
      INPUTS:
            x:   the argument of the function. All calculations are done 
                 in doble precision.
             i:    the exponent. If this optional argument is not entered, it 
             is set to 2.
     
      KEYWORD PARAMETERS:
         POLARIZATION: 0 Total 
                   1 Parallel       (l2=1, l3=0, in Sokolov&Ternov notation)
                   2 Perpendicular  (l2=0, l3=1)
     
      OUTPUTS:
            returns the value  of the sync_hi function
     
      PROCEDURE:
            Uses the relation ship Hi(x) =  x^i * sync_f(0,x)
     
      MODIFICATION HISTORY:
            Written by:     M. Sanchez del Rio, srio@esrf.fr, 2002-05-23
            20120208 srio@esrf.eu: python version
      
    """
    x=numpy.array(x)
    x.shape=-1
    y1 = numpy.power(x,i) * sync_f(0,x,polarization=polarization)
    return y1


def sync_ang(flag,angle_mrad,polarization=0, \
    e_gev=1.0,i_a=0.001,hdiv_mrad=1.0,r_m=1.0,energy=1.0,ec_ev=1.0):
    r""" Calculates the synchrotron radiation angular distribution

      NAME:
            sync_ang
     
      PURPOSE:
            Calculates the synchrotron radiation angular distribution
     
      CATEGORY:
            Mathematics.
     
      CALLING SEQUENCE:
            Result = sync_ang(flag, angle )
     
      INPUTS:
         flag:     0 Flux fully integrated in photon energy
             1 Flux at a given photon energy
        angle:  the angle array [in mrad]
     
      KEYWORD PARAMETERS:
         polarization: 0 Total 
                   1 Parallel       (l2=1, l3=0, in Sokolov&Ternov notation)
                   2 Perpendicular  (l2=0, l3=1)
     
         IF flag=0 THE FOLLOWING KEYWORDS MUST BE ENTERED
             e_geV= The electron energy [in GeV]  (default=1.0)
             i_a= the electron beam intensity [in A] (default=1.0D-3)
             hdiv_mrad= the horizontal divergence [in mrad] (default=1)
             r_m= the bending magnet radius [in m] (default=1.0)
     
         IF flag=1 THE FOLLOWING KEYWORDS MUST BE ENTERED
             All keywords for FLAG=0, except r_m,  plus:
             energy = the energy value [in eV] (default=1)
             ec_ev= The critical energy [eV] (default=1)
     
      OUTPUTS:
            returns the array with the angular distribution 
             IF flag ==1 power density [Watts/mrad(Psi)]
     
      PROCEDURE:
     
         References: 
             G K Green, "Spectra and optics of synchrotron radiation" 
                 BNL 50522 report (1976)
             A A Sokolov and I M Ternov, Synchrotron Radiation, 
                 Akademik-Verlag, Berlin, 1968
     
      MODIFICATION HISTORY:
            Written by:     M. Sanchez del Rio, srio@esrf.fr, 2002-06-03
            20120208 srio@esrf.eu: python version
      
    """

    angle_mrad = numpy.array(angle_mrad)
    angle_mrad.shape = -1
    
    if flag == 0:
        # fully integrated in photon energy
        a8 = 3e10*codata_c*codata_ec/numpy.power(codata_mee,5) # 41.357
        gamma = e_gev*1e3/codata_mee
        a5 = sync_f(angle_mrad*gamma/1e3,polarization=polarization)* \
             a8*i_a*hdiv_mrad/r_m*numpy.power(e_gev,5)
        return a5
    
    if flag == 1:
        #a8 = 1.3264d13
        a8 = codata_ec/numpy.power(codata_mee,2)/codata_h*(9e-2/2/numpy.pi) 
        energy = numpy.array(energy)
        energy.shape=-1
        eene = energy/ec_ev
        gamma = e_gev*1e3/codata_mee
        a5=sync_f(angle_mrad*gamma/1e3,eene,polarization=polarization)* \
            numpy.power(eene,2)*a8*i_a*hdiv_mrad*numpy.power(e_gev,2)
        return a5


def sync_ene(f_psi,energy_ev,ec_ev=1.0,polarization=0,  \
             e_gev=1.0,i_a=0.001,hdiv_mrad=1.0, \
             psi_min=0.0, psi_max=0.0, psi_npoints=1): 
    r""" Calculates the synchrotron radiation energy spectrum
      NAME:
            sync_ene
     
      PURPOSE:
            Calculates the synchrotron radiation energy spectrum
     
      CATEGORY:
            Mathematics.
     
      CALLING SEQUENCE:
            Result = sync_ene(flag, Energy )
     
      INPUTS:
         flag:     0 Flux fully integrated in angle (Psi)
             1 Flux at Psi=0
             2 Flux integrated in the angular interval [Psi_Min,Psi_Max]
             3 Flux at Psi=Psi_Min
             4 Flux versus angle in [Psi_Min,Psi_Max] and energy

            energy:  the energy array [in eV]
     
      KEYWORD PARAMETERS:
         polarization: 0 Total 
                   1 Parallel       (l2=1, l3=0, in Sokolov&Ternov notation)
                   2 Perpendicular  (l2=0, l3=1)
     
         If flag=0 or flag=1 the following keywords MUST BE ENTERED
     
             ec_ev= The critical energy [eV]
             e_geV= The electron energy [in GeV] 
             i_a= the electron beam intensity [in A]
             hdiv_mrad= the horizontal divergence [in mrad]
     
         If flag=2, in addition to the mentioned keywords, the following 
             ones must be present:
     
             Psi_Min the minimum integration angle [in mrad]
             Psi_Max the maximum integration angle [in mrad]
             Psi_NPoints the number of points in psi for integration
     
         If flag=3, in addition to the mentioned keywords for flag=0 OR 
             flag=1, the following kewford must be defined: 
     
             psi_min the Psi angular value [in mrad]
     
      KEYWORD PARAMETERS (OUTPUT):
     
         IF flag=2, the following keywords can be used to obtain additional info:
     
             fmatrix=a two dimensional variable containing the matrix of 
                 flux as a function of angle [first index] and energy 
                 [second index]
             angle_mrad= a: one-dim array with the angular points [in mrad]
     
      OUTPUTS:
            returns the array with the flux [photons/sec/0.1%bw] for FLAG=0,2
            and the flux [photons/sec/0.1%bw/mrad] for FLAG=1,3,4
     
      PROCEDURE:
     
         References: 
             G K Green, "Spectra and optics of synchrotron radiation" 
                 BNL 50522 report (1976)
             A A Sokolov and I M Ternov, Synchrotron Radiation, 
                 Akademik-Verlag, Berlin, 1968
     
      EXAMPLE:
         The following program was used for testing sync_ene
         
         
        #create 10-points energy array in [20,30] keV
        e=numpy.linspace(20000.0,30000.0,10)
         
         ;
         ; test of spectra at Psi=0
         ;
         ; at psi=0 (i.e., flag=1)
        In [274]: srfunc.sync_ene(1,e,ec_ev=19166.0,e_gev=6,i_a=0.1,hdiv_mrad=1)
        Out[274]: 
              array([[  6.89307648e+13,   6.81126315e+13,   6.71581119e+13,
          6.60866137e+13,   6.49155481e+13,   6.36605395e+13,
          6.23356084e+13,   6.09533305e+13,   5.95249788e+13,
          5.80606485e+13]])


         ; at psi_min (FLAG=3)
        In [279]: srfunc.sync_ene(3,e,ec_ev=19166.0,e_gev=6,i_a=0.1, \
                  hdiv_mrad=1,psi_min=0.0)
        Out[279]: 
        array([[  6.89307648e+13,   6.81126315e+13,   6.71581119e+13,
          6.60866137e+13,   6.49155481e+13,   6.36605395e+13,
          6.23356084e+13,   6.09533305e+13,   5.95249788e+13,
          5.80606485e+13]])

         ;
         ; test of integrated spectra 
         ;
         
         ; Integrating (by hand) using flag=3
        # a is large enough to cover the full radiation fan
         a = numpy.linspace(-0.2,0.2,50) 
        #create 10-points energy array in [20,30] keV
        e=numpy.linspace(20000.0,30000.0,10)
         
        y3=e*0.0
         for i in range(a.size): 
             y2=srfunc.sync_ene(3,e,ec_ev=19166.0,e_gev=6,i_a=0.1,hdiv_mrad=1,psi_min=a[i])
             y3[i] = y3[i] + y2
         y3=y3*(a[1]-a[0])
         
         ; Integrating (automatically) using FLAG=2
         y4 = srfunc.sync_ene(2,e,ec_ev=19166.0,e_gev=6,i_a=0.1,hdiv_mrad=1,\
        psi_min=-0.2,psi_max=0.2,psi_npoints=50)
         
         ; Integrated (over all angles) using FLAG=0
         y5 = srfunc.sync_ene(0,e,ec_ev=19166.0,e_gev=6,i_a=0.1,hdiv_mrad=1)
         
        In [475]: for i in range(y3.size):
            print e[i],y3[i],y4[i],y5[i]
           .....:     
           .....:     

         The results obtained are: 
        energy        int_by_hand       int_num           int
        20000.0       9.32554203564e+12 9.32554203564e+12 9.33199803948e+12
        21111.1111111 8.95286605221e+12 8.95286605221e+12 8.9590640148e+12
        22222.2222222 8.58856640727e+12 8.58856640727e+12 8.59451215453e+12
        23333.3333333 8.2334342483e+12 8.2334342483e+12 8.2391341364e+12
        24444.4444444 7.88805461639e+12 7.88805461639e+12 7.89351540031e+12
        25555.5555556 7.55284456882e+12 7.55284456882e+12 7.55807329003e+12
        26666.6666667 7.22808379127e+12 7.22808379127e+12 7.23308768405e+12
        27777.7777778 6.91393939677e+12 6.91393939677e+12 6.91872581084e+12
        28888.8888889 6.61048616971e+12 6.61048616971e+12 6.61506250643e+12
        30000.0       6.31772320182e+12 6.31772320182e+12 6.32209686189e+12

         EXAMPLE 2
             Surface plot of flux versus angle ane energy
             e = numpy.linspace(20000,30000,20)
             tmp1,fm,a = srfunc.sync_ene(2,e,ec_ev=19166.0,e_gev=6,i_a=0.1,\
                      hdiv_mrad=1,psi_min=-0.2,psi_max=0.2,psi_npoints=50)
             surface,fm,a,e, ztitle='Flux[phot/sec/0.1%bw/mradPsi', $
                 xtitle='Angle [mrad]',ytitle='Energy [eV]'

      MODIFICATION HISTORY:
            Written by:     M. Sanchez del Rio, srio@esrf.fr, 2002-06-03
         2007-05-14 srio@esrf.fr debug with FLAG=2. The bandwith in 
             angle depends on the number of points. Now it is 1mrad
             Bug reported by flori@n-nolz.de
             Added default values. 
        2007-12-13  srio@esrf.eu fixes bug reported by Gernot.Buth@iss.fzk.de
                    concerning the normalization of the angular integral.
     
        20120208 srio@esrf.eu: python version
      
     -
    """

    energy_ev = numpy.array(energy_ev)
    oldshape = energy_ev.shape
    energy_ev.shape = -1

    if f_psi == 0: # fully integrated in Psi
    # numerical cte for integrated flux
        a8 = numpy.sqrt(3e0)*9e6*codata_ec/codata_h/codata_c/codata_mee
        a5 = a8*e_gev*i_a*hdiv_mrad* \
             sync_g1(energy_ev/ec_ev,polarization=polarization)
        #TODO: check this 
        if len(energy_ev) == len(a5): 
            a5.shape = oldshape
    elif f_psi == 1: #at Psi = 0
        #a8 =  1.3264d13
        a8 = codata_ec/numpy.power(codata_mee,2)/codata_h*(9e-2/2/numpy.pi) 
        a5 = a8*numpy.power(e_gev,2)*i_a*hdiv_mrad* \
             sync_hi(energy_ev/ec_ev,polarization=polarization)
        a5.shape = oldshape
    elif f_psi == 2: #between PsiMin and PsiMax
        # a8 = 1.3264d13
        a8 = codata_ec/numpy.power(codata_mee,2)/codata_h*(9e-2/2/numpy.pi) 
        eene = energy_ev/ec_ev
        gamma = e_gev*1e3/codata_mee
        angle_mrad = numpy.linspace(psi_min,psi_max,psi_npoints)
        eene2 = numpy.outer(angle_mrad*0.0e0+1,eene)
        a5=sync_f(angle_mrad*gamma/1e3,eene,polarization=polarization)

        
        a5 = a5*numpy.power(eene2,2)*a8*i_a*hdiv_mrad*numpy.power(e_gev,2)
        fMatrix = a5
        #a5 = Total(fMatrix,1) 
        # corrected srio@esrf.eu 2007/12/13 
        # bug reported by Gernot.Buth@iss.fzk.de
        angle_step = (float(psi_max)-psi_min)/(psi_npoints-1.0)
        a5 = fMatrix.sum(axis=0) * angle_step
    elif f_psi == 3: #at PsiMin
        a8 = codata_ec/numpy.power(codata_mee,2)/codata_h*(9e-2/2/numpy.pi) 
        #a8 = 1.3264d13
        eene = energy_ev/ec_ev
        gamma = e_gev*1e3/codata_mee
        angle_mrad = psi_min
        a5=sync_f(angle_mrad*gamma/1e3,eene,polarization=polarization)
        a5 = a5*numpy.power(eene,2)*a8*i_a*hdiv_mrad*numpy.power(e_gev,2)
        a5.shape = oldshape
    elif f_psi == 4: #at every point in [PsyMin,PsiMax]
        # a8 = 1.3264d13
        a8 = codata_ec/numpy.power(codata_mee,2)/codata_h*(9e-2/2/numpy.pi)
        eene = energy_ev/ec_ev
        gamma = e_gev*1e3/codata_mee
        angle_mrad = numpy.linspace(psi_min,psi_max,psi_npoints)
        eene2 = numpy.outer(angle_mrad*0.0e0+1,eene)
        a5=sync_f(angle_mrad*gamma/1e3,eene,polarization=polarization)

        a5 = a5*numpy.power(eene2,2)*a8*i_a*hdiv_mrad*numpy.power(e_gev,2)

    else:
        raise Exception("Invalid value f_psi=%d"%f_psi)

    return a5

#
#------------------------- WIGGLER FUNCTIONS -----------------------------------
#
def wiggler_spectrum(traj, enerMin=1000.0, enerMax=100000.0, nPoints=100, \
                     per=0.2, electronCurrent=0.2, outFile="", elliptical=False):
    r"""
     NAME:
           wiggler_spectrum
    
     PURPOSE:
           Calculates the spectrum of a wiggler using a trajectory file
    	as input. 
    
     CATEGORY:
           Synchrotron radiation
    
     CALLING SEQUENCE:
    	wiggler_spectrum(traj)
    
     INPUTS:
    
           traj:      The array with the electron trajectory
    		       (created by wiggler_trajectory)
           enerMin:    Minimum photon energy [eV]
           enerMax:    Maximum photon energy [eV]
           nPoints:     Number of energy points
    	   electronCurrent:     The electron beam current in mA
    	   per:         The ID period in m
    
           outFile:    The name of the file with results
           elliptical: False (for elliptical wigglers, not yet implemented)
    
     OUTPUTS:
           An array with the resulting spectrum
    
     PROCEDURE:
    	Based on SHADOW's wiggler_spectrum. Uses wiggler_nphoton
    
     MODIFICATION HISTORY:
           Written by:     M. Sanchez del Rio, srio@esrf.fr, 2002-07-15
    	2002-07-18 srio@esrf.fr adds doc. Use "current" value 
    	2006-06-18 srio@esrf.fr uses hc from Physical_Constants()
        2012-10-08 srio@esrf.eu python version
    
    """

    x = traj[0,:]
    y = traj[1,:]
    z = traj[2,:]
    betax = traj[3,:]
    betay = traj[4,:]
    betaz = traj[5,:]
    curv =  traj[6,:]
    b_t =  traj[7,:]

    step = numpy.sqrt(numpy.power(y[2]-y[1],2) + numpy.power(x[2]-x[1],2) +  \
           numpy.power(z[2]-z[1],2))
    #;
    #; Compute gamma and the beam energy
    #;
    gamma = 1.0/numpy.sqrt(1 - numpy.power(betay[1],2) - \
                               numpy.power(betax[1],2) - \
                               numpy.power(betaz[1],2))
    bener = gamma*(9.109e-31)*numpy.power(2.998e8,2)/(1.602e-19)*1.0e-9
    print("\nElectron beam energy (from velocities) = %f GeV "%(bener))
    print("\ngamma (from velocities) = %f "%(gamma))


    #;
    #; Figure out the limit of photon energy.
    #;

    curv_max = numpy.abs(curv).max()
    curv_min = numpy.abs(curv).min()


    print("Radius of curvature (max) = %f m "%(1.0/curv_min))
    print("                    (min) = %f m "%(1.0/curv_max))

    TOANGS  =  m2ev*1e10 
    phot_min = TOANGS*3.0*numpy.power(gamma,3)/4.0/numpy.pi/1.0e10*curv_min
    phot_max = TOANGS*3.0*numpy.power(gamma,3)/4.0/numpy.pi/1.0e10*curv_max

    print("Critical Energy (max.) = %f eV"%(phot_max))
    print("                (min.) = %f eV"%(phot_min))

    out = numpy.zeros((3,nPoints))
    #;
    #; starts the loop in energy
    #;
    mul_fac=numpy.abs(curv)*numpy.sqrt(1+numpy.power(betax/betay,2)+ \
                                         numpy.power(betaz/betay,2))*1.0e3
    energy_array = numpy.linspace(enerMin,enerMax,nPoints) 
    for i in range(len(energy_array)):
        energy = energy_array[i]

        #
        #;
        #; wiggler_nphoton computes the no. of photons per mrad (ANG_NUM) at 
        #; each point. It is then used to generate the no. of photons per axial 
        #; length (PHOT_NUM) along the trajectory S.
        phot_num = numpy.zeros(len(curv))
        rad= numpy.abs(1.0/curv)
        ang_num = wiggler_nphoton(rad,electronEnergy=bener,photonEnergy=energy)
        phot_num=ang_num*mul_fac
        #;
        #; Computes CDF of the no. of photon along the trajectory S.
        #; In the elliptical case, the entire traversed path length (DS) is 
        #; computed. In the normal case, only the component (Y) in the direction 
        #; of propagation computed.
        #;
        #

        i_wig = 1 # 1=normalWiggler, 2=ellipticalWiggler (not implemented)

        if i_wig == 2:
            dx = x-numpy.roll(x,1)
            dy = y-numpy.roll(y,1)
            dz = z-numpy.roll(z,1)
            ds1 = numpy.sqrt(dx*dx+dy*dy+dz*dz) 
            ds1[0]=0.0
            ds = numpy.zeros(np)
            for j in range(len(curv)): 
                ds[j]=numpy.sum(ds1[0:j])
            phot_cdf=(numpy.roll(phot_num,1)+phot_num)*0.5e0*(ds-numpy.roll(ds,1))
            phot_cdf[0]=0.0
            tot_num = numpy.sum(phot_cdf)
        else:
            phot_cdf = (numpy.roll(phot_num,1)+phot_num)*0.5e0*(y-numpy.roll(y,1))
            phot_cdf[0]=0.0
            tot_num = numpy.sum(phot_cdf)

        out[0,i] = energy
        out[1,i] = tot_num*(electronCurrent*1e3)
        out[2,i] = tot_num*(electronCurrent*1e3)*codata_ec*1e3

    if outFile != "":
        f = open(outFile,"w")
        f.write("#F "+outFile+"\n")
        f.write("\n#S 1 Wiggler spectrum\n")
        f.write("#N 3\n")
        f.write("#L PhotonEnergy[eV]  Flux[phot/s/0.1%bw]  PowerDensity[W/eV]\n")  
        for j in range(out.shape[1]): 
            f.write(("%19.12e  "*out.shape[0]+"\n")%tuple(out[i,j] for i in range(out.shape[0])))
        f.close()
        print("File with wiggler spectrum written to file: "+outFile)

    return out[0,:].copy(), out[1,:].copy(), out[2,:].copy()


def wiggler_cdf(traj, enerMin=10000.0, enerMax=10010.0, enerPoints=101, \
                outFile="", elliptical=False):
    r"""
     NAME:
           wiggler_cdf
    
     PURPOSE:
           Calculates the cdf (cumulative density function) of a wiggler 
           using a trajectory as input. 
    
     CATEGORY:
           Synchrotron radiation. Shadow preprocessors
    
     CALLING SEQUENCE:
    	wiggler_spectrum(traj)
    
     INPUTS:
    
           traj:      The array with the electron trajectory
    		       (created by wiggler_trajectory)
           enerMin:    Minimum photon energy [eV]
           enerMax:    Maximum photon energy [eV]
           enerPoints: Number of points in photon energy, for internal 
                       integration.
           outFile:    If != "", the name of the file where results
                       are stored in ASCII, in the format: 

                       np step bener 1.0/curv_max 1.0/curv_min enerMin enerMax
                       np fields with: 
                           x[i] y[i] cdf[i] angle[i] curv[i]

                       where
                           np: number of points in trajectory
                           step: the step in y in m  from trajectory
                           bener: the electron energy in GeV
                           curv_max: trajectory maximum curvature in m^-1
                           curv_min: trajectory maximum curvature in m^-1

                           x,y: electron coordinates in m
                           cdf: cumulative distribution function
                           angle: the deviation angle of the electron in rad
                           curv: the trajectory curvature in m^-1

           elliptical: False (for elliptical wigglers, not yet implemented)
    
     OUTPUTS:
           None. If wanted (as usual), a file with the resulting info to be 
           used by Shadow.
    
     PROCEDURE:
    	   Based on SHADOW's one. Uses wiggler_nphoton 
    
     MODIFICATION HISTORY:
           Written by:     M. Sanchez del Rio, srio@esrf.eu, 2014-10-21
    
    """

    x = traj[0,:]
    y = traj[1,:]
    z = traj[2,:]
    betax = traj[3,:]
    betay = traj[4,:]
    betaz = traj[5,:]
    curv =  traj[6,:]
    b_t =  traj[7,:]

    np = len(x)

    step = numpy.sqrt(numpy.power(y[2]-y[1],2) +  \
                      numpy.power(x[2]-x[1],2) +  \
                      numpy.power(z[2]-z[1],2))
    #;
    #; Compute gamma and the beam energy
    #;
    gamma = 1.0/numpy.sqrt(1.0 - numpy.power(betay[1],2) - \
                                 numpy.power(betax[1],2) - \
                                 numpy.power(betaz[1],2))
    bener = gamma*(9.109e-31)*numpy.power(2.998e8,2)/(1.602e-19)*1.0e-9
    print("\nwiggler_cdf: Electron beam energy (from velocities) = %f GeV "%(bener))
    print("\nwiggler_cdf: gamma (from velocities) = %f GeV "%(gamma))


    #;
    #; Figure out the limit of photon energy.
    #;
    curv_max = 0.0
    curv_min = 1.0e20

    curv_max = numpy.abs(curv).max()
    curv_min = numpy.abs(curv).min()

    if curv_min == 0.0: curv_min = 1e-15

    print("wiggler_cdf: Curvature (min)) = %f m^-1 "%(curv_min))
    print("wiggler_cdf:           (max)    %f m^-1 "%(curv_max))
    print("wiggler_cdf: Radius of curvature (max) = %f m "%(1.0/curv_min))
    print("wiggler_cdf:                     (min) = %f m "%(1.0/curv_max))


    TOANGS  =  m2ev*1e10 
    phot_min = TOANGS*3.0*numpy.power(gamma,3)/4.0/numpy.pi/1.0e10*curv_min
    phot_max = TOANGS*3.0*numpy.power(gamma,3)/4.0/numpy.pi/1.0e10*curv_max

    print("wiggler_cdf: Critical Energy (max.) = %f eV"%(phot_max))
    print("wiggler_cdf:                 (min.) = %f eV"%(phot_min))


    #TODO: here it is necessary to define an array in energy for
    # performing the energy integration via wiggler_nphoton.
    # Originally, in Shadow, nphoton gives already the integrated 
    # values. It could be possible to make more precise integration
    # by finding the parametrized integral of sync_ene()
    # (basically the integrat of G1). May be it can be parametrized using
    # Mathematica? 
    phot_num = numpy.zeros(np) 
    energy_array = numpy.linspace(enerMin,enerMax,enerPoints) 
    energy_step = energy_array[1] - energy_array[0]


    rad = numpy.abs(1.0/curv)
    ang_num = numpy.zeros(len(curv))

    for j in range(len(energy_array)):
        tmp = wiggler_nphoton(rad,electronEnergy=bener,\
                      photonEnergy=energy_array[j],polarization=0)
        ang_num += tmp / (0.001 * energy_array[j]) * energy_step

    phot_num = ang_num*numpy.abs(curv)*numpy.sqrt(1.0+\
               numpy.power((betax/betay),2) + \
               numpy.power((betaz/betay),2))*1.0e3

    #!C
    #!C Computes CDF of the no. of photon along the trajectory S.
    #!C In the elliptical case, the entire traversed path length (DS) 
    #!C is computed.
    #!C In the normal case, only the component (Y) in the direction of 
    #!C propagation is computed.
    #!C
    if False:  # loop version
        phot_cdf = numpy.zeros(np) 
        for i in range(1,np):
            if elliptical: 
                ds[1] = 0.0
                dx[i] = x[i] - x[i-1]
                dy[i] = y[i] - y[i-1]
                dz[i] = z[i] - z[i-1]
                ds[i] = numpy.sqrt(numpy.power(dx[i],2) + numpy.power(dy[i],2) + numpy.power(dz[i],2)) + ds[i-1]
                phot_cdf[i]   = phot_cdf[i-1] +  \
                    (phot_num[i-1] + phot_num[i]) * 0.5 * (ds[i] - ds[i-1])
            else:
                phot_cdf[i]   = phot_cdf[i-1] +  \
                    (phot_num[i-1] + phot_num[i]) * 0.5 * (y[i] - y[i-1])
    else: # vector version
        if elliptical: 
            pass  # TODO: fill this part
        else:
            y0  = numpy.roll(y,1) 
            y0[0] = y0[1]
            phot_num0 = numpy.roll(phot_num,1)
            phot_num0[0] = phot_num0[1]
            phot_cdf = numpy.cumsum(phot_num0+phot_num) * 0.5 * ( y - y0 ) 

    tot_num   = phot_cdf[-1]
    print("wiggler_cdf: Total no.of photons = %e "%(tot_num))

    cdf = phot_cdf / tot_num
    angle = numpy.arctan2( betax, betay )

    if outFile != "":
        f = open(outFile,"w")
        #NP_TRAJ,PATH_STEP,BENER,RAD_MIN,RAD_MAX,PH1,PH2
        f.write("%d %e %e %e %e %e %e \n"%(np,step,bener,1.0/curv_max, 1.0/curv_min,energy_array[0],energy_array[-1]))



        if elliptical: 
            raise Exception(NotImplemented)
        else:
            for i in range(np):
                f.write("%e %e %e %e %e \n"%(x[i],y[i],cdf[i],angle[i],curv[i]))

        f.close()
        print("wiggler_cdf: File with wiggler cdf written to file: %s "%(outFile))

    return {"np":np,"step":step,"bener":bener,"inv_curv_max":(1.0/curv_max),"inv_curv_min":(1.0/curv_min),
                    "emin":energy_array[0],"emax":energy_array[-1],
                    "x":x,"y":y,"cdf":cdf,"angle":angle,"curv":curv}


def wiggler_trajectory(b_from=0, inData="", nPer=12, nTrajPoints=100, \
                       ener_gev=6.04, per=0.125, kValue=14.0, trajFile="",
                       shift_x_flag=0,shift_x_value=0.0,
                       shift_betax_flag=0,shift_betax_value=0.0):
    r"""
    Calculates the trajectory of the electrons under a magnetic field in Z direction.

     PROCEDURE:
        Based on btraj.f, a utility written in 10/91 by M. Sanchez del Rio and C. Vettier to input
        asymmetric wigglers in SHADOW.

        See formulae in ESRF red book    pag CIV-297
    :param b_from: A Flag for the type of inpyt magnetic field:
                   0: kValue (deflecting parameters) is given
                   1: A file with the magnetic field (y[m] B[T]) is given
                   2: A file with the magnetic field harmonics (n Bn[T]) is given
    :param inData: A string with the file with the file name containing the field information (for b_from:1,2),
                    or a [2,npoint] numpy array with the field information.
    :param nPer: Number of periods
    :param nTrajPoints: Number of trajectory points (per period)
    :param ener_gev: The electron energy in GeV
    :param per:  Wiggler period in meters  (for b_from:1,3)
    :param kValue: The K (deflecting parameter) value
    :param trajFile: The name of a file where the resut is written. (Default: "" no written file).
    :param shift_x_flag: a flag to indicate how to shift x (electron transversal coordinate):
                        0 = No shift
                        1 = Half excursion
                        2 = Minimum
                        3 = Maximum
                        4 = Value at zero
                        5 = User value
    :param shift_x_value: The user value for shift x (if shift_x_flag=5)
    :param shift_betax_flag: a flag to indicate how to shift the betax coordinate (electron transversal velocity):
                        0 = No shift
                        1 = Half excursion
                        2 = Minimum
                        3 = Maximum
                        4 = Value at zero
                        5 = User value
    :param shift_betax_value: The user value for shift betax (if shift_betax_flag=5)
    :return:           (traj,pars)
           traj: a variable to store the output matrix (8 colums with:
           x[m]  y[m]  z[m]  BetaX  BetaY  BetaZ  Curvature  B[T] )
           pars: a variable with text info


     MODIFICATION HISTORY:
           Written by:     M. Sanchez del Rio, srio@esrf.fr, 2002-07-17
        2002-07-17 srio@esrf.fr
        2012-10-08 srio@esrf.eu python version

     TODO: implement the elliptical wiggler (Bx!=0)
    """


    if b_from == 0:
        nharm = 1
        n=[1.0]
        cte = codata_ec/2.0/numpy.pi/ codata_me/codata_c # 93.372904
        bh=[kValue/(cte*per)]
        ystep = per/(nTrajPoints-1)

    # get period [m], number of periods, harmonics [T]
    if b_from >= 1:
        if type(inData) == type(""):   # file input
            if os.path.isfile(inData) == False:
                sys.exit('File nor found: '+inData)
            a = numpy.loadtxt(inData)
        else:  # numpy array input
            a = inData

    if b_from == 1:
        yy = a.T[0,:]
        bz = a.T[1,:]
        nmax = len(bz) + 1
        nTrajPoints = len(bz)
        ystep = yy[1]-yy[0]
        per  = yy[-1]-yy[0]
        print("%d points read "%(nTrajPoints))
        print("Period of the ID is %f m "%(per))


    if b_from == 2:
        n = a.T[0,:]
        bh = a.T[1,:]
        nharm = len(n)
        print("%d harmonics read"%(nharm))
        print("Period of the ID is %f m "%(per))
        print("Number of periods is %d "%(nPer))

    gamma = 1.0e3/codata_mee*ener_gev
    beta02 = 1.0 - numpy.power(gamma,-2)
    beta0 = numpy.sqrt(beta02)
    print(" gamma = %f"%(gamma))
    print(" beta = v/c = %20.10f"%(beta0))
    print(" beta^2 = %20.12f"%(beta02))
    start_len = per * (nPer-1) / 2.

    #;
    #;    calculates the integral of the magnetic field bz (proportional
    #;    to speeds)
    #;    

    #from scipy.integrate import simps
    #import int_tabulated
    #ystep = yy[1]-yy[0] 
    #mystep = yy  - numpy.roll(yy,1)
    #mystep[0] = 0.0

    if b_from == 1:
        betax = numpy.zeros(nTrajPoints)
        for i in range(1,nTrajPoints):
            # be careful [0:i] goes from 0 to i-1!!!!!!
            #betax[i] = scipy.integrate.simps(bz[0:i+1],yy[0:i+1])
            #betax[i] = numpy.sum(bz[0:i+1]*mystep[0:i+1])

            #betax[i] = numpy.trapz(bz[0:i+1],x=yy[0:i+1])
            #srio@esrf.eu changed sign of magnetic field to get the right curvature for ELECTRONS!
            betax[i] = numpy.trapz(-bz[0:i+1],x=yy[0:i+1])
        yInt = betax[-1]
    else:
        phase0 = numpy.zeros(nTrajPoints) - numpy.pi
        yy = numpy.linspace(-per/2.,per/2.,nTrajPoints)
        bz = numpy.zeros(nTrajPoints)
        betax = numpy.zeros(nTrajPoints)
        phase = 2.0*numpy.pi*(yy/per)
        for n in range(nharm): 
            bz = bz + bh[n] * numpy.cos(phase*(n+1))
            #srio@esrf.eu changed sign of magnetic field to get the right curvature for ELECTRONS!
            #tmp = (numpy.sin(phase*(n+1))-numpy.sin(phase0*(n+1)))*(bh[n]/(n+1))
            tmp = (numpy.sin(phase*(n+1))-numpy.sin(phase0*(n+1)))*(-bh[n]/(n+1))
            betax = betax + tmp
        betax = betax*(per/2.0/numpy.pi)

    print(" integral of b = %20.12e T.m "%(betax[nTrajPoints-1]))
    #;
    #;    rescale b to speeds v/c = 0.3/e[gev] * integral (bz [t] ds[m])
    #;
    betax = -codata_c*1e-9/ener_gev*betax
    betay = numpy.sqrt(beta0*beta0 - betax*betax)
    emc = codata_ec/gamma/codata_me/codata_c
    #;
    #;    calculates positions as the integral of velocities
    #;
    yx = numpy.zeros(nTrajPoints)
    if b_from == 1:
        for i in range(1,nTrajPoints):
            # be carefil [0:i] goes from 0 to i-1!!!!!!
            yx[i] = numpy.trapz(betax[0:i+1],x=yy[0:i+1])
            #yx[i] = scipy.integrate.simps(betax[0:i+1],yy[0:i+1])
            #yx[i] = numpy.sum(betax[0:i+1]*mystep[0:i+1])
    else:
        for n in range(nharm):
            phase = yy * (2.0*numpy.pi/per)
            yx = yx - (numpy.cos(phase*(n+1)) - numpy.cos(phase0*(n+1))) * \
                      (-bh[n]/numpy.power(n+1,2))

        yx = yx * (-3.e-1/ener_gev) * numpy.power(per/2.0/numpy.pi,2)
    #;
    #;    creates parameters text
    #;

    # full power per electronCurrent = 0.1A 
    #see slide 21 in http://www.cockcroft.ac.uk/education/PG_courses_2009-10/Spring_2010/CLarke%20Lecture%202.pdf 
    tmpP = (bz*bz).sum()*(yy[1]-yy[0])
    tmpP = tmpP*1265.5*ener_gev*ener_gev*0.1

    pars = ""
    pars += "\nPeriod of the ID is [m] :      "+repr(per)
    pars += "\nNumber of periods :            "+repr(nPer)
    pars += "\nNumber of points per periods : "+repr(nTrajPoints)
    pars += "\nBeta = v/c :                   %30.15e"%(beta0)
    #pars += "\nB = Cte*curv; Cte:             "+repr(-1.0*beta0/emc)
    pars += "\nElectron energy [GeV] :        "+repr(ener_gev)
    pars += "\nTotal emitted power [W] per 100 mA electron current:  "+repr(tmpP)
    pars += "\nMaximum of B = %f T "%(bz.max())
    #pars += "\nMaximum of K = %f T "%(93.36*bz.max()*per)
    pars += "\nMaximum of K = %f "%(codata_ec/codata_me/codata_c/2/numpy.pi*bz.max()*per)

    if b_from == 0:
        pars += "\nK value :                      "+repr(kValue)
        #pars += "\nMaximum of B [T] :             "+repr(bh[0])

    if b_from == 1:
        pars += "\nB[T] profile from file or array."
        pars += "\nIntegral of B = %f T.m "%(betax[1])

    if b_from == 2:
        pars += "\nB harmonics from file or array."
        pars += "\nNumber of harmonics = %d "%(nharm)

    #
    # calculate curvature
    #

    curv = emc*bz/beta0

    #
    # correct orbits
    #
    if shift_betax_flag == 0: # no shift
        shift_betax_local = 0.0
    elif shift_betax_flag == 1:  # shift half value
        shift_betax_local = -0.5 * (betax.max() + betax.min())
    elif shift_betax_flag == 2:  # shift min
        shift_betax_local = -betax.min()
    elif shift_betax_flag == 3:  # shift max
        shift_betax_local = -betax.max()
    elif shift_betax_flag == 4: # valus at zero
        shift_betax_local = -numpy.interp(0.0,yy,betax)
    elif shift_betax_flag == 5:  # user value
        shift_betax_local = shift_betax_value
    else:
        raise "Unknown value for shift_betax_flag"

    if shift_betax_local != 0.0:
        betax += shift_betax_local
        betay = numpy.sqrt(beta02 - betax**2)
        yx += shift_betax_local * yy


    if shift_x_flag == 0: # no shift
        shift_x_local = 0.0
    elif shift_x_flag == 1: # half excursion
        shift_x_local = -0.5 * ( yx.max() + yx.min() )
    elif shift_x_flag == 2: # min
        shift_x_local = -yx.min()
    elif shift_x_flag == 3: # max
        shift_x_local = -yx.max()
    elif shift_x_flag == 4: # valus at zero
        shift_x_local = -numpy.interp(0.0,yy,yx)
    elif shift_x_flag == 5: # user value
        shift_x_local =  shift_x_value
    else:
        raise "Unknown value for shift_x_flag"

    if shift_x_local != 0.0:
        yx +=  shift_x_local

    #;
    #;    creates trajectory and file 
    #;

    nPointsTot = nTrajPoints+(nPer-1)*(nTrajPoints-1)
    traj = numpy.zeros((8,nPointsTot))


    ii = -1
    for j in range(nPer):
        nn = 0
        if (j > 0): 
           nn = 1     # to avoid overlap
        #look here!!
        for i in range(nn,nTrajPoints):
            ii = ii + 1
            traj[0,ii] = yx[i]
            traj[1,ii] = yy[i]+j * per - start_len
            traj[2,ii] = 0.0
            traj[3,ii] = betax[i]
            traj[4,ii] = betay[i]
            traj[5,ii] = 0.0
            traj[6,ii] = curv[i]
            traj[7,ii] = bz[i]

    if trajFile != "":
        f = open(trajFile,"w")
        f.write("#F "+trajFile+"\n")

        f.write("\n#S 1 Electron trajectory and velocity")
        f.write(pars.replace('\n','\n#UD '))
        f.write("\n#N 8\n")
        f.write("#L x[m]  y[m]  z[m]  BetaX  BetaY  BetaZ  Curvature  B[T]\n")
        for j in range(nPointsTot): 
            tmp = traj[:,j]
            f.write(("%19.12e  "*8+"\n")%tuple(tmp[i] for i in range(len(tmp))) )
        f.close()
        print("File with trajectory written to file: "+trajFile)

    return (traj,pars)

def wiggler_nphoton(r_m,electronEnergy=1.0,photonEnergy=1000.0,polarization=0):
    r""" 
     NAME:
           wiggler_nphoton
    
     PURPOSE:
           Calculates the synchrotron radiation spectrum versus bending radius
        Assumptions:
            Electron current = 1 mA
            Horizontal divergence = 1mrad
            Energy bandwidth = 1 eV
    
     CATEGORY:
           Synchrotron radiation
    
     CALLING SEQUENCE:
           result = wiggler_nPhoton(r_m,e_geV,energy_ev)
    
     INPUTS:
        r_m: the array with the bending radii in m

     KEYWORD PARAMETERS:
        electronEnergy: electrons energy in GeV
        photonEnergy:  photon energy in eV
        polarization: 0 Total 
                      1 Parallel       (l2=1, l3=0, in Sokolov&Ternov notation)
                      2 Perpendicular  (l2=0, l3=1)
    
    
     OUTPUTS:
           returns the array with the flux [photons/sec/1eV/mrad/mA)
    
     PROCEDURE:
    
        It uses sync_ene
    
        References: 
            G K Green, "Spectra and optics of synchrotron radiation" 
                BNL 50522 report (1976)
            A A Sokolov and I M Ternov, Synchrotron Radiation, 
                Akademik-Verlag, Berlin, 1968
    
     EXAMPLE:
        The following program was used for testing nphoton
        
        
        IDL> r_m = makearray1(100,1,500)
        IDL> xplot,r_m,nphoton(r_m,6.0,10000.0)    
    
    
     MODIFICATION HISTORY:
        Written by:     M. Sanchez del Rio, srio@esrf.fr, 2002-06-24
        2012-10-08 srio@esrf.eu python version
    
    """

    #cte = (3.0d0/4/!dpi)*physical_constants('h')*physical_constants('c')* $
    #  (1d3/physical_constants('mee'))^3/physical_constants('ec')
    #cte = 2218.2873 

    cte = (3.0e0/4/numpy.pi) * codata_h * codata_c * numpy.power(1e3 / codata_mee,3) / codata_ec


    ec_ev = cte * numpy.power(electronEnergy,3) / r_m


    nn = sync_ene(0,photonEnergy,ec_ev=ec_ev,e_gev=electronEnergy,i_a=1e-3, hdiv_mrad=1.0, polarization=polarization)

    return nn 

def wiggler_harmonics(Bs,Nh=41,fileOutH=""):

    r"""
     NAME:
           wiggler_harmonics
    
     PURPOSE:
           Calculates the harmonic decomposition of the magnetic field map B[s]
    
     CATEGORY:
           Synchrotron radiation
    
     CALLING SEQUENCE:
            hh = wiggler_harmonics(Bs,Nh=41,fileOutH="tmp.h")
    
     INPUTS:
           Bs:  An array [npoints,2] with B [T] vs s[m]
    
     KEYWORD PARAMETERS:
           Nh:     The number of harmonics for the decomposition
           fileOutH:    A file name where to write the  resulting decomposition.
                        If set to "" (default) no file is written. 
     OUTPUT: 
           hh an array [2,nH] with the number of harmonic and the coefficient. 
    
     PROCEDURE:
        Based on Fourier filtering. See 
        http://stackoverflow.com/questions/5843085/fourier-series-in-numpy-question-about-previous-answer

     EXAMPLE

     MODIFICATION HISTORY:
           Written by:     M. Sanchez del Rio, srio@esrf.fr, 2014-10-08
        2012-10-08 srio@esrf.eu python version
    
    """


    #see http://stackoverflow.com/questions/5843085/fourier-series-in-numpy-question-about-previous-answer
    cn = lambda n: ((  2.0*bb*numpy.exp(-1j*2*n*numpy.pi*yy/(yy[-1]-yy[0]))).sum()/bb.size).real
    
    #
    # define magnetic field [T] vs s [m]
    #
    yy = Bs[:,0]
    bb = Bs[:,1]

    #get coeffs
    hh = numpy.zeros((Nh-1,2))
    for i in range(1,Nh):
        hh[i-1,0] = i
        hh[i-1,1] = cn(i)
    
    if fileOutH != "":
        f = open(fileOutH,'w')
        for i in range(hh.shape[0]):
            f.write('%d %f \n'%(hh[i,0],hh[i,1]))
        f.close()
        print("file "+fileOutH+" written to disk (harmonics).")
    
    return hh


#
# tests
#
def test_xraybooklet_fig2_1(pltOk=False):
    print("# ")
    print("# example 1, Fig 2-1 in http://xdb.lbl.gov/Section2/Sec_2-1.html ")
    print("# ")
    y = numpy.logspace(-3,1,100)  #  from 0.001 to 10, 100 points
    g1 = sync_g1(y,polarization=0)
    h2 = sync_hi(y,i=2,polarization=0)
    # TODO: check transpose
    h2 = h2.T

    toptitle = "Synchrotron Universal Functions $G_1$ and $H_2$"
    xtitle = "$y=E/E_c$"
    ytitle = "$G_1(y)$ and $H_2(y)$"

    #pltOk = 0
    if pltOk:
        plt.figure(1)
        plt.loglog(y,g1,'b',label="$G_1(y)$")
        plt.loglog(y,h2,'r',label="$H_2(y)$")
        plt.title(toptitle)
        plt.xlabel(xtitle)
        plt.ylabel(ytitle)
        plt.ylim((1e-3,10))
        ax = plt.subplot(111)
        ax.legend(bbox_to_anchor=(1.1, 1.05))
    else:
        print("\n\n\n\n\n#########  %s ######### "%(toptitle))
        print("\n  %s  %s "%(xtitle,ytitle))
        for i in range(len(y)):
            print(" %f  %e %e "%(y[i],g1[i],h2[i]))

def test_xraybooklet_fig2_2(pltOk=False):
    print("# ")
    print("# example  2, Fig 2-2 in http://xdb.lbl.gov/Section2/Sec_2-1.html ")
    print("# ")

    y = numpy.linspace(0,8,100)  #  from 0.001 to 10, 100 points
    f3   = sync_f(y,3.0,polarization=1)
    f3pi = sync_f(y,3.0,polarization=2)
    f1   = sync_f(y,1.0,polarization=1)
    f1pi = sync_f(y,1.0,polarization=2)
    fp1   = sync_f(y,0.1,polarization=1)
    fp1pi = sync_f(y,0.1,polarization=2)
    fp01   = sync_f(y,0.01,polarization=1)
    fp01pi = sync_f(y,0.01,polarization=2)

    toptitle = "Synchrotron Angular Emission"
    xtitle = "$\gamma \Psi$"
    ytitle = "$N(\%)$"

    f3.shape = -1
    f3pi.shape = -1
    f3max = f3.max()*1e-2 # to get %

    f1.shape = -1
    f1pi.shape = -1
    f1max = f1.max()*1e-2 # to get %

    fp01.shape = -1
    fp01pi.shape = -1
    fp01max = fp01.max()*1e-2 # to get %

    fp1.shape = -1
    fp1pi.shape = -1
    fp1max = fp1.max()*1e-2 # to get %

    if pltOk:
        plt.figure(2)
        plt.plot(y,f3/f3max,'b',label="E/Ec=3 $\sigma$-pol")
        plt.plot(y,f3pi/f3max,'b--',label="E/Ec=3 $\pi$-pol")

        plt.plot(y,f1/f1max,'g',label="E/Ec=1 $\sigma$-pol")
        plt.plot(y,f1pi/f1max,'g--',label="E/Ec=1 $\pi$-pol")

        plt.plot(y,fp1/fp1max,'k',label="E/Ec=0.1 $\sigma$-pol")
        plt.plot(y,fp1pi/fp1max,'k--',label="E/Ec=0.1 $\pi$-pol")

        plt.plot(y,fp01/fp01max,'r',label="E/Ec=0.01 $\sigma$-pol")
        plt.plot(y,fp01pi/fp01max,'r--',label="E/Ec=0.01 $\pi$-pol")

        plt.title(toptitle)
        plt.xlabel(xtitle)
        plt.ylabel(ytitle)
        ax = plt.subplot(111)
        ax.legend(bbox_to_anchor=(1.1, 1.05))
    else:
        print("\n\n\n\n\n#########  %s ######### "%(toptitle))
        print("\n  %s  %s "%(xtitle,ytitle))
        for i in range(len(y)):
            print(" %f  %e %e "%(y[i],f3[i],f3pi[i]))

def test_esrf_bm_spectrum(pltOk=False):
    print("#")
    print("# example 3, ESRF BM spectrum")
    print("#")
    # input for ESRF
    e_gev = 6.04    # electron energy in GeV
    r_m = 25.0      # magnetic radius in m
    i_a = 0.2       # electron current in A
    # calculate critical energy in eV
    gamma = e_gev*1e3/codata_mee
    ec_m = 4.0*numpy.pi*r_m/3.0/numpy.power(gamma,3) # wavelength in m
    ec_ev = m2ev/ec_m

    energy_ev = numpy.linspace(100.0,100000.0,99) # photon energy grid
    f_psi = 0    # flag: full angular integration
    flux = sync_ene(f_psi,energy_ev,ec_ev=ec_ev,polarization=0,  \
           e_gev=e_gev,i_a=i_a,hdiv_mrad=1.0, \
           psi_min=0.0, psi_max=0.0, psi_npoints=1)

    toptitle = "ESRF Bending Magnet emission"
    xtitle = "Photon energy [eV]"
    ytitle = "Photon flux [Ph/s/mrad/0.1%bw]"

    if pltOk:
        plt.figure(3)
        plt.loglog(energy_ev,flux)
        plt.title(toptitle)
        plt.xlabel(xtitle)
        plt.ylabel(ytitle)
    else:
        print("\n\n\n\n\n#########  %s ######### "%(toptitle))
        print("\n  %s  %s "%(xtitle,ytitle))
        for i in range(len(flux)):
            print(" %f  %12.3e"%(energy_ev[i],flux[i]))

def test_esrf_bm_angle_power(pltOk=False):
    print("#")
    print("# example 4: ESRF BM angular emission of power")
    print("#")
    # input for ESRF
    e_gev = 6.04    # electron energy in GeV
    r_m = 25.0      # magnetic radius in m
    i_a = 0.2       # electron current in A

    angle_mrad = numpy.linspace(-1.0,1.0,201) # angle grid
    flag = 0 # full energy integration
    flux = sync_ang(flag,angle_mrad,polarization=0, \
           e_gev=e_gev,i_a=i_a,hdiv_mrad=1.0,r_m=r_m)

    #TODO: integrate curve and compare with total power
    toptitle = "ESRF Bending Magnet angular emission (all energies)"
    xtitle   = "Psi[mrad]"
    ytitle   = "Photon Power [Watts/mrad(Psi)]"
    if pltOk:
        plt.figure(4)
        plt.plot(angle_mrad,flux)
        plt.title(toptitle)
        plt.xlabel(xtitle)
        plt.ylabel(ytitle)
    else:
        print("\n\n\n\n\n#########  %s ######### "%(toptitle))
        print("\n  %s  %s "%(xtitle,ytitle))
        for i in range(len(flux)):
            print("  %f  %f"%(angle_mrad[i],flux[i]))

def test_esrf_bm_angle_flux(pltOk=False):
    print("#")
    print("# example 5: ESRF BM angular emission of flux")
    print("#")

    # input for ESRF
    e_gev = 6.04    # electron energy in GeV
    r_m = 25.0      # magnetic radius in m
    i_a = 0.2       # electron current in A
    # calculate critical energy in eV
    gamma = e_gev*1e3/codata_mee
    ec_m = 4.0*numpy.pi*r_m/3.0/numpy.power(gamma,3) # wavelength in m
    ec_ev = m2ev/ec_m

    angle_mrad = numpy.linspace(-1.0,1.0,201) # angle grid
    flag = 1 # at at given energy
    fluxEc = sync_ang(flag,angle_mrad,polarization=0, \
           e_gev=e_gev,i_a=i_a,hdiv_mrad=1.0,energy=ec_ev, ec_ev=ec_ev)
    flux10Ec = sync_ang(flag,angle_mrad,polarization=0, \
           e_gev=e_gev,i_a=i_a,hdiv_mrad=1.0,energy=10*ec_ev, ec_ev=ec_ev)
    fluxp1Ec = sync_ang(flag,angle_mrad,polarization=0, \
           e_gev=e_gev,i_a=i_a,hdiv_mrad=1.0,energy=0.1*ec_ev, ec_ev=ec_ev)

    toptitle = "ESRF Bending Magnet angular emission (all energies)"
    xtitle   = "Psi[mrad]"
    ytitle   = "Flux [phot/s/0.1%bw/mrad(Psi)]"
    if pltOk:
        plt.figure(5)
        plt.plot(angle_mrad,fluxp1Ec,'g',label="E=0.1Ec=%.3f keV"%(ec_ev*.1*1e-3))
        plt.plot(angle_mrad,fluxEc,'r',label="E=Ec=%.3f keV"%(ec_ev*1e-3))
        plt.plot(angle_mrad,flux10Ec,'b',label="E=10 Ec=%.3f keV"%(ec_ev*10*1e-3))
        plt.title(toptitle)
        plt.xlabel(xtitle)
        plt.ylabel(ytitle)
        ax = plt.subplot(111)
        ax.legend(bbox_to_anchor=(1.1, 1.05))
    else:
        print("\n\n\n\n\n#########  %s ######### "%(toptitle))
        print("\n  %s  %s "%(xtitle,ytitle))
        for i in range(len(fluxEc)):
            print("  %f  %f"%(angle_mrad[i],fluxEc[i]))

def test_clarke_43(pltOk=False):
    print("#")
    print("# Example 6 Slide 35 of")
    print("# http:https://www.cockcroft.ac.uk/wp-content/uploads/2014/12/Lecture-1.pdf")
    print("#")


    def calcFWHM(h,binSize):
        t = numpy.where(h>=max(h)*0.5)
        return binSize*(t[0][-1]-t[0][0]+1), t[0][-1], t[0][0]
    #
    e_gev = 3.0    # electron energy in GeV
    b_t = 1.4      # magnetic radius in m
    i_a = 0.3      # electron current in A


    gamma = e_gev*1e3/codata_mee

    #calculates Magnetic radius
    #cte = codata_me*codata_c/codata_ec*(1/(codata_mee*1e-3)) # 0.3
    #r_m = cte*e_gev/b_t
    #more exactly
    r_m = codata_me*codata_c/codata_ec/b_t*numpy.sqrt( gamma*gamma - 1)

    # calculate critical energy in eV
    ec_m = 4.0*numpy.pi*r_m/3.0/numpy.power(gamma,3) # wavelength in m
    ec_ev = m2ev/ec_m


    print("Gamma: %f \n"%(gamma))
    print("Critical wavelength [A]: %f \n"%(1e10*ec_m))
    print("Critical photon energy [eV]: %f \n"%(ec_ev))

    e = numpy.array([0.1*ec_ev,ec_ev,10*ec_ev])
    a = numpy.linspace(-0.6,0.6,150)

    fm = sync_ene(4,e,ec_ev=ec_ev,e_gev=e_gev,i_a=i_a,\
        hdiv_mrad=1,psi_min=-0.6,psi_max=0.6,psi_npoints=150)

    fmPar = sync_ene(4,e,ec_ev=ec_ev,e_gev=e_gev,i_a=i_a,\
        hdiv_mrad=1,psi_min=-0.6,psi_max=0.6,psi_npoints=150,polarization=1)
    fmPer = sync_ene(4,e,ec_ev=ec_ev,e_gev=e_gev,i_a=i_a,\
        hdiv_mrad=1,psi_min=-0.6,psi_max=0.6,psi_npoints=150,polarization=2)
    toptitle='Flux vs vertical angle '
    xtitle  ='angle [mrad]'
    ytitle = "Photon flux [Ph/s/mrad/0.1%bw]"

    print("for E = 0.1 Ec FWHM=%f mrad "%( calcFWHM(fm[:,0],a[1]-a[0])[0]))
    print("for E =     Ec FWHM=%f mrad "%( calcFWHM(fm[:,1],a[1]-a[0])[0]))
    print("for E = 10  Ec FWHM=%f mrad "%( calcFWHM(fm[:,2],a[1]-a[0])[0]))

    print("Using approximated formula: ")
    print("for E = 0.1 Ec FWHM=%f mrad "%( 0.682 / e_gev * numpy.power(10.0,0.425) ))
    print("for E =     Ec FWHM=%f mrad "%( 0.682 / e_gev * numpy.power(1.0,0.425) ))
    print("for E = 10  Ec FWHM=%f mrad "%( 0.682 / e_gev * numpy.power(0.1,0.425) ))

    if pltOk:
        plt.figure(61)
        plt.plot(a,fm[:,0],'b',label="0.1*$\omega_c$")
        plt.plot(a,fmPar[:,0],"b--")
        plt.plot(a,fmPer[:,0],"b-.")
        plt.title(toptitle+"E=0.1Ec")
        plt.xlabel(xtitle)
        plt.ylabel(ytitle)

        plt.figure(61)
        plt.plot(a,fm[:,1],'red',label="$\omega_c$")
        plt.plot(a,fmPar[:,1],"r--")
        plt.plot(a,fmPer[:,1],"r-.")
        plt.title(toptitle+"E=Ec")
        plt.xlabel(xtitle)
        plt.ylabel(ytitle)

        plt.figure(62)
        plt.plot(a,fm[:,2],'green',label="10*$\omega_c$")
        plt.plot(a,fmPar[:,2],"g--")
        plt.plot(a,fmPer[:,2],"g-.")
        plt.title(toptitle+"E=10Ec")
        plt.xlabel(xtitle)
        plt.ylabel(ytitle)

    else:
        print("\n\n\n\n\n#########  %s ######### "%(toptitle))
        print("\n  %s  %s %s "%(ytitle,xtitle,"Flux"))
        for j in range(len(e)):
            for i in range(len(a)):
                print("  %f  %f  %e   "%(e[j],a[i],fm[i,j]))

def test_esrf_bm_2d(pltOk=False):
    print("#")
    print("# Example 7, ESRF flux vs energy and angle")
    print("#")
    # input for ESRF
    e_gev = 6.04    # electron energy in GeV
    r_m = 25.0      # magnetic radius in m
    i_a = 0.2       # electron current in A
    # calculate critical energy in eV
    gamma = e_gev*1e3/codata_mee
    ec_m = 4.0*numpy.pi*r_m/3.0/numpy.power(gamma,3) # wavelength in m
    ec_ev = m2ev/ec_m

    a = numpy.linspace(-0.2,0.2,50)
    e = numpy.linspace(20000,80000,80)

    fm = sync_ene(4,e,ec_ev=ec_ev,e_gev=e_gev,i_a=i_a,\
        hdiv_mrad=1,psi_min=a.min(),psi_max=a.max(),psi_npoints=a.size)
    toptitle='Flux vs vertical angle and photon energy'
    xtitle  ='angle [mrad]'
    ytitle  ='energy [eV]'
    ztitle = "Photon flux [Ph/s/mrad/0.1%bw]"

    if pltOk:
        fig = plt.figure(7)
        ax = fig.add_subplot(111, projection='3d')
        fa, fe = numpy.meshgrid(a, e)
        surf = ax.plot_surface(fa, fe, fm.T, \
            rstride=1, cstride=1, \
            linewidth=0, antialiased=False)

        plt.title(toptitle)
        ax.set_xlabel(xtitle)
        ax.set_ylabel(ytitle)
        ax.set_zlabel(ztitle)

    else:
        print("\n\n\n\n\n#########  %s ######### "%(toptitle))
        print("\n  %s  %s %s "%(xtitle,ytitle,ztitle))
        for i in range(len(a)):
            for j in range(len(e)):
                print("  %f  %f  %e   "%(a[i],e[j],fm[i,j]))

def test_wiggler_flux_vs_r(pltOk=False):
    print("#")
    print("# Example 8 (Wiggler flux vs bending radius at a given photon energy)")
    print("#")
    r_m = numpy.linspace(1.0,500.0,100)
    flux = wiggler_nphoton(r_m,electronEnergy=6.04,photonEnergy=10000.0)
    toptitle = "Wiggler flux vs bending radius at photon energy E=10 keV"
    xtitle   = "Magneric radius R [m]"
    ytitle   = "Flux [phot/s/eV/1mradH/mA]"

    if pltOk:
        plt.figure(8)
        plt.plot(r_m,flux)
        plt.title(toptitle)
        plt.xlabel(xtitle)
        plt.ylabel(ytitle)
    else:
        print("\n\n\n\n\n#########  %s ######### "%(toptitle))
        print("\n  %s  %s "%(xtitle,ytitle))
        for i in range(len(r_m)):
            print("  %f  %e"%(r_m[i],flux[i]))


def test_wiggler_external_b(pltOk=False):
    print("#")
    print("# Example 9 (Wiggler trajectory and flux for a 3pole wiggler ")
    print("#")

    # this is the B(s) map (T, m)
    b_t = numpy.array([[ -1.00000000e-01,   1.08200000e-03],
   [ -9.80000000e-02,   8.23000000e-04],
   [ -9.60000000e-02,   4.45000000e-04],
   [ -9.40000000e-02,   8.60000000e-05],
   [ -9.20000000e-02,  -4.93000000e-04],
   [ -9.00000000e-02,  -1.20800000e-03],
   [ -8.80000000e-02,  -2.16100000e-03],
   [ -8.60000000e-02,  -3.44500000e-03],
   [ -8.40000000e-02,  -5.10500000e-03],
   [ -8.20000000e-02,  -7.34500000e-03],
   [ -8.00000000e-02,  -1.03050000e-02],
   [ -7.80000000e-02,  -1.42800000e-02],
   [ -7.60000000e-02,  -1.96770000e-02],
   [ -7.40000000e-02,  -2.70560000e-02],
   [ -7.20000000e-02,  -3.73750000e-02],
   [ -7.00000000e-02,  -5.20600000e-02],
   [ -6.80000000e-02,  -7.35170000e-02],
   [ -6.60000000e-02,  -1.05680000e-01],
   [ -6.40000000e-02,  -1.54678000e-01],
   [ -6.20000000e-02,  -2.28784000e-01],
   [ -6.00000000e-02,  -3.34838000e-01],
   [ -5.80000000e-02,  -4.70272000e-01],
   [ -5.60000000e-02,  -6.16678000e-01],
   [ -5.40000000e-02,  -7.46308000e-01],
   [ -5.20000000e-02,  -8.39919000e-01],
   [ -5.00000000e-02,  -8.96470000e-01],
   [ -4.80000000e-02,  -9.26065000e-01],
   [ -4.60000000e-02,  -9.38915000e-01],
   [ -4.40000000e-02,  -9.40738000e-01],
   [ -4.20000000e-02,  -9.32236000e-01],
   [ -4.00000000e-02,  -9.08918000e-01],
   [ -3.80000000e-02,  -8.60733000e-01],
   [ -3.60000000e-02,  -7.73534000e-01],
   [ -3.40000000e-02,  -6.36577000e-01],
   [ -3.20000000e-02,  -4.52611000e-01],
   [ -3.00000000e-02,  -2.37233000e-01],
   [ -2.80000000e-02,  -7.09700000e-03],
   [ -2.60000000e-02,   2.26731000e-01],
   [ -2.40000000e-02,   4.54558000e-01],
   [ -2.20000000e-02,   6.61571000e-01],
   [ -2.00000000e-02,   8.29058000e-01],
   [ -1.80000000e-02,   9.45984000e-01],
   [ -1.60000000e-02,   1.01683300e+00],
   [ -1.40000000e-02,   1.05536200e+00],
   [ -1.20000000e-02,   1.07490000e+00],
   [ -1.00000000e-02,   1.08444200e+00],
   [ -8.00000000e-03,   1.08898000e+00],
   [ -6.00000000e-03,   1.09111200e+00],
   [ -4.00000000e-03,   1.09208300e+00],
   [ -2.00000000e-03,   1.09249400e+00],
   [  0.00000000e+00,   1.09262000e+00],
   [  2.00000000e-03,   1.09249400e+00],
   [  4.00000000e-03,   1.09208300e+00],
   [  6.00000000e-03,   1.09111200e+00],
   [  8.00000000e-03,   1.08898000e+00],
   [  1.00000000e-02,   1.08444200e+00],
   [  1.20000000e-02,   1.07490000e+00],
   [  1.40000000e-02,   1.05536200e+00],
   [  1.60000000e-02,   1.01683300e+00],
   [  1.80000000e-02,   9.45984000e-01],
   [  2.00000000e-02,   8.29058000e-01],
   [  2.20000000e-02,   6.61571000e-01],
   [  2.40000000e-02,   4.54558000e-01],
   [  2.60000000e-02,   2.26731000e-01],
   [  2.80000000e-02,  -7.09700000e-03],
   [  3.00000000e-02,  -2.37233000e-01],
   [  3.20000000e-02,  -4.52611000e-01],
   [  3.40000000e-02,  -6.36577000e-01],
   [  3.60000000e-02,  -7.73534000e-01],
   [  3.80000000e-02,  -8.60733000e-01],
   [  4.00000000e-02,  -9.08918000e-01],
   [  4.20000000e-02,  -9.32236000e-01],
   [  4.40000000e-02,  -9.40738000e-01],
   [  4.60000000e-02,  -9.38915000e-01],
   [  4.80000000e-02,  -9.26065000e-01],
   [  5.00000000e-02,  -8.96470000e-01],
   [  5.20000000e-02,  -8.39919000e-01],
   [  5.40000000e-02,  -7.46308000e-01],
   [  5.60000000e-02,  -6.16678000e-01],
   [  5.80000000e-02,  -4.70272000e-01],
   [  6.00000000e-02,  -3.34838000e-01],
   [  6.20000000e-02,  -2.28784000e-01],
   [  6.40000000e-02,  -1.54678000e-01],
   [  6.60000000e-02,  -1.05680000e-01],
   [  6.80000000e-02,  -7.35170000e-02],
   [  7.00000000e-02,  -5.20600000e-02],
   [  7.20000000e-02,  -3.73750000e-02],
   [  7.40000000e-02,  -2.70560000e-02],
   [  7.60000000e-02,  -1.96770000e-02],
   [  7.80000000e-02,  -1.42800000e-02],
   [  8.00000000e-02,  -1.03050000e-02],
   [  8.20000000e-02,  -7.34500000e-03],
   [  8.40000000e-02,  -5.10500000e-03],
   [  8.60000000e-02,  -3.44500000e-03],
   [  8.80000000e-02,  -2.16100000e-03],
   [  9.00000000e-02,  -1.20800000e-03],
   [  9.20000000e-02,  -4.93000000e-04],
   [  9.40000000e-02,   8.60000000e-05],
   [  9.60000000e-02,   4.45000000e-04],
   [  9.80000000e-02,   8.23000000e-04],
   [  1.00000000e-01,   1.08200000e-03]])



    # normal (sinusoidal) wiggler
    t0,p = wiggler_trajectory(b_from=0, nPer=1, nTrajPoints=100,  \
                             ener_gev=6.04, per=0.2, kValue=7.75, \
                             trajFile="tmpS.traj")

    # magnetic field from B(s) map
    t1,p = wiggler_trajectory(b_from=1, nPer=1, nTrajPoints=100,  \
                   ener_gev=6.04, inData=b_t,trajFile="tmpB.traj")
    # magnetic field from harmonics
    hh = wiggler_harmonics(b_t,Nh=41,fileOutH="tmp.h")
    t2,p = wiggler_trajectory(b_from=2, nPer=1, nTrajPoints=100,  \
                   ener_gev=6.04, per=0.2, inData=hh,trajFile="tmpH.traj")


    toptitle = "3-pole ESRF wiggler trajectory"
    xtitle   = "y[m]"
    ytitle   = "x[um]"

    if pltOk:
        plt.figure(91)
        plt.plot(t0[1,:],1e6*t0[0,:],'b',label="Sinusoidal")
        plt.plot(t1[1,:],1e6*t1[0,:],'r',label="Numeric tmp.b")
        plt.plot(t2[1,:],1e6*t2[0,:],'g',label="Harmonics tmp.h")
        plt.title(toptitle)
        plt.xlabel(xtitle)
        plt.ylabel(ytitle)
        ax = plt.subplot(111)
        ax.legend(bbox_to_anchor=(1.1, 1.05))
    else:
        print("\n\n\n\n\n#########  %s ######### "%(toptitle))
        print(" x[m]  y[m]  z[m]  BetaX  BetaY  BetaZ  Curvature  B[T] ")
        for i in range(t2.shape[1]):
            print(("%.2e "*8+"\n")%( tuple(t2[0,i] for i in range(t2.shape[0])   )))
    if True:
        #
        # now spectra
        #
        e, f0, tmp = wiggler_spectrum(t0,enerMin=100.0,enerMax=100000.0,nPoints=100, \
                     electronCurrent=0.2, outFile="tmp.dat", elliptical=False)
        e, f1, tmp = wiggler_spectrum(t1,enerMin=100.0,enerMax=100000.0,nPoints=100, \
                     electronCurrent=0.2, outFile="tmp.dat", elliptical=False)
        e, f2, tmp = wiggler_spectrum(t2,enerMin=100.0,enerMax=100000.0,nPoints=100, \
                     electronCurrent=0.2, outFile="tmp.dat", elliptical=False)

        toptitle = "3-pole ESRF wiggler spectrum"
        xtitle   = "E [eV]"
        ytitle   = "Flux[phot/s/0.1%bw]"

        if pltOk:
            plt.figure(92)
            plt.plot(e,f0,'b',label="Sinusoidal")
            plt.plot(e,f1,'r',label="Numeric 3PW_1.1T.b")
            plt.plot(e,f2,'g',label="Harmonics 3PW_1.1T.h")
            plt.title(toptitle)
            plt.xlabel(xtitle)
            plt.ylabel(ytitle)
            ax = plt.subplot(111)
            ax.legend(bbox_to_anchor=(1.1, 1.05))
        else:
            print("\n\n\n\n\n#########  %s ######### "%(toptitle))
            print(" energy[eV] flux_sinusoidal  flux_fromB  flux_fromHarmonics  ")
            for i in range(t2.shape[1]):
                print(("%.2e "*4+"\n")%( e[i],f0[i], f1[i], f2[i] ))


#
#------------------------- MAIN ------------------------------------------------
#
if __name__ == '__main__':

    pltOk = False
    try:
        import matplotlib.pylab as plt
        from mpl_toolkits.mplot3d import Axes3D  # need for example 6
        pltOk = True
    except ImportError:
        print("failed to import matplotlib. No on-line plots.")

    test_xraybooklet_fig2_1(pltOk=pltOk)
    test_xraybooklet_fig2_2(pltOk=pltOk)
    test_esrf_bm_spectrum(pltOk=pltOk)
    test_esrf_bm_angle_power(pltOk=pltOk)
    test_esrf_bm_angle_flux(pltOk=pltOk)
    test_clarke_43(pltOk=pltOk)
    test_esrf_bm_2d(pltOk=pltOk)
    test_wiggler_flux_vs_r(pltOk=pltOk)
    test_wiggler_external_b(pltOk=pltOk)

    if pltOk: plt.show()

