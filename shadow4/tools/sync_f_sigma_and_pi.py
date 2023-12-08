# This is similar to sync_f in srxraylib but faster
import numpy
import scipy
# todo: move to sr-xraylib
def sync_f_sigma_and_pi(rAngle, rEnergy):
    r""" angular dependency of synchrotron radiation emission

      NAME:
            sync_f_sigma_and_pi

      PURPOSE:
            Calculates the function used for calculating the angular
         dependence of synchrotron radiation.

      CATEGORY:
            Mathematics.

      CALLING SEQUENCE:
            Result = sync_f_sigma_and_pi(rAngle,rEnergy)

      INPUTS:
            rAngle:  (array) the reduced angle, i.e., angle[rads]*Gamma. It can be a
             scalar or a vector.
            rEnergy:  (scalar) a value for the reduced photon energy, i.e.,
             energy/critical_energy.

      KEYWORD PARAMETERS:


      OUTPUTS:
            returns the value  of the sync_f for sigma and pi polarizations
             The result is an array of the same dimension as rAngle.

      PROCEDURE:
            The number of emitted photons versus vertical angle Psi is
         proportional to sync_f, which value is given by the formulas
         in the references.


         References:
             G K Green, "Spectra and optics of synchrotron radiation"
                 BNL 50522 report (1976)
             A A Sokolov and I M Ternov, Synchrotron Radiation,
                 Akademik-Verlag, Berlin, 1968

      OUTPUTS:
            returns the value  of the sync_f function

      PROCEDURE:
            Uses BeselK() function

      MODIFICATION HISTORY:
            Written by:     M. Sanchez del Rio, srio@esrf.fr, 2002-05-23
         2002-07-12 srio@esrf.fr adds circular polarization term for
             wavelength integrated spectrum (S&T formula 5.25)
         2012-02-08 srio@esrf.eu: python version
         2019-10-31 srio@lbl.gov  speed-up changes for shadow4

    """

    #
    # ; Eq 11 in Pag 6 in Green 1975
    #
    ji = numpy.sqrt((1.0 + rAngle**2)**3) * rEnergy / 2.0
    efe_sigma = scipy.special.kv(2.0 / 3.0, ji) * (1.0 + rAngle**2)
    efe_pi = rAngle * scipy.special.kv(1.0 / 3.0, ji) / numpy.sqrt(1.0 + rAngle ** 2) * (1.0 + rAngle ** 2)
    return efe_sigma**2, efe_pi**2