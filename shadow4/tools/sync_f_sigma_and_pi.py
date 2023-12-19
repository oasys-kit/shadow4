# This is similar to sync_f in srxraylib but faster
import numpy
from scipy.special import kv, gamma
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
    efe_sigma = kv(2.0 / 3.0, ji) * (1.0 + rAngle**2)
    efe_pi = rAngle * kv(1.0 / 3.0, ji) / numpy.sqrt(1.0 + rAngle ** 2) * (1.0 + rAngle ** 2)
    return efe_sigma**2, efe_pi**2

def sync_f_sigma_and_pi_approx(rAngle, rEnergy):
    r""" angular dependency of synchrotron radiation emission. Using approximated Modified Bessel functions.

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
    efe_sigma = kv_approx(2.0 / 3.0, ji) * (1.0 + rAngle**2)
    efe_pi = rAngle * kv_approx(1.0 / 3.0, ji) / numpy.sqrt(1.0 + rAngle ** 2) * (1.0 + rAngle ** 2)
    return efe_sigma**2, efe_pi**2

def kv_approx(nu, x):
    # Approximated expressions for the modified Bessel functions K1/3, K2/3 and K5/3
    # Coefficients have been fitted using the expression in:
    # https://goi.org/10.1088/1674-4527/13/6/007
    # See file shadow4-tests/shadow4tests/devel/fit_bessel_kv.py
    if numpy.abs(nu - 1/3) < 1e-10:
        coeffs = [-0.31902416, -0.81577317, -0.78202672,  0.30405889,  0.70028439, -1.16431336,
  0.24015406, -0.0261485 ]
    elif numpy.abs(nu - 2/3) < 1e-10:
        coeffs = [-0.37896593, -0.34095854, -0.62947205,  0.05467015,  0.62890735, -1.07260337,
  1.66367831, -1.78893917]
    elif numpy.abs(nu - 5/3) < 1e-10:
        coeffs = [-2.35033577e-01,  2.17241138e-02, -7.04622366e-03,  9.65554026e-04,
  7.64819524e-01, -4.54068899e+00,  1.11791188e+01, -7.25096908e+00]
    else:
        raise Exception("Fit coefficients not available for nu=%f" % nu)

    gammanu = gamma(nu)

    H1 = coeffs[0] * x ** (1 / 1) + coeffs[1] * x ** (1 / 2) + coeffs[2] * x ** (1 / 3) + coeffs[3] * x ** (1 / 4)
    H2 = coeffs[4] * x ** (1 / 1) + coeffs[5] * x ** (1 / 2) + coeffs[6] * x ** (1 / 3) + coeffs[7] * x ** (1 / 4)
    delta1 = numpy.exp(H1)
    delta2 = 1 - numpy.exp(H2)

    A1 = 0.5 * gammanu * (x / 2)**(-nu)
    A2 = numpy.sqrt(numpy.pi / (2 * x)) * numpy.exp(-x)
    out =  A1 * delta1 + A2 * delta2
    return out

if __name__ == "__main__":
    import time

    if False:
        # test speed
        psi1 = 0.0
        psi_interval_number_of_points = 1001
        angle_array_reduced = numpy.linspace(-0.5 * psi1, 0.5 * psi1, psi_interval_number_of_points)

        t0 = time.time()
        for i in range(15000):
            tmp = sync_f_sigma_and_pi(angle_array_reduced, 100.)
        print("Time (accurate): ", time.time()-t0)

        t0 = time.time()
        for i in range(15000):
            tmp = sync_f_sigma_and_pi_approx(angle_array_reduced, 100.)
        print("Time (approximated): ", time.time()-t0)

    if False:
        # display/compare Kv approximated
        ji_interval_number_of_points = 1001
        ji_array = numpy.linspace(0, 5, ji_interval_number_of_points)

        k13 = kv(1.0 / 3.0, ji_array)
        k23 = kv(2.0 / 3.0, ji_array)
        k53 = kv(5.0 / 3.0, ji_array)

        k13approx = kv_approx(1.0 / 3.0, ji_array)
        k23approx = kv_approx(2.0 / 3.0, ji_array)
        k53approx = kv_approx(5.0 / 3.0, ji_array)

        from srxraylib.plot.gol import plot

        plot(ji_array, k13,
             ji_array, k23,
             ji_array, k53,
             ji_array, k13approx,
             ji_array, k23approx,
             ji_array, k53approx,
             legend=['k1/3','k2/3','k5/3', 'k1/3approx','k2/3approx','k5/3approx'],
             marker=[None, None, None, '+', '+', '+'], linestyle=[None,None,None,'','',''], xlog=0, ylog=0, yrange=[1e-18,10])


    if True:
        # test flatness (wiggler interpolation error).
        psi1 = 50.0
        psi_interval_number_of_points = 1001
        angle_array_reduced = numpy.linspace(-0.5 * psi1, 0.5 * psi1, psi_interval_number_of_points)

        eene = numpy.linspace(1, 1e4, 100)
        t0 = time.time()
        for i in range(eene.size):
            fm_s, fm_p = sync_f_sigma_and_pi(angle_array_reduced, eene[i])
            print(eene[i], fm_s.min(), fm_s.max(), fm_s.max() - fm_s.min())

