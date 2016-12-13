#
# Working version for Undulator preprocessors in python
#
# this code replaces SHADOW's undul_phot.
#
# It calculates the undulator radiation as a function of energy, theta and phi. Phi is the polar angle.
#
# Inputs: from xshundul.json (written by ShadowVui)
# Output: file uphot.dat, like in shadow
#
#
# TODO:
#    Follow changes in Sophie's code
#    Calculate degree of polarization (set to one now)
#
#
import numpy as np
import scipy.constants as codata
import json
import sys

angstroms_to_eV = codata.h*codata.c/codata.e*1e10


#
# COPY OF SOPHIE'S CODE
#

def trajectory_undulator_reference(K=1.87 , gamma=2544.03131115, lambda_u=0.020, Nb_period=10, Nb_point=10, Beta_et=0.99993):
    N = Nb_period * Nb_point + 1
    ku= 2.0*np.pi/lambda_u
    omega_u = Beta_et*codata.c *ku
    # trajectory =
    #   [t........]
    # 	[ X/c......]
    # 	[ Y/c ......]
    #	[ Z/c ......]
    # 	[ Vx/c .....]
    # 	[ Vy/c .....]
    # 	[ Vz/c .....]
    # 	[ Ax/c .....]
    # 	[ Ay/c .....]
    # 	[ Az/c .....]
    trajectory = np.zeros((10, N))
    # t
    trajectory[0] = np.linspace(-(lambda_u / (codata.c * Beta_et)) * (Nb_period / 2), (lambda_u / (codata.c * Beta_et)) * (Nb_period / 2), N)
    #trajectory[0] = np.linspace(0.0,(lambda_u / (c * Beta_et)) * (Nb_period), N)
    # X et Z en fct de t
    trajectory[3] = Beta_et*trajectory[0] - ((K/gamma)**2) * (1.0/(8.0*ku*codata.c))* np.sin(2.0*omega_u * trajectory[0])
    trajectory[1] = (K/(gamma*ku*codata.c))* np.sin(omega_u * trajectory[0])
    # Vx et Vz en fct de t
    trajectory[6] = Beta_et - ((K/gamma)**2) * (2.0*omega_u/(8.0*ku*codata.c ))*np.cos(2.0*omega_u * trajectory[0])
    trajectory[4] = (K/(gamma*ku*codata.c))*omega_u* np.cos(omega_u * trajectory[0])
    # Ax et Az en fct de t
    trajectory[9] = ((2.0*omega_u*K/gamma)**2) * (1.0/(8.0*ku*codata.c))*np.sin(2.0*omega_u * trajectory[0])
    trajectory[7] = -(K/(gamma*ku*codata.c))*(omega_u**2)* np.sin(omega_u * trajectory[0])
    # trajectory *= codata.c
    # trajectory[0] *= (1.0/codata.c)
    return trajectory

def energy_radiated(omega=2.53465927101*10**17,trajectory=np.zeros((11,10)) , x=0.00 , y=0.0, D=None):
    N = trajectory.shape[1]

    if D == None:
        # in radian :
        n_chap = np.array([x, y, 1.0 - 0.5 * (x ** 2 + y ** 2)])
    else:
        # in meters :
        R = np.sqrt(x ** 2 + y ** 2 + D ** 2)
        n_chap = np.array([x, y, D]) / R

    E = np.full((3,), 0. + 1j * 0., dtype=np.complex)
    integrand = np.full((3,N), 0. + 1j * 0., dtype=np.complex)
    Alpha=trajectory[7]*(n_chap[2]-trajectory[6]) - (n_chap[0]-trajectory[4])*trajectory[9]
    Alpha2=np.exp(0. + 1j * omega * (trajectory[0] - n_chap[0]*trajectory[1]-n_chap[2]*trajectory[3]))
    integrand[0] += (-(n_chap[1]**2)*trajectory[7]-n_chap[2]*Alpha)*Alpha2
    integrand[1] += n_chap[1]*(n_chap[0]*trajectory[7]+n_chap[2]*trajectory[9])*Alpha2
    integrand[2] += (-(n_chap[1]**2)*trajectory[9]+n_chap[0]*Alpha)*Alpha2
    integrand *= (1.0 / (1.0 - n_chap[0]*trajectory[4]-n_chap[2]*trajectory[6])) ** 2
    for k in range(3):
        E[k] = np.trapz(integrand[k], trajectory[0])
        #E[k] = integrate.simps(integrand[k], trajectory[0])
    # np.linalg.norm
    return (np.abs(E[0]) ** 2 + np.abs(E[1])** 2 + np.abs(E[2])** 2)
#
# END COPY OF SOPHIE'S CODE
#

def undul_phot(myinput):
    #
    # read inputs from a file created by ShadowVUI ----------------------------
    #
    if isinstance(myinput,str):
        inFileTxt = myinput # "xshundul.json"
        with open(inFileTxt, mode='r') as f1:
            h = json.load(f1)

    elif isinstance(myinput,dict):
        h = myinput
    else:
        raise Exception("Unknown input")


    # list all non-empty keywords
    print ("-----------------------------------------------------")
    for i,j in h.items():
        if (j != None):
            print ("%s = %s" % (i,j))
    print ("-----------------------------------------------------")
    print ("k: ",h['K'])

    #
    # calculate trajectory
    #
    gamma = h["E_ENERGY"] * 1e9 / 0.511e6
    print("GAMMA:",gamma)
    Beta = np.sqrt(1.0 - (1.0 / gamma ** 2))
    Beta_et = Beta * (1.0 - (h['K'] / (2.0 * gamma)) ** 2)
    T = trajectory_undulator_reference(K=h['K'], gamma=gamma, lambda_u=h["LAMBDAU"], Nb_period=h["NPERIODS"],
                                        Nb_point=20,Beta_et=Beta_et)



    E = np.linspace(h["EMIN"],h["EMAX"],h["NG_E"])
    wavelength_array_in_A = angstroms_to_eV / E
    omega_array = 2*np.pi * codata.c / (wavelength_array_in_A * 1e-10)

    #
    # cartesian grid
    #
    gridding = 1 # 0=cartesian, 1=polar

    if gridding == 0:
        D = None
        X = np.linspace(0.0,h["MAXANGLE"]*1e-3,h["NG_T"])
        Y = np.linspace(0.0,h["MAXANGLE"]*1e-3,h["NG_T"])

    else:
        D = 100.0 # placed far away (100 m)
        theta = np.linspace(0,h["MAXANGLE"]*1e-3,h["NG_T"])
        phi = np.linspace(0,np.pi/2,h["NG_P"])

    c6= codata.e*1e-10/(8.0*np.pi**2*codata.epsilon_0*codata.c*codata.h)

    if gridding == 0:
        Z2 = np.zeros((E.size,X.size,Y.size))
        for o in range(omega_array.size):
            print("Calculating energy %g eV (%d of %d)"%(E[o],o+1,omega_array.size))
            for i in range(X.size):
                for j in range(Y.size):
                    Z2[o,i,j] = c6*energy_radiated(omega=omega_array[o],trajectory=T , x=X[i] , y=Y[j], D=D )
    elif gridding == 1:
        Z2 = np.zeros((omega_array.size,theta.size,phi.size))
        for o in range(omega_array.size):
            print("Calculating energy %g eV (%d of %d)"%(E[o],o+1,omega_array.size))
            for t in range(theta.size):
                for p in range(phi.size):
                    R = D / np.cos(theta[t])
                    r = R * np.sin(theta[t])
                    X = r * np.cos(phi[p])
                    Y = r * np.sin(phi[p])
                    Z2[o,t,p] = c6*energy_radiated(omega=omega_array[o],trajectory=T , x=X , y=Y, D=D )


    #
    # create uphot.dat file (like in SHADOW undul_phot)
    #
    file_out = "uphot.dat"
    f = open(file_out,'w')
    f.write("%d  %d  %d \n"%(h["NG_E"],h["NG_T"],h["NG_P"]))
    for e in E:
        f.write("%g \n"%(e))

    for e in E:
        for t in theta:
            f.write("%g \n"%t)

    for e in E:
        for t in theta:
            for p in phi:
                f.write("%g \n"%p)


    for e in range(E.size):
        for t in range(theta.size):
            for p in range(phi.size):
                f.write("%g \n"%Z2[e,t,p])

    for e in range(E.size):
        for t in range(theta.size):
            for p in range(phi.size):
                f.write("1.0 \n")

    f.close()
    print("File written to disk: %s"%file_out)


def test_undul_phot():
    tmp = \
        """
        {
        "LAMBDAU":     0.0320000015,
        "K":      0.250000000,
        "E_ENERGY":       6.03999996,
        "E_ENERGY_SPREAD":    0.00100000005,
        "NPERIODS": 50,
        "EMIN":       10500.0000,
        "EMAX":       10550.0000,
        "INTENSITY":      0.200000003,
        "MAXANGLE":     0.0149999997,
        "NG_E": 11,
        "NG_T": 51,
        "NG_P": 11,
        "NG_PLOT(1)":"1",
        "NG_PLOT(2)":"No",
        "NG_PLOT(3)":"Yes",
        "UNDUL_PHOT_FLAG(1)":"4",
        "UNDUL_PHOT_FLAG(2)":"Shadow code",
        "UNDUL_PHOT_FLAG(3)":"Urgent code",
        "UNDUL_PHOT_FLAG(4)":"SRW code",
        "UNDUL_PHOT_FLAG(5)":"Gaussian Approx",
        "UNDUL_PHOT_FLAG(6)":"python code by Sophie",
        "SEED": 36255,
        "SX":     0.0399999991,
        "SZ":    0.00100000005,
        "EX":   4.00000005E-07,
        "EZ":   3.99999989E-09,
        "FLAG_EMITTANCE(1)":"1",
        "FLAG_EMITTANCE(2)":"No",
        "FLAG_EMITTANCE(3)":"Yes",
        "NRAYS": 15000,
        "F_BOUND_SOUR": 0,
        "FILE_BOUND":"NONESPECIFIED",
        "SLIT_DISTANCE":       1000.00000,
        "SLIT_XMIN":      -1.00000000,
        "SLIT_XMAX":       1.00000000,
        "SLIT_ZMIN":      -1.00000000,
        "SLIT_ZMAX":       1.00000000,
        "NTOTALPOINT": 10000000,
        "JUNK4JSON":0
        }
        """
    h = json.loads(tmp)
    undul_phot(h)
    tmp = np.loadtxt("uphot.dat",skiprows=1)
    print("Obtained result[700]: %g (comparing to 6.09766e+16)"%tmp[7000])
    assert( np.abs(tmp[7000] - 6.09766e+16) < 1e13)


if __name__ == "__main__":
    undul_phot("xshundul.json")
    # test_undul_phot()
