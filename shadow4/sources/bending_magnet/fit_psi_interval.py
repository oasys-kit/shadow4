"""
TODO: delete this file
"""

import numpy
from srxraylib.sources.srfunc import  sync_f
from srxraylib.plot.gol import plot

#
# fits the dependence of the width of the emission vs energy to get an estimation of the sampling interval
#


# RENERGY = numpy.array([0.01,0.1,1.0])
RENERGY = numpy.logspace(-6,2,100)  # energies from 1e-6 to 100 times the critical energy
OUT = RENERGY * 0.0

for i,rEnergy in enumerate(RENERGY):


    # rEnergy = 0.01

    rAngle = numpy.linspace(-1000,1000,2000)
    f_s   = sync_f(rAngle, rEnergy=rEnergy, polarization=1, gauss=0, l2=1, l3=0)

    #CALCULATE fwhm
    tt = numpy.where(f_s>=max(f_s)*0.01)  # tjis is the width at 0.01 of the height
    if f_s[tt].size > 1:
        binSize = rAngle[1]-rAngle[0]
        fwhm = binSize*(tt[0][-1]-tt[0][0])
        fwhm_coordinates = (rAngle[tt[0][0]],rAngle[tt[0][-1]])
    print("FWHM: ",fwhm)
    print("FWHM coordinates: ",fwhm_coordinates)
    OUT[i] = fwhm_coordinates[1]

    if False:
        f_tot = sync_f(rAngle, rEnergy=rEnergy, polarization=0, gauss=0, l2=1, l3=0)
        f_p   = sync_f(rAngle, rEnergy=rEnergy, polarization=2, gauss=0, l2=1, l3=0)
        plot(rAngle, f_s/f_tot[rAngle.size//2], rAngle, f_tot/f_s[rAngle.size//2],title="E/Ec=%f (%f)"%(rEnergy,fwhm_coordinates[1]))





# plot(RENERGY,OUT,xlog=1)



x = numpy.log10(RENERGY)
y = numpy.log10(OUT)

c = numpy.polyfit(x,y,1)

y_fit = c[1]+c[0]*x
print(c)

plot(x,y, x, y_fit, xtitle="log(E/Ec)",ytitle="log of width in units of gamma",legend=["data","fit"] )


# for i in range(RENERGY.size):
#     print("%f  %f"%(RENERGY[i],OUT[i]))
#
# plot(x,numpy.log10(y))
#
plot(numpy.log10(RENERGY),OUT,
        numpy.log10(RENERGY),10**y_fit,
        xlog=0,xtitle="log(E/Ec)",ytitle="width in units of gamma",legend=["data","fit"])


#example

rEnergy = 0.4 / 3000.0
x = numpy.log10(rEnergy)
y_fit = c[1] + c[0] * x
psi_interval_in_units_one_over_gamma = int(10 ** y_fit)  # this is the semiinterval
psi_interval_in_units_one_over_gamma *= 4
if psi_interval_in_units_one_over_gamma < 4:
    psi_interval_in_units_one_over_gamma = 4
print("psi_interval_in_units_one_over_gamma: ",psi_interval_in_units_one_over_gamma)

xx = numpy.linspace(-0.5*psi_interval_in_units_one_over_gamma,0.5*psi_interval_in_units_one_over_gamma,1001)
f_s   = sync_f( xx , rEnergy=rEnergy, polarization=1, gauss=0, l2=1, l3=0)

plot(xx,f_s,title="Example Ec=%f"%rEnergy)
