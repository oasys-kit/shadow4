#
# Tests of the python implementation of the shadow3/undulator preprocessors (SourceUndulatorFactory)
#


#
# switch on/off plots
#
DO_PLOT = False


import unittest
import numpy
import json

# CODE TO TEST
from SourceUndulator import SourceUndulator
from SourceUndulatorFactory import undul_phot, undul_cdf
from SourceUndulatorFactory import undul_phot_pysru,undul_phot_srw
# input/output
from SourceUndulatorInputOutput import load_uphot_dot_dat,write_uphot_dot_dat
from SourceUndulatorInputOutput import load_xshundul_dot_sha,write_xshundul_dot_sha

if DO_PLOT:
    try:
        from srxraylib.plot.gol import plot_image,plot, plot_show
    except:
        print("srxraylib not available (for plots). Plots switched off.")
        DO_PLOT = False
#
# auxiliary functions
#

# for test purposes only
# def run_shadow3_using_preprocessors(jsn):
#     from SourceUndulator import SourceUndulator
#     u = SourceUndulator()
#     u.load_json_shadowvui_dictionary(jsn)
#     u.run_using_preprocessors()




#
# Tests
#

class TestSourceUndulatorFactory(unittest.TestCase):

    def test_undul_phot(self):

        print("\n#                                                            ")
        print("# test_undul_phot  ")
        print("#                                                              ")

        h = {}
        h["E_ENERGY"] = 6.04
        h["INTENSITY"] = 0.2
        h["LAMBDAU"] = 0.032
        h["NPERIODS"] = 50
        h["K"] = 0.25
        h["EMIN"] = 10200.0
        h["EMAX"] = 10650.0
        h["NG_E"] = 11
        h["MAXANGLE"] = 0.015
        h["NG_T"] = 51
        h["NG_P"] = 11

        # internal code
        udict = undul_phot(E_ENERGY = h["E_ENERGY"],INTENSITY = h["INTENSITY"],
                                        LAMBDAU = h["LAMBDAU"],NPERIODS = h["NPERIODS"],K = h["K"],
                                        EMIN = h["EMIN"],EMAX = h["EMAX"],NG_E = h["NG_E"],
                                        MAXANGLE = h["MAXANGLE"],NG_T = h["NG_T"],
                                        NG_P = h["NG_P"])


        photon_energy = numpy.linspace(h["EMIN"],h["EMAX"],h["NG_E"],dtype=float)
        theta = numpy.linspace(0,h["MAXANGLE"]*1e-3,h["NG_T"],dtype=float)
        phi = numpy.linspace(0,numpy.pi/2,h["NG_P"],dtype=float)

        numpy.testing.assert_almost_equal(udict["photon_energy"],photon_energy)
        numpy.testing.assert_almost_equal(udict["theta"],theta)
        numpy.testing.assert_almost_equal(udict["phi"],phi)

        rad = udict["radiation"]
        pol = udict["polarization"]

        print("   radiation[1,1,2]", rad[1,1,2] )
        print("   radiation[1,5,7]", rad[1,5,7] )
        print("polarization[1,1,2]", pol[1,1,2] )
        print("polarization[1,5,7]", pol[1,5,7] )

        diff1 = (rad[1,1,2] - 4.42001096822e+20) / 4.42001096822e+20
        diff2 = (rad[1,5,7] - 3.99227535348e+20) / 3.99227535348e+20
        diff3 = (pol[1,1,2] - 0.999999776021) / 0.999999776021
        diff4 = (pol[1,5,7] - 0.99999794449) / 0.99999794449

        print("Relative difference    radiation[1,1,2]", diff1)
        print("Relative difference    radiation[1,5,7]", diff2)
        print("Relative difference polarization[1,1,2]", diff3)
        print("Relative difference polarization[1,5,7]", diff4)

        self.assertAlmostEqual(diff1,0.00,delta=1e-4)
        self.assertAlmostEqual(diff2,0.00,delta=1e-4)
        self.assertAlmostEqual(diff3,0.00,delta=1e-4)
        self.assertAlmostEqual(diff4,0.00,delta=1e-4)

    def test_undul_phot_NG_E_one(self):

        print("\n#                                                            ")
        print("# test_undul_phot  ")
        print("#                                                              ")

        h = {}
        h["E_ENERGY"] = 6.04
        h["INTENSITY"] = 0.2
        h["LAMBDAU"] = 0.032
        h["NPERIODS"] = 50
        h["K"] = 0.25
        h["EMIN"] = 10200.0
        h["EMAX"] = 10650.0
        h["NG_E"] = 1
        h["MAXANGLE"] = 0.015
        h["NG_T"] = 51
        h["NG_P"] = 11

        # internal code
        udict = undul_phot(E_ENERGY = h["E_ENERGY"],INTENSITY = h["INTENSITY"],
                                        LAMBDAU = h["LAMBDAU"],NPERIODS = h["NPERIODS"],K = h["K"],
                                        EMIN = h["EMIN"],EMAX = h["EMAX"],NG_E = h["NG_E"],
                                        MAXANGLE = h["MAXANGLE"],NG_T = h["NG_T"],
                                        NG_P = h["NG_P"])


        photon_energy = numpy.linspace(h["EMIN"],h["EMAX"],h["NG_E"],dtype=float)
        theta = numpy.linspace(0,h["MAXANGLE"]*1e-3,h["NG_T"],dtype=float)
        phi = numpy.linspace(0,numpy.pi/2,h["NG_P"],dtype=float)

        numpy.testing.assert_almost_equal(udict["photon_energy"],photon_energy)
        numpy.testing.assert_almost_equal(udict["theta"],theta)
        numpy.testing.assert_almost_equal(udict["phi"],phi)

        rad = udict["radiation"]
        pol = udict["polarization"]

        print("   radiation[0,1,2]", rad[0,1,2] )
        print("   radiation[0,5,7]", rad[0,5,7] )
        print("polarization[0,1,2]", pol[0,1,2] )
        print("polarization[0,5,7]", pol[0,5,7] )

        diff1 = (rad[0,1,2] - 8.26623377929e+20) / 8.26623377929e+20
        diff2 = (rad[0,5,7] - 8.22039896005e+20) / 8.22039896005e+20
        diff3 = (pol[0,1,2] - 0.999999688687) / 0.999999688687
        diff4 = (pol[0,5,7] - 0.999997183752) / 0.999997183752

        print("Relative difference    radiation[0,1,2]", diff1)
        print("Relative difference    radiation[0,5,7]", diff2)
        print("Relative difference polarization[0,1,2]", diff3)
        print("Relative difference polarization[0,5,7]", diff4)

        self.assertAlmostEqual(diff1,0.00,delta=1e-4)
        self.assertAlmostEqual(diff2,0.00,delta=1e-4)
        self.assertAlmostEqual(diff3,0.00,delta=1e-4)
        self.assertAlmostEqual(diff4,0.00,delta=1e-4)



    def test_undul_phot_pysru(self):
        print("\n#                                                            ")
        print("# test_undul_phot_pysru  ")
        print("#                                                              ")

        try:
            import pySRU
        except:
            print("......................Skipping: pySRU not available...")
            return


        hh = {}
        hh["E_ENERGY"] = 6.04
        hh["INTENSITY"] = 0.2
        hh["LAMBDAU"] = 0.032
        hh["NPERIODS"] = 50
        hh["K"] = 0.25
        hh["EMIN"] = 10200.0
        hh["EMAX"] = 10650.0
        hh["NG_E"] = 11
        hh["MAXANGLE"] = 0.015
        hh["NG_T"] = 51
        hh["NG_P"] = 11

        # internal code
        udict = undul_phot_pysru(E_ENERGY = hh["E_ENERGY"],INTENSITY = hh["INTENSITY"],
                                        LAMBDAU = hh["LAMBDAU"],NPERIODS = hh["NPERIODS"],K = hh["K"],
                                        EMIN = hh["EMIN"],EMAX = hh["EMAX"],NG_E = hh["NG_E"],
                                        MAXANGLE = hh["MAXANGLE"],NG_T = hh["NG_T"],
                                        NG_P = hh["NG_P"])


        photon_energy = numpy.linspace(hh["EMIN"],hh["EMAX"],hh["NG_E"],dtype=float)
        theta = numpy.linspace(0,hh["MAXANGLE"]*1e-3,hh["NG_T"],dtype=float)
        phi = numpy.linspace(0,numpy.pi/2,hh["NG_P"],dtype=float)

        numpy.testing.assert_almost_equal(udict["photon_energy"],photon_energy)
        numpy.testing.assert_almost_equal(udict["theta"],theta)
        numpy.testing.assert_almost_equal(udict["phi"],phi)

        rad = udict["radiation"]
        pol = udict["polarization"]

        print("   radiation[1,1,2]", rad[1,1,2] )
        print("   radiation[1,5,7]", rad[1,5,7] )
        print("polarization[1,1,2]", pol[1,1,2] )
        print("polarization[1,5,7]", pol[1,5,7] )

        diff1 = (rad[1,1,2] - 4.42891350998e+20) / 4.42891350998e+20
        diff2 = (rad[1,5,7] - 4.00147694552e+20) / 4.00147694552e+20
        diff3 = (pol[1,1,2] - 0.999999778472) / 0.999999778472
        diff4 = (pol[1,5,7] - 0.999997902917) / 0.999997902917

        print("Relative difference    radiation[1,1,2]", diff1)
        print("Relative difference    radiation[1,5,7]", diff2)
        print("Relative difference polarization[1,1,2]", diff3)
        print("Relative difference polarization[1,5,7]", diff4)

        self.assertAlmostEqual(diff1,0.00,delta=1e-4)
        self.assertAlmostEqual(diff2,0.00,delta=1e-4)
        self.assertAlmostEqual(diff3,0.00,delta=1e-4)
        self.assertAlmostEqual(diff4,0.00,delta=1e-4)


    def test_comparison_undul_phot(self,do_plot_intensity=DO_PLOT,do_plot_polarization=DO_PLOT,do_plot_trajectory=DO_PLOT):

        print("\n#                                                            ")
        print("# test_comparison_undul_phot  ")
        print("#                                                              ")

        #
        # test undul_phot (undulator radiation)
        #

        try:
            import pySRU
            is_available_pysru = True
        except:
            is_available_pysru = False

        try:
            import srwlib
            is_available_srw = True
        except:
            is_available_srw = False


        # "EMIN":       10500.0000,
        # "EMAX":       10550.0000,

        tmp = \
            """
            {
            "LAMBDAU":     0.0320000015,
            "K":      0.250000000,
            "E_ENERGY":       6.03999996,
            "E_ENERGY_SPREAD":    0.00100000005,
            "NPERIODS": 50,
            "EMIN":       10200.0000,
            "EMAX":       10650.0000,
            "INTENSITY":      0.2,
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

        # SHADOW3 preprocessor
        # run_shadow3_using_preprocessors(h)

        u = SourceUndulator()
        u.load_json_shadowvui_dictionary(h)
        u.run_using_preprocessors()



        undul_phot_preprocessor_dict = load_uphot_dot_dat("uphot.dat")

        if do_plot_intensity: plot_image(undul_phot_preprocessor_dict['radiation'][0,:,:],undul_phot_preprocessor_dict['theta']*1e6,undul_phot_preprocessor_dict['phi']*180/numpy.pi,
                   title="INTENS UNDUL_PHOT_PREPROCESSOR: RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)

        if do_plot_polarization: plot_image(undul_phot_preprocessor_dict['polarization'][0,:,:],undul_phot_preprocessor_dict['theta']*1e6,undul_phot_preprocessor_dict['phi']*180/numpy.pi,
                   title="POL_DEG UNDUL_PHOT_PREPROCESSOR: RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)


        # internal code
        undul_phot_dict = undul_phot(E_ENERGY = h["E_ENERGY"],INTENSITY = h["INTENSITY"],
                                        LAMBDAU = h["LAMBDAU"],NPERIODS = h["NPERIODS"],K = h["K"],
                                        EMIN = h["EMIN"],EMAX = h["EMAX"],NG_E = h["NG_E"],
                                        MAXANGLE = h["MAXANGLE"],NG_T = h["NG_T"],
                                        NG_P = h["NG_P"])

        if do_plot_intensity: plot_image(undul_phot_dict['radiation'][0,:,:],undul_phot_dict['theta']*1e6,undul_phot_dict['phi']*180/numpy.pi,
                   title="INTENS UNDUL_PHOT: RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)
        if do_plot_polarization: plot_image(undul_phot_dict['polarization'][0,:,:],undul_phot_dict['theta']*1e6,undul_phot_dict['phi']*180/numpy.pi,
                   title="POL_DEG UNDUL_PHOT: RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)

        # pySRU
        if is_available_pysru:
            undul_phot_pysru_dict = undul_phot_pysru(E_ENERGY = h["E_ENERGY"],INTENSITY = h["INTENSITY"],
                                            LAMBDAU = h["LAMBDAU"],NPERIODS = h["NPERIODS"],K = h["K"],
                                            EMIN = h["EMIN"],EMAX = h["EMAX"],NG_E = h["NG_E"],
                                            MAXANGLE = h["MAXANGLE"],NG_T = h["NG_T"],
                                            NG_P = h["NG_P"])
            if do_plot_intensity: plot_image(undul_phot_pysru_dict['radiation'][0,:,:],undul_phot_pysru_dict['theta']*1e6,undul_phot_pysru_dict['phi']*180/numpy.pi,
                       title="INTENS UNDUL_PHOT_PYSRU: RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)
            if do_plot_polarization: plot_image(undul_phot_pysru_dict['polarization'][0,:,:],undul_phot_pysru_dict['theta']*1e6,undul_phot_pysru_dict['phi']*180/numpy.pi,
                       title="POL_DEG UNDUL_PHOT_PYSRU: RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)

        # srw
        if is_available_srw:
            undul_phot_srw_dict = undul_phot_srw(E_ENERGY = h["E_ENERGY"],INTENSITY = h["INTENSITY"],
                                            LAMBDAU = h["LAMBDAU"],NPERIODS = h["NPERIODS"],K = h["K"],
                                            EMIN = h["EMIN"],EMAX = h["EMAX"],NG_E = h["NG_E"],
                                            MAXANGLE = h["MAXANGLE"],NG_T = h["NG_T"],
                                            NG_P = h["NG_P"])
            if do_plot_intensity: plot_image(undul_phot_srw_dict['radiation'][0,:,:],undul_phot_srw_dict['theta']*1e6,undul_phot_srw_dict['phi']*180/numpy.pi,
                       title="INTENS UNDUL_PHOT_SRW: RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)
            if do_plot_polarization: plot_image(undul_phot_srw_dict['polarization'][0,:,:],undul_phot_srw_dict['theta']*1e6,undul_phot_srw_dict['phi']*180/numpy.pi,
                       title="POL_DEG UNDUL_PHOT_SRW: RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)


        x = undul_phot_dict["photon_energy"]
        y0 = (undul_phot_preprocessor_dict["radiation"]).sum(axis=2).sum(axis=1)
        y1 =              (undul_phot_dict["radiation"]).sum(axis=2).sum(axis=1)
        if is_available_pysru: y2 = (undul_phot_pysru_dict["radiation"]).sum(axis=2).sum(axis=1)
        if is_available_srw:   y3 = (undul_phot_srw_dict["radiation"]).sum(axis=2).sum(axis=1)

        if do_plot_intensity:
            if is_available_pysru and is_available_srw:
                plot(x,y0,x,y1,x,y2,x,y3,xtitle="Photon energy [eV]",ytitle="Flux[photons/s/eV/rad^2]",legend=["preprocessor","internal","pySRU","SRW"])
            else:
                if is_available_pysru:
                    plot(x,y0,x,y1,x,y2,xtitle="Photon energy [eV]",ytitle="Flux[photons/s/eV/rad^2]",legend=["preprocessor","internal","pySRU"])
                if is_available_srw:
                    plot(x,y0,x,y1,x,y3,xtitle="Photon energy [eV]",ytitle="Flux[photons/s/eV/rad^2]",legend=["preprocessor","internal","SRW"])
        tmp = numpy.where(y0 > 0.1)

        print("\n")
        print(">>> test_undul_phot: preprocessor/internal: %4.2f %% "%(numpy.average( 100*numpy.abs((y0[tmp]-y1[tmp])/y1[tmp]) )))
        self.assertLess( numpy.average( 100*numpy.abs((y0[tmp]-y1[tmp])/y1[tmp]) ), 5 )
        if is_available_pysru:
            print(">>> test_undul_phot:        pySRU/internal: %4.2f %% "%(numpy.average( 100*numpy.abs((y2[tmp]-y1[tmp])/y1[tmp]) )))
            self.assertLess( numpy.average( 100*numpy.abs((y2[tmp]-y1[tmp])/y1[tmp]) ), 1 )
        if is_available_srw:
            print(">>> test_undul_phot:          SRW/internal: %4.2f %% "%(numpy.average( 100*numpy.abs((y3[tmp]-y1[tmp])/y1[tmp]) )))
            self.assertLess( numpy.average( 100*numpy.abs((y3[tmp]-y1[tmp])/y1[tmp]) ), 5 )


        #
        # trajectory
        #
        if do_plot_trajectory:
            # Trajectory is only in undul_phot (internal) and und_phot_pysru
            # t0 = (undul_phot_preprocessor_dict["trajectory"])
            t1 = (             undul_phot_dict["trajectory"])
            if is_available_pysru:
                t2 = (       undul_phot_pysru_dict["trajectory"])
                plot(t1[3],1e6*t1[1],t2[3],1e6*t2[1],title='Trajectory',xtitle='z [m]',ytitle='x[um]',legend=['internal','pysru'],show=False)
            else:
                plot(t1[3],1e6*t1[1],title='Trajectory',xtitle='z [m]',ytitle='x[um]',legend=['internal'],show=False)


        if do_plot_polarization or do_plot_intensity or do_plot_trajectory: plot_show()



    def test_undul_cdf(self,do_plot=DO_PLOT):

        print("\n#                                                            ")
        print("# test_undul_cdf  ")
        print("#                                                              ")
        tmp = \
            """
            {
            "LAMBDAU":     0.0320000015,
            "K":      0.250000000,
            "E_ENERGY":       6.03999996,
            "E_ENERGY_SPREAD":    0.00100000005,
            "NPERIODS": 50,
            "EMIN":       10200.0000,
            "EMAX":       10650.0000,
            "INTENSITY":      0.2,
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
        #
        # run_shadow3_using_preprocessors(h) # uphot.dat must exist

        u = SourceUndulator()
        u.load_json_shadowvui_dictionary(h)
        u.run_using_preprocessors()

        #
        #
        #
        radiation = load_uphot_dot_dat(file_in="uphot.dat")

        cdf2 = undul_cdf(radiation,method='sum',do_plot=False)
        write_xshundul_dot_sha(cdf2,file_out="xshundul2.sha")


        cdf3 = undul_cdf(radiation,method='trapz',do_plot=False)
        write_xshundul_dot_sha(cdf3,file_out="xshundul3.sha")

        cdf1 = load_xshundul_dot_sha(file_in="xshundul.sha", do_plot=False,show=False)
        cdf2 = load_xshundul_dot_sha(file_in="xshundul2.sha",do_plot=False,show=False)
        cdf3 = load_xshundul_dot_sha(file_in="xshundul3.sha",do_plot=False,show=False)




        ZERO1 = cdf1['cdf_Energy']
        ONE1 = cdf1['cdf_EnergyTheta']
        TWO1 = cdf1['cdf_EnergyThetaPhi']

        ZERO2 = cdf2['cdf_Energy']
        ONE2 = cdf2['cdf_EnergyTheta']
        TWO2 = cdf2['cdf_EnergyThetaPhi']

        ZERO3 = cdf3['cdf_Energy']
        ONE3 = cdf3['cdf_EnergyTheta']
        TWO3 = cdf3['cdf_EnergyThetaPhi']

        tmp = numpy.where(ZERO1 > 0.1*ZERO1.max())
        print("test_undul_cdf: ZERO:   sum/shadow3 %4.2f %%: "%(numpy.average( 100*numpy.abs((ZERO2[tmp]-ZERO1[tmp])/ZERO1[tmp]) )))
        print("test_undul_cdf: ZERO: trapz/shadow3 %4.2f %%: "%(numpy.average( 100*numpy.abs((ZERO3[tmp]-ZERO1[tmp])/ZERO1[tmp]) )))

        tmp = numpy.where(ONE1 > 0.1*ONE1.max())
        print(r"test_undul_cdf: ONE:   sum/shadow3 %4.2f %%: "%(numpy.average( 100*numpy.abs((ONE2[tmp]-ONE1[tmp])/ONE1[tmp]) )))
        print(r"test_undul_cdf: ONE: trapz/shadow3 %4.2f %%: "%(numpy.average( 100*numpy.abs((ONE3[tmp]-ONE1[tmp])/ONE1[tmp]) )))

        tmp = numpy.where(TWO1 > 0.1*TWO1.max())
        print("test_undul_cdf: TWO:   sum/shadow3 %4.2f %%: "%(numpy.average( 100*numpy.abs((TWO2[tmp]-TWO1[tmp])/TWO1[tmp]) )))
        print("test_undul_cdf: TWO: trapz/shadow3 %4.2f %%: "%(numpy.average( 100*numpy.abs((TWO3[tmp]-TWO1[tmp])/TWO1[tmp]) )))


        if do_plot:

            plot(cdf1["energy"],cdf1["cdf_EnergyThetaPhi"],cdf2["energy"],cdf2["cdf_EnergyThetaPhi"],cdf3["energy"],cdf3["cdf_EnergyThetaPhi"],
                title="cdf vs energy ",xtitle="photon energy [eV]",ytitle="cdf (integrated in theta,phi)",
                legend=["preprocessor","internal by sumation","internal by trapezoidal integration"],show=False)


            plot_image(cdf1['cdf_Energy'][0,:,:],1e6*radiation['theta'],radiation['phi'],title="PREPROCESSORS cdf_Energy[0]",
                       xtitle="Theta [urad]",ytitle="Phi",aspect='auto',show=False)
            plot_image(cdf3['cdf_Energy'][0,:,:],1e6*radiation['theta'],radiation['phi'],title="internal-trapezoidal cdf_Energy[0]",
                       xtitle="Theta [urad]",ytitle="Phi",aspect='auto',show=False)

            plot_image(cdf1['cdf_EnergyTheta'],1e6*radiation['theta'],radiation['phi'],title="PREPROCESSORS cdf_EnergyTheta",
                       xtitle="Theta [urad]",ytitle="Phi",aspect='auto',show=False)
            plot_image(cdf3['cdf_EnergyTheta'],1e6*radiation['theta'],radiation['phi'],title="internal-trapezoidal cdf_EnergyTheta",
                       xtitle="Theta [urad]",ytitle="Phi",aspect='auto',show=False)
            plot_show()

    def test_undul_cdf_NG_E_one(self,do_plot=DO_PLOT):

        print("\n#                                                            ")
        print("# test_undul_cdf_NG_E_one  ")
        print("#                                                              ")

        case = 0 # 2 = vui file

        if case == 0:

            tmp = \
                """
                {
                "LAMBDAU":     0.0320000015,
                "K":      0.250000000,
                "E_ENERGY":       6.03999996,
                "E_ENERGY_SPREAD":    0.00100000005,
                "NPERIODS": 50,
                "EMIN":       10498.0000,
                "EMAX":       10499.0000,
                "INTENSITY":      0.200000003,
                "MAXANGLE":      0.100000001,
                "NG_E": 1,
                "NG_T": 51,
                "NG_P": 11,
                "NG_PLOT(1)":"0",
                "NG_PLOT(2)":"No",
                "NG_PLOT(3)":"Yes",
                "UNDUL_PHOT_FLAG(1)":"0",
                "UNDUL_PHOT_FLAG(2)":"Shadow code",
                "UNDUL_PHOT_FLAG(3)":"Urgent code",
                "UNDUL_PHOT_FLAG(4)":"SRW code",
                "UNDUL_PHOT_FLAG(5)":"Gaussian Approximation",
                "UNDUL_PHOT_FLAG(6)":"ESRF python code",
                "SEED": 36255,
                "SX":     0.0399999991,
                "SZ":    0.00100000005,
                "EX":   4.00000005E-07,
                "EZ":   3.99999989E-09,
                "FLAG_EMITTANCE(1)":"0",
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

            u = SourceUndulator()
            u.load_json_shadowvui_dictionary(h)
            u.run_using_preprocessors()

        elif case == 1:

            # some inputs
            E_ENERGY=6.04
            E_ENERGY_SPREAD=0.001
            INTENSITY=0.2
            SX=0.04
            SZ=0.001
            SXP=10e-6
            SZP=4e-6
            FLAG_EMITTANCE=1
            LAMBDAU=0.032
            NPERIODS=50
            K=0.25
            EMIN= 10498.0000
            EMAX= 10499.0000
            NG_E=101
            MAXANGLE=0.1
            NG_T=51
            NG_P=11
            SEED=36255
            NRAYS=15000

            u = SourceUndulator()

            u.set_from_keywords(
                E_ENERGY = E_ENERGY,
                E_ENERGY_SPREAD = E_ENERGY_SPREAD,
                INTENSITY = INTENSITY,
                SX = SX,
                SZ = SZ,
                SXP = SXP,
                SZP = SZP,
                FLAG_EMITTANCE = FLAG_EMITTANCE,
                LAMBDAU = LAMBDAU,
                NPERIODS = NPERIODS,
                K = K,
                EMIN = EMIN,
                EMAX = EMAX,
                NG_E = NG_E,
                MAXANGLE = MAXANGLE,
                NG_T = NG_T,
                NG_P = NG_P,
                SEED = SEED,
                NRAYS = NRAYS,
                )

            u.set_energy_monochromatic_at_resonance(harmonic_number=1)
        elif case == 2:

            u = SourceUndulator()
            u.load_json_shadowvui_file("xshundul.json")
            u.run_using_preprocessors()

        else:
            raise Exception("Undefined")

        #
        #
        #
        radiation = load_uphot_dot_dat(file_in="uphot.dat")

        cdf2 = undul_cdf(radiation,method='sum',do_plot=False)
        write_xshundul_dot_sha(cdf2,file_out="xshundul2.sha")


        cdf3 = undul_cdf(radiation,method='trapz',do_plot=False)
        write_xshundul_dot_sha(cdf3,file_out="xshundul3.sha")

        cdf1 = load_xshundul_dot_sha(file_in="xshundul.sha", do_plot=False,show=False)
        cdf2 = load_xshundul_dot_sha(file_in="xshundul2.sha",do_plot=False,show=False)
        cdf3 = load_xshundul_dot_sha(file_in="xshundul3.sha",do_plot=False,show=False)


        ZERO1 = cdf1['cdf_Energy']
        ONE1 = cdf1['cdf_EnergyTheta']
        TWO1 = cdf1['cdf_EnergyThetaPhi']

        ZERO2 = cdf2['cdf_Energy']
        ONE2 = cdf2['cdf_EnergyTheta']
        TWO2 = cdf2['cdf_EnergyThetaPhi']

        ZERO3 = cdf3['cdf_Energy']
        ONE3 = cdf3['cdf_EnergyTheta']
        TWO3 = cdf3['cdf_EnergyThetaPhi']

        #
        #
        NG_E = (radiation["photon_energy"]).size
        NG_T = (radiation["theta"]).size
        NG_P = (radiation["phi"]).size
        #
        print("-----> shape comparison",NG_E,NG_T,NG_P,ZERO1.shape,ZERO2.shape,ZERO3.shape)
        print("-----> shape comparison",NG_E,NG_T,NG_P,ONE1.shape,ONE2.shape,ONE3.shape)
        print("-----> shape comparison",NG_E,NG_T,NG_P,TWO1.shape,TWO2.shape,TWO3.shape)
        for ie in range(NG_E):
            for it in range(NG_T):
                for ip in range(NG_P):
                    # print(">>>>",TWO1[ie],TWO2[ie],TWO3[ie])
                    print(">>>>[%d %d %d]>>"%(ie,it,ip),ONE1[ie,it],ONE2[ie,it],ONE3[ie,it])
                    # print(">>>>[%d %d %d]>>"%(ie,it,ip),ZERO1[ie,it,ip],ZERO2[ie,it,ip],ZERO3[ie,it,ip])

        tmp = numpy.where(ZERO1 > 0.1*ZERO1.max())
        print("test_undul_cdf: ZERO:   sum/shadow3 %4.2f %%: "%(numpy.average( 100*numpy.abs((ZERO2[tmp]-ZERO1[tmp])/ZERO1[tmp]) )))
        print("test_undul_cdf: ZERO: trapz/shadow3 %4.2f %%: "%(numpy.average( 100*numpy.abs((ZERO3[tmp]-ZERO1[tmp])/ZERO1[tmp]) )))

        tmp = numpy.where(ONE1 > 0.1*ONE1.max())
        print(r"test_undul_cdf: ONE:   sum/shadow3 %4.2f %%: "%(numpy.average( 100*numpy.abs((ONE2[tmp]-ONE1[tmp])/ONE1[tmp]) )))
        print(r"test_undul_cdf: ONE: trapz/shadow3 %4.2f %%: "%(numpy.average( 100*numpy.abs((ONE3[tmp]-ONE1[tmp])/ONE1[tmp]) )))


        if do_plot:


            plot_image(cdf1['cdf_Energy'][0,:,:],1e6*radiation['theta'],radiation['phi'],title="PREPROCESSORS cdf_Energy[0]",
                       xtitle="Theta [urad]",ytitle="Phi",aspect='auto',show=False)
            plot_image(cdf3['cdf_Energy'][0,:,:],1e6*radiation['theta'],radiation['phi'],title="internal-trapezoidal cdf_Energy[0]",
                       xtitle="Theta [urad]",ytitle="Phi",aspect='auto',show=False)

            plot_image(cdf1['cdf_EnergyTheta'],1e6*radiation['theta'],radiation['phi'],title="PREPROCESSORS cdf_EnergyTheta",
                       xtitle="Theta [urad]",ytitle="Phi",aspect='auto',show=False)
            plot_image(cdf3['cdf_EnergyTheta'],1e6*radiation['theta'],radiation['phi'],title="internal-trapezoidal cdf_EnergyTheta",
                       xtitle="Theta [urad]",ytitle="Phi",aspect='auto',show=False)
            plot_show()
