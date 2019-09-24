__authors__ = ["M Sanchez del Rio - ESRF ISDD Advanced Analysis and Modelling"]
__license__ = "MIT"
__date__ = "12/01/2017"
#
# Tests of the python implementation of the shadow3/undulator preprocessors (SourceUndulatorFactory)
#


#
#
#
DO_PLOT = False  # switch on/off plots
SHADOW3_BINARY = "/users/srio/OASYS1.1/shadow3/shadow3"


import unittest
import numpy
import json
import os

# CODE TO TEST
# from SourceUndulator import SourceUndulator
from orangecontrib.shadow.util.undulator.source_undulator_factory import SourceUndulatorFactory
from orangecontrib.shadow.util.undulator.source_undulator_factory_pysru import SourceUndulatorFactoryPysru
from orangecontrib.shadow.util.undulator.source_undulator_factory_srw import SourceUndulatorFactorySrw
# input/output
from orangecontrib.shadow.util.undulator.source_undulator_input_output import SourceUndulatorInputOutput
#load_file_undul_phot,write_file_undul_phot
#load_file_undul_cdf,write_file_undul_cdf
#plot_undul_cdf,plot_undul_phot

import Shadow

if DO_PLOT:
    try:
        from srxraylib.plot.gol import plot_image,plot, plot_show
    except:
        print("srxraylib not available (for plots). Plots switched off.")
        DO_PLOT = False

#
# Auxiliary functions
#


def _shadow3_commands(commands="exit\n",input_file="shadow3_tmp.inp"):
    # for internal use
    f = open(input_file,'w')
    f.write(commands)
    f.close()
    os.system(SHADOW3_BINARY+" < "+input_file)


def _calculate_shadow3_beam_using_preprocessors(jsn):

    os.system("rm -f start.00 systemfile.dat begin.dat xshundul.plt xshundul.par xshundul.traj xshundul.info xshundul.sha")

    # epath
    commands = "epath\n2\n%f \n%f \n%f \n%d \n1.\nxshundul.par\nxshundul.traj\n1\nxshundul.plt\nexit\n"% \
    (jsn["LAMBDAU"],jsn["K"],jsn["E_ENERGY"],101)
    _shadow3_commands(commands=commands,input_file="shadow3_epath.inp")

    # undul_set
    NG_E = jsn["_NG_E"]
    # TODO: is seems a bug in shadow3: undul_set must be _NG_E>1 (otherwise nosense)
    # but if emin=emaxthe  resulting uphot.nml has _NG_E=1
    if NG_E == 1: NG_E = 2

    commands = "undul_set\n0\n0\n%d \n%d \n%d \nxshundul.traj\n%d\n%f\n%f\n%f\n%f\n0\n1000\nexit\n"% \
        (NG_E,jsn["_NG_T"],jsn["_NG_P"],
         jsn["NPERIODS"],jsn["_EMIN"],jsn["_EMAX"],jsn["INTENSITY"],1e3*jsn["_MAXANGLE"])
    _shadow3_commands(commands=commands,input_file="shadow3_undul_set.inp")

    # undul_phot
    _shadow3_commands(commands="undul_phot\nexit\n",input_file="shadow3_undul_phot.inp")
    _shadow3_commands(commands="undul_phot_dump\nexit\n",input_file="shadow3_undul_phot_dump.inp")

    # undul_cdf
    _shadow3_commands(commands="undul_cdf\n0\n1\nxshundul.sha\nxshundul.info\nexit\n",
                    input_file="shadow3_undul_cdf.inp")


    # input source
    if int(jsn["_FLAG_EMITTANCE(1)"]):
        commands = "input_source\n1\n0\n%d \n%d \n0 \n2 \nxshundul.sha\n%g\n%g\n%g\n%d\n%g\n%d\n%d\n%d\n%d\nexit\n"% \
        (jsn["NRAYS"],jsn["SEED"],jsn["SX"],jsn["SZ"],jsn["EX"],0,jsn["EZ"],0,3,1,1)
    else:
        commands = "input_source\n1\n0\n%d \n%d \n0 \n2 \nxshundul.sha\n%g\n%g\n%g\n%d\n%g\n%d\n%d\n%d\n%d\nexit\n"% \
        (jsn["NRAYS"],jsn["SEED"],0,0,0,0,0,0,3,1,1)

    _shadow3_commands(commands=commands,input_file="shadow3_input_source.inp")

    # run source
    commands = "source\nsystemfile\nexit\n"
    _shadow3_commands(commands=commands,input_file="shadow3_source.inp")


    # return shadow3 beam
    beam = Shadow.Beam()
    beam.load("begin.dat")
    return beam



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
        h["_EMIN"] = 10200.0
        h["_EMAX"] = 10650.0
        h["_NG_E"] = 11
        h["_MAXANGLE"] = 15e-6
        h["_NG_T"] = 51
        h["_NG_P"] = 11

        # internal code
        udict = SourceUndulatorFactory.undul_phot(E_ENERGY = h["E_ENERGY"],INTENSITY = h["INTENSITY"],
                                        LAMBDAU = h["LAMBDAU"],NPERIODS = h["NPERIODS"],K = h["K"],
                                        EMIN = h["_EMIN"],EMAX = h["_EMAX"],NG_E = h["_NG_E"],
                                        MAXANGLE = h["_MAXANGLE"],NG_T = h["_NG_T"],
                                        NG_P = h["_NG_P"])


        photon_energy = numpy.linspace(h["_EMIN"],h["_EMAX"],h["_NG_E"],dtype=float)
        theta = numpy.linspace(0,h["_MAXANGLE"],h["_NG_T"],dtype=float)
        phi = numpy.linspace(0,numpy.pi/2,h["_NG_P"],dtype=float)

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
        self.assertAlmostEqual(diff3,0.00,delta=5e-3)
        self.assertAlmostEqual(diff4,0.00,delta=5e-3)

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
        h["_EMIN"] = 10200.0
        h["_EMAX"] = 10650.0
        h["_NG_E"] = 1
        h["_MAXANGLE"] = 15e-6
        h["_NG_T"] = 51
        h["_NG_P"] = 11

        # internal code
        udict = SourceUndulatorFactory.undul_phot(E_ENERGY = h["E_ENERGY"],INTENSITY = h["INTENSITY"],
                                        LAMBDAU = h["LAMBDAU"],NPERIODS = h["NPERIODS"],K = h["K"],
                                        EMIN = h["_EMIN"],EMAX = h["_EMAX"],NG_E = h["_NG_E"],
                                        MAXANGLE = h["_MAXANGLE"],NG_T = h["_NG_T"],
                                        NG_P = h["_NG_P"])


        photon_energy = numpy.linspace(h["_EMIN"],h["_EMAX"],h["_NG_E"],dtype=float)
        theta = numpy.linspace(0,h["_MAXANGLE"],h["_NG_T"],dtype=float)
        phi = numpy.linspace(0,numpy.pi/2,h["_NG_P"],dtype=float)

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
        self.assertAlmostEqual(diff3,0.00,delta=5e-3)
        self.assertAlmostEqual(diff4,0.00,delta=5e-3)



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
        hh["_EMIN"] = 10200.0
        hh["_EMAX"] = 10650.0
        hh["_NG_E"] = 11
        hh["_MAXANGLE"] = 15e-6
        hh["_NG_T"] = 51
        hh["_NG_P"] = 11

        # internal code
        udict = SourceUndulatorFactoryPysru.undul_phot(E_ENERGY = hh["E_ENERGY"],INTENSITY = hh["INTENSITY"],
                                        LAMBDAU = hh["LAMBDAU"],NPERIODS = hh["NPERIODS"],K = hh["K"],
                                        EMIN = hh["_EMIN"],EMAX = hh["_EMAX"],NG_E = hh["_NG_E"],
                                        MAXANGLE = hh["_MAXANGLE"],NG_T = hh["_NG_T"],
                                        NG_P = hh["_NG_P"])


        photon_energy = numpy.linspace(hh["_EMIN"],hh["_EMAX"],hh["_NG_E"],dtype=float)
        theta = numpy.linspace(0,hh["_MAXANGLE"],hh["_NG_T"],dtype=float)
        phi = numpy.linspace(0,numpy.pi/2,hh["_NG_P"],dtype=float)

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
        self.assertAlmostEqual(diff3,0.00,delta=2e-3)
        self.assertAlmostEqual(diff4,0.00,delta=2e-3)


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

        tmp = \
            """
            {
            "LAMBDAU":     0.0320000015,
            "K":      0.250000000,
            "E_ENERGY":       6.03999996,
            "E_ENERGY_SPREAD":    0.00100000005,
            "NPERIODS": 50,
            "_EMIN":       10200.0000,
            "_EMAX":       10650.0000,
            "INTENSITY":      0.2,
            "_MAXANGLE":     0.000015,
            "_NG_E": 11,
            "_NG_T": 51,
            "_NG_P": 11,
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
            "_FLAG_EMITTANCE(1)":"1",
            "_FLAG_EMITTANCE(2)":"No",
            "_FLAG_EMITTANCE(3)":"Yes",
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

        # u = SourceUndulator()
        # u.load_json_shadowvui_dictionary(h)
        _calculate_shadow3_beam_using_preprocessors(h)



        undul_phot_preprocessor_dict = SourceUndulatorInputOutput.load_file_undul_phot("uphot.dat")

        # if do_plot_intensity: plot_image(undul_phot_preprocessor_dict['radiation'][0,:,:],undul_phot_preprocessor_dict['theta']*1e6,undul_phot_preprocessor_dict['phi']*180/numpy.pi,
        #            title="INTENS UNDUL_PHOT_PREPROCESSOR: RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)
        #
        # if do_plot_polarization: plot_image(undul_phot_preprocessor_dict['polarization'][0,:,:],undul_phot_preprocessor_dict['theta']*1e6,undul_phot_preprocessor_dict['phi']*180/numpy.pi,
        #            title="POL_DEG UNDUL_PHOT_PREPROCESSOR: RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)


        SourceUndulatorInputOutput.plot_undul_phot(undul_phot_preprocessor_dict,title="UNDUL_PHOT_PREPROCESSOR",
                        do_plot_intensity=do_plot_intensity,do_plot_polarization=do_plot_polarization,do_show=False)


        # internal code
        undul_phot_dict = SourceUndulatorFactory.undul_phot(E_ENERGY = h["E_ENERGY"],INTENSITY = h["INTENSITY"],
                                        LAMBDAU = h["LAMBDAU"],NPERIODS = h["NPERIODS"],K = h["K"],
                                        EMIN = h["_EMIN"],EMAX = h["_EMAX"],NG_E = h["_NG_E"],
                                        MAXANGLE = h["_MAXANGLE"],NG_T = h["_NG_T"],
                                        NG_P = h["_NG_P"])

        # if do_plot_intensity: plot_image(undul_phot_dict['radiation'][0,:,:],undul_phot_dict['theta']*1e6,undul_phot_dict['phi']*180/numpy.pi,
        #            title="INTENS UNDUL_PHOT: RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)
        # if do_plot_polarization: plot_image(undul_phot_dict['polarization'][0,:,:],undul_phot_dict['theta']*1e6,undul_phot_dict['phi']*180/numpy.pi,
        #            title="POL_DEG UNDUL_PHOT: RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)

        SourceUndulatorInputOutput.plot_undul_phot(undul_phot_dict,title="UNDUL_PHOT",
                        do_plot_intensity=do_plot_intensity,do_plot_polarization=do_plot_polarization,do_show=False)

        # pySRU
        if is_available_pysru:
            undul_phot_pysru_dict = SourceUndulatorFactoryPysru.undul_phot(E_ENERGY = h["E_ENERGY"],INTENSITY = h["INTENSITY"],
                                            LAMBDAU = h["LAMBDAU"],NPERIODS = h["NPERIODS"],K = h["K"],
                                            EMIN = h["_EMIN"],EMAX = h["_EMAX"],NG_E = h["_NG_E"],
                                            MAXANGLE = h["_MAXANGLE"],NG_T = h["_NG_T"],
                                            NG_P = h["_NG_P"])
            # if do_plot_intensity: plot_image(undul_phot_pysru_dict['radiation'][0,:,:],undul_phot_pysru_dict['theta']*1e6,undul_phot_pysru_dict['phi']*180/numpy.pi,
            #            title="INTENS UNDUL_PHOT_PYSRU: RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)
            # if do_plot_polarization: plot_image(undul_phot_pysru_dict['polarization'][0,:,:],undul_phot_pysru_dict['theta']*1e6,undul_phot_pysru_dict['phi']*180/numpy.pi,
            #            title="POL_DEG UNDUL_PHOT_PYSRU: RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)

            SourceUndulatorInputOutput.plot_undul_phot(undul_phot_pysru_dict,title="UNDUL_PHOT_PYSRU",
                                    do_plot_intensity=do_plot_intensity,do_plot_polarization=do_plot_polarization,do_show=False)

        # srw
        if is_available_srw:
            undul_phot_srw_dict = SourceUndulatorFactorySrw.undul_phot(E_ENERGY = h["E_ENERGY"],INTENSITY = h["INTENSITY"],
                                            LAMBDAU = h["LAMBDAU"],NPERIODS = h["NPERIODS"],K = h["K"],
                                            EMIN = h["_EMIN"],EMAX = h["_EMAX"],NG_E = h["_NG_E"],
                                            MAXANGLE = h["_MAXANGLE"],NG_T = h["_NG_T"],
                                            NG_P = h["_NG_P"])
            # if do_plot_intensity: plot_image(undul_phot_srw_dict['radiation'][0,:,:],undul_phot_srw_dict['theta']*1e6,undul_phot_srw_dict['phi']*180/numpy.pi,
            #            title="INTENS UNDUL_PHOT_SRW: RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)
            # if do_plot_polarization: plot_image(undul_phot_srw_dict['polarization'][0,:,:],undul_phot_srw_dict['theta']*1e6,undul_phot_srw_dict['phi']*180/numpy.pi,
            #            title="POL_DEG UNDUL_PHOT_SRW: RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)

            SourceUndulatorInputOutput.plot_undul_phot(undul_phot_srw_dict,title="UNDUL_PHOT_SRW",
                                    do_plot_intensity=do_plot_intensity,do_plot_polarization=do_plot_polarization,do_show=False)

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
            "_EMIN":       10200.0000,
            "_EMAX":       10650.0000,
            "INTENSITY":      0.2,
            "_MAXANGLE":     0.000015,
            "_NG_E": 11,
            "_NG_T": 51,
            "_NG_P": 11,
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
            "_FLAG_EMITTANCE(1)":"1",
            "_FLAG_EMITTANCE(2)":"No",
            "_FLAG_EMITTANCE(3)":"Yes",
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

        # u = SourceUndulator()
        # u.load_json_shadowvui_dictionary(h)
        _calculate_shadow3_beam_using_preprocessors(h)

        #
        #
        #
        radiation = SourceUndulatorInputOutput.load_file_undul_phot(file_in="uphot.dat")

        cdf2 = SourceUndulatorFactory.undul_cdf(radiation,method='sum')
        SourceUndulatorInputOutput.write_file_undul_cdf(cdf2,file_out="xshundul2.sha")


        cdf3 = SourceUndulatorFactory.undul_cdf(radiation,method='trapz')
        SourceUndulatorInputOutput.write_file_undul_cdf(cdf3,file_out="xshundul3.sha")

        cdf1 = SourceUndulatorInputOutput.load_file_undul_cdf(file_in="xshundul.sha")
        cdf2 = SourceUndulatorInputOutput.load_file_undul_cdf(file_in="xshundul2.sha")
        cdf3 = SourceUndulatorInputOutput.load_file_undul_cdf(file_in="xshundul3.sha")




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
                "_EMIN":       10498.0000,
                "_EMAX":       10499.0000,
                "INTENSITY":      0.200000003,
                "_MAXANGLE":      0.000100,
                "_NG_E": 1,
                "_NG_T": 51,
                "_NG_P": 11,
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
                "_FLAG_EMITTANCE(1)":"0",
                "_FLAG_EMITTANCE(2)":"No",
                "_FLAG_EMITTANCE(3)":"Yes",
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

            # u = SourceUndulator()
            # u.load_json_shadowvui_dictionary(h)
            _calculate_shadow3_beam_using_preprocessors(h)

        elif case == 1:

            # some inputs
            E_ENERGY=6.04
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
            MAXANGLE=100e-6
            NG_T=51
            NG_P=11
            N_J=20
            SEED=36255
            NRAYS=15000

            u = SourceUndulator()

            u.set_from_keywords(
                E_ENERGY = E_ENERGY,
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
                N_J = N_J,
                SEED = SEED,
                NRAYS = NRAYS,
                )

            u.set_energy_monochromatic_at_resonance(harmonic_number=1)
        elif case == 2:

            # u = SourceUndulator()
            # u.load_json_shadowvui_file("xshundul.json")
            # calculate_shadow3_beam_using_preprocessors(json.loads(tmp))

            pass
        else:
            raise Exception("Undefined")

        #
        #
        #
        radiation = SourceUndulatorInputOutput.load_file_undul_phot(file_in="uphot.dat")

        cdf2 = SourceUndulatorFactory.undul_cdf(radiation,method='sum')
        SourceUndulatorInputOutput.write_file_undul_cdf(cdf2,file_out="xshundul2.sha")


        cdf3 = SourceUndulatorFactory.undul_cdf(radiation,method='trapz')
        SourceUndulatorInputOutput.write_file_undul_cdf(cdf3,file_out="xshundul3.sha")

        cdf1 = SourceUndulatorInputOutput.load_file_undul_cdf(file_in="xshundul.sha",)
        cdf2 = SourceUndulatorInputOutput.load_file_undul_cdf(file_in="xshundul2.sha")
        cdf3 = SourceUndulatorInputOutput.load_file_undul_cdf(file_in="xshundul3.sha")


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
