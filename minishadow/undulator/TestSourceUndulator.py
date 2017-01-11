#
# tests for SourceUndulator
#

import json
import os
import unittest
import numpy

from numpy.testing import assert_almost_equal
from SourceUndulator import SourceUndulator
from SourceUndulatorInputOutput import load_file_undul_phot

import Shadow
from srxraylib.plot.gol import plot,plot_image,plot_show

#
# switch on/off plots
#
DO_PLOT = False

#
# Tests
#

class TestSourceUndulator(unittest.TestCase):
    #
    # auxiliary functions
    #
    def compare_shadow3_files(self,file1,file2,do_plot=DO_PLOT,do_assert=True):
        print("Comparing shadow3 binary files: %s %s"%(file1,file2))

        if do_plot:
            Shadow.ShadowTools.plotxy(file1,4,6,nbins=101,nolost=1,title=file1)
            Shadow.ShadowTools.plotxy(file2,4,6,nbins=101,nolost=1,title=file2)

        if do_assert:
            begin1 = Shadow.Beam()
            begin1.load(file1)
            begin2     = Shadow.Beam()
            begin2.load(file2)
            assert_almost_equal(begin1.rays[:,0:6],begin2.rays[:,0:6],3)

    def compare_undul_phot_files(self,file1,file2,do_plot=DO_PLOT,do_assert=True):
        print("Comparing undul_phot output files: %s %s"%(file1,file2))

        dict1 = load_file_undul_phot(file_in=file1)
        dict2 = load_file_undul_phot(file_in=file2)

        rad1 = dict1["radiation"]
        # Do not compare polarizartion, I believe the preprocessor one is wrong
        # pol1 = dict1["polarization"]
        e1   = dict1["photon_energy"]
        t1   = dict1["theta"]
        p1   = dict1["phi"]

        rad2 = dict2["radiation"]
        # pol2 = dict2["polarization"]
        e2   = dict2["photon_energy"]
        t2   = dict2["theta"]
        p2   = dict2["phi"]

        print(r"---> Max diff E array %f "%( (e2-e1).max() ))
        print(r"---> Max diff T array %f "%( (t2-t1).max() ))
        print(r"---> Max diff P array %f "%( (p2-p1).max() ))


        rad_max = numpy.max( (rad1,rad2) )
        diff_max = numpy.max( (rad1-rad2) )

        print(r"---> diff_rad_max/rad_max = %f %%"%(100*diff_max/rad_max))

        if do_plot:
            plot_image(dict1['radiation'][0,:,:],dict1['theta']*1e6,dict1['phi']*180/numpy.pi,
                       title="INTENS UNDUL_PHOT_PREPROCESSOR: RN0[0]"+file1,xtitle="Theta [urad]",ytitle="Phi [deg]",
                       aspect='auto',show=False)

            plot_image(dict2['radiation'][0,:,:],dict1['theta']*1e6,dict1['phi']*180/numpy.pi,
                       title="INTENS UNDUL_PHOT_PREPROCESSOR: RN0[0]"+file2,xtitle="Theta [urad]",ytitle="Phi [deg]",
                       aspect='auto',show=False)

            plot_show()

        if do_assert:
            assert_almost_equal(e1,e2)
            assert_almost_equal(t1,t2)
            assert_almost_equal(p1,p2)
            # compare only points with appreciable intensity
            # accept if differences are less that 15%
            for ie,e in enumerate(e1):
                for it,t in enumerate(t1):
                    for ip,p in enumerate(p1):
                        if rad1[ie,it,ip] > 0.1*rad_max:
                            mydiff =  100*numpy.abs(rad1[ie,it,ip]-rad2[ie,it,ip])/rad1[ie,it,ip]
                            print(r"--> intensity first:%g second:%g  diff:%g %%"%(rad1[ie,it,ip],rad2[ie,it,ip],mydiff))
                            self.assertLess( mydiff, 15. )


    def test_compare_preprocessor_and_internal_from_shadowvui_json_file(self,shadowvui_json_file=None):

        if shadowvui_json_file == None:
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
            "NG_E": 101,
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
            shadowvui_json_file = "tmp.json"
            f = open(shadowvui_json_file,'w')
            f.write(tmp)
            f.close()
            print("File %s written to disk."%shadowvui_json_file)



        #
        # clean
        #
        os.system("rm start.00 begin*.dat uphot*.dat xshundul*.sha")

        #
        # run
        #

        methods = ['preprocessor','internal']

        for method in methods:
            u = SourceUndulator()
            u.load_json_shadowvui_file(shadowvui_json_file)
            print(u.info())

            if method == 'preprocessor':
                # run using binary shadow3 (with preprocessors)
                u.run_using_preprocessors()
            else:
                beam = u.run(code_undul_phot=method,dump_uphot_dot_dat=True,dump_start_files=True)
                beam.write("begin.dat")

            os.system("cp begin.dat begin_%s.dat"%method)
            os.system("cp uphot.dat uphot_%s.dat"%method)
            os.system("cp xshundul.sha xshundul_%s.sha"%method)

            # make script
            # oe0 = Shadow.Source()
            # oe0.load("start.00")
            # Shadow.ShadowTools.make_python_script_from_list([oe0],script_file="script_undulator.py")


        self.compare_undul_phot_files("uphot_%s.dat"%(methods[0]),"uphot_%s.dat"%(methods[1]),do_plot=DO_PLOT,do_assert=True)
        self.compare_shadow3_files("begin_%s.dat"%(methods[0]),"begin_%s.dat"%(methods[1]),do_plot=DO_PLOT,do_assert=True)


    def test_setters_and_getters(self):

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

        self.assertEqual( u.E_ENERGY, E_ENERGY)
        self.assertEqual( u.E_ENERGY_SPREAD, E_ENERGY_SPREAD)
        self.assertEqual( u.INTENSITY, INTENSITY)
        self.assertEqual( u.SX, SX)
        self.assertEqual( u.SZ, SZ)
        self.assertEqual( u.SXP, SXP)
        self.assertEqual( u.SZP, SZP)
        self.assertEqual( u.FLAG_EMITTANCE, FLAG_EMITTANCE)
        self.assertEqual( u.LAMBDAU, LAMBDAU)
        self.assertEqual( u.NPERIODS, NPERIODS)
        self.assertEqual( u.K, K)
        self.assertEqual( u.EMIN, EMIN)
        self.assertEqual( u.EMAX, EMAX)
        self.assertEqual( u.NG_E, NG_E)
        self.assertEqual( u.MAXANGLE, MAXANGLE)
        self.assertEqual( u.NG_T, NG_T)
        self.assertEqual( u.NG_P, NG_P)
        self.assertEqual( u.SEED, SEED)
        self.assertEqual( u.NRAYS, NRAYS)

        tmp = u.to_dictionary()
        self.assertEqual(tmp["E_ENERGY"] , E_ENERGY)
        self.assertEqual(tmp["E_ENERGY_SPREAD"] , E_ENERGY_SPREAD)
        self.assertEqual(tmp["INTENSITY"] , INTENSITY)
        self.assertEqual(tmp["SX"] , SX)
        self.assertEqual(tmp["SZ"] , SZ)
        self.assertEqual(tmp["SXP"] , SXP)
        self.assertEqual(tmp["SZP"] , SZP)
        self.assertEqual(tmp["FLAG_EMITTANCE"] , FLAG_EMITTANCE)
        self.assertEqual(tmp["LAMBDAU"] , LAMBDAU)
        self.assertEqual(tmp["NPERIODS"] , NPERIODS)
        self.assertEqual(tmp["K"] , K)
        self.assertEqual(tmp["EMIN"] , EMIN)
        self.assertEqual(tmp["EMAX"] , EMAX)
        self.assertEqual(tmp["NG_E"] , NG_E)
        self.assertEqual(tmp["MAXANGLE"] , MAXANGLE)
        self.assertEqual(tmp["NG_T"] , NG_T)
        self.assertEqual(tmp["NG_P"] , NG_P)
        self.assertEqual(tmp["SEED"] , SEED)
        self.assertEqual(tmp["NRAYS"] , NRAYS)

        for key in tmp:
            tmp_old = tmp[key]
            tmp_new = tmp_old * 5
            tmp[key] = tmp_new
            #print("<><> %s changed from %f to %f"%(key,tmp_old,tmp[key]))


        u.set_from_dictionary(tmp)
        u.info()

        tmp_new = u.to_dictionary()
        self.assertEqual(tmp_new["E_ENERGY"] , E_ENERGY * 5)
        self.assertEqual(tmp_new["E_ENERGY_SPREAD"] , E_ENERGY_SPREAD * 5)
        self.assertEqual(tmp_new["INTENSITY"] , INTENSITY * 5)
        self.assertEqual(tmp_new["SX"] , SX * 5)
        self.assertEqual(tmp_new["SZ"] , SZ * 5)
        self.assertEqual(tmp_new["SXP"] , SXP * 5)
        self.assertEqual(tmp_new["SZP"] , SZP * 5)
        self.assertEqual(tmp_new["FLAG_EMITTANCE"] , FLAG_EMITTANCE * 5)
        self.assertEqual(tmp_new["LAMBDAU"] , LAMBDAU * 5)
        self.assertEqual(tmp_new["NPERIODS"] , NPERIODS * 5)
        self.assertEqual(tmp_new["K"] , K * 5)
        self.assertEqual(tmp_new["EMIN"] , EMIN * 5)
        self.assertEqual(tmp_new["EMAX"] , EMAX * 5)
        self.assertEqual(tmp_new["NG_E"] , NG_E * 5)
        self.assertEqual(tmp_new["MAXANGLE"] , MAXANGLE * 5)
        self.assertEqual(tmp_new["NG_T"] , NG_T * 5)
        self.assertEqual(tmp_new["NG_P"] , NG_P * 5)
        self.assertEqual(tmp_new["SEED"] , SEED * 5)
        self.assertEqual(tmp_new["NRAYS"] , NRAYS * 5)


    def test_file_dump(self):

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

        # dump file
        u.write(file_out='startj.00')

        uDict = u.to_dictionary()

        with open('startj.00') as data_file:
            data = json.load(data_file)

        for key in data:
            self.assertEqual(data[key],uDict[key])
            data[key] *= 5

        # read file

        with open('tmp.json', 'w') as outfile:
            json.dump(data, outfile, indent=4, sort_keys=True, separators=(',', ':'))

        u.load('tmp.json')

        self.assertEqual( u.E_ENERGY , data["E_ENERGY"])
        self.assertEqual( u.E_ENERGY_SPREAD , data["E_ENERGY_SPREAD"])
        self.assertEqual( u.INTENSITY , data["INTENSITY"])
        self.assertEqual( u.SX , data["SX"])
        self.assertEqual( u.SZ , data["SZ"])
        self.assertEqual( u.SXP , data["SXP"])
        self.assertEqual( u.SZP , data["SZP"])
        self.assertEqual( u.FLAG_EMITTANCE , data["FLAG_EMITTANCE"])
        self.assertEqual( u.LAMBDAU , data["LAMBDAU"])
        self.assertEqual( u.NPERIODS , data["NPERIODS"])
        self.assertEqual( u.K , data["K"])
        self.assertEqual( u.EMIN , data["EMIN"])
        self.assertEqual( u.EMAX , data["EMAX"])
        self.assertEqual( u.NG_E , data["NG_E"])
        self.assertEqual( u.MAXANGLE , data["MAXANGLE"])
        self.assertEqual( u.NG_T , data["NG_T"])
        self.assertEqual( u.NG_P , data["NG_P"])
        self.assertEqual( u.SEED , data["SEED"])
        self.assertEqual( u.NRAYS , data["NRAYS"])

    def do_info(self):

        # some inputs
        E_ENERGY=6.04
        E_ENERGY_SPREAD=0.001
        INTENSITY=0.2
        SX=0.04
        SZ=0.001
        SXP=10e-6
        SZP=4e-6
        FLAG_EMITTANCE=0
        LAMBDAU=0.032
        NPERIODS=50
        K=0.25
        EMIN= 10498.0000
        EMAX= 10499.0000
        NG_E=11
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

        # print(u.info(debug=True))

        # u.set_energy_monochromatic_at_resonance(harmonic_number=1)
        # # print(u.info)
        # # u.EMAX = u.EMIN  + 0.01
        # # u.NG_E = 1
        #
        # print(u.info(debug=True))
        # beam = u.run(code_undul_phot='internal',dump_uphot_dot_dat=True,dump_start_files=True)
        # beam.write("begin.dat")

