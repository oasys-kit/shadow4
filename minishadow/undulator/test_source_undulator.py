__authors__ = ["M Sanchez del Rio - ESRF ISDD Advanced Analysis and Modelling"]
__license__ = "MIT"
__date__ = "12/01/2017"

#
# tests for SourceUndulator
#

import json
import os
import unittest
import numpy

from numpy.testing import assert_almost_equal
from orangecontrib.shadow.util.undulator.source_undulator import SourceUndulator
from orangecontrib.shadow.util.undulator.source_undulator_input_output import SourceUndulatorInputOutput

import Shadow
from srxraylib.plot.gol import plot,plot_image,plot_show

from orangecontrib.shadow.util.undulator.test_source_undulator_factory import _calculate_shadow3_beam_using_preprocessors, _shadow3_commands, SHADOW3_BINARY

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
            for i in [1,3,4,6]:
                print("asserting column: ",i)
                assert_almost_equal(begin1.rays[:,i-1].std(),begin2.rays[:,i-1].std(),3)


    def compare_undul_phot_files(self,file1,file2,do_plot=DO_PLOT,do_assert=True):
        print("Comparing undul_phot output files: %s %s"%(file1,file2))

        dict1 = SourceUndulatorInputOutput.load_file_undul_phot(file_in=file1)
        dict2 = SourceUndulatorInputOutput.load_file_undul_phot(file_in=file2)

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
            "_EMIN":       10498.0000,
            "_EMAX":       10499.0000,
            "INTENSITY":      0.200000003,
            "_MAXANGLE":      0.000100,
            "_NG_E": 101,
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
        # clean
        #
        os.system("rm start.00 begin*.dat uphot*.dat xshundul*.sha")

        #
        # run
        #

        methods = ['preprocessor','internal']

        for method in methods:


            if method == 'preprocessor':
                # run using binary shadow3 (with preprocessors)
                _calculate_shadow3_beam_using_preprocessors(h)
            else:
                from syned.storage_ring.electron_beam import ElectronBeam
                from syned.storage_ring.magnetic_structures.undulator import Undulator


                #
                # syned
                #

                su = Undulator.initialize_as_vertical_undulator(K=h["K"],period_length=h["LAMBDAU"],periods_number=h["NPERIODS"])

                ebeam = ElectronBeam(energy_in_GeV = h["E_ENERGY"],
                             energy_spread         = h["E_ENERGY_SPREAD"],
                             current               = h["INTENSITY"],
                             number_of_bunches     = 1,
                             moment_xx             = (1e-2*h["SX"])**2,
                             moment_xxp            = 0.0,
                             moment_xpxp           = (h["EX"]/h["SX"])**2,
                             moment_yy             = (1e-2*h["SZ"])**2,
                             moment_yyp            = 0.0,
                             moment_ypyp           = (h["EZ"]/h["SZ"])**2,
                                     )

                u = SourceUndulator(name="test",syned_electron_beam=ebeam,syned_undulator=su,
                                flag_emittance=int(h["_FLAG_EMITTANCE(1)"]),flag_size=0,
                                emin=h["_EMIN"],emax=h["_EMAX"],ng_e=h["_NG_E"],
                                maxangle=h["_MAXANGLE"],ng_t=h["_NG_T"],ng_p=h["_NG_P"],
                                code_undul_phot="internal")

                print(u.info())
                # beam = u.calculate_shadow3_beam(user_unit_to_m=1e-2,SEED=36255,NRAYS=h["NRAYS"],)

                rays = u.calculate_rays(user_unit_to_m=1e-2,SEED=36255,NRAYS=h["NRAYS"])
                beam = Shadow.Beam(N=rays.shape[0])
                beam.rays = rays
                beam.write("begin.dat")


            os.system("cp begin.dat begin_%s.dat"%method)
            os.system("cp uphot.dat uphot_%s.dat"%method)


        self.compare_undul_phot_files("uphot_%s.dat"%(methods[0]),"uphot_%s.dat"%(methods[1]),do_plot=DO_PLOT,do_assert=True)
        self.compare_shadow3_files("begin_%s.dat"%(methods[0]),"begin_%s.dat"%(methods[1]),do_plot=DO_PLOT,do_assert=True)




