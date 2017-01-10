#
# tests for SourceUndulator
#

import json
import os
import unittest

from numpy.testing import assert_equal, assert_almost_equal
from SourceUndulator import SourceUndulator

import Shadow


#
# switch on/off plots
#
DO_PLOT = False






#
# auxiliary functions
#
def compare_shadow3_files(file1,file2,do_plot=DO_PLOT,do_assert=True):
    print("Comparing shadow3 binary files files: %s %s"%(file1,file2))
    if do_plot:
        Shadow.ShadowTools.plotxy(file1,4,6,nbins=101,nolost=1,title=file1)
        Shadow.ShadowTools.plotxy(file2,4,6,nbins=101,nolost=1,title=file2)
    if do_assert:
        begin1 = Shadow.Beam()
        begin1.load(file1)
        begin2     = Shadow.Beam()
        begin2.load(file2)
        assert_almost_equal(begin1.rays[:,0:6],begin2.rays[:,0:6],3)

def run_from_shadowvui_json_file(method='preprocessor'):

    u = SourceUndulator()
    u.load_json_shadowvui_file("xshundul.json")

    # print(u.to_dictionary())
    print(u.info())

    if method == 'preprocessor':
        # run using binary shadow3 (with preprocessors)
        u.run_using_preprocessors()
    else:
        beam = u.run(code_undul_phot=method,dump_uphot_dot_dat=True,dump_start_files=True)
        beam.write("begin.dat")

#
# Tests
#

class TestSourceUndulator(unittest.TestCase):

    def test_compare_preprocessor_and_internal_from_shadowvui_json_file(self):

        #
        # clean
        #
        os.system("rm begin*.dat")
        os.system("rm uphot*.dat")
        os.system("rm xshundul*.sha")

        #
        # run
        #

        methods = ['preprocessor','internal']

        for method in methods:
            run_from_shadowvui_json_file(method=method)

            os.system("cp begin.dat begin_%s.dat"%method)
            os.system("cp uphot.dat uphot_%s.dat"%method)
            os.system("cp xshundul.sha xshundul_%s.sha"%method)

            # make script
            # oe0 = Shadow.Source()
            # oe0.load("start.00")
            # Shadow.ShadowTools.make_python_script_from_list([oe0],script_file="script_undulator.py")

        compare_shadow3_files("begin_%s.dat"%(methods[0]),"begin_%s.dat"%(methods[1]),do_plot=DO_PLOT,do_assert=True)
