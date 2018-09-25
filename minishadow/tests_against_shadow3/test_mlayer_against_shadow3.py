import os
import numpy
from numpy.testing import assert_array_almost_equal
from minishadow.mlayer.mlayer import MLayer


__author__ = 'srio'


pre_mlayer_scan_dot_inp = """pre_mlayer_scan
mlayer.dat
501
1
5000
20000
0.75
exit
"""



if __name__ == "__main__":

    from srxraylib.plot.gol import plot

    #
    # shadow3
    #

    f = open("pre_mlayer_scan.inp",'w')
    f.write(pre_mlayer_scan_dot_inp)
    f.close()
    print("File written to disk: pre_mlayer_scan.inp")

    os.system("/Users/srio/OASYS1.1/shadow3/shadow3 < pre_mlayer_scan.inp")

    shadow3_data = numpy.loadtxt("pre_mlayer_scan.dat")

    print(shadow3_data.shape)

    x0 = shadow3_data[:,0].copy()
    y0 = shadow3_data[:,2].copy()

    #
    # python implementation
    #

    a = MLayer()
    a.read_preprocessor_file("/Users/srio/Oasys/mlayer.dat")

    rs, rp, e, t = a.scan(fileOut=None, #"pre_mlayer_scan.dat",
            energyN=501,energy1=5000.0,energy2=20000.0,
            thetaN=1,theta1=0.75,theta2=0.75)

    x1 = e
    y1 = rs[:,0]
    plot(x0,y0,x1,y1,legend=["shadow3","python MLayer"])

    assert_array_almost_equal(x0,x1,3)
    assert_array_almost_equal(y0,y1,3)

