#
# Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
#
import Shadow
import numpy

# write (1) or not (0) SHADOW files start.xx end.xx star.xx
iwrite = 1

#
# initialize shadow3 source (oe0) and beam
#
beam = Shadow.Beam()
oe0 = Shadow.Source()
oe1 = Shadow.OE()

#
# Define variables. See meaning of variables in: 
#  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml 
#  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
#

oe0.FSOUR = 0
oe0.F_PHOT = 0
oe0.HDIV1 = 5e-06
oe0.HDIV2 = 5e-06
oe0.NPOINT = 50000
oe0.PH1 = 1000.0

oe1.DUMMY = 1.0
oe1.FHIT_C = 1
oe1.FMIRR = 1
oe1.F_ANGLE = 1
oe1.RLEN1 = 100.0
oe1.RLEN2 = 100.0
oe1.RWIDX1 = 100.0
oe1.RWIDX2 = 100.0
oe1.T_IMAGE = 300.0
oe1.T_INCIDENCE = 89.714
oe1.T_REFLECTION = 89.714
oe1.T_SOURCE = 1000.0



#Run SHADOW to create the source

if iwrite:
    oe0.write("start.00")

beam.genSource(oe0)

if iwrite:
    oe0.write("end.00")
    beam.write("begin.dat")


#
#run optical element 1
#
print("    Running optical element: %d"%(1))
if iwrite:
    oe1.write("start.01")
beam.traceOE(oe1,1)
if iwrite:
    oe1.write("end.01")
    beam.write("star.01")


Shadow.ShadowTools.plotxy(beam,1,3,nbins=101,nolost=1,title="Real space")
# Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
# Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")
    