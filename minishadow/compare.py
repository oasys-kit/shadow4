#
# Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
#
import Shadow
import numpy



Shadow.ShadowTools.plotxy("minimirr.01",2,1,nbins=101,nolost=1,title="Mirror (Python)")
Shadow.ShadowTools.plotxy("../idl/mirr.01",2,1,nbins=101,nolost=1,title="Mirror (IDL)")

# Shadow.ShadowTools.plotxy("shrot.01",       4,6,nbins=101,nolost=1,title="Rotated (Python)")
# Shadow.ShadowTools.plotxy("../idl/shrot.01",4,6,nbins=101,nolost=1,title="Mirror (IDL)")

Shadow.ShadowTools.plotxy("ministar.01",1,3,nbins=101,nolost=1,title="Image (Python)")
Shadow.ShadowTools.plotxy("../idl/star.01",1,3,nbins=101,nolost=1,title="Image (IDL)")


