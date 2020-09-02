#
# this test corresponds to the first example in screen_pattern.ows in the tutorial
#


import numpy

from syned.beamline.beamline import BeamlineElement
from syned.beamline.element_coordinates import ElementCoordinates
from shadow4.syned.absorbers.slit import Slit as SySlit                         # TODO: syned.beamline.optical_elements.
from shadow4.syned.absorbers.beam_stopper import BeamStopper as SyBeamStopper   # TODO: syned.beamline.optical_elements.

from shadow4.beam.beam import Beam
from shadow4.optical_elements.screen import Screen

from Shadow.ShadowTools import plotxy
from shadow4.compatibility.beam3 import Beam3

from shadow4.syned.shape import MultiplePatch

def write_grid_pattern(factor=1e-2):
    #
    # creates grid.pol
    #


    patch = MultiplePatch()

    corners = numpy.array([-3.0,-1.5,3,1.5]) # x_leftbottom,y_leftbottom,x_rightup,y_roghtup
    t = numpy.array([9.0,4.5])               # translation vector (i.e., horiz. and V preiods)
    n = numpy.array([2,3])                   # number of translation (H,V)


    file_out = "grid.pol"
    f = open(file_out,'w')
    nn = (2*n[0]+1)*(2*n[1]+1)
    f.write("%d\n"%nn)

    #pay attention that the last element is not included...
    n0 = numpy.arange(-n[0],n[0]+1)
    n1 = numpy.arange(-n[1],n[1]+1)
    for i in n0:
        for j in n1:
            f.write("%d\n"%4)
            f.write(  "%f  %f  \n"%(corners[0]+i*t[0],  corners[1]+j*t[1]) )
            f.write(  "%f  %f  \n"%(corners[0]+i*t[0],  corners[3]+j*t[1]) )
            f.write(  "%f  %f  \n"%(corners[2]+i*t[0],  corners[3]+j*t[1]) )
            f.write(  "%f  %f  \n"%(corners[2]+i*t[0],  corners[1]+j*t[1]) )
            patch.append_polygon(x = factor*numpy.array([corners[0]+i*t[0],corners[0]+i*t[0],corners[2]+i*t[0],corners[2]+i*t[0]]),
                                 y = factor*numpy.array([corners[1]+j*t[1],corners[3]+j*t[1],corners[3]+j*t[1],corners[1]+j*t[1]]))

    f.close()
    print('file %s written to disk.'%file_out)

    return patch


def run_shadow3():
    #
    # Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
    #
    import Shadow
    import numpy

    # write (1) or not (0) SHADOW files start.xx end.xx star.xx
    iwrite = 0

    #
    # initialize shadow3 source (oe0) and beam
    #
    source = Shadow.Beam()
    beam = Shadow.Beam()
    oe0 = Shadow.Source()
    oe1 = Shadow.OE()

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.CONE_MAX = 0.05
    oe0.FDISTR = 5
    oe0.FSOUR = 0
    oe0.IDO_VX = 0
    oe0.IDO_VZ = 0
    oe0.IDO_X_S = 0
    oe0.IDO_Y_S = 0
    oe0.IDO_Z_S = 0
    oe0.NPOINT = 25000
    oe0.PH1 = 1.47642

    oe1.DUMMY = 1.0
    oe1.FILE_SCR_EXT = numpy.array([b'grid.pol', b'', b'', b'', b'', b'', b'', b'', b'', b''])
    oe1.FWRITE = 3
    oe1.F_REFRAC = 2
    oe1.F_SCREEN = 1
    oe1.I_SLIT = numpy.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    oe1.I_STOP = numpy.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    oe1.K_SLIT = numpy.array([2, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    oe1.N_SCREEN = 1
    oe1.T_IMAGE = 5.0
    oe1.T_INCIDENCE = 0.0
    oe1.T_REFLECTION = 180.0
    oe1.T_SOURCE = 322.971

    # Run SHADOW to create the source

    if iwrite:
        oe0.write("start.00")

    source.genSource(oe0)
    beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")

    #
    # run optical element 1
    #
    print("    Running optical element: %d" % (1))
    if iwrite:
        oe1.write("start.01")

    beam.traceOE(oe1, 1)

    if iwrite:
        oe1.write("end.01")
        beam.write("star.01")

    # Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, title="Real space")
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")

    return source, beam

if __name__ == "__main__":

    from srxraylib.plot.gol import set_qt
    from numpy.testing import assert_almost_equal

    set_qt()
    do_plot = True
    do_assert = True

    patches = write_grid_pattern()

    print(patches.info())

    for i,patch in enumerate(patches.get_patches()):
        print(">> patch ", i, patch.info())

    # patches = MultiplePatch()
    # patches.append_polygon([-1,-1,1,1],[-1,1,1,-1])

    source3, beam3 = run_shadow3()

    #
    if do_plot:
        plotxy(beam3,1,3,nbins=100,title="SHADOW3 BEAMSTOPPER", nolost=True)



    beam = Beam.initialize_from_array(source3.rays)


    sy1 = SyBeamStopper(name="Undefined",boundary_shape=patches)   # this is beamstopper
    # sy1 = SySlit(name="Undefined", boundary_shape=patches)         # this is slit (negative)

    coordinates_syned = ElementCoordinates(p=322.971*1e-2, q=5.0*1e-2)

    beamline_element_syned = BeamlineElement(optical_element=sy1, coordinates=coordinates_syned)

    slit1 = Screen(beamline_element_syned=beamline_element_syned)

    print(slit1.info())


    #
    # trace
    #

    beam2 = slit1.trace_beam(beam, flag_lost_value=-101)

    #
    if do_plot:
        beam2_3 = Beam3.initialize_from_shadow4_beam(beam2)
        plotxy(beam2_3,1,3,nbins=100,title="SHADOW4 BEAMSTOPPER", nolost=True)

    print("col#   shadow4  shadow3")
    for i in range(18):
        if i in [0,1,2,12]:
            factor = 1e2
        else:
            factor = 1.0
        print("col%d   %f  %f  " % (i+1, factor * beam2.rays[10,i], beam3.rays[10,i]))
        if do_assert:
            assert_almost_equal (factor * beam2.rays[:,i], beam3.rays[:,i], 4)