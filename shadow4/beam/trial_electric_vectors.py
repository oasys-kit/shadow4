import numpy

from syned.beamline.beamline import BeamlineElement
from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.optical_elements.mirrors.mirror import Mirror as SyMirror

from shadow4.beam.beam import Beam
from shadow4.sources.source_geometrical.gaussian import SourceGaussian

from shadow4.syned.shape import Rectangle, Ellipse, TwoEllipses # TODO from syned.beamline.shape
from shadow4.syned.shape import Toroidal, Conic, NumericalMesh # TODO from syned.beamline.shape


from shadow4.optical_elements.mirror import Mirror

from Shadow.ShadowTools import plotxy
from shadow4.compatibility.beam3 import Beam3



if __name__ == "__main__":

    """
     Exit from SOURCE
     Call to RESET
     Exit from RESET
     Call to SETSOUR
     Exit from SETSOUR
     Call to IMREF
     Exit from IMREF
     Call to OPTAXIS
     Exit from OPTAXIS
     Call to MSETUP
     Exit from MSETUP
     Call to RESTART
     Exit from RESTART
     Call to MIRROR
     Exit from MIRROR
     Call to IMAGE
     Exit from IMAGE
     Call to DEALLOC
     Exit from DEALLOC
    
    """
    #
    # source
    #
    do_plot=0
    source = SourceGaussian.initialize_from_keywords(number_of_rays=100000,
                                                     sigmaX=0.0,
                                                     sigmaY=0.0,
                                                     sigmaZ=0.0,
                                                     sigmaXprime=0e-6,
                                                     sigmaZprime=0e-6, )
    beam0 = Beam()
    beam0.genSource(source)
    print(beam0.info())

    if do_plot:
        beam0s3 = Beam3.initialize_from_shadow4_beam(beam0)
        plotxy(beam0s3, 4, 6, title="Image 0", nbins=201)

    # #
    # # syned definitopns
    # #
    #
    # # surface shape
    #
    # surface_shape = Conic(conic_coefficients=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0])
    # # surface_shape = Toroidal(min_radius=5.0, maj_radius=2.0) # Plane() # SurfaceShape()
    #
    #
    # # boundaries
    # boundary_shape = None
    #
    # symirror1 = SyMirror(
    #     name="M1",
    #     surface_shape=surface_shape,
    #     boundary_shape=boundary_shape,
    #     coating=None,
    #     coating_thickness=None)
    #
    # coordinates_syned = ElementCoordinates(p=10.0,
    #                                        q=14.0,
    #                                        angle_radial=45 * numpy.pi / 180,
    #                                        angle_azimuthal=  numpy.pi / 2)
    #
    # beamline_element_syned = BeamlineElement(optical_element=symirror1, coordinates=coordinates_syned)
    #
    # #
    # # shadow definitions
    # #
    # mirror1 = Mirror(beamline_element_syned=beamline_element_syned)
    # print(mirror1.info())
    #
    # #
    # # run
    # #
    # beam1, mirr1 = mirror1.trace_beam(beam0)
    # print(mirr1.info())
    #
    # #
    # # check
    # #
    #
    # if do_plot:
    #     beam1s3 = Beam3.initialize_from_shadow4_beam(beam1)
    #     plotxy(beam1s3, 1, 3, title="Image 1", nbins=101, nolost=1)
    #     mirr1s3 = Beam3.initialize_from_shadow4_beam(mirr1)
    #     plotxy(mirr1s3, 2, 1, title="Footprint 1", nbins=101, nolost=1)


    print("Source")
    Esx,Esy,Esz = beam0.get_columns([7,8,9])
    Epx,Epy,Epz = beam0.get_columns([16,17,18])
    print(Esx,Esy,Esz)
    print(Epx,Epy,Epz)

    beam1 = beam0.duplicate()
    beam1.rotate(numpy.pi/2,axis=2,rad=1)
    print("Image")
    Esx,Esy,Esz = beam1.get_columns([7,8,9])
    Epx,Epy,Epz = beam1.get_columns([16,17,18])
    print(Esx,Esy,Esz)
    print(Epx,Epy,Epz)
