
from shadow4.sources.source_geometrical.grid_cartesian import SourceGridCartesian
from shadow4.compatibility.beam3 import Beam3
import unittest

DO_PLOT = True

# TODO: add assert

class TestGridCartesian(unittest.TestCase):
    def do_point_source(self):

        nx = 20
        ny = 30
        a = SourceGridCartesian.initialize_point_source(
                    direction_space        = [2e-3,2e-3],
                    direction_space_points = [nx, ny],
                    direction_space_center = [0.0, 0.0] )
        print(a.info())

        x,y,z = a.get_arrays_real_space()
        print("x:",x.shape)
        print("y:",y.shape)
        print("z:",z.shape)
        self.assertEqual(1, x.size*y.size*z.size)

        xp,zp = a.get_arrays_direction_space()
        # print("xp:",xp)
        # print("zp:",zp)

        XP,ZP = a.get_mesh_divergences()
        print("XP.shape ZP.shape",XP.shape,ZP.shape)
        self.assertEqual(nx*ny, XP.size)
        self.assertEqual(nx*ny, ZP.size)

        VP = a.get_volume_divergences()
        print("VP",VP.shape,VP.size)
        self.assertEqual(3*nx*ny, VP.size)

        Vx = a.get_volume_real_space()
        print("Vx: ",Vx.shape)

        V = a.get_volume()
        print("V: ",V.shape)


        beam_shadow3 = a.get_beam_shadow3()
        beam_shadow3.write("begin.dat")

        if DO_PLOT:
            import Shadow
            Shadow.ShadowTools.plotxy(beam_shadow3,4,6)


    def test_collimated_source(self):

        #
        #
        #
        a = SourceGridCartesian.initialize_collimated_source(real_space=[1.,0.0,1.0],real_space_points=[100,1,100])
        print(a.info())
        beam_shadow3 = Beam3.initialize_from_shadow4_beam(a.get_beam())


        beam = a.get_beam()

        if DO_PLOT:

            from srxraylib.plot.gol import plot_scatter, set_qt
            set_qt()
            plot_scatter(beam.get_column(1),beam.get_column(3))
            import Shadow
            Shadow.ShadowTools.plotxy(beam_shadow3,1,3)

        print(beam.info())