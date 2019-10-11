

import numpy
from numpy.testing import assert_almost_equal

from Shadow import Beam as Beam_shadow3
from shadow4.beam.beam import Beam as Beam_shadow4

class Beam3(Beam_shadow3):

    def __init__(self,N):
        super().__init__(N)
        self.beam_shadow4 = None

    @classmethod
    def initialize_from_shadow4_beam(cls,beam):
        if not isinstance(beam,Beam_shadow4):
            raise Exception("Bad input. It must be a shadow4 beam")

        b3 = Beam3(N=beam.get_number_of_rays())

        rays4 = beam.get_rays()

        # rays4[:,10] = beam.get_column(26) * 1e-8

        b3.rays = rays4

        b3.beam_shadow4 = beam

        return b3

    @classmethod
    def load_h5(cls,filename,simulation_name="run001",beam_name="begin"):

        import h5py

        f = h5py.File(filename, 'r')

        column_names = Beam_shadow4.column_names()

        try:
            x = (f["%s/%s/x"%(simulation_name,beam_name)])[:]

            rays = numpy.zeros( (x.size,18))
            for i,column_name in enumerate(column_names):
                rays[:,i] = (f["%s/%s/%s"%(simulation_name,beam_name,column_name)])[:]

        except:
            f.close()
            raise Exception("Cannot find data in %s:/%s/%s"%(filename,simulation_name,beam_name))

        f.close()

        b3 = Beam3(N=x.size)
        b3.rays = rays

        return b3

    @classmethod
    def initialize_from_array(cls,array):
        return Beam3.initialize_from_shadow4_beam(Beam_shadow4.initialize_from_array(array))

    def identical(self,beam2):

        try:
            assert_almost_equal(self.rays,beam2.rays)
            return True
        except:
            return False

    def difference(self,beam2):
        raysnew = beam2.rays
        fact  = 1.0
        for i in range(18):
            m0 = (raysnew[:, i] * fact).mean()
            m1 = self.rays[:, i].mean()
            if numpy.abs(m1) > 1e-10:
                print("\ncol %d, mean: beam_tocheck %g, beam %g , diff/beam %g: " % (i + 1, m0, m1, numpy.abs(m0 - m1) / numpy.abs(m1)))
            else:
                print("\ncol %d, mean: beam_tocheck %g, beam %g " % (i + 1, m0, m1))

            std0 = (raysnew[:, i] * fact).std()
            std1 = self.rays[:, i].std()
            if numpy.abs(std1) > 1e-10:
                print("col %d, std: beam_tocheck %g, beam %g , diff/beam %g" % (i + 1, std0, std1, numpy.abs(std0 - std1) / numpy.abs(std1)))
            else:
                print("col %d, std: beam_tocheck %g, beam %g " % (i + 1, std0, std1))

if __name__ == "__main__":

    from numpy.testing import assert_almost_equal
    import Shadow
    from shadow4.sources.source_geometrical.grid_cartesian import  SourceGridCartesian
    from srxraylib.plot.gol import set_qt

    set_qt()

    source = SourceGridCartesian.initialize_collimated_source(real_space=[10., 0.0, 10.0], real_space_points=[100, 1, 100])

    b4 = source.get_beam()

    print(b4)

    b3 = Beam3.initialize_from_shadow4_beam(b4)


    Shadow.ShadowTools.plotxy(b3,1,3)

    assert(isinstance(b3, Beam_shadow3))

    assert_almost_equal(b3.rays[:,10],b4.rays[:,10])

    assert_almost_equal(b3.getshonecol(11), b4.get_photon_energy_eV(),3)


    #
    # test reader
    #

    a = Beam3.load_h5("../examples/begin.h5")

    Shadow.ShadowTools.plotxy(a,1,3,title="imported from h5 (shadow4)")


