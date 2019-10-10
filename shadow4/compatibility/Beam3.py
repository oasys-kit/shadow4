

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

    assert(b3, Beam_shadow3)

    assert_almost_equal(b3.rays[:,10],b4.rays[:,10])

    assert_almost_equal(b3.getshonecol(11), b4.get_photon_energy_eV(),3)


