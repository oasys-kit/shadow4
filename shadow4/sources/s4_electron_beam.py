
from syned.storage_ring.electron_beam import ElectronBeam

class S4ElectronBeam(ElectronBeam):
    def __init__(self,
                 energy_in_GeV     = 1.0,
                 energy_spread     = 0.0,
                 current           = 0.1,
                 number_of_bunches = 400,
                 moment_xx         = 0.0,
                 moment_xxp        = 0.0,
                 moment_xpxp       = 0.0,
                 moment_yy         = 0.0,
                 moment_yyp        = 0.0,
                 moment_ypyp       = 0.0):

        super().__init__(
            energy_in_GeV     = energy_in_GeV    ,
            energy_spread     = energy_spread    ,
            current           = current          ,
            number_of_bunches = number_of_bunches,
            moment_xx         = moment_xx        ,
            moment_xxp        = moment_xxp       ,
            moment_xpxp       = moment_xpxp      ,
            moment_yy         = moment_yy        ,
            moment_yyp        = moment_yyp       ,
            moment_ypyp       = moment_ypyp      ,
            )

    def to_python_code(self):

        script = "\n# electron beam"
        script += "\nfrom shadow4.sources.s4_electron_beam import S4ElectronBeam"
        script += "\nelectron_beam = S4ElectronBeam(energy_in_GeV=%g,energy_spread=%g,current=%g)" % \
                  (self.energy(), self._energy_spread, self.current())


        moment_xx, moment_xxp, moment_xpxp = self.get_moments_horizontal()
        moment_yy, moment_yyp, moment_ypyp = self.get_moments_vertical()

        if moment_yyp == 0 and moment_xxp == 0:
            sigma_x, sigma_y, sigma_xp, sigma_yp = self.get_sigmas_all()
            script += "\nelectron_beam.set_sigmas_all(sigma_x=%g,sigma_y=%g,sigma_xp=%g,sigma_yp=%g)" % \
                      (sigma_x, sigma_y, sigma_xp, sigma_yp)
        else:
            script += "\nelectron_beam.set_moments_horizontal(%g,%g,%g)" % (
            moment_xx, moment_xxp, moment_xpxp)
            script += "\nelectron_beam.set_moments_vertical(%g,%g,%g)" % (
            moment_yy, moment_yyp, moment_ypyp)

        return script

if __name__ == "__main__":
    a = S4ElectronBeam()
    print(a.to_python_code())