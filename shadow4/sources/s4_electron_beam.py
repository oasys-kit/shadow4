"""
Defines the shadow4 electron beam.
"""
from syned.storage_ring.electron_beam import ElectronBeam

class S4ElectronBeam(ElectronBeam):
    """
    Defines an electron beam at a given point of the storage ring.

    Parameters
    ----------
    energy_in_GeV : float, optional
        The electron energy in GeV.
    energy_spread : float, optional
        The electron energy spread (in a fraction of the energy_in_GeV).
    current : float, optional
        The electron beam current intensity in A.
    moment_xx : float, optional
        The <x^2> moment.
    moment_xxp : float, optional
        The <x x'> moment.
    moment_xpxp : float, optional
        The <x'^2> moment.
    moment_yy : float, optional
        The <y^2> moment.
    moment_yyp : float, optional
        The <y y'> moment.
    moment_ypyp : float, optional
        The <y'^2> moment.
    dispersion_x : float, optional
        The eta_x parameter, spatial dispersion
    dispersion_y : float, optional
        The eta_y parameter, spatial dispersion
    dispersionp_x : float, optional
        The eta'_x parameter, angular dispersion
    dispersionp_y : float, optional
        The eta'_y parameter, angular dispersion
    """
    def __init__(self,
                 energy_in_GeV     = 1.0,
                 energy_spread     = 0.0,  # not used, but defined for the future.
                 current           = 0.1,
                 number_of_bunches = 0,    # todo: delete
                 moment_xx         = 0.0,
                 moment_xxp        = 0.0,
                 moment_xpxp       = 0.0,
                 moment_yy         = 0.0,
                 moment_yyp        = 0.0,
                 moment_ypyp       = 0.0,
                 dispersion_x      = 0.0,
                 dispersion_y      = 0.0,
                 dispersionp_x     = 0.0,
                 dispersionp_y     = 0.0):

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
            dispersion_x      = dispersion_x,
            dispersion_y      = dispersion_y,
            dispersionp_x     = dispersionp_x,
            dispersionp_y     = dispersionp_y,
            )

    def to_python_code(self):
        """
        Returns the python code to create the electron beam.

        Returns
        -------
        str
            The python code.

        """
        script = "\n# electron beam"
        script += "\nfrom shadow4.sources.s4_electron_beam import S4ElectronBeam"
        script += "\nelectron_beam = S4ElectronBeam(energy_in_GeV=%g,energy_spread=%g,current=%g)" % \
                  (self.energy(), self._energy_spread, self.current())

        moment_xx, moment_xxp, moment_xpxp = self.get_moments_horizontal(dispersion=False)
        moment_yy, moment_yyp, moment_ypyp = self.get_moments_vertical(dispersion=False)

        if moment_yyp == 0 and moment_xxp == 0:
            sigma_x, sigma_xp, sigma_y, sigma_yp = self.get_sigmas_all(dispersion=False)
            script += "\nelectron_beam.set_sigmas_all(sigma_x=%g, sigma_y=%g, sigma_xp=%g, sigma_yp=%g)" % \
                      (sigma_x, sigma_y, sigma_xp, sigma_yp)
        else:
            script += "\nelectron_beam.set_moments_horizontal(%g,%g,%g)" % (moment_xx, moment_xxp, moment_xpxp)
            script += "\nelectron_beam.set_moments_vertical(%g,%g,%g)" % (moment_yy, moment_yyp, moment_ypyp)

        dispersion_x, dispersionp_x, dispersion_y, dispersionp_y = self.get_dispersion_all()

        script += "\nelectron_beam.set_dispersion_all(%g,%g,%g, %g)" % (dispersion_x, dispersionp_x, dispersion_y, dispersionp_y)


        return script

if __name__ == "__main__":
    a = S4ElectronBeam()
    print(a.to_python_code())