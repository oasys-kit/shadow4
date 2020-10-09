"""
Represents a diffraction setup using a shadow3 preprocessor (bragg) file
Except for energy, all units are in SI.
"""

import numpy
import scipy.constants as codata

from crystalpy.diffraction.DiffractionSetupAbstract import DiffractionSetupAbstract
from shadow4.physical_models.bragg.bragg import Bragg
from crystalpy.diffraction.GeometryType import BraggDiffraction

class S4DiffractionSetup(DiffractionSetupAbstract, Bragg):

    def __init__(self, geometry_type=BraggDiffraction, crystal_name="", thickness=1e-6,
                 miller_h=1, miller_k=1, miller_l=1,
                 asymmetry_angle=0.0,
                 azimuthal_angle=0.0,
                 preprocessor_file=""):
        """
        Constructor.
        :param geometry_type: GeometryType (BraggDiffraction,...).
        :param crystal_name: The name of the crystal, e.g. Si.
        :param thickness: The crystal thickness.
        :param miller_h: Miller index H.
        :param miller_k: Miller index K.
        :param miller_l: Miller index L.
        :param asymmetry_angle: The asymmetry angle between surface normal and Bragg normal (radians).
        :param azimuthal_angle: The angle between the projection of the Bragg normal
                                on the crystal surface plane and the x axis (radians).
        """

        DiffractionSetupAbstract.__init__(self, geometry_type=geometry_type, crystal_name=crystal_name,
                                          thickness=thickness,
                                          miller_h=miller_h, miller_k=miller_k, miller_l=miller_l,
                                          asymmetry_angle=asymmetry_angle, azimuthal_angle=azimuthal_angle)

        Bragg.__init__(self, preprocessor_file=preprocessor_file, preprocessor_dictionary=None)
        if preprocessor_file != "":
            self.load_preprocessor_file()

    #
    # implementation of abstract methods uin DiffractionSetupAbstract
    #
    def F0(self, energy, ratio=None):
        """
        Calculate F0 from Zachariasen.
        :param energy: photon energy in eV.
        :return: F0
        """
        F_0, FH, FH_BAR, STRUCT, FA, FB = self.structure_factor(energy, ratio)
        return F_0

    def FH(self, energy, ratio=None):
        """
        Calculate FH from Zachariasen.
        :param energy: photon energy in eV.
        :return: FH
        """
        F_0, FH, FH_BAR, STRUCT, FA, FB = self.structure_factor(energy, ratio)
        return FH

    def FH_bar(self, energy, ratio=None):
        """
        Calculate FH_bar from Zachariasen.
        :param energy: photon energy in eV.
        :return: FH_bar
        """
        F_0, FH, FH_BAR, STRUCT, FA, FB = self.structure_factor(energy, ratio)
        return FH_BAR

    def dSpacing(self):
        """
        Returns the lattice spacing d in A
        :return: Lattice spacing.
        """
        return self._preprocessor_dictionary["dspacing_in_cm"] * 1e8

    def angleBragg(self, energy):
        """
        Returns the Bragg angle for a given energy.
        :param energy: Energy to calculate the Bragg angle for.
        :return: Bragg angle.
        """
        return numpy.arcsin( self.wavelength_in_A(energy) / 2 / self.dSpacing())

    def unitcellVolume(self):
        """
        Returns the unit cell volume.

        :return: Unit cell volume
        """
        # Retrieve unit cell volume from xraylib.
        one_over_volume_times_electron_radius_in_cm = self.get_preprocessor_dictionary()["one_over_volume_times_electron_radius_in_cm"]
        # codata_e2_mc2 = 2.81794032e-15 = Classical electron radius in S.I.
        codata_e2_mc2 = codata.hbar * codata.alpha / codata.m_e / codata.c
        one_over_volume_in_cm = one_over_volume_times_electron_radius_in_cm / (codata_e2_mc2 * 1e2)
        unit_cell_volume = 1.0 / one_over_volume_in_cm * (1e-2)**3
        return unit_cell_volume * (1e10)**3

if __name__ == "__main__":
    from crystalpy.diffraction.GeometryType import BraggDiffraction

    s4diffraction_setup = S4DiffractionSetup(geometry_type          = BraggDiffraction,
                                               crystal_name           = "Si",
                                               thickness              = 100e-6,
                                               miller_h               = 1,
                                               miller_k               = 1,
                                               miller_l               = 1,
                                               asymmetry_angle        = 0.0,
                                               azimuthal_angle        = 0.0,
                                               preprocessor_file="bragg_xop.dat")

    print(s4diffraction_setup.info(8000))