from syned.beamline.optical_element import OpticalElement
from shadow4.beamline.s4_optical_element import S4OpticalElement
from shadow4.beamline.s4_beamline_element import S4BeamlineElement
from syned.beamline.element_coordinates import ElementCoordinates


class S4BeamMovement(OpticalElement, S4OpticalElement):
    def __init__(self, name="Undefined",
                 apply_flag=0, # 0=movement off, 1=movement on
                 translation_x=0.0,
                 translation_y=0.0,
                 translation_z=0.0,
                 rotation_x=0.0,
                 rotation_y=0.0,
                 rotation_z=0.0,
                 ):
        super().__init__(name=name)
        self.apply_flag    = apply_flag
        self.translation_x = translation_x
        self.translation_y = translation_y
        self.translation_z = translation_z
        self.rotation_x    = rotation_x
        self.rotation_y    = rotation_y
        self.rotation_z    = rotation_z


class S4BeamMovementElement(S4BeamlineElement):

    def __init__(self, optical_element=None, coordinates=None):
        super().__init__(optical_element if optical_element is not None else S4BeamMovement(),
                         coordinates if coordinates is not None else ElementCoordinates())

    def trace_beam(self,beam_in=None):

        verbose = 0
        #
        beam = beam_in.duplicate()

        #
        # apply movement to beam
        #
        oe = self.get_optical_element()

        if oe.apply_flag:
            if verbose:
                print("Beam translated [%f,%f,%f] um" % (1e6 * oe.translation_x,
                                                         1e6 * oe.translation_y,
                                                         1e6 * oe.translation_z, ))
            beam.translation([oe.translation_x, oe.translation_y, oe.translation_z])
            if oe.rotation_x != 0.0:
                if verbose: print("Beam rotated %f urad along X" % (1e6 * self.rotation_x))
                beam.rotate(oe.rotation_x, axis=1, rad=True)
            if oe.rotation_y != 0.0:
                if verbose: print("Beam rotated %f urad along Y" % (1e6 * oe.rotation_y))
                beam.rotate(oe.rotation_y, axis=2, rad=True)
            if oe.rotation_z != 0.0:
                if verbose: print("Beam rotated %f urad along Z" % (1e6 * oe.rotation_z))
                beam.rotate(oe.rotation_z, axis=3, rad=True)

        #
        # from oe reference system to image plane  (nothing to do...)
        #
        beam_out = beam.duplicate()
        return beam_out, beam

    def to_python_code(self):
        return '\n\n#Error retrieving BeamMovement code\n'

