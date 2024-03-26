"""
Container for the variables describing the movements of a S4 beamline element.

"""
class S4BeamlineElementMovements(object):
    """
    Constructor.

    Parameters
    ----------
    f_move : int, optional
        flag to apply movements: 0=No, 1=yes.
    offset_x : float, optional
        The offset in X direction in m.
    offset_y : float, optional
        The offset in Y direction in m.
    offset_z : float, optional
        The offset in Z direction in m.
    rotation_x : float, optional
        The rotation angle around the X axis in rad.
    rotation_y : float, optional
        The rotation angle around the Y axis in rad.
    rotation_z : float, optional
        The rotation angle around the Z axis in rad.
    """
    def __init__(self,
                 f_move     =0,  # flag  0=Inactive, 1=Active
                 offset_x   =0.0,
                 offset_y   =0.0,
                 offset_z   =0.0,
                 rotation_x =0.0,  # in rads
                 rotation_y =0.0,
                 rotation_z =0.0,
                 ):
        self.f_move = f_move
        self.offset_x = offset_x
        self.offset_y = offset_y
        self.offset_z = offset_z
        self.rotation_x = rotation_x
        self.rotation_y = rotation_y
        self.rotation_z = rotation_z

    def get_offsets(self):
        """
        Returns the offsets.

        Returns
        -------
        tuple
            (offset_x, offset_y, offset_z) the offsets along X, Y, and Z.
        """
        return self.offset_x, self.offset_y, self.offset_z

    def get_rotations(self):
        """
        Returns the rotation angles.

        Returns
        -------
        tuple
            (rotation_x, rotation_y, rotation_z) the rotations around X, Y, and Z axes.
        """
        return self.rotation_x, self.rotation_y, self.rotation_z
