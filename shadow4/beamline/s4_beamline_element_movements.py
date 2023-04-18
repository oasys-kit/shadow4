class S4BeamlineElementMovements(object):

    def __init__(self,
                 f_move=0,  # flag  0=Inactive, 1=Active
                 offset_x=0,
                 offset_y=0,
                 offset_z=0,
                 rotation_x=0,  # in rads
                 rotation_y=0,
                 rotation_z=0,
                 ):
        self.f_move = f_move
        self.offset_x = offset_x
        self.offset_y = offset_y
        self.offset_z = offset_z
        self.rotation_x = rotation_x
        self.rotation_y = rotation_y
        self.rotation_z = rotation_z

    def get_offsets(self):
        return self.offset_x, self.offset_y, self.offset_z

    def get_rotations(self):
        return self.rotation_x, self.rotation_y, self.rotation_z
