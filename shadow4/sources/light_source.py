class LightSource(object):

    def __init__(self):
        pass

    def to_python_code(self, data=None):
        raise NotImplementedError()

    def info(self):
        raise NotImplementedError()

    def get_beam(self, beam=None):
        raise NotImplementedError()