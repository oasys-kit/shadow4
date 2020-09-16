

class S4OpticalElement(object):

    def __init__(self):
        pass

    def info(self):
        raise NotImplementedError()

    def to_python_code(self, data=None):
        raise NotImplementedError()
