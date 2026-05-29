import numpy

def repr_legacy(x):
    """
    Mimics repr() behavior from Python 3.8 / NumPy < 1.14
    where numpy scalars printed like plain Python floats/ints.
    """
    if isinstance(x, numpy.floating):
        return repr(float(x))
    elif isinstance(x, numpy.integer):
        return repr(int(x))
    elif isinstance(x, numpy.complexfloating):
        return repr(complex(x))
    else:
        return repr(x)