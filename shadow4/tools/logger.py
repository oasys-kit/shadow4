"""
Logging verbosity helpers for shadow4.

The module wraps Python's standard :mod:`logging` root logger with a simple
two-level interface (verbose / debug) so that the rest of shadow4 does not
need to import :mod:`logging` directly.

Default level is WARNING, meaning neither :func:`is_verbose` nor
:func:`is_debug` returns ``True`` until explicitly enabled.

Example usage::

    from shadow4.tools.logger import set_verbose, printlog

    set_verbose(1)
    printlog("ray tracing started")   # emits at INFO level
"""

import logging

def is_verbose():
    """Return ``True`` if the root logger level is INFO or lower."""
    return logging.root.level <= logging.INFO

def is_debug():
    """Return ``True`` if the root logger level is DEBUG or lower."""
    return logging.root.level <= logging.DEBUG

def set_verbose(status=1):
    """Set the root logger to INFO (verbose) or back to WARNING.

    Parameters
    ----------
    status : int or bool
        Non-zero enables INFO-level output; zero restores WARNING.
    """
    if status:
        logging.getLogger().setLevel(logging.INFO)
    else:
        logging.getLogger().setLevel(logging.WARNING)

def set_debug(status):
    """Set the root logger to DEBUG or back to WARNING.

    Parameters
    ----------
    status : int or bool
        Non-zero enables DEBUG-level output; zero restores WARNING.
    """
    if status:
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        logging.getLogger().setLevel(logging.WARNING)

def printlog(*args):
    """Emit a message at INFO level via the root logger.

    Parameters
    ----------
    *args :
        Passed directly to :func:`logging.info`.
    """
    logging.info(*args)

if __name__ == "__main__":
    from shadow4.tools.logger import is_verbose, is_debug

    print("is_verbose:", is_verbose())
    logging.info("This is not seen")
    set_verbose(1)
    print("is_verbose:", is_verbose())
    logging.info("This is \n seen")
    printlog("via pringlog...")
    print("is_debug; ", is_debug())

