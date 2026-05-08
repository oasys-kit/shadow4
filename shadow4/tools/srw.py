"""
Optional SRW (Synchrotron Radiation Workshop) import helper.

Attempts to import SRW in the following order:

1. ``srwpy`` (OASYS 2.0, SRW 4)
2. ``srwlib`` (SRW 3, standalone)
3. ``oasys_srw`` (OASYS 1.x)

Sets the module-level flag ``SRW_INSTALLED = False`` if none of the above
is available, so callers can guard SRW-dependent code with::

    from shadow4.tools.srw import SRW_INSTALLED
    if SRW_INSTALLED:
        ...
"""

SRW_INSTALLED = True
try:
    try:
        from srwpy.srwlib import *   # OASYS 2.0 - SRW 4
        from srwpy.uti_plot import *
    except:
        from srwlib import *     # SRW 3
        from uti_plot import *
except:
    try:
        from oasys_srw.srwlib import *   # OASYS 1.X
        from oasys_srw.uti_plot import *
    except:
        SRW_INSTALLED = False