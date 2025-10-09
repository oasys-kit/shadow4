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