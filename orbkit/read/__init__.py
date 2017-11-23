'''
Build new ORBKIT io interface
'''

from .high_level import main_read
from .molden import read_molden
from .gamess import read_gamess
from .gaussian_fchk import read_gaussian_fchk
from .gaussian_log import read_gaussian_log
from .aomix import read_aomix
from .wfx import read_wfx
from .wfn import read_wfn
from .native import read_native
from .cclib_parser import read_with_cclib
