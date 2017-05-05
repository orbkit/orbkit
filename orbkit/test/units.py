'''
Test orbkit units module
'''
from orbkit.units import *
from orbkit.test.tools import equal

#Reference values for conversion
ref_ha2ev = 27.211386
ref_a02aa = 0.529177210
ref_debye2ae = 0.393456

equal(ha2ev*ev2ha, 1)
equal(a02aa*aa2a0, 1)
equal(debye2ae*ae2debye, 1)

equal(ha2ev, ref_ha2ev)
equal(a02aa, ref_a02aa)
equal(debye2ae, ref_debye2ae)
