'''
Test orbkit units module
'''
from orbkit.units import *
from orbkit.test.tools import equal

#Reference values for conversion
ref_ha2ev = 27.211386
ref_a02aa = 0.529177210
ref_debye2ea0 = 0.393456
ref_u2me = 1822.888486192

equal(ha2ev * ev2ha, 1)
equal(a02aa * aa2a0, 1)
equal(debye2ea0 * ea02debye, 1)
equal(u2me * me2u, 1)

equal(ha2ev, ref_ha2ev)
equal(a02aa, ref_a02aa)
equal(u2me, ref_u2me)

#Debye vaue is not precise enouth - test less stringently
equal(debye2ea0, ref_debye2ea0 , 1e-3)
