'''
Test orbkit units module
'''
from orbkit.units import *
from orbkit.test.tools import equal

#Reference values for conversion
ref_ha_to_ev = 27.211386
ref_a0_to_aa = 0.529177210
ref_debye_to_ea0 = 0.393456
ref_u_to_me = 1822.888486192

equal(ha_to_ev * ev_to_ha, 1)
equal(a0_to_aa * aa_to_a0, 1)
equal(debye_to_ea0 * ea0_to_debye, 1)
equal(u_to_me * me_to_u, 1)

equal(ha_to_ev, ref_ha_to_ev)
equal(a0_to_aa, ref_a0_to_aa)
equal(u_to_me, ref_u_to_me)

#Debye vaue is not precise enouth - test less stringently
equal(debye_to_ea0, ref_debye_to_ea0 , 1e-3)
