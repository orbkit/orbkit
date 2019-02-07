import numpy


'''
CODATA recommended values of the fundamental physical constants: 2014
doi: 10.1103/RevModPhys.88.035009
'''

pi = numpy.pi

_h = 6.626070040*1e-34 #J s
_hbar = _h / (2*pi)                                                   
_e = 1.6021766208*1e-19 #C                                              
_me = 9.10938356*1e-31 #kg                                               
_e0 = 8.854187817*1e-12 #F m**(-1)                                      
_u = 1.660539040*1e-27 #kg
_c = 299792458. # m s**-1
                                                                        
#Derived units                                                          
_a0 = 4*pi*_e0*_hbar**2 / (_me*_e**2) #Bohr radius                   
_aa = 1e-10 #Angstrom                                                   
_ha = _hbar**2 / (_me*_a0**2) #Hartree
_statC = 1. / 2997924580 #Statcoulomb
_debye = 1e-20*_statC # Debye
                                                                        
#Conversion factors                                                     
a0_to_aa = _a0 / _aa                                                       
aa_to_a0 = _aa / _a0                                                       
                                                                        
u_to_me = _u / _me
me_to_u = _me / _u                                                        
                                                                        
ev_to_ha = _e / _ha
ha_to_ev = _ha / _e 

debye_to_ea0 = _debye / (_e * _a0)
ea0_to_debye = _e * _a0 / _debye

hbarha_to_fs = _ha/_hbar*1e-15
fs_to_hbarha = _hbar/_ha*1e15

ha_to_kJ_mol = 2625.50
