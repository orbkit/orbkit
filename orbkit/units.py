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
sec2atu = _ha / _hbar # Seconds to atomic time units
atu2sec = 1. / sec2atu
                                                                                                                           
a02aa = _a0 / _aa                                                       
aa2a0 = _aa / _a0                                                       
                                                                        
u2me = _u / _me
me2u = _me / _u                                                        
                                                                        
ev2ha = _e / _ha
ha2ev = _ha / _e 

debye2ea0 = _debye / (_e * _a0)
ea02debye = _e * _a0 / _debye

c_in_au = _c*(aa2a0*10**10)/(sec2atu) # Speed of light in atomic units