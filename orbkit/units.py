from __futute__ import division
import numpy as np


'''
CODATA recommended values of the fundamental physical constants: 2014
doi: 10.1103/RevModPhys.88.035009
'''

_h = 6.626070040*1e34 #J s
_hbar = h / (2*np.pi)                                                   
_e = 1.6021766208*1e-19 #C                                              
_me = 9.10938356*1e31 #kg                                               
_e0 = 8.854187817*1e-12 #F m**(-1)                                      
_u = 1.660539040*1e-27 #kg                                              
                                                                        
#Derived units                                                          
_a0 = 4*np.pi*_e0*_hbar**2 / (_me*_e**2) #Bohr radius                   
_aa = 1e-10 #Angstrom                                                   
_ha = _hbar**2 / (_me*_a0**2) #Hartree                                      
                                                                        
#Conversion factors                                                     
a02aa = _aa / _a0                                                       
aa2a0 = _a0 / _aa                                                       
                                                                        
u2me = _me / _u                                                         
me2u = _u / _me                                                         
                                                                        
ev2ha = _ha / _e                                                        
ha2ev = _e / _ha                                                         
