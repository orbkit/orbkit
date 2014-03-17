# -*- coding: iso-8859-1 -*-

'''
orbkit
Gunter Hermann, Vincent Pohl, and Axel Schild

Institut fuer Chemie und Biochemie, Freie Universitaet Berlin, 14195 Berlin, Germany

This file is part of orbkit.

orbkit is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or any later version.

orbkit is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with orbkit.  If not, see <http://www.gnu.org/licenses/>.
'''

import numpy
from orbkit import grid

def integrate_sector(p0,p1,rho):
  # --- FUNCTION integrate_sector --- integrate rho from p0 to p1 --------------------------------
  # r  = numpy.sqrt(x**2+y**2)         # vector of all radii
  # x1 = r * numpy.cos(p1) # x-values of the limits
  # y1 = r * numpy.sin(p1) # y-values of the limits
  I_all = numpy.zeros((len(p1)));
  rho = numpy.sum(rho,2)
    
  pt = numpy.arctan2(grid.x,grid.y)
  for kk in range(len(p1)):
    I = 0.;
    for ii in range(len(grid.x)):
      if (pt[ii] <= p1[kk]) and (pt[ii] > p0):
        for jj in range(len(grid.y)):
           I = I + rho[ii,jj]
    I_all[kk]= I * grid.d3r
  
  return I_all

