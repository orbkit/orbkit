# -*- coding: iso-8859-1 -*-
'''Module for organizing multiprocessing tasks. 
'''
'''
orbkit
Gunter Hermann, Vincent Pohl, Lukas Eugen Marsoner Steinkasserer, Axel Schild, and Jean Christophe Tremblay

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
import sys
import numpy
from time import time
from multiprocessing import Pool

quiet = False
def slicer(N,slice_length=1e4,numproc=1):
  i = 0
  slice_length = 1 if int(slice_length) <= 0.0 else int(slice_length)
  sNum = int((N/(slice_length)))
  xx = []
  if numproc > 1:
    for s in range(sNum):
      if i == N:
        N -= 1
        break
      elif (i + slice_length) >= N:
        xx.append((numpy.array([i,N],dtype=int)))      
      else:
        xx.append((numpy.array([i,i + slice_length],dtype=int)))
      i += slice_length
  else:
    xx.append((numpy.array([0,N],dtype=int))) 
  return xx

def initializer(gargs):
  global global_args
  global_args = gargs

def display(string):
  print(string)

def run(f,x=numpy.arange(10).reshape((-1,1)),numproc=1,display=display,
        initializer=lambda x: x, global_args=None):
  #--- Start the worker processes --
  if numproc > 1:
    pool = Pool(processes=numproc, initializer=initializer, initargs=(global_args,))
    it = pool.imap(f, x)
  else:
    initializer(global_args)
  
  #--- Initialize some additional user information ---
  status_old = 0
  s_old = 0
  t = [time()]
  return_val = [None for i in x]
  try:
    for l,ifid in enumerate(x):
      return_val[l] = it.next() if numproc > 1 else f(ifid)  
      #--- Print out the progress of the computation ---
      status = numpy.floor(l*10/float(len(x)))*10
      if not quiet and display is not None and not status % 10 and status != status_old:
        t.append(time())
        display("\tFinished %(f)d%% (%(s)d slices in %(t).3fs)"
                        % {'f': status,
                            's': l + 1 - s_old,
                            't': t[-1]-t[-2]})
        status_old = status
        s_old = l + 1
  
  finally:
    if numproc > 1:
      pool.close()
      pool.join()
      pool.terminate()
  
  if not quiet and display is not None:
    display("\t" + 40*"-")
    display("\tComputation required %(t).3fs" % {'t': t[-1]-t[0]})
  
  return return_val
