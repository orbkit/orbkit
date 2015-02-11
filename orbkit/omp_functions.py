import sys
import numpy
from time import time
from multiprocessing import Pool

def slicer(N,vector=1e4,numproc=1):
  i = 0
  vector = 1 if int(vector) <= 0.0 else int(vector)
  sNum = int((N/(vector))+1)
  xx = []
  if numproc > 1:
    for s in range(sNum):
      if i == N:
        N -= 1
        break
      elif (i + vector) >= N:
        xx.append((numpy.array([i,N],dtype=int)))      
      else:
        xx.append((numpy.array([i,i + vector],dtype=int)))
      i += vector
  else:
    xx.append((numpy.array([0,N],dtype=int))) 
  return xx

def run(f,x=numpy.arange(10).reshape((-1,1)),numproc=1,display=sys.stdout.write):
  #--- Start the worker processes --
  if numproc > 1:
    pool = Pool(processes=numproc)
    it = pool.imap(f, x)
  
  #--- Initialize some additional user information ---
  status_old = 0
  s_old = 0
  t = [time()]
  return_val = [None for i in x]
  for l,ifid in enumerate(x):
    return_val[l] = it.next() if numproc > 1 else f(ifid)  
    #--- Print out the progress of the computation ---
    status = numpy.floor(l*10/float(len(x)))*10
    if not status % 10 and status != status_old:
      t.append(time())
      display("\tFinished %(f)d%% (%(s)d slices in %(t).3fs)"
                      % {'f': status,
                          's': l + 1 - s_old,
                          't': t[-1]-t[-2]})
      status_old = status
      s_old = l + 1
  
  if numproc > 1:
    pool.close()
    pool.join()
    pool.terminate()
  
  return return_val