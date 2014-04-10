# -*- coding: iso-8859-1 -*-
'''Module for writing the .oklog files and printing the terminal output.'''
from orbkit import options
from os import remove,path

is_initiated = False #: If True, logfile is initialized.

def init(name=None):
  '''Sets the name of the .oklog file and removes the old .oklog file.'''
  global log_fid,is_initiated
  if name is None:
    log_fid = '%s.oklog' % options.outputname
  else:
    log_fid = name
  if not options.no_log:
    try:
      remove(log_fid)
    except OSError:
      pass
    if not options.quiet:
      print 'Writing log to %s\n' % path.relpath(log_fid)
  is_initiated = True

def display(string):
  '''Prints :literal:`string` to the terminal output and to the .oklog file.'''
  try:
    if not options.quiet:
      print(string)
    if not options.no_log:
      if not is_initiated:
        init(name=options.outputname)
      f = open(log_fid, 'a')
      f.write('%s\n' % string)
      f.close()
  except AttributeError:
    print(string)
