# -*- coding: iso-8859-1 -*-
'''Module for writing the .oklog files and printing the terminal output.'''
from orbkit import options
from os import remove,path

is_initiated = False #: If True, logfile is initialized.
log_fid = None #: Specifies the filename of the oklog file

def init_display(name=None):
  '''Sets the name of the .oklog file and removes the old .oklog file.'''
  global log_fid,is_initiated
  if name is not None:
    log_fid = '%s.oklog' % name
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
  global log_fid
  if not options.quiet:
    print(string)
  if not options.no_log:
    if log_fid is None:
      if options.outputname is None:
        return
      else:
        log_fid = '%s.oklog' % options.outputname
    if not is_initiated:
      init_display()
    f = open(log_fid, 'a')
    f.write('%s\n' % string)
    f.close()
