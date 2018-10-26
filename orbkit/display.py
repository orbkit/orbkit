# -*- coding: iso-8859-1 -*-
'''Module for writing the .oklog files and printing the terminal output.'''
from orbkit import options
from os import remove,path
import platform

is_linux = (platform.system() == 'Linux')
if is_linux:
  import resource

is_initiated = False #: If True, logfile is initialized.
log_fid = None #: Specifies the filename of the oklog file

def init_display(name=None):
  '''Sets the name of the .oklog file and removes the old .oklog file.'''
  global log_fid,is_initiated
  if name is not None or name != '':
    log_fid = '%s.oklog' % name if not name.endswith('.oklog') else name
  if (not options.no_log) or log_fid is not None:
    try:
      remove(log_fid)
    except OSError:
      pass
    if not options.quiet:
      print('Writing log to %s\n' % path.relpath(log_fid))
  is_initiated = True

def display(string):
  '''Prints :literal:`string` to the terminal output and to the .oklog file.'''
  global log_fid
  if not options.quiet:
    print(string)
  if not options.no_log:
    if log_fid is None:
      if options.outputname is None or options.outputname == '':
        return
      else:
        log_fid = '%s.oklog' % options.outputname.split('@')[0]
    if not is_initiated:
      init_display(log_fid)
    f = open(log_fid, 'a')
    f.write('%s\n' % string)
    f.close()

def tForm(string,T,extra=''):
  t_diff = int(round(T))
  tF = {}
  tF['str'] = string
  tF['extra'] = extra
  tF['min'] = t_diff/60
  tF['sec'] = t_diff%60
  tF['h']   = tF['min']/60
  tF['min'] = tF['min']%60
  if tF['h'] == 0:
    if tF['min'] == 0: 
      return ('\n%(str)s took %(sec).3fs%(extra)s.\n' % 
            {'str':string, 'sec': T, 'extra':extra})
    else: 
      return ('\n%(str)s took %(min)dmin and %(sec)ds%(extra)s.\n' % tF)
  else: return ('\n%(str)s took %(h)dh, %(min)dmin and %(sec)ds%(extra)s.\n' % tF)
  # tForm 

def good_bye_message(t):
  if is_linux:
    ram_requirement = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    ram_req = '\nand required %.2f MB of RAM' % (ram_requirement/1000.)
  else:
    ram_req = ''
  msg = tForm('The calculation',t[-1]-t[0],extra=ram_req)
  msg += '\nThank you. Good bye.'
  display(msg)
  return msg