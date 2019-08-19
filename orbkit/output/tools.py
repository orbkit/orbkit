import numpy

def get_amira_colormap(colors,filename,opacity=1.0):
  
  f = open(filename,'w')
  f.write('''# AmiraMesh 3D ASCII 2.0


  define Lattice %d

  Parameters {
      ContentType "Colormap",
      MinMax 0 1,
      Interpolate 1,
      LabelField 0
  }

  Lattice { float[4] Data } @1

  # Data section follows
  @1
  '''%len(colors))

  for i,j,k in colors:
    f.write('%.15e %.15e %.15e %.15e\n' % (i,j,k,opacity))

  f.close()


def colormap_creator(rho,filename,n_peaks=5,start=0.01,stop=0.999,peak_width=0.1):
  '''Creates a .cmap colormap for ZIBAmira adjusted to the density.

  **Default:** Isosurface values between 1% and 99.9% of the total density.
  '''
  
  # Where do we have the start and stop percentage? 
  rho_min, rho_max = determine_rho_range(rho,start=start,stop=stop)
  
  # Compute the distance between two isosurface values 
  delta_peak =(rho_max-rho_min)/(n_peaks-1)
  
  # Open a cmap file 
  fid = open('%(f)s.cmap' % {'f': filename}, 'w')
  
  # Write the header 
  fid.write('<!DOCTYPE Colormap>\n')
  fid.write('<ColormapVisage2.0 Name="%(f)s">\n' % {'f': filename})
  fid.write('  <Graph Active="1" Type="0" Name="">\n')
  
  # Initialize a counter for the contour values 
  counter = 0
  
  # Initialize a string for the contours 
  c_str = ('    <Control Opacity="%(o)f" Number="%(c)d" Blue="%(v)f"' + 
            ' Red="%(v)f" Green="%(v)f" Value="%(p)f"/>\n' )
  
  # Write the initial value at zero with zero opacity 
  fid.write(c_str % {'o': 0, 'c': counter, 'v': 0, 'p': 0})
  counter += 1
  
  # Loop over the contour values 
  for ii in range(n_peaks):
    # Calculate the value for the isosurface peak and two values 
    # next to the peak with zero opacity 
    peak = rho_min+delta_peak*(ii+1)
    peak_minus = peak * (1 - peak_width/2.)
    peak_plus = peak * (1 + peak_width/2.)
    
    # Calculate a value for the opacity and the color
    value = 1-(float(ii+1)/(n_peaks+1)*0.9)
    
    # Write the peak 
    # Left edge of the peak 
    fid.write(c_str % {'o': 0, 'c': counter, 'v': value, 'p': peak_minus})
    counter += 1
    # The peak 
    fid.write(c_str % {'o': 1-value, 'c': counter, 'v': value, 'p': peak})
    counter += 1
    # Right edge of the peak 
    fid.write(c_str % {'o': 0, 'c': counter, 'v': value, 'p': peak_plus})
    counter += 1
  
  # Write a final value at 1.5 * (final peak value) with zero opacity 
  fid.write(c_str % {'o': 0, 'c': counter, 'v': 0, 'p': peak_plus*1.5})
  
  # Finalize the file 
  fid.write('  </Graph>\n')
  fid.write('</ColormapVisage2.0>')
  
  # Close the file
  fid.close()

def determine_rho_range(rho,start=0.01,stop=0.999):
  '''Get a range for the isosurface values for the colormap_creator.'''
  # Sort the density values 
  a=numpy.reshape(rho,-1)
  a=a[numpy.argsort(a)]
  
  # Where do we have the start and stop percentage? 
  n1=int(start*len(a))
  n2=int(stop*len(a))
  rho_min=a[n1]
  rho_max=a[n2]  
  return rho_min, rho_max

def colormap_creator_peaks(filename,peaks,peak_width=0.02,peak_minus=None,
  peak_plus=None,alpha=0.2,rgb=0.2):  
  '''Creates a ZIBAmira colomap for selected data values.
  
  **Parameters:**
  
  filename : str
    Specifies the filename of the colormap.
  peaks : list
    Determines the values for peaks in the colormap.
  peak_width : float
    Specifies the width of of the peaks in the colormap.
  peak_min : None or float, optional
    Specifies the lower boundary of the colomap. (Peak with no hight.)
    If None, is set to the smallest data value minus 2*peak_width.
  peak_min : None or float, optional
    Specifies the upper boundary of the colomap. (Peak with no hight.)
    If None, is set to the larges data value plus 2*peak_width.
  alpha : float, data_range={0..1}
    Determines the opacity of the peak. (alpha = 1 - opacity)
  rgb : float or list or numpy.ndarray, data_range={0..1}
    If float or shape=(len(peaks),), specifies a the grey tone.
    Else, specifies the color in [red,green,blue].
  
  '''
  peaks = numpy.sort(peaks)
  if peak_minus is None:
    peak_minus = peaks[0]-(2*peak_width)
  if peak_plus is None:
    peak_plus = peaks[-1]+(2*peak_width)
  
  if isinstance(rgb,(float,int)):
    rgb = numpy.zeros((len(peaks),3)) + rgb
  else:
    rgb = numpy.array(rgb, dtype=float)
    if rgb.ndim == 1:
      if len(rgb) == 3:
        rgb = numpy.zeros((len(peaks),3)) + rgb[numpy.newaxis,:]
      elif len(rgb) == len(peaks):
        rgb = numpy.zeros((len(peaks),3)) + rgb[:,numpy.newaxis]
      else:
        raise ValueError("Wrong shape of 'rgb'")
    elif not (rgb.ndim == 2 and rgb.shape == (len(peaks),3)):
        raise ValueError("Wrong shape of 'rgb'")
  
  # Open a cmap file 
  fid = open('%(f)s.cmap' % {'f': filename}, 'w')
  
  # Write the header 
  fid.write('<!DOCTYPE Colormap>\n')
  fid.write('<ColormapVisage2.0 Name="%(f)s">\n' % {'f': filename})
  fid.write('  <Graph Active="1" Type="0" Name="">\n')
  
  # Initialize a counter for the contour values 
  counter = 0
  
  # Initialize a string for the contours 
  c_str = ('    <Control Opacity="%(o)f" Number="%(c)d" Blue="%(b)f"' + 
            ' Red="%(r)f" Green="%(g)f" Value="%(p)f"/>\n' )
  
  # Write the initial value at zero with zero opacity 
  fid.write(c_str % {'o': 0, 'c': counter, 'r': 0, 'g': 0, 'b': 0, 
                     'p': peak_minus})
  counter += 1
  
  for i,p in enumerate(peaks):
    # Write the peak 
    # Left edge of the peak
    fid.write(c_str % {'o': 0, 'c': counter, 'p': p-peak_width/2., 
                       'r': rgb[i,0], 'g': rgb[i,1], 'b': rgb[i,2]})
    counter += 1
    # The peak 
    fid.write(c_str % {'o': 1-alpha, 'c': counter, 'p': p, 
                       'r': rgb[i,0], 'g': rgb[i,1], 'b': rgb[i,2]})
    counter += 1
    # Right edge of the peak
    fid.write(c_str % {'o': 0, 'c': counter, 'p': p+peak_width/2., 
                       'r': rgb[i,0], 'g': rgb[i,1], 'b': rgb[i,2]})
    counter += 1
  
  # Write a final value at 1.5 * (final peak value) with zero opacity 
  fid.write(c_str % {'o': 0, 'c': counter, 'r': 0, 'g': 0, 'b': 0, 
                     'p': peak_plus})
  
  # Finalize the file 
  fid.write('  </Graph>\n')
  fid.write('</ColormapVisage2.0>')
  
  # Close the file
  fid.close()

def meshgrid2(*arrs):
  '''adapted from:
  http://stackoverflow.com/a/1830192
  '''
  arrs = tuple(reversed(arrs)) 
  lens = list(map(len, arrs))
  dim = len(arrs)
  
  sz = 1
  for s in lens:
      sz*=s
  
  ans = []    
  for i, arr in enumerate(arrs):
    slc = [1]*dim
    slc[i] = lens[i]
    arr2 = numpy.asarray(arr).reshape(slc)
    for j, sz in enumerate(lens):
      if j!=i:
        arr2 = arr2.repeat(sz, axis=j) 
    ans.append(arr2)
  
  return tuple(ans[::-1])
