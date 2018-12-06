from __future__ import division, print_function
import numpy
import io
class Bibliography:
  '''Class handling parsing, plotting and printing of ORBKIT bibliographic data for the website'''
  def __init__(self):
    self.unique_entries = None
    self.entries = []
    self.filenames = []
    self.keys = [key for key in self.get_template()]
    self.months = {'jan': 0, 'feb': 1, 'mar': 2,  'apr': 3,
              'may': 4, 'jun': 5, 'jul': 6,  'aug': 7,
              'sep': 8, 'oct': 9, 'nov': 10, 'dec': 11,}
    self.journals = {'Molecular Physics': 'Mol. Phys.',
                'Physical Review A': 'Phys. Rev. A',
                'Physical Review B': 'Phys. Rev. B',
                'The Journal of Physical Chemistry A': 'J. Phys. Chem. A',
                'The Journal of Physical Chemistry C': 'J. Phys. Chem. C',
                'Molecules': 'Molecules',
                'Journal of Computational Chemistry': 'J. Comput. Chem.',
                'Chemical Physics': 'Chem. Phys.',
                'Chemical Physics Letters': 'Chem. Phys. Lett.',
                'Journal of Molecular Graphics and Modelling': 'J. Mol. Graph. Model.',
                'The Journal of Chemical Physics': 'J. Chem. Phys.',
                'Journal of Computer-Aided Molecular Design': 'J. Comput. Aided Mol. Des.',
                'The Journal of Physical Chemistry Letters': 'J. Phys. Chem. Lett.',
                'Physica B: Condensed Matter': 'Physica B',
                'Physical Chemistry Chemical Physics': 'Phys. Chem. Chem. Phys.',
                'International Journal of Quantum Chemistry':   'Int. J. Quantum Chem.',
                'Inorganic Chemistry': 'Inorg. Chem.',
                'Journal of the American Chemical Society': 'J. Am. Chem. Soc.',
                'Journal of Chemical Theory and Computation': 'J. Chem. Theory Comput.',
               }

  def get_template(self):
    return {'year': None,
            'month': None,
            'ym': None,
            'title': None,
            'author': None,
            'volume': None,
            'number': None,
            'journal': None,
            'pages': None,
            'doi': None,
            'url': None
           }

  def read_bib(self, filename):
    '''Reads a bib file and parses the data.
    '''

    with io.open(filename + '.bib', 'r',encoding='ISO-8859-1') as fd:
      entry = []
      for line in fd.readlines():
        if '@article' in line:
          entry.append(self.get_template())
        else:
          splitline = line.split('=')
          for key in self.keys:
            if key in splitline[0]:
              for i in range(10):
                if splitline[1][0] in [',', '{', ' ', '\t']:
                  splitline[1] = splitline[1][1:]
                if splitline[1][-1] in [',', '}', ' ', '\t', '\n']:
                  splitline[1] = splitline[1][:-1]
              entry[-1][key] = splitline[1]
    for i in range(len(entry)):
      name = entry[i]['author'].split('and')[0]
      if ',' in name:
        name = name.split(',')
        name = ' '.join([name[1].strip(),name[0].strip()])
      if len(entry[i]['author'].split('and')) == 1:
        entry[i]['author'] = name + ' '
      else:
        entry[i]['author'] = name + ' *et al.* '
      entry[i]['author'] = entry[i]['author']
      if 'arXiv' in entry[i]['journal']:
        entry[i]['pages'] = entry[i]['journal'].split()[-1]
        entry[i]['journal'] = 'arXiv preprint'
      else:
        entry[i]['journal'] = self.journals[entry[i]['journal']]
      
      entry[i]['ym'] = int(entry[i]['year']) + self.months[entry[i]['month'].lower()] / 12.
    self.entries.append(entry)
    self.filenames.append(filename)

  def sort(self):
    '''Sorts bib entries by year and month of publication and writes the results to a .csv file.
    '''
    for ie, entry in enumerate(self.entries):
      entry_sort = []
      ym = numpy.zeros(len(entry))
      for i, item in enumerate(entry):
        ym[i] = item['ym']

      for i in numpy.argsort(ym):
        entry_sort.append(self.entries[ie][i])
  
      self.entries[ie] = entry_sort

  def url_exists(self, url):
    exists = False
    for entry in self.unique_entries:
      if entry['url'] == url:
        return True

  def clean(self):
    '''Uses url to determine unique citations'''
    self.unique_entries = []
    for entry in self.entries:
      for item in entry:
        if not self.url_exists(item['url']):
          self.unique_entries.append(item)
  
  def wirte_csv(self):
    with io.open('citations.csv', 'w',encoding='ISO-8859-1') as fd:
      for i, entry in enumerate(self.unique_entries):
        if 'arxiv' in entry['journal'].lower():
          print(u'{i},"{author} `{journal} <{url}>`__ {pages} ({year})."'.format(i=i+1,**entry).replace(u'  ', u' '), file=fd)
        elif entry['volume'] is None or entry['pages'] is None:
          print(u'{i},"{author} `{journal} <{url}>`__ {year}."'.format(i=i+1,**entry).replace(u'  ', u' '), file=fd)
        else:
          print(u'{i},"{author} `{journal} <{url}>`__ **{volume}**, {pages} ({year})."'.format(i=i+1,**entry).replace(u'  ', u' '), file=fd)
          
  def plot(self):
    '''Creates a plot from the .csv files containing the bibliographic data.
    '''
    # These are the "Tableau 20" colors as RGB.    
    tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
                 (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
                 (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
                 (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
                 (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    

    # Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
    for i in range(len(tableau20)):    
      r, g, b = tableau20[i]    
      tableau20[i] = (r / 255., g / 255., b / 255.)   

    labels = {'orbkit': 'ORBKIT', 'detci': 'detCI@ORBKIT'}
    import matplotlib.pyplot as plt
    plt.figure(figsize=(8,6))
    xticks = []
    max_y = float('-inf')
    for ie, entry in enumerate(self.entries):
      ydata = []
      xdata = []
      for i, item in enumerate(entry):
        xdata.append(item['ym'])
        ydata.append(i+1)

      xdata = numpy.array(xdata)
      ydata = numpy.array(ydata)

      if numpy.amax(ydata) > max_y:
        max_y = int(round(numpy.amax(ydata),0))
      for x in xdata:
        if int(round(x,0)) not in xticks:
         xticks.append(int(round(x,0)))

      plt.plot(xdata, ydata, linewidth=2.5, color=tableau20[2*ie], label=labels[self.filenames[ie]])
      plt.plot(xdata, ydata,'ro',marker='o', ms=8 , color=tableau20[2*ie])

    #plt.tight_layout()
    plt.rc('font',family='Serif')
    plt.legend(bbox_to_anchor=(1.0, 1.0), loc=1, borderaxespad=0.0)

    dx = int(round((xticks[-1] - xticks[0]) / 6, 0))
    dy = int(round(max_y / 6, 0))
    plt.xticks(xticks, [str(x) for x in xticks], fontsize=14)    
    plt.axis([xticks[0]-.17, xticks[-1] + .17, -dy/5, max_y + 1.2*dy])
    plt.xlabel('date', fontsize=14)
    plt.ylabel('total citations', fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.savefig('citations.png', dpi=100)
    #plt.show()



