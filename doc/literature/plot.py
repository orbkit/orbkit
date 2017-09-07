from parse_bib import Bibliography
files = ['orbkit', 'detci']
bibliography = Bibliography()

for fname in files:
  bibliography.read_bib(fname)

bibliography.sort()
bibliography.clean()
bibliography.plot()
bibliography.wirte_csv()
