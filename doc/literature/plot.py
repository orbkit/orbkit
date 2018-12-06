from parse_bib import Bibliography
'''google scholar ClusterIDs

orbkit:  12919164687750559136
detciI: 2491925015538118396
detciII: 17273187972992586419

python2 scholar.py -C 12919164687750559136  --citations-only > orbkit.txt
python2 scholar.py -C  2491925015538118396  --citations-only > detciI.txt
python2 scholar.py -C  17273187972992586419  --citations-only > detciII.txt
cat detci*.txt | grep Title  | sort | uniq > tmp.txt

Enter titles in https://search.crossref.org and download the bibtex files

Do not forget citations before first ORBKIT Paper:
  doi = {10.1080/00268976.2015.1122843},
  doi = {10.1103/PhysRevA.89.052504},
  doi = {10.1021/acs.jpca.5b00907},
  doi = {10.1021/acs.jpcc.5b08606},
  doi = {10.3390/molecules200813830},
  doi = {10.1103/PhysRevA.93.012504},
'''

files = ['orbkit', 'detci']
bibliography = Bibliography()

for fname in files:
  bibliography.read_bib(fname)

bibliography.sort()
bibliography.clean()
bibliography.plot()
bibliography.wirte_csv()
