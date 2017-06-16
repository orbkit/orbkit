import os, inspect
import numpy

from orbkit import multiple_files as mult
from orbkit import options
from orbkit.test.tools import equal

options.quiet = True

tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))
folder = os.path.join(tests_home, 'outputs_for_testing')
filepath = os.path.join(folder, 'NaCl_molden_files.tar.gz')

mult.read(filepath, all_mo=True, nosym=False)


index_list, mo_overlap = mult.order_using_analytical_overlap(None)

ref_index_list = numpy.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51])

ref_mo_overlap = numpy.array([
   4.13745151e-01,   2.15209531e-03,   2.40850414e-01,  -2.07815132e-01,
   3.90468693e-04,   4.71516377e-07,  -6.04865743e-02,  -9.27526663e-02,
   1.62857101e-03,  -1.00590028e-02,   2.73144332e-04,   1.63299615e-03,
  -5.18414926e-12,   1.13089955e-02,   1.47904182e-03,   3.64802085e-03,
   1.27676532e-02,   4.83763397e-03,   4.16172762e-13,   2.08800853e-02,
  -2.39703107e-13,  -7.54349302e-03,   1.33143892e-02,   1.45756792e-02,
   6.67309109e-13,  -6.82625732e-03,  -3.67718123e-03,   3.15176712e-02,
   3.93694019e-13,   2.66813819e-02,   4.84897659e-13,   2.77389964e-02,
  -2.56704784e-03,   1.08086793e-02,   1.31446440e-13,  -7.66053565e-03,
   3.45599160e-02,  -5.70839661e-14,   2.12447532e-02,   7.70203986e-02,
   6.44096989e-03,  -3.66681252e-02,   1.53864579e-13,   3.26282367e-02,
   1.43079111e-14,   1.94622960e-03,  -1.21481229e-01,   6.01625651e-02,
  -5.28532595e-02,  -8.08826853e-02,  -9.41869139e-02,  -1.09207539e-01])

equal(numpy.array(index_list[0][0]), ref_index_list)
equal(numpy.array(mo_overlap[0][0][0]), ref_mo_overlap)
