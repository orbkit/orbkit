from . import AOIntegrals, FCIDUMP
from . import symmetry

def generate_fcidump(qc, filename='FCIDUMP', core=0, occ=None, nosymm=False, max_dims=0):
  '''Creates a FCIDUMP file from a QCInfo object.

    **Parameters:**

      filename : string
      core, occ : int or list of ints
        active space (core and occupied) orbitals per IRREP
      nosymm : boolean
        Do not exploit symmetry.
      max_dims : int
        If > 0, calculate AO Integrals in slices containing no more than max_dims AOs (but at least one shell). Slower, but requires less memory.
  '''

  ao = AOIntegrals(qc)
  Nelec = sum([mo['occ_num'] for mo in qc.mo_spec])

  if not nosymm:
    pass
  else:
    if occ is None:
      occ = ao.Norb
    nmopi = [occ]
    MOranges = [range(occ)]

  fcidump = FCIDUMP(Norb=occ, Nelec=Nelec, spin=1)
  fcidump.set_OrbSym(nmopi)
  fcidump.nuclear_repulsion = qc.nuclear_repulsion

  # Hcore
  for ii, MOrange in enumerate(MOranges):
    ao.add_MO_block_1e(MOrange, MOrange)

  blocks = ao.Hcore(asMO=1, max_dims=max_dims)

  if len(MOranges) == 1:
    fcidump.H = blocks
  else:
    for ii, block in enumerate(blocks):
      fcidump.set_H(block, irrep=ii+1)

  # ERI
  ind = []
  for i in range(len(nmopi)):
    for j in range(i, len(nmopi)):       # hermitian
      for k in range(i, len(nmopi)):     # exchange of electronic coordinates
        for l in range(k, len(nmopi)):   # hermitian

          _i, _j, _k, _l = i+1, j+1, k+1, l+1
          if symmetry.irrep_mul(_i, _j, _k, _l) != 1:
            continue

          ao.add_MO_block_2e(
            MOrangei=MOranges[i],
            MOrangej=MOranges[j],
            MOrangek=MOranges[k],
            MOrangel=MOranges[l],
          )
          ind.append((_i, _j, _k, _l))

  blocks = ao.Vee(asMO=1, max_dims=max_dims)

  if len(MOranges) == 1:
    fcidump.G = blocks
  else:
    for i, block in enumerate(blocks):
      i, j, k, l = ind[i]
      fcidump.set_G(block, i, j, k, l)

  # incorporate core orbitals
  fcidump.reduce_active_space(core=core)

  # write to file
  fcidump.store(filename)
