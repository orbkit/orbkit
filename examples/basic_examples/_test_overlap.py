from orbkit import read,tools, analytical_integrals
from orbkit import grid,cy_grid,cy_core,cy_overlap
from orbkit.tools import *

drv=2
qc = read.main_read('h2o.molden',itype='molden',all_mo=True)
a=analytical_integrals.get_ao_overlap2(qc.geo_spec,qc.geo_spec,qc.ao_spec)
b=analytical_integrals.get_ao_overlap(qc.geo_spec,qc.geo_spec,qc.ao_spec)
#geo_spec = qc.geo_spec
#ao_spec = qc.ao_spec
#coord_a=geo_spec
#coord_b=geo_spec
#lxlylz_a = ao_spec.get_lxlylz()
#lxlylz_b=lxlylz_a
#ra = numpy.zeros((0,3))
#rb = numpy.array(ra, copy=True)
#for ao in ao_spec:
  #ra = numpy.append(ra,coord_a[ao['atom']][numpy.newaxis,:],axis=0)
  #rb = numpy.append(rb,coord_b[ao['atom']][numpy.newaxis,:],axis=0)
#coeffs = ao_spec.get_lmpao()
#index = ao_spec.get_lmprim2cont(return_l=True)

#ra = require(ra,dtype='f')
#rb = require(rb,dtype='f')

##lxlylz_a comes from an AOClass instance so its type is checked there
##lxlylz_b might come from anywhere so we neet to make sure it has the right type
#lxlylz_b = require(lxlylz_b,dtype='i')

#a = cy_overlap.aooverlap2(ra,rb,lxlylz_a,lxlylz_b,coeffs,index,drv,int(ao_spec.normalized))

##geo_spec = require(geo_spec, dtype='f')
##lxlylz,assign = ao_spec.get_lxlylz(get_assign=True,bincount=True)
##ao_coeffs,pnum_list,atom_indices = prepare_ao_calc(ao_spec)
##is_normalized = each_ao_is_normalized(ao_spec)
##lxlylz = require(lxlylz,dtype='i')
##assign = require(assign,dtype='i')
##b=cy_overlap.aooverlap(geo_spec,geo_spec,lxlylz,lxlylz,assign,ao_coeffs,pnum_list,atom_indices,drv,0)


assert all(a==b)

#c = 0
#c_p=c_ao=0
#c_i=0
#for i in range(assign.shape[0]):
      #print c_ao,c_p,c_i
      #for i_ao in range(assign[i]):   
            #for p_ao in range(pnum_list[i]):  
              ##print c,sum(assign[:i]*pnum_list[:i])+i_ao+p_ao, c_i+i_ao+p_ao
              ##print sum(assign[:i])+i_ao, sum(pnum_list[:i])+p_ao,c_p+p_ao,c_ao+i_ao
              #c+=1
      #c_ao += assign[i]
      #c_p += pnum_list[i]
      #c_i += assign[i]*pnum_list[i]
      
    
