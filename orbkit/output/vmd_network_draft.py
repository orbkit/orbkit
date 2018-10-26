'''Drafts for vmd network files'''

#: Draft for the depiction of a single molecular orbital
mo_string= '''mol new %(n1)s type cube first 0 last -1 step 1 fidbonds 1 autobonds 1 waitfor all
mol delrep 0 top
mol representation CPK 1.000000 0.300000 100.000000 100.000000
mol color Name
mol selection {all}
mol material Opaque
mol addrep top
mol selupdate 0 top 0
mol colupdate 0 top 0
mol scaleminmax top 0 0.000000 0.000000
mol smoothrep top 0 0
mol drawframes top 0 {now}
mol representation Isosurface %(isoblue)0.6f 0 0 0 1 1
mol color ColorID 0
mol selection {all}
mol material Opaque
mol addrep top
mol selupdate 1 top 0
mol colupdate 1 top 0
mol scaleminmax top 1 0.000000 0.000000
mol smoothrep top 1 0
mol drawframes top 1 {now}
mol representation Isosurface %(isored)0.6f 0 0 0 1 1
mol color ColorID 1
mol selection {all}
mol material Opaque
mol addrep top
mol selupdate 2 top 0
mol colupdate 2 top 0
mol scaleminmax top 2 0.000000 0.000000
mol smoothrep top 2 0
mol drawframes top 2 {now}
mol rename top "%(n2)s"

%(mo_options)s

%(render)srender TachyonInternal %(n1)s.tga

molinfo top set drawn 0

# done with molecule %(c)03d

'''

#: Draft for the vmd network
vmd_string = '''#!/usr/local/bin/vmd
display depthcue   off

proc vmdrestoremycolors {} {
color scale colors RWB {1.0 0.0 0.0} {1.0 1.0 1.0} {0.0 0.0 1.0}
color scale colors BWR {0.0 0.0 1.0} {1.0 1.0 1.0} {1.0 0.0 0.0}
color scale colors RGryB {1.0 0.0 0.0} {0.5 0.5 0.5} {0.0 0.0 1.0}
color scale colors BGryR {0.0 0.0 1.0} {0.5 0.5 0.5} {1.0 0.0 0.0}
color scale colors RGB {1.0 0.0 0.0} {0.0 1.0 0.0} {0.0 0.0 1.0}
color scale colors BGR {0.0 0.0 1.0} {0.0 1.0 0.0} {1.0 0.0 0.0}
color scale colors RWG {1.0 0.0 0.0} {1.0 1.0 1.0} {0.0 1.0 0.0}
color scale colors GWR {0.0 1.0 0.0} {1.0 1.0 1.0} {1.0 0.0 0.0}
color scale colors GWB {0.0 1.0 0.0} {1.0 1.0 1.0} {0.0 0.0 1.0}
color scale colors BWG {0.0 0.0 1.0} {1.0 1.0 1.0} {0.0 1.0 0.0}
color scale colors BlkW {0.0 0.0 0.0} {0.5 0.5 0.5} {1.0 1.0 1.0}
color scale colors WBlk {1.0 1.0 1.0} {0.5 0.5 0.5} {0.0 0.0 0.0}
  color scale method RWB
  set colorcmds {
    {color Display {Background} white}
    {color Display {BackgroundTop} black}
    {color Display {BackgroundBot} blue2}
    {color Display {FPS} white}
    {color Name {C} black}
    {color Name {LPA} green}
    {color Name {LPB} green}
    {color Type {LP} green}
    {color Type {DRUD} pink}
    {color Element {X} cyan}
    {color Element {Ac} ochre}
    {color Element {Ag} ochre}
    {color Element {Al} ochre}
    {color Element {Am} ochre}
    {color Element {Ar} ochre}
    {color Element {As} ochre}
    {color Element {At} ochre}
    {color Element {Au} ochre}
    {color Element {B} ochre}
    {color Element {Ba} ochre}
    {color Element {Be} ochre}
    {color Element {Bh} ochre}
    {color Element {Bi} ochre}
    {color Element {Bk} ochre}
    {color Element {Br} ochre}
    {color Element {Ca} ochre}
    {color Element {Cd} ochre}
    {color Element {Ce} ochre}
    {color Element {Cf} ochre}
    {color Element {Cl} ochre}
    {color Element {Cm} ochre}
    {color Element {Co} ochre}
    {color Element {Cr} ochre}
    {color Element {Cs} ochre}
    {color Element {Cu} ochre}
    {color Element {Db} ochre}
    {color Element {Ds} ochre}
    {color Element {Dy} ochre}
    {color Element {Er} ochre}
    {color Element {Es} ochre}
    {color Element {Eu} ochre}
    {color Element {F} ochre}
    {color Element {Fe} ochre}
    {color Element {Fm} ochre}
    {color Element {Fr} ochre}
    {color Element {Ga} ochre}
    {color Element {Gd} ochre}
    {color Element {Ge} ochre}
    {color Element {He} ochre}
    {color Element {Hf} ochre}
    {color Element {Hg} ochre}
    {color Element {Ho} ochre}
    {color Element {Hs} ochre}
    {color Element {I} ochre}
    {color Element {In} ochre}
    {color Element {Ir} ochre}
    {color Element {K} ochre}
    {color Element {Kr} ochre}
    {color Element {La} ochre}
    {color Element {Li} ochre}
    {color Element {Lr} ochre}
    {color Element {Lu} ochre}
    {color Element {Md} ochre}
    {color Element {Mg} ochre}
    {color Element {Mn} ochre}
    {color Element {Mo} ochre}
    {color Element {Mt} ochre}
    {color Element {Na} ochre}
    {color Element {Nb} ochre}
    {color Element {Nd} ochre}
    {color Element {Ne} ochre}
    {color Element {Ni} ochre}
    {color Element {No} ochre}
    {color Element {Np} ochre}
    {color Element {Os} ochre}
    {color Element {Pa} ochre}
    {color Element {Pb} ochre}
    {color Element {Pd} ochre}
    {color Element {Pm} ochre}
    {color Element {Po} ochre}
    {color Element {Pr} ochre}
    {color Element {Pt} ochre}
    {color Element {Pu} ochre}
    {color Element {Ra} ochre}
    {color Element {Rb} ochre}
    {color Element {Re} ochre}
    {color Element {Rf} ochre}
    {color Element {Rg} ochre}
    {color Element {Rh} ochre}
    {color Element {Rn} ochre}
    {color Element {Ru} ochre}
    {color Element {Sb} ochre}
    {color Element {Sc} ochre}
    {color Element {Se} ochre}
    {color Element {Sg} ochre}
    {color Element {Si} ochre}
    {color Element {Sm} ochre}
    {color Element {Sn} ochre}
    {color Element {Sr} ochre}
    {color Element {Ta} ochre}
    {color Element {Tb} ochre}
    {color Element {Tc} ochre}
    {color Element {Te} ochre}
    {color Element {Th} ochre}
    {color Element {Ti} ochre}
    {color Element {Tl} ochre}
    {color Element {Tm} ochre}
    {color Element {U} ochre}
    {color Element {V} ochre}
    {color Element {W} ochre}
    {color Element {Xe} ochre}
    {color Element {Y} ochre}
    {color Element {Yb} ochre}
    {color Element {Zr} ochre}
    {color Resname {} silver}
    {color Chain {X} gray}
    {color Segname {} gray}
    {color Conformation {all} blue}
    {color Molecule {0} blue}
    {color Molecule {1_hf.000_MO_69.1.cb} blue}
    {color Molecule {1} gray}
    {color Structure {3_10_Helix} blue}
    {color Surface {Grasp} gray}
    {color Labels {Atoms} black}
    {color Labels {Bonds} black}
    {color Labels {Springs} orange}
    {color Stage {Even} gray}
    {color Stage {Odd} silver}
  }
  foreach colcmd $colorcmds {
    set val [catch {eval $colcmd}]
  }
  color change rgb 0 0.0 0.0 1.0
  color change rgb 2 0.3499999940395355 0.3499999940395355 0.3499999940395355
  color change rgb 3 1.0 0.5 0.0
  color change rgb 4 1.0 1.0 0.0
  color change rgb 5 0.5 0.5 0.20000000298023224
  color change rgb 6 0.6000000238418579 0.6000000238418579 0.6000000238418579
  color change rgb 7 0.0 1.0 0.0
  color change rgb 9 1.0 0.6000000238418579 0.6000000238418579
  color change rgb 11 0.6499999761581421 0.0 0.6499999761581421
  color change rgb 12 0.5 0.8999999761581421 0.4000000059604645
  color change rgb 13 0.8999999761581421 0.4000000059604645 0.699999988079071
  color change rgb 14 0.5 0.30000001192092896 0.0
  color change rgb 15 0.5 0.5 0.75
  color change rgb 17 0.8799999952316284 0.9700000286102295 0.019999999552965164
  color change rgb 18 0.550000011920929 0.8999999761581421 0.019999999552965164
  color change rgb 19 0.0 0.8999999761581421 0.03999999910593033
  color change rgb 20 0.0 0.8999999761581421 0.5
  color change rgb 21 0.0 0.8799999952316284 1.0
  color change rgb 22 0.0 0.7599999904632568 1.0
  color change rgb 23 0.019999999552965164 0.3799999952316284 0.6700000166893005
  color change rgb 24 0.009999999776482582 0.03999999910593033 0.9300000071525574
  color change rgb 25 0.27000001072883606 0.0 0.9800000190734863
  color change rgb 26 0.44999998807907104 0.0 0.8999999761581421
  color change rgb 27 0.8999999761581421 0.0 0.8999999761581421
  color change rgb 28 1.0 0.0 0.6600000262260437
  color change rgb 29 0.9800000190734863 0.0 0.23000000417232513
  color change rgb 30 0.8100000023841858 0.0 0.0
  color change rgb 31 0.8899999856948853 0.3499999940395355 0.0
  color change rgb 32 0.9599999785423279 0.7200000286102295 0.0
}
vmdrestoremycolors

%(mo)s

'''
