import os
from ase.io import read, write
from pathlib import Path
from ase.calculators.abacus.abacus import Abacus

cs_dir = './CsPbBr3'
cs_cif = Path(cs_dir, 'ICSD_CollCode14610.cif')
cs_atoms = read(cs_cif, format='cif')
#cs_vasp = Path(cs_dir, 'POSCAR')
#cs_atoms = read(cs_vasp, format='vasp')
cs_stru = Path(cs_dir, 'STRU')
pp = {'Cs':'Cs_ONCV_PBE-1.0.upf', 'Br':'Br_ONCV_PBE-1.0.upf', 'Pb':'Pb_ONCV_PBE-1.0.upf'}
basis = {'Cs':'dpsi_Cs.dat', 'Br':'dpsi_Br.dat', 'Pb':'dpsi_Pb.dat'}
#write(cs_stru, cs_atoms, format='abacus', pp=pp, basis=basis)
kpts = [8, 8, 8]
os.chdir(cs_dir)
calc = Abacus(atoms=cs_atoms, ntype=5, ecutwfc=100, scf_nmax=50, smearing_method='gaussian', smearing_sigma=0.001, basis_type='lcao', ks_solver='genelpa', 
            mixing_type='pulay', mixing_beta='0.4', scf_thr=1e-8, vdw_method='d3_bj', out_chg=1, calculation='cell-relax', force_thr=0.001, stress_thr=5,
            cal_force=1, cal_stress=1, out_stru=1, pp=pp, basis=basis, kpts=kpts)
calc.write_input(cs_atoms)
os.chdir('../')