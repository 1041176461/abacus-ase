import os
from ase.io import read, write
from pathlib import Path
from ase.calculators.abacus.abacus import Abacus
from ase.optimize import BFGS

cs_dir = './'
cs_stru = Path(cs_dir, 'STRU')
pp = {'Al':'Al.PD04.PBE.UPF'}
basis = {'Al':'Al_gga_10au_100Ry_3s3p2d.orb'}
cs_atoms = read(cs_stru,format='abacus')
kpts = [4,4,4]
calc = Abacus(atoms=cs_atoms, ntype=1, ecutwfc=20, scf_nmax=50, smearing_method='gaussian', smearing_sigma=0.01, basis_type='lcao', ks_solver='genelpa', 
            mixing_type='pulay', mixing_beta='0.7', scf_thr=1e-8, out_chg=1, calculation='scf', force_thr=0.001, stress_thr=5,
            cal_force=1, cal_stress=1, out_stru=1, pp=pp, basis=basis, kpts=kpts,command='mpirun -n 8 abacus')
calc.write_input(cs_atoms)

opt = BFGS(cs_atoms,trajectory='opt-abacus.traj')
opt.run(fmax=0.05)
