import os
from ase.io import read, write
from pathlib import Path
from ase.calculators.abacus.abacus import Abacus
from ase.calculators.espresso import Espresso


cs_dir = './'
cs_stru = Path(cs_dir, 'STRU')
cs_atoms = read(cs_stru, format='abacus')
pp = {'Al':'Al.PD04.PBE.UPF'}
input_data = {
        'control':{
            'calculation':'scf',
            'restart_mode':'from_scratch',
            'pseudo_dir':'./',
            'wf_collect':True,
            'etot_conv_thr':1.0e-6,
            'verbosity':'high',
            'tstress':True,
            'tprnfor':True
            },
        'system':{
            'ibrav':0,
            'nat':4,
            'ntyp':1,
            'nbnd':16,
            'ecutwfc':50,
            'occupations':'smearing',
            'smearing':'gauss',
            'degauss':0.01
            },
        'electrons':{
            'diagonalization':'cg',
            'mixing_mode':'plain',
            'mixing_beta':0.7,
            'conv_thr':1.0e-8
            }
        }
kpts = [21,21,21]
calc = Espresso(pseudopotentials=pp,input_data=input_data,kpts=kpts)
calc.write_input(cs_atoms)
