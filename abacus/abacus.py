from __future__ import print_function
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 16:33:38 2018

Modified on Wed Jun 20 15:00:00 2018
@author: Shen Zhen-Xiong

Modified on Wed Jun 03 23:00:00 2022
@author: Ji Yu-yang
"""

import subprocess
from os.path import join
import numpy as np
from ase.calculators.abacus.create_input import AbacusInput
from ase.calculators.calculator import FileIOCalculator, all_changes  #Calculator


class Abacus(AbacusInput, FileIOCalculator):
    # Initialize parameters and get some information -START-
    name = 'abacus'
    implemented_properties = ['energy', 'forces', 'fermi', 'stress']

    default_parameters = dict(calculation='scf',
                              ecutwfc=50,
                              smearing_method='gaussian',
                              mixing_type='pulay-kerker',
                              basis_type='lcao',
                              gamma_only=1,
                              ks_solver="genelpa",
                              stru_file='STRU',
                              )

    def __init__(self,
                 restart=None,
                 ignore_bad_restart_file=False,
                 directory='.',
                 label='abacus',
                 atoms=None,
                 command=None,
                 txt='abacus.out',
                 **kwargs):

        self.results = {}

        # Initialize parameter dictionaries
        AbacusInput.__init__(self, restart)

        # Set directory and label
        self.directory = directory
        self.label = label

        FileIOCalculator.__init__(self,
                                  restart,
                                  ignore_bad_restart_file,
                                  label,
                                  atoms,
                                  **kwargs)

        self.restart = restart
        self.pseudo_dir = pseudo_dir
        self.potential_name = potential_name
        self.basis_dir = basis_dir
        self.basis_name = basis_name
        self.offsite_basis_dir = offsite_basis_dir
        self.offsite_basis_name = offsite_basis_name
        self.fix = fix
        self.stru_filename = stru_filename
        self.coordinates_type = coordinates_type

        self.out_path = ''


        if log_file is not None:
            self.log_file = log_file

        AbacusInput.set(self, **self.parameters)
        AbacusInput.set(self, **kwargs)

    # Initialize parameters and get some information -END-

    def check_state(self, atoms):
        system_changes = FileIOCalculator.check_state(self, atoms)
        # Ignore boundary conditions:
        if 'pbc' in system_changes:
            system_changes.remove('pbc')
        return system_changes

    def initialize(self, atoms):
        numbers = atoms.get_atomic_numbers().copy()
        self.species = []
        for a, Z in enumerate(numbers):
            if Z not in self.species:
                self.species.append(Z)
        self.general_params["ntype"] = len(self.species)

    # Run abacus
    def calculate(self,
                  atoms=None,
                  properties=None,
                  system_changes=all_changes):
        FileIOCalculator.calculate(self,
                                   atoms,
                                   properties,
                                   system_changes)

    def run(self):
        with open(self.log_file, 'a') as f:
            run = subprocess.Popen(self.command,
                                   stderr=f,
                                   stdin=f,
                                   stdout=f,
                                   cwd=self.directory,
                                   shell=True)
            return run.communicate()

    def get_fermi_level(self):
        return self.results['fermi']

    """
    def get_potential_energy(self, atoms):
        return self.get_property('energy', atoms)

    def get_forces(self, atoms):
        return self.get_property('forces', atoms)

    def get_property(self, name, atoms = None, allow_calculation = True):
        if atoms is None:
            atoms = self.atoms
            system_changes = []
        else:
            system_changes = self.check_state(atoms)
            if system_changes:
                self.reset()
        if name not in self.results:
            if not allow_calculation:
                return None
            self.calculate(atoms)
        result = self.results[name]
        return result
    """


if __name__ == "__main__":
    pass
