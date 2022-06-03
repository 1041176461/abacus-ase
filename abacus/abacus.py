from __future__ import print_function
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 16:33:38 2018

Modified on Wed Jun 20 15:00:00 2018
@author: Shen Zhen-Xiong

Modified on Wed Jun 03 23:00:00 2022
@author: Ji Yu-yang
"""

import os
import subprocess
import numpy as np

from ase.io import write
from ase.calculators.abacus.create_input import AbacusInput
from ase.calculators.calculator import FileIOCalculator  # Calculator


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

        self.command = command
        self.txt = txt

    # Initialize parameters and get some information -END-

    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)

        if changed_parameters:
            self.reset()
        return changed_parameters

    def set_atoms(self, atoms):
        self.atoms = atoms

    def check_state(self, atoms):
        system_changes = FileIOCalculator.check_state(self, atoms)
        # Ignore boundary conditions:
        if 'pbc' in system_changes:
            system_changes.remove('pbc')
        return system_changes

    def initialize(self, atoms):
        numbers = np.unique(atoms.get_atomic_numbers())
        self.system_params["ntype"] = len(numbers)

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)

        if scaled is None:
            scaled = np.all(atoms.get_pbc())

        write(os.path.join(self.directory, 'STRU'), atoms, **self.parameters)

        self.initialize(atoms)
        AbacusInput.write_input(self, directory=self.directory)
        AbacusInput.write_kpt(self, directory=self.directory)
        AbacusInput.write_pp(
            self, pp=self.parameters['pp'], directory=self.directory, pseudo_dir=self.parameters.pop('pseudo_dir', None))
        if 'basis' in self.parameters.keys():
            AbacusInput.write_orb(
                self, basis=self.parameters['basis'], directory=self.directory, pseudo_dir=self.parameters.pop('basis_dir', None))
        if 'offsite_basis' in self.parameters.keys():
            AbacusInput.write_abfs(self, offsite_basis=self.parameters['offsite_basis'], directory=self.directory, pseudo_dir=self.parameters.pop(
                'offsite_basis_dir', None))

    def run(self):
        with open(self.txt, 'a') as f:
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
