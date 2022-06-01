# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 10:31:30 2018

@author: shenzx
"""

from os.path import join, basename, exists
from ase import Atoms
import numpy as np
import shutil
from ase.calculators.abacus.potential import PotentialDict
from ase.calculators.abacus.basis import BasisDict
from ase.utils import reader, writer
# from abacus.potential import PotentialDict
# from abacus.basis import BasisDict
# import sys
# sys.path.append("E:\Git\project\ase-abacus\ase-abacus")


def potential_list():
    return list(PotentialDict.keys())


def basis_list():
    return list(BasisDict.keys())


def judge_exist_stru(stru=None):
    if stru is None:
        return False
    else:
        return True


def read_ase_stru(stru=None, coordinates_type="Cartesian"):
    if judge_exist_stru(stru):
        atoms_list = []
        atoms_position = []
        atoms_masses = []
        atoms_magnetism = []
        atoms_all = stru.get_chemical_symbols()

        # sort atoms according to atoms
        for atoms_all_name in atoms_all:
            temp = True
            for atoms_list_name in atoms_list:
                if atoms_all_name == atoms_list_name:
                    temp = False
                    break

            if temp:
                atoms_list.append(atoms_all_name)

        for atoms_list_name in atoms_list:
            atoms_position.append([])
            atoms_masses.append([])
            atoms_magnetism.append([])

        # get position, masses, magnetism from ase atoms
        if coordinates_type == 'Cartesian':
            for i in range(len(atoms_list)):
                for j in range(len(atoms_all)):
                    if atoms_all[j] == atoms_list[i]:
                        atoms_position[i].append(list(
                            stru.get_positions()[j]))
                        atoms_masses[i] = stru.get_masses()[j]
                        atoms_magnetism[i] = list(
                            stru.get_initial_magnetic_moments())

        elif coordinates_type == 'Direct':
            for i in range(len(atoms_list)):
                for j in range(len(atoms_all)):
                    if atoms_all[j] == atoms_list[i]:
                        atoms_position[i].append(list(
                            stru.get_scaled_positions()[j]))
                        atoms_masses[i] = stru.get_masses()[j]
                        atoms_magnetism[i] = list(
                            stru.get_initial_magnetic_moments())
        else:
            raise ValueError("'coordinates_type' is ERROR,"
                             "please set to 'Cartesian' or 'Direct'")

        return atoms_list, atoms_masses, atoms_position, atoms_magnetism


def write_input_stru_core(fd,
                          stru=None,
                          potential=None,
                          pseudo_dir='./',
                          basis=None,
                          orbital_dir='./',
                          offsite_basis=None,
                          offsite_orbital_dir='./',
                          coordinates_type="Cartesian",
                          atoms_list=None,
                          atoms_position=None,
                          atoms_masses=None,
                          atoms_magnetism=None,
                          fix=1):
    if not judge_exist_stru(stru):
        return "No input structure!"

    elif (atoms_list is None):
        return "Please set right atoms list"
    elif(atoms_position is None):
        return "Please set right atoms position"
    elif(atoms_masses is None):
        return "Please set right atoms masses"
    elif(atoms_magnetism is None):
        return "Please set right atoms magnetism"
    else:
        fd.write('ATOMIC_SPECIES\n')
        for i, elem in enumerate(atoms_list):
            if not exists(potential[elem]):
                pseudofile = join(pseudo_dir, basename(potential[elem]))
            else:
                pseudofile = potential[elem]
            temp1 = ' ' * (4-len(atoms_list[i]))
            temp2 = ' ' * (14-len(str(atoms_masses[i])))
            atomic_species = (atoms_list[i] + temp1
                              + str(atoms_masses[i]) + temp2
                              + pseudofile)

            fd.write(atomic_species)
            fd.write('\n')

        if basis is not None:
            fd.write('\n')
            fd.write('NUMERICAL_ORBITAL\n')
            for i, elem in enumerate(atoms_list):
                if not exists(basis[elem]):
                    orbitalfile = join(orbital_dir, basename(basis[elem]))
                else:
                    orbitalfile = basis[elem]
                fd.write(orbitalfile)
                fd.write('\n')

        if offsite_basis is not None:
            fd.write('\n')
            fd.write('ABFS_ORBITAL\n')
            for i, elem in enumerate(atoms_list):
                if not exists(offsite_basis[elem]):
                    orbitalfile = join(offsite_orbital_dir,
                                       basename(offsite_basis[elem]))
                else:
                    orbitalfile = offsite_basis[elem]
            fd.write(orbitalfile)
            fd.write('\n')

        fd.write('\n')
        fd.write('LATTICE_CONSTANT\n')
        fd.write('1.889726125 \n')
        fd.write('\n')

        fd.write('LATTICE_VECTORS\n')
        for i in range(3):
            for j in range(3):
                temp3 = str("{:0<12f}".format(
                    stru.get_cell()[i][j])) + ' ' * 3
                fd.write(temp3)
                fd.write('   ')
            fd.write('\n')
        fd.write('\n')

        fd.write('ATOMIC_POSITIONS\n')
        fd.write(coordinates_type)
        fd.write('\n')
        fd.write('\n')
        for i in range(len(atoms_list)):
            fd.write(atoms_list[i])
            fd.write('\n')
            fd.write(str("{:0<12f}".format(atoms_magnetism[i][0])))
            fd.write('\n')
            fd.write(str(len(atoms_position[i])))
            fd.write('\n')

            for j in range(len(atoms_position[i])):
                temp4 = str("{:0<12f}".format(
                    atoms_position[i][j][0])) + ' ' * 3
                temp5 = str("{:0<12f}".format(
                    atoms_position[i][j][1])) + ' ' * 3
                temp6 = str("{:0<12f}".format(
                    atoms_position[i][j][2])) + ' ' * 3
                sym_pos = (temp4 + temp5 + temp6 +
                           (str(fix) + '   ') * 3)
                fd.write(sym_pos)
                fd.write('\n')
            fd.write('\n')


@writer
def write_abacus(fd,
                 stru=None,
                 pseudo_dir='./',
                 pps=None,
                 basis_dir='./',
                 basis=None,
                 offsite_basis_dir=None,
                 offsite_basis=None,
                 fix=1,
                 scaled=False,
                 **kwargs):

    if scaled:
        coordinates_type = 'Direct'
    else:
        coordinates_type = 'Cartesian'

    if not judge_exist_stru(stru):
        return "No input structure!"

    else:
        (atoms_list,
         atoms_masses,
         atoms_position,
         atoms_magnetism) = read_ase_stru(stru, coordinates_type)

        write_input_stru_core(fd,
                              stru,
                              pps,
                              pseudo_dir,
                              basis,
                              basis_dir,
                              offsite_basis,
                              offsite_basis_dir,
                              coordinates_type,
                              atoms_list,
                              atoms_position,
                              atoms_masses,
                              atoms_magnetism,
                              fix)


@reader
def read_abacus(fd,
                ase=True,
                **kwargs):
    # Read structure information from abacus structure file
    # try:
    #     f = open(join(directory, filename), 'r')
    # except Exception:
    #     return "Failed to open 'STRU', Please Check!"
    # else:
    #     lines = f.readlines()
    #     f.close()

    lines = fd.readlines()
    # initialize reading information
    temp = []
    for line in lines:
        line = line.strip()
        line = line.replace('\n', ' ')
        line = line.replace('\t', ' ')
        line = line.replace('//', ' ')
        line = line.replace('#', ' ')

        if len(line) != 0:
            temp.append(line)

    # print(temp)
    atom_species = 0
    for i in range(len(temp)):
        if temp[i] == 'NUMERICAL_ORBITAL':
            atom_species = i - 1
            break

    atom_symbol = []
    atom_mass = []
    atom_potential = []
    atom_number = []
    atom_magnetism = []
    atom_positions = []
    atom_appendix = []

    # get symbol, mass, potential
    for i in range(1, atom_species+1):
        atom_symbol.append(temp[i].split()[0])
        atom_mass.append(float(temp[i].split()[1]))
        atom_potential.append(temp[i].split()[2])
        atom_number.append(0)
        atom_magnetism.append(0)
        atom_positions.append([])
        atom_appendix.append([])

    # get abfs basis
    atom_basis = []
    atom_offsite_basis = []
    for i in range(atom_species+2, (atom_species+1) * 2):
        atom_basis.append(temp[i].split()[0])
    if 'ABFS_ORBITAL' in temp:
        scale = 3
        for i in range((atom_species+1) * 2 + 1, (atom_species+1) * 3):
            atom_offsite_basis.append(temp[i].split()[0])
    else:
        scale = 2

    # get lattice
    atom_lattice_scale = float(temp[(atom_species+1) * scale + 1].split()[0])
    atom_lattice = np.array(
        [[float(temp[(atom_species+1) * scale + 3 + i].split()[:3][j])
          for j in range(3)] for i in range(3)])

    # get coordinates type
    atom_coor = temp[(atom_species + 1) * scale + 7].split()[0]

    # get position,  atoms number, magnetism, fix
    for i in range(atom_species):
        pos_start = (atom_species + 1) * scale + 8 + 3 * i
        for j in range(i):
            pos_start += atom_number[j]
        atom_it = atom_symbol.index(temp[pos_start].split()[0])
        atom_magnetism[atom_it] = float(temp[pos_start + 1].split()[0])
        atom_number[atom_it] = int(temp[pos_start + 2].split()[0])

        atom_positions[atom_it] = np.array(
            [[float(temp[pos_start + 3 + i].split()[:3][j])
              for j in range(3)] for i in range(atom_number[atom_it])])

        atom_appendix[atom_it] = np.array(
            [[temp[pos_start + 3 + i].split()[3:6][j]
              for j in range(3)]for i in range(atom_number[atom_it])])

    # Reset structure information and return results
    formula_symbol = ''
    formula_positions = []
    for i in range(atom_species):
        if atom_number[i] == 1:
            formula_symbol += atom_symbol[i]

        else:
            formula_symbol += atom_symbol[i] + str(atom_number[i])

        for j in range(atom_number[i]):
            formula_positions.append(atom_positions[i][j])

    formula_cell = atom_lattice * atom_lattice_scale * 0.529177210903

    if ase is True:
        if atom_coor == 'Direct':
            return Atoms(symbols=formula_symbol,
                         cell=formula_cell,
                         scaled_positions=formula_positions)

        elif atom_coor == 'Cartesian':
            return Atoms(symbols=formula_symbol,
                         cell=formula_cell,
                         positions=formula_positions)

        else:
            raise ValueError("atomic coordinate type is ERROR")

    else:
        return (formula_symbol,
                formula_cell,
                formula_positions,
                atom_potential,
                atom_basis,
                atom_offsite_basis)
