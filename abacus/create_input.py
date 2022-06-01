# -*- coding: utf-8 -*-
"""
Created on Sun Dec 16 17:02:37 2018

@author: shenzx

Modified on Wed Jun 01 15:00:00 2022
@author: Ji Yu-yang
"""
from __future__ import print_function
import warnings
import shutil
from os.path import join
import numpy as np
from ase.io.abacus import write_input_stru
from ase.calculators.calculator import FileIOCalculator
# copyright © Key Lab of Quantum Information, CAS, China
"""This module defines an ASE interface to ABACUS

Developed on the basis of modules by Zhen-Xiong Shen.
 The path of the directory containing the
 pseudopotential and basis directories (LDA, PBE, SG15, ORBITAL, ...)
 should be set by the enviromental flag $ABACUS_PP_PATH, $ABACUS_ORBITAL_PATH.

The user should also set the enviroment flag
 $ABACUS_SCRIPT pointing to a python script looking

like::
    import os
    exitcode = os.system('abacus')
http://abacus.ustc.edu.cn/
"""

# Parameters list that can be set in INPUT.  -START-
# 1
general_keys = [
        'suffix',              # the name of main output directory
        'latname',             # the name of lattice name
        'stru_file',           # the filename of file containing atom positions
        'kpoint_file',         # the name of file containing k points
        'pseudo_dir',          # the directory containing pseudo files
        'orbital_dir',         # the directory containing orbital files
        'pseudo_type',         # the type pseudo files
        'pseudo_rcut',         # cut-off radius for radial integration
        'pseudo_mesh',         # 0: use our own mesh to do radial renormalization; 1: use mesh as in QE
        'lmaxmax',             # maximum of l channels used
        'dft_functional',      # exchange correlation functional
        'calculation',         # test; scf; relax; nscf; ienvelope; istate;
        'ntype',               # atom species number
        'nspin',               # 1: single spin; 2: up and down spin; 4: noncollinear spin
        'nbands',              # number of bands
        'nbands_sto',          # number of stochastic bands
        'nbands_istate',       # number of bands around Fermi level for istate calulation
        'nche_sto',            # number of orders for Chebyshev expansion in stochastic DFT
        'symmetry',            # turn symmetry on or off
        'init_vel',            # read velocity from STRU or not
        'symmetry_prec',       # accuracy for symmetry
        'nelec',               # input number of electrons
        'tot_magnetization',   # total magnetization of the system
        'out_mul',             # mulliken  charge or not
        'noncolin',            # using non-collinear-spin
        'lspinorb',            # consider the spin-orbit interaction
        'gamma_only'           # It is an important parameter only to be used in localized orbitals set. 
        ]
# 2
pw_keys = [
        'ecutwfc',             # energy cutoff for wave functions
        'pw_diag_thr',         # threshold for eigenvalues is cg electron iterations
        'scf_thr',             # charge density error
        'init_wfc',            # start wave functions are from 'atomic' or 'file'
        'init_chg',            # start charge is from 'atomic' or file
        'chg_extrap',          # atomic; first-order; second-order; dm:coefficients of SIA
        'out_chg',             # >0 output charge density for selected electron steps
        'out_pot',             # output realspace potential
        'out_wfc_pw',          # output wave functions
        'out_wfc_r',           # output wave functions in realspace
        'out_dos',             # output energy and dos
        'out_band',            # output k-index and band structure
        'out_proj_band',       # output projected band structure
        'restart_save',        # print to disk every step for restart
        'restart_load',        # restart from disk
        'read_file_dir',       # directory of files for reading
        'nx',                  # number of points along x axis for FFT grid
        'ny',                  # number of points along y axis for FFT grid
        'nz',                  # number of points along z axis for FFT grid
        'cell_factor'          # used in the construction of the pseudopotential tables
        ] 
# 3
relaxation_keys = [
        'ks_solver',                # cg; dav; lapack; genelpa; hpseps; scalapack_gvx
        'scf_nmax',                 # #number of electron iterations
        'out_force',                # output the out_force or not
        'relax_nmax',               # number of ion iteration steps
        'out_stru',                 # output the structure files after each ion step
        'force_thr',                # force threshold, unit: Ry/Bohr
        'force_thr_ev',             # force threshold, unit: eV/Angstrom
        'force_thr_ev2',            # force invalid threshold, unit: eV/Angstrom
        'relax_cg_thr',             # threshold for switching from cg to bfgs, unit: eV/Angstrom
        'stress_thr',               # stress threshold
        'press1',                   # target pressure, unit: KBar
        'press2',                   # target pressure, unit: KBar
        'press3',                   # target pressure, unit: KBar
        'relax_bfgs_w1',            # wolfe condition 1 for bfgs
        'relax_bfgs_w2',            # wolfe condition 2 for bfgs
        'relax_bfgs_rmax',          # maximal trust radius, unit: Bohr
        'relax_bfgs_rmin',          # minimal trust radius, unit: Bohr
        'relax_bfgs_init',          # initial trust radius, unit: Bohr
        'cal_stress',               # calculate the stress or not
        'fixed_axes',               # which axes are fixed
        'relax_method',             # bfgs; sd; cg; cg_bfgs;
        'out_level',                # ie(for electrons); i(for ions);
        'out_dm',                   # >0 output density matrix
        'deepks_out_labels',        # >0 compute descriptor for deepks
        'deepks_scf',               # >0 add V_delta to Hamiltonian
        'deepks_bandgap',           # >0 for bandgap label
        'deepks_out_unittest',      # if set 1, prints intermediate quantities that shall be used for making unit test
        'deepks_model',             # file dir of traced pytorch model: 'model.ptg
        'deepks_descriptor_lmax2',  # lmax used in generating descriptor
        'deepks_descriptor_rcut0',  # rcut used in generating descriptor
        'deepks_descriptor_ecut0'   # ecut used in generating descriptor
        ]
# 4
lcao_keys = [
        'basis_type',          # PW; LCAO in pw; LCAO
        'nb2d',                # 2d distribution of atoms
        'search_radius',       # input search radius (Bohr)
        'search_pbc',          # input periodic boundary condition
        'lcao_ecut',           # energy cutoff for LCAO
        'lcao_dk',             # delta k for 1D integration in LCAO
        'lcao_dr',             # delta r for 1D integration in LCAO
        'lcao_rmax',           # max R for 1D two-center integration table
        'out_mat_hs',          # output H and S matrix
        'out_mat_hs2',         # output H(R) and S(R) matrix
        'out_mat_r',           # output r(R) matrix
        'out_wfc_lcao',        # ouput LCAO wave functions
        'bx',                  # division of an element grid in FFT grid along x
        'by',                  # division of an element grid in FFT grid along y
        'bz'                   # division of an element grid in FFT grid along z
        ]
# 5
smearing_keys = [
        'smearing_method',     # type of smearing_method: gauss; fd; fixed; mp; mp2; mv
        'smearing_sigma'       # energy range for smearing
        ]
# 6
charge_mixing_keys = [
        'mixing_type',      # plain; kerker; pulay; pulay-kerker
        'mixing_beta',      # mixing parameter: 0 means no new charge
        'mixing_ndim',      # mixing dimension in pulay
        'mixing_gg0'        # mixing parameter in kerker
        ]
# 7
dos_keys = [
        'dos_emin_ev',      # minimal range for dos
        'dos_emax_ev',      # maximal range for dos
        'dos_edelta_ev',    # delta energy for dos
        'dos_scale',        # scale dos range by
        'dos_sigma'         # gauss b coefficeinet(default=0.07)        
        ]
# 9
molecular_dynamics_keys = [
        'md_type',            # choose ensemble
        'md_nstep',           # md steps
        'md_ensolver',        # choose potential
        'md_dt',              # time step
        'md_mnhc',            # number of Nose-Hoover chains
        'md_tfirst',          # temperature first
        'md_tlast',           # temperature last
        'md_dumpfreq',        # The period to dump MD information
        'md_restartfreq',     # The period to output MD restart information
        'md_restart',         # whether restart
        'lj_rcut',            # cutoff radius of LJ potential
        'lj_epsilon',         # the value of epsilon for LJ potential
        'lj_sigma',           # the value of sigma for LJ potential
        'msst_direction',     # the direction of shock wave
        'msst_vel',           # the velocity of shock wave
        'msst_vis',           # artificial viscosity
        'msst_tscale',        # reduction in initial temperature
        'msst_qmass',         # mass of thermostat
        'md_tfreq',           # oscillation frequency, used to determine qmass of NHC
        'md_damp'             # damping parameter (time units) used to add force in Langevin method
        ]

# 10
test_keys = [
        'out_alllog',       # output information for each processor, when parallel
        'nurse',            # for coders
        'colour',           # for coders, make their live colourful
        't_in_h',           # calculate the kinetic energy or not
        'vl_in_h',          # calculate the local potential or not
        'vnl_in_h',         # calculate the nonlocal potential or not
        'vh_in_h',          # calculate the hartree potential or not
        'vion_in_h',        # calculate the local ionic potential or not
        'test_force',       # test the force
        'test_stress',      # test the force
        ]
# 11
vdw_keys = [
        'vdw_method',       # the method of calculating vdw (none ; d2 ; d3_0 ; d3_bj
        'vdw_s6',           # scale parameter of d2/d3_0/d3_bj
        'vdw_s8',           # scale parameter of d3_0/d3_bj
        'vdw_a1',           # damping parameter of d3_0/d3_bj
        'vdw_a2',           # damping parameter of d3_bj
        'vdw_d',            # damping parameter of d2
        'vdw_abc',          # third-order term?
        'vdw_C6_file',      # filename of C6
        'vdw_C6_unit',      # unit of C6, Jnm6/mol or eVA6
        'vdw_R0_file',      # filename of R0
        'vdw_R0_unit',      # unit of R0, A or Bohr
        'vdw_model',        # expression model of periodic structure, radius or period
        'vdw_radius',       # radius cutoff for periodic structure
        'vdw_radius_unit',  # unit of radius cutoff for periodic structure
        'vdw_cn_thr',       # radius cutoff for cn
        'vdw_cn_thr_unit',  # unit of cn_thr, Bohr or Angstrom
        'vdw_period'        # periods of periodic structure
        ]
# 12
tddft_keys = [
        'tddft',               # calculate tddft or not
        'td_scf_thr',          # threshold for electronic iteration of tddft
        'td_dt',               # time of ion step
        'td_force_dt',         # time of force change
        'td_val_elec_01',      # td_val_elec_01
        'td_val_elec_02',      # td_val_elec_02
        'td_val_elec_03',      # td_val_elec_03
        'td_vext',             # add extern potential or not
        'td_vext_dire',        # extern potential direction
        'td_timescale',        # extern potential td_timescale
        'td_vexttype',         # extern potential type
        'td_vextout',          # output extern potential or not
        'td_dipoleout',        # output dipole or not
        'ocp',                 # change occupation or not
        'ocp_set'              # set occupation
] 
# 13
exx_keys = [
        'exx_hybrid_alpha',    
        'exx_hse_omega',       
        'exx_separate_loop',   
        'exx_hybrid_step',     
        'exx_lambda',          
        'exx_pca_threshold',   
        'exx_c_threshold',     
        'exx_v_threshold',     
        'exx_dm_threshold',    
        'exx_schwarz_threshold', 
        'exx_cauchy_threshold',
        'exx_ccp_threshold',   
        'exx_ccp_rmesh_times', 
        'exx_distribute_type', 
        'exx_opt_orb_lmax',    
        'exx_opt_orb_ecut',    
        'exx_opt_orb_tolerence'
]
# 14
berry_keys = [
    'berry_phase',        # calculate berry phase or not
    'gdir',               # calculate the polarization in the direction of the lattice vector
    'towannier90',        # use wannier90 code interface or not
    'nnkpfile',           # the wannier90 code nnkp file name
    'wannier_spin'        # calculate spin in wannier90 code interface
]

# 15
isol_keys = [
    'imp_sol',             # calculate implicit solvation correction or not
    'eb_k',                # the relative permittivity of the bulk solvent
    'tau',                 # the effective surface tension parameter
    'sigma_k',             # the width of the diffuse cavity
    'nc_k'                 # the cut-off charge density
]

# Parameters list that can be set in INPUT.  -END-

class AbacusInput(object):
    
    # Initialize internal dictionary of input parameters to None  -START-
    def __init__(self, restart=None):
        """
        self.directory = './'        # shenzx v20200724
        self.stru_filename = 'STRU'  # shenzx v20200724
        self.pseudo_dir = './'   # shenzx v20200724
        self.potential_name = None   # shenzx v20200724
        self.basis_dir = './'        # shenzx v20200724
        self.basis_name = None       # shenzx v20200724
        self.fix = 1                 # shenzx v20200724
        self.coordinates_type = 'Cartesian'      # shenzx v20200724
        """
        self.general_params = {}
        self.pw_params = {}
        self.relaxation_params = {}
        self.lcao_params = {}
        self.smearing_params = {}
        self.charge_mixing_params = {}
        self.dos_params = {}
        self.molecular_dynamics_params = {}
        self.test_params = {}
        self.vdw_params = {}
        self.tddft_params = {}
        self.exx_params = {}
        self.berry_params = {}
        self.isol_params = {}
        
        for key in general_keys:
            self.general_params[key] = None
        for key in pw_keys:
            self.pw_params[key] = None
        for key in relaxation_keys:
            self.relaxation_params[key] = None
        for key in lcao_keys:
            self.lcao_params[key] = None
        for key in smearing_keys:
            self.smearing_params[key] = None
        for key in charge_mixing_keys:
            self.charge_mixing_params[key] = None
        for key in dos_keys:
            self.dos_params[key] = None
        for key in molecular_dynamics_keys:
            self.molecular_dynamics_params[key] = None
        for key in test_keys:
            self.test_params[key] = None
        for key in vdw_keys:
            self.vdw_params[key] = None
        for key in tddft_keys:
            self.tddft_params[key] = None
        for key in exx_keys:
            self.exx_params[key] = None
        for key in berry_keys:
            self.berry_params[key] = None
        for key in isol_keys:
            self.isol_params[key] = None
        # Initialize internal dictionary of input parameters to None  -END-

        # Appoint the KPT parameters which are not INPUT parameters  -START-
        self.kpt_params = {
                'knumber': 0,           # The number of K points
                'kmode': 'Gamma',       # Mode of K points, can be Gamma, MP, Line, Direct, Cartesian
                'kpts': [1, 1, 1, 0, 0, 0]  # Give the K points
                }
        # Appoint the KPT parameters which are not INPUT parameters  -END-

    # Set the INPUT and KPT parameters  -START-
    def set(self, **kwargs):
        for key in kwargs:
            if key in self.general_params:
                self.general_params[key] = kwargs[key]
            elif key in self.pw_params:
                self.pw_params[key] = kwargs[key]
            elif key in self.relaxation_params:
                self.relaxation_params[key] = kwargs[key]
            elif key in self.lcao_params:
                self.lcao_params[key] = kwargs[key]
            elif key in self.smearing_params:
                self.smearing_params[key] = kwargs[key]
            elif key in self.charge_mixing_params:
                self.charge_mixing_params[key] = kwargs[key]
            elif key in self.dos_params:
                self.dos_params[key] = kwargs[key]
            elif key in self.molecular_dynamics_params:
                self.molecular_dynamics_params[key] = kwargs[key]
            elif key in self.test_params:
                self.test_params[key] = kwargs[key]
            elif key in self.vdw_params:
                self.vdw_params[key] = kwargs[key]
            elif key in self.tddft_params:
                self.spectrum_params[key] = kwargs[key]
            elif key in self.kpt_params:
                self.kpt_params[key] = kwargs[key]
            elif key in self.exx_params:
                self.exx_params[key] = kwargs[key]
            elif key in self.berry_params:
                self.berry_params[key] = kwargs[key]
            elif key in self.isol_params:
                self.isol_params[key] = kwargs[key]
            else:
                raise TypeError('Parameter not defined:  ' + key)
    # Set the INPUT and KPT parameters  -END-

    # Write INPUT file  -START-
    def write_input_input(self, directory='./', **kwargs):
        with open(join(directory, 'INPUT'), 'w') as input_file:
            input_file.write('INPUT_PARAMETERS\n')
            input_file.write('# Created by Atomic Simulation Enviroment\n')
            for key, val in self.general_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
            
            for key, val in self.pw_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
                    
            for key, val in self.relaxation_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
                    
            for key, val in self.lcao_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
                    
            for key, val in self.smearing_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
                    
            for key, val in self.charge_mixing_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
                    
            for key, val in self.dos_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
                    
            for key, val in self.molecular_dynamics_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
                    
            for key, val in self.test_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
                    
            for key, val in self.other_method_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
                    
            for key, val in self.vdw_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
                    
            for key, val in self.tddft_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
            
            for key, val in self.exx_params.items():
                if val is not None:
                    params = str(key) + ' ' * (30 - len(key)) + str(val)
                    input_file.write('%s\n' % params)

            for key, val in self.berry_params.items():
                if val is not None:
                    params = str(key) + ' ' * (30 - len(key)) + str(val)
                    input_file.write('%s\n' % params)

            for key, val in self.isol_params.items():
                if val is not None:
                    params = str(key) + ' ' * (30 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
    # Write INPUT file  -END-
    # Read  INPUT  file  --START-

    def read_input_input(self,
                         filename='INPUT',
                         directory='./',
                         **kwargs):
        with open(join(directory, filename), 'r') as file:
            file.readline()
            lines = file.readlines()

        for line in lines: 
            try:
                line = line.replace("# ", "#  ")
                data = line.split()
                if len(data) == 0:
                    continue
                elif data[0][0] == "# ":
                    continue
                
                key = data[0]
                if key in general_keys:
                    self.general_params[key] = data[1]
                elif key in pw_keys:
                    self.pw_params[key] = data[1]
                elif key in relaxation_keys:
                    self.relaxation_params[key] = data[1]
                elif key in lcao_keys:
                    self.lcao_params[key] = data[1]
                elif key in smearing_keys:
                    self.smearing_params[key] = data[1]
                elif key in charge_mixing_keys:
                    self.charge_mixing_params[key] = data[1]
                elif key in dos_keys:
                    self.dos_params[key] = data[1]
                elif key in molecular_dynamics_keys:
                    self.molecular_dynamics_params[key] = data[1]
                elif key in test_keys:
                    self.test_params[key] = data[1]
                elif key in vdw_keys:
                    if key == 'vdw_period':
                        self.vdw_params[key] = (data[1] + '  '
                        + data[2] + '  ' + data[3])
                    else:
                        self.vdw_params[key] = data[1]
                elif key in tddft_keys:
                    self.tddft_params[key] = data[1]
                elif key in exx_keys:
                    self.exx_params[key] = data[1]
                elif key in berry_keys:
                    self.berry_params[key] = data[1]
                elif key in isol_keys:
                    self.isol_params[key] = data[1]

                return 'ok'  
            
            except  KeyError:
                raise IOError('keyword "%s" in INPUT is'
                                  'not know by calculator.' % key)
                    
            except IndexError:
                raise IOError('Value missing for keyword "%s" .' % key)
    # Read  INPUT  file  --END- 

    # Write KPT  -START-
    def write_input_kpt(self,
                        directory='./',
                        filename='KPT',
                        **kwargs):
        k = self.kpt_params
        if self.general_params['gamma_only'] is None:
            return warnings.warn(" 'gamma_only' parameter has not been set, "
                                 "please set it to 0 or 1")

        elif self.general_params['gamma_only'] == 1:
            with open(join(directory, filename), 'w') as kpoint:
                kpoint.write('K_POINTS\n')
                kpoint.write('0\n')
                kpoint.write('Gamma\n')
                kpoint.write('1 1 1 0 0 0')

        elif self.general_params['gamma_only'] == 0:
            with open(join(directory, filename), 'w') as kpoint:
                kpoint.write('K_POINTS\n')
                kpoint.write('%s\n' % str(k['knumber']))
                kpoint.write('%s\n' % str(k['kmode']))
                if k['kmode'] in ['Gamma', 'MP']:
                    for n in range(len(k['kpts'])):
                        kpoint.write('%s  ' % str(k['kpts'][n]))
                            
                elif k['kmode'] in ['Direct', 'Cartesian', 'Line']:
                    for n in range(len(k['kpts'])):
                        for i in range(len(k['kpts'][n])):
                            kpoint.write('%s  ' % str(k['kpts'][n][i]))
                        kpoint.write('\n')

                else:
                    raise ValueError("The value of kmode is not right, set to "
                                     "Gamma, MP, Direct, Cartesian, or Line.")
        else:
            return warnings.warn("The value of gamma_only is not right, "
                                 "please set it to 0 or 1")
    # Write KPT  -END-

    # Read KPT file  -START-

    def read_kpt(self,
                 filename='KPT',
                 directory='./',
                 **kwargs):
        with open(filename, 'r') as file:
            lines = file.readlines()

        if lines[2][-1] == '\n':
            kmode = lines[2][:-1]
        else:
            kmode = lines[2]

        if kmode in ['Gamma', 'MP']:
            self.kpt_params['kmode'] = lines[2][:-1]
            self.kpt_params['knumber'] = lines[1].split()[0]
            self.kpt_params['kpts'] = np.array(lines[3].split())

        elif kmode in ['Cartesian', 'Direct', 'Line']:
            self.kpt_params['kmode'] = lines[2][:-1]
            self.kpt_params['knumber'] = lines[1].split()[0]
            self.kpt_params['kpts'] = np.array([list(map(float, line.split())) 
                                     for line in lines[3:]])

        else:
            raise ValueError("The value of kmode is not right, set to "
                                     "Gamma, MP, Direct, Cartesian, or Line.")
    # Read KPT file  -END-

    # Write and read potential  -START-
    def write_potential(self,
                        pseudo_dir='./',
                        potential_name=None,
                        directory='./',
                        **kwargs):
        if pseudo_dir == directory:
            return 'It is ok,  pseudo_dir is in work directory'
        else:
            if self.potential_name == None:
                raise ValueError('The value of "potential_name" is not right, '
                                 'please set it to a list')
            else:
                self.pseudo_dir = pseudo_dir
                self.potential_name = potential_name
                for name in self.potential_name:
                    shutil.copyfile(join(self.pseudo_dir, name),
                                    join(directory, name))

    def read_potential(self,
                       pseudo_dir='./',
                       potential_name=None,
                       **kwargs):
        if self.potential_name is None:
            raise ValueError('The value of "potential_name" is not right, '
                             'please set it to a list')
        else:
            self.pseudo_dir = pseudo_dir
            self.potential_name = potential_name
    # Write and read potential  -END-

    # Write and read orbital basis  -START-
    def write_basis(self,
                    basis_dir='./',
                    basis_name=None,
                    directory='./',
                    **kwargs):
        if basis_dir == directory:
            print('It is ok,  basis_dir is in work directory')
        else:
            if self.basis_name is None:
                raise ValueError('The value of "basis_name" is not right, '
                                 'please set it to a list')
            else:
                self.basis_dir = basis_dir
                self.basis_name = basis_name
                for name in self.basis_name:
                    shutil.copyfile(join(self.basis_dir, name),
                                    join(directory, name))

        print("basis_dir = %s, basis_name = %s"%(basis_dir, basis_name))

    def read_basis(self,
                   basis_dir='./',
                   basis_name=None,
                   directory='./',
                   **kwargs):
        if self.basis_name is None:
            raise ValueError('The value of "basis_name" is not right, '
                             'please set it to a list')
        else:
            self.basis_dir = basis_dir
            self.basis_name = basis_name
            
    # Write and read orbital basis  -START-
    def write_offsite_basis(self,
                            offsite_basis_dir=None,
                            offsite_basis_name=None,
                            directory='./',
                            **kwargs):
        if offsite_basis_dir == directory:
            print('It is ok,  offsite_basis_dir is in work directory')
        else:
            if self.offsite_basis_name is not None:
                self.offsite_basis_dir = offsite_basis_dir
                self.offsite_basis_name = offsite_basis_name
                for name in self.offsite_basis_name:
                    shutil.copyfile(join(self.offsite_basis_dir, name),
                                    join(directory, name))

        print("offsite_basis_dir = %s, offsite_basis_name = %s"%(offsite_basis_dir, offsite_basis_name))

    def read_offsite_basis(self,
                           offsite_basis_dir=None,
                           offsite_basis_name=None,
                           directory='./',
                           **kwargs):
        if self.offsite_basis_name is not None:
            self.offsite_basis_dir = offsite_basis_dir
            self.offsite_basis_name = offsite_basis_name

    # Write and read orbital basis  -START-
    def write_input(self,
                    atoms,
                    properties=None,
                    system_changes=None):
        """Write input parameters to files-file."""

        FileIOCalculator.write_input(self,
                                     atoms,
                                     properties,
                                     system_changes)

        if(system_changes is None):  # shenzx v20200724
            system_changes = '  '    # shenzx v20200724
        if ('numbers' in system_changes or 
                'initial_magmoms' in system_changes):
            self.initialize(atoms)

        write_input_stru(stru=atoms,
                         filename=self.stru_filename,
                         pseudo_dir=self.pseudo_dir,
                         potential_name=self.potential_name,
                         basis_dir=self.basis_dir,
                         basis_name=self.basis_name,
                         offsite_basis_dir=self.offsite_basis_dir,
                         offsite_basis_name=self.offsite_basis_name,
                         fix=self.fix,
                         directory=self.directory,
                         coordinates_type=self.coordinates_type)
        self.write_input_input(directory=self.directory)
        self.write_input_kpt(directory=self.directory)
    # Write all input file -END-


if __name__ == "__main__":
    # Test a writing...
    import os
    print(os.getcwd())
    from ase import Atoms  # just use when test this module
    test = AbacusInput()
    test.set(gamma_only = 1)
    test.write_input(atoms = Atoms("CO"))
