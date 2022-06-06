# abacus-ase
[ABACUS](https://github.com/abacusmodeling/abacus-develop) calculator for [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/index.html) (ASE) which is developed based on codes written by [ZhenXiong Shen](https://gitee.com/wszhang/ase_calculator_abacus). This updated interface optimizes support for ABACUS v2.2.0 onwards and for many ASE functions, such as `ase.io.read` and `ase.io.write`.

## Features
1. Only versions after ABACUS v2.2.0 are supported, if you want to use ASE with previous ABACUS, please see [ZhenXiong Shen](https://gitee.com/wszhang/ase_calculator_abacus)
2. ABACUS 'STRU' and 'running_*log' files can be easily parsed by `read(filename='STRU', format='abacus')` and `read(filename='running_md.log', format='abacus-out')`
3. One can use `write('STRU', images, format='abacus')` to write ASE `Atoms` objects to ABACUS 'STRU' files 
4. Energy, Force, Stress, Fermi level, K-points, Eigenvalues, Occupations and MD information can be easily get by `atoms.calc.get_*` method.
5. Pseudopotential and orbital settings have been changed to two parameters with `dict` type:   `pp` and `basis`, e.g.
    ```python
    pp = {'Na': 'Na.UPF',
    'Cl': 'Cl.UPF'}
    basis = {'Na': 'Na.orb',
    'Cl': 'Cl.orb'}
    ```  
    The directory for the pseudopotential and orbitals files can be set with the `pseudo_dir` and `basis_dir` parameters, respectively.
6. K-point settings have been changed if one want to generate a Monkhorst-Pack grid in ABACUS calculations, e.g. 
    ```
    knumber = 0,           # The number of K points
    kmode = Gamma,        # Mode of K points, can be Gamma, MP, Line, Direct, Cartesian
    kpts = [1, 1, 1],       # Give the K points
    koffset = [0, 0, 0]    # Give the displacement of K points
    ```
For more tutorials for ASE, please see [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/index.html)

## Installation
1. Copy `abacus` folder under `ase.calculator` folder
2. Copy `abacus.py` folder under `ase.io` folder
3. Add abacus-related codes to 'calculator.py' under `ase.calculator` folder and to 'formats.py' under `ase.io` folder, you'd better not override these two files directly using 'calculator.py' and 'formats.py' provided here, due to version differences

## Issues
1. Because of my ability is limited, only partial input parameters of ABACUS can be parsed, if you want to use ASE to generate ABACUS 'INPUT' file, there will be something wrong.