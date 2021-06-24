# MIP high-throuput workflow scripts
scripts for high_throuput workflow run_pbe.lsf for pbe. all correspoding INCAR files, such as INCAR.SCF, for every steps are needed.
There are all programs or scripts which used in high_throuput workflow from test_spin, relax, scf, to band
## 1. Three linux scripts needed, inculding `GETPOTCARPAW`, `GETINCAR` and `GETBANDK.sh`, 

* [GETPOTCARPAW](./program/GETPOTCARPAW): need two input files (PAW and POSCAR), element row in POSCAR is needed, pay attention to write a correct path in GETPOTCARPAW for reading potentials and PAW file.

* [GETINCAR](./program/GETINCAR) is for LDA+U calculations. Two files (LDAU POSCAR) is needed. Pay attention to the correct path.

* [GETBANDK.sh](./program/GETBANDK.sh) need to install phonopy, and different versions of phonopy will lead to different "$spacegroup and $symbol reading rows", GETBANDK.sh in this folder is for phonopy-1.12.

## 2. Some programs are needed

* [lattice_parameters](./program/lattice_parameters) for read lattice parameters, input file is POSCAR. Compiler: ifort -o lattice_parameters lattice_parameters.f90, you can put lattice_parameters in your bin folder.

* [getgap](./program/getgap/getgap) for reading band gap. Compiler: ifort -o getgap getgap.f90

* [chgsum.pl](./program/chgsum.pl) and [bader](./program/bader) files for bader analysis.

* [dospar](./program/dospar/dospar) for density of states analysis. Compiler:make dospar

* [getfulleig](./program/getfulleig.f90) for fermisurface which is in getbxsf folder: Compiler: ifort -o getfulleig getfulleig.f90

* [TransOpt](https://github.com/yangjio4849/TransOpt.git) for  calculate electrical transport properties, details in https://github.com/yangjio4849/TransOpt.git

