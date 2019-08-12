# coildes
Stellarator Coil Design Thesis Code


## Installation Requirements

The commands below are given for Ubuntu.  They have been verified using Ubuntu 18.04.

### Fortran compiler

The first step is to install the [GNU Fortran compiler gfortran](https://gcc.gnu.org/wiki/GFortran)

```
sudo apt-get install gfortran
```

### OPAM
```
sudo apt install opam
opam init
```
### ODEPACK
```
opam depext odepack.0.6.8
opam install odepack
```
The library to link is ```~/.opam/system/lib/odepack/odepack.a```
[ODEPACK](https://computing.llnl.gov/casc/odepack/)


## matricks library

