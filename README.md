# coildes
Stellarator Coil Design Thesis Code


## Installation Requirements


## The LSODE solver library from ODEPACK

[ODEPACK](https://computing.llnl.gov/casc/odepack/)

### A Fortran compiler
Ane of these fortran compilers: gfortran, g95, g77, f77

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

## matricks library

