# coildes

```diff
- CURRENTLY NOT READY FOR USE. CHECK BACK IN SEPTEMBER 2019 -
```

[Stellarator](https://en.wikipedia.org/wiki/Stellarator) Coil Design Thesis Code


## Installation Requirements

The commands below are given for [Ubuntu Linux](https://en.wikipedia.org/wiki/Ubuntu).  They have been verified using Ubuntu 18.04.

### Fortran compiler

The first step is to install the GNU Fortran compiler [gfortran](https://en.wikipedia.org/wiki/GNU_Fortran)

```
sudo apt update
sudo apt install gfortran
```

### Numerical packages BLAS and LAPACK

The following installs [BLAS](https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms), [LAPACK](https://en.wikipedia.org/wiki/LAPACK), as well as [LAPACKE](https://www.netlib.org/lapack/lapacke.html).

```
sudo apt update
sudo apt install liblapack3
sudo apt install liblapack-dev
sudo apt install libopenblas-base
sudo apt install libopenblas-dev
sudo apt install liblapacke-dev
sudo apt install liblapack-dev
```

### OPAM
```
sudo apt update
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

