# coildes

```diff
- CURRENTLY NOT READY FOR USE. CHECK BACK IN SEPTEMBER 2019 -
```

[Stellarator](https://en.wikipedia.org/wiki/Stellarator) Coil Design Thesis Code


## Installation Requirements

The commands below are given for [Ubuntu Linux](https://en.wikipedia.org/wiki/Ubuntu).  They have been verified using Ubuntu 18.04.

The first step is to update your Ubuntu installation

```
sudo apt update
```


### Fortran compiler

Install the GNU Fortran compiler [gfortran](https://en.wikipedia.org/wiki/GNU_Fortran)

```
sudo apt install gfortran
```

### Numerical packages BLAS and LAPACK

Install [BLAS](https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms), [LAPACK](https://en.wikipedia.org/wiki/LAPACK), as well as [LAPACKE](https://www.netlib.org/lapack/lapacke.html).

```
sudo apt install liblapack3
sudo apt install liblapack-dev
sudo apt install libopenblas-base
sudo apt install libopenblas-dev
sudo apt install liblapacke-dev
```

### Fourier Transform Package FFTW

Install [FFTW](https://en.wikipedia.org/wiki/FFTW).

```
sudo apt install libfftw3-dev libfftw3-doc
```

## matricks library

The [matricks](https://github.com/ron2015schmitt/matricks) linear algebra library is used throughout the code.

The matricks library is included as a submodule.

```
cd matricks
./configure
cd ..
```

### ODEPACK

The package [ODEPACK](https://computing.llnl.gov/casc/odepack/) is used to follow the field lines of the magnetic field configuration.

A [ODEPACK library mirror](https://github.com/jacobwilliams/odepack) is included as a submodule.




