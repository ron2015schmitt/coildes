
#include "coils.hpp"

// Standard C++ libraries

#include <iostream>
#include <complex>

using namespace std;






void omegasubspace(LAvector<p3vector<double> >& grad_r, Matrix<complex<double> >& fs, 
		   LAvector<double>& nn, LAvector<double>& mm,
		   Matrix<complex<double> >& fsR2,
		   LAvector<double>& nnR2, LAvector<double>& mmR2);
