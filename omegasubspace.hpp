
#include "coils.hpp"

// Standard C++ libraries

#include <iostream>
#include <complex>

using namespace std;






void omegasubspace(Vector<p3vector<double> >& grad_r, Matrix<complex<double> >& fs, 
		   Vector<double>& nn, Vector<double>& mm,
		   Matrix<complex<double> >& fsR2,
		   Vector<double>& nnR2, Vector<double>& mmR2);
