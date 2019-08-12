/************************************************************************* 
 * 
 *   File Name    :  gammamatrix.hpp
 *   Platform     :  gcc compiler (v3.2.2)
 *   Author       :  Ron Schmitt
 *   Date         :  16 May 2003
 * 
 *
 *   SYNOPSIS     
 *     Generates gamma matrix (shape parameter displacement)
 *
 **************************************************************************/


void gammamatrix( Matrix<complex<double> > & gamma,
		  const Matrix<complex<double> >& fs,
		  const FourierSurface& fsurface, 
		  const LAvector<p3vector<double> > A,
		  const LAvector<double>& thetas, const LAvector<double>& phis,
		  const unsigned int Ntheta, const unsigned int Nphi);


void fullgammamatrix( Matrix<complex<double> > & gamma,
		      const LAvector<double>& nn, const LAvector<double>& mm, 
		      const Matrix<complex<double> >& fsIN,
		      const Matrix<complex<double> >& fsOUT,
		      const LAvector<p3vector<double> > A,
		      const LAvector<double>& thetas, const LAvector<double>& phis,
		      const unsigned int Ntheta, const unsigned int Nphi);

