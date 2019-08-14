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
		  const Vector<p3vector<double> > A,
		  const Vector<double>& thetas, const Vector<double>& phis,
		  const unsigned int Ntheta, const unsigned int Nphi);


void fullgammamatrix( Matrix<complex<double> > & gamma,
		      const Vector<double>& nn, const Vector<double>& mm, 
		      const Matrix<complex<double> >& fsIN,
		      const Matrix<complex<double> >& fsOUT,
		      const Vector<p3vector<double> > A,
		      const Vector<double>& thetas, const Vector<double>& phis,
		      const unsigned int Ntheta, const unsigned int Nphi);

