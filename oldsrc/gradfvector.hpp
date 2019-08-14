/************************************************************************* 
 * 
 *   File Name    :  gradfvector.hpp
 *   Platform     :  gcc compiler (v3.2.2)
 *   Author       :  Ron Schmitt
 *   Date         :  7 May 2004
 * 
 *
 *   SYNOPSIS     
 *     Generates gradient of the orthonormal fourier series
 *
 **************************************************************************/

void gradfvector( const Vector<double>& thetas, const Vector<double>& phis,
		  const Vector<double>& mm, const Vector<double>& nn, 
		  const Vector<p3vector<double> >& grad_theta, 
		  const Vector<p3vector<double> >& grad_phi,
		  Matrix<p3vector<complex<double> > >& grad_f ) ;


void gradlambda( const Vector<double>& thetas, const Vector<double>& phis,
		 const Vector<double>& mm, const Vector<double>& nn, 
		 const Vector<complex<double> > lambdaF,
		 Vector<double>& dlambda_dtheta, 
		 Vector<double>& dlambda_dphi);

