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

void gradfvector( const LAvector<double>& thetas, const LAvector<double>& phis,
		  const LAvector<double>& mm, const LAvector<double>& nn, 
		  const LAvector<p3vector<double> >& grad_theta, 
		  const LAvector<p3vector<double> >& grad_phi,
		  Matrix<p3vector<complex<double> > >& grad_f ) ;


void gradlambda( const LAvector<double>& thetas, const LAvector<double>& phis,
		 const LAvector<double>& mm, const LAvector<double>& nn, 
		 const LAvector<complex<double> > lambdaF,
		 LAvector<double>& dlambda_dtheta, 
		 LAvector<double>& dlambda_dphi);

