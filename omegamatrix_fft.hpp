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






// use weight[i] = J[i] * norm(grad_r[i]) to get same results as omegamatrix();
// can *NOT*  use reduced set grad_fdel

void omegamatrix_fft( Matrix<complex<double> > & omega,
		      const Vector<p3vector<double> > & B,
		      const Vector<double> & weight,
		      const Matrix<p3vector<complex<double> > > & grad_f,
		      const unsigned int Nphi, const unsigned int Ntheta,
		      const unsigned int Nnn, const unsigned int Nmm,
		      const unsigned int Nfund, const unsigned int Mfund);


// use weight[i] = J[i] * norm(grad_r[i]) to get same results as new_omegamatrix();
// can use reduced set grad_fdel

void new_omegamatrix_fft( Matrix<complex<double> > & omega,
		      const Vector<p3vector<double> > & B,
		      const Vector<double> & weight,
		      const Matrix<p3vector<complex<double> > > & grad_fdel,
		      const unsigned int Nphi, const unsigned int Ntheta,
		      const unsigned int Nnn, const unsigned int Nmm,
			  const unsigned int Nfund, const unsigned int Mfund);
 

void new_omega_from_lambda_fft( Matrix<complex<double> > & omega,
				const Vector<double> & dlambda_dtheta,
				const Vector<double> &  dlambda_dphi,
				const Matrix<complex<double> > & fsR,
				const Vector<double> & nnR, 
				const Vector<double> & mmR,
				const unsigned int Nphi, const unsigned int Ntheta,
				const unsigned int Nnn, const unsigned int Nmm,
				const unsigned int Nfund, const unsigned int Mfund,
				const double iota, const double fluxshear);
  
void modified_omega_from_lambda_fft( Matrix<complex<double> > & omega,
				const Vector<double> & dlambda_dtheta,
				const Vector<double> &  dlambda_dphi,
				const Matrix<complex<double> > & fsR,
				const Vector<double> & nnR, 
				const Vector<double> & mmR,
				const unsigned int Nphi, const unsigned int Ntheta,
				const unsigned int Nnn, const unsigned int Nmm,
				const unsigned int Nfund, const unsigned int Mfund,
				const double iota, const double fluxshear);




void OmegaN_PF_from_lambda_fft( Matrix<complex<double> > & PF,
			    const  Vector<double> & lambda,
			    const Vector<double> & nn, 
			    const Vector<double> & mm,
			    const Vector<double>& thetas, 
			    const Vector<double>& phis,
			    const unsigned int Nphi, const unsigned int Ntheta,
		            const unsigned int Nnn, const unsigned int Nmm,
			    const unsigned int Nfund, const unsigned int Mfund);

void NormalizedOmegamatrix_fft( Matrix<complex<double> > & omegaN,
				const Vector<p3vector<double> > & B,
				const Matrix<p3vector<complex<double> > > & grad_f,
				const Vector<double>& Bp,
				const unsigned int Nphi, const unsigned int Ntheta,
				const unsigned int Nnn, const unsigned int Nmm,
				const unsigned int Nfund, const unsigned int Mfund);

void NormalizedOmegamatrix_fft_fast( Matrix<complex<double> > & omegaN,
				const Vector<double> & Bt,
				const Vector<double>& Bp,
				const Matrix<complex<double> > & fsR,
				const Vector<double>& nnR,  const Vector<double>& mmR,
				const unsigned int Nphi, const unsigned int Ntheta,
				const unsigned int Nnn, const unsigned int Nmm,
				     const unsigned int Nfund, const unsigned int Mfund);
