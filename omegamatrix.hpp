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




void omegamatrix( Matrix<complex<double> > & omega,
		  const LAvector<p3vector<double> > & B,
		  const Matrix<p3vector<complex<double> > > & grad_f,
		  const Matrix<complex<double> > & f_del,
		  const LAvector<p3vector<double> > dx_dtheta,
		  const LAvector<p3vector<double> > dx_dphi);

// this func has same output, but uses a different computational method. 
// useful for checking correctness of program but *much* (7 times) slower!
void omegamatrixB( Matrix<complex<double> > & omega,
		   const LAvector<p3vector<double> > & B,
		   const Matrix<p3vector<complex<double> > > & grad_f,
		   const Matrix<complex<double> > & f_del,
		   const LAvector<p3vector<double> > dx_dr,
		   const LAvector<p3vector<double> > dx_dtheta,
		   const LAvector<p3vector<double> > dx_dphi,
		   const LAvector<p3vector<double> > grad_r);

void new_omegamatrix( Matrix<complex<double> > & omega,
		      const LAvector<p3vector<double> > & B,
		      const LAvector<double> & weight,
		      const Matrix<p3vector<complex<double> > > & grad_fdel,
		      const Matrix<complex<double> > & f);

void other_omegamatrix( Matrix<complex<double> > & omega,
			const LAvector<p3vector<double> > & B,
			const Matrix<p3vector<complex<double> > > & grad_fdel,
			const Matrix<complex<double> > & f,
			const LAvector<p3vector<double> > dx_dtheta,
			const LAvector<p3vector<double> > dx_dphi);



void NormalizedOmegamatrix( Matrix<complex<double> > & omegaN,
			    const LAvector<p3vector<double> > & B,
			    const Matrix<complex<double> > & fs,
			    const Matrix<p3vector<complex<double> > > & grad_f,
			    const LAvector<double> Bphi);


void NormalizedOmegamatrixB( Matrix<complex<double> > & omegaN,
			    const LAvector<p3vector<double> > & B,
			    const Matrix<complex<double> > & fs,
			    const Matrix<p3vector<complex<double> > > & grad_f,
			    const LAvector<p3vector<double> > grad_phi);


void OmegaN_PF_from_lambda( Matrix<complex<double> > & PF,
			       const  LAvector<double> & lambda,
			       const Matrix<complex<double> > & fs,
			       const LAvector<double> & nn, 
			       const LAvector<double> & mm,
			       const LAvector<double>& thetas, 
			       const LAvector<double>& phis
			    );

void OmegaN_PFinv_from_lambda( Matrix<complex<double> > & PFinv,
			       const  LAvector<double> & lambda,
			       const  LAvector<double> & dlambda_dtheta,
			       const Matrix<complex<double> > & fs,
			       const LAvector<double> & nn, 
			       const LAvector<double> & mm,
			       const LAvector<double>& thetas, 
			       const LAvector<double>& phis
			       );


void divfree_omegamatrix( Matrix<complex<double> > & omega,
			  const LAvector<complex<double> > & JBtF,
			  const LAvector<complex<double> > & JBpF,
			  const LAvector<double> & nn,
			  const LAvector<double> & mm);


void lambda_omegamatrix( Matrix<complex<double> > & omega,
			 const LAvector<complex<double> > & lambdaF,
			 const LAvector<double> & nn,
			 const LAvector<double> & mm);


void lambda_omegamatrix( Matrix<complex<double> > & omega,
			 const LAvector<complex<double> > & lambdaF,
			 const LAvector<double> & nnR,
			 const LAvector<double> & mmR,
			 const double fluxshear,
			 const double iota);



// calculates Fourier to Magnetic Fourier transformation using lambda

inline void F2MAGF_from_lambda( Matrix<complex<double> > & PF,
			    const  LAvector<double> & lambda,
			    const Matrix<complex<double> > & fs,
			    const LAvector<double> & nn, 
			    const LAvector<double> & mm,
			    const LAvector<double>& thetas, 
			    const LAvector<double>& phis
			    )
{
   OmegaN_PF_from_lambda( PF, lambda, fs, nn, mm, thetas, phis);

}
