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
		  const Vector<p3vector<double> > & B,
		  const Matrix<p3vector<complex<double> > > & grad_f,
		  const Matrix<complex<double> > & f_del,
		  const Vector<p3vector<double> > dx_dtheta,
		  const Vector<p3vector<double> > dx_dphi);

// this func has same output, but uses a different computational method. 
// useful for checking correctness of program but *much* (7 times) slower!
void omegamatrixB( Matrix<complex<double> > & omega,
		   const Vector<p3vector<double> > & B,
		   const Matrix<p3vector<complex<double> > > & grad_f,
		   const Matrix<complex<double> > & f_del,
		   const Vector<p3vector<double> > dx_dr,
		   const Vector<p3vector<double> > dx_dtheta,
		   const Vector<p3vector<double> > dx_dphi,
		   const Vector<p3vector<double> > grad_r);

void new_omegamatrix( Matrix<complex<double> > & omega,
		      const Vector<p3vector<double> > & B,
		      const Vector<double> & weight,
		      const Matrix<p3vector<complex<double> > > & grad_fdel,
		      const Matrix<complex<double> > & f);

void other_omegamatrix( Matrix<complex<double> > & omega,
			const Vector<p3vector<double> > & B,
			const Matrix<p3vector<complex<double> > > & grad_fdel,
			const Matrix<complex<double> > & f,
			const Vector<p3vector<double> > dx_dtheta,
			const Vector<p3vector<double> > dx_dphi);



void NormalizedOmegamatrix( Matrix<complex<double> > & omegaN,
			    const Vector<p3vector<double> > & B,
			    const Matrix<complex<double> > & fs,
			    const Matrix<p3vector<complex<double> > > & grad_f,
			    const Vector<double> Bphi);


void NormalizedOmegamatrixB( Matrix<complex<double> > & omegaN,
			    const Vector<p3vector<double> > & B,
			    const Matrix<complex<double> > & fs,
			    const Matrix<p3vector<complex<double> > > & grad_f,
			    const Vector<p3vector<double> > grad_phi);


void OmegaN_PF_from_lambda( Matrix<complex<double> > & PF,
			       const  Vector<double> & lambda,
			       const Matrix<complex<double> > & fs,
			       const Vector<double> & nn, 
			       const Vector<double> & mm,
			       const Vector<double>& thetas, 
			       const Vector<double>& phis
			    );

void OmegaN_PFinv_from_lambda( Matrix<complex<double> > & PFinv,
			       const  Vector<double> & lambda,
			       const  Vector<double> & dlambda_dtheta,
			       const Matrix<complex<double> > & fs,
			       const Vector<double> & nn, 
			       const Vector<double> & mm,
			       const Vector<double>& thetas, 
			       const Vector<double>& phis
			       );


void divfree_omegamatrix( Matrix<complex<double> > & omega,
			  const Vector<complex<double> > & JBtF,
			  const Vector<complex<double> > & JBpF,
			  const Vector<double> & nn,
			  const Vector<double> & mm);


void lambda_omegamatrix( Matrix<complex<double> > & omega,
			 const Vector<complex<double> > & lambdaF,
			 const Vector<double> & nn,
			 const Vector<double> & mm);


void lambda_omegamatrix( Matrix<complex<double> > & omega,
			 const Vector<complex<double> > & lambdaF,
			 const Vector<double> & nnR,
			 const Vector<double> & mmR,
			 const double fluxshear,
			 const double iota);



// calculates Fourier to Magnetic Fourier transformation using lambda

inline void F2MAGF_from_lambda( Matrix<complex<double> > & PF,
			    const  Vector<double> & lambda,
			    const Matrix<complex<double> > & fs,
			    const Vector<double> & nn, 
			    const Vector<double> & mm,
			    const Vector<double>& thetas, 
			    const Vector<double>& phis
			    )
{
   OmegaN_PF_from_lambda( PF, lambda, fs, nn, mm, thetas, phis);

}
