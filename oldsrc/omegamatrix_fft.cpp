/************************************************************************* 
 * 
 *   File Name    :  omegagamma.cpp
 *   Platform     :  gcc compiler (v3.2.2)
 *   Author       :  Ron Schmitt
 *   Date         :  9 jul 2004
 * 
 *
 *   SYNOPSIS     
 *     Generates omega matrix (tangent b-field matrix)
 *
 **************************************************************************/




// Standard C libraries
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>

// Standard C++ libraries

#include <iostream>

using namespace std;


// Coils libraries

#include "coils.hpp"
#include "surface.hpp"
#include "coilfft.hpp"

#include "omegamatrix_fft.hpp"

inline double dummy(const double a)
{
return a;
}


// use weight[i] = J[i] * norm(grad_r[i]) to get same results as omegamatrix();
// can *NOT*  use reduced set grad_fdel

void omegamatrix_fft( Matrix<complex<double> > & omega,
		      const Vector<p3vector<double> > & B,
		      const Vector<double> & weight,
		      const Matrix<p3vector<complex<double> > > & grad_f,
		      const unsigned int Nphi, const unsigned int Ntheta,
		      const unsigned int Nnn, const unsigned int Nmm,
		      const unsigned int Nfund, const unsigned int Mfund)
{ 
  
   const unsigned int Npts = B.size();
   const unsigned int NFgrad = grad_f.Ncols();
  
   //  const double dphi_dtheta= (2*PI)*(2*PI)/Npts;

   Matrix<complex<double> >  cBdotgradf(NFgrad,Npts,"omegamatrixW::cBdotgradf");
  
   for (unsigned int i = 0; i<Npts ; i++) {
      const p3vector<double> vtemp = weight[i]*B[i];
      for (unsigned int k = 0; k<NFgrad ; k++) {
	 cBdotgradf(k,i) = -conj(vtemp|grad_f(i,k));
      }
   }

   // perform the integration over (theta,phi) space, using FFT
   omega_half_fft(cBdotgradf, omega, Nphi, Ntheta,Nnn,Nmm, Nfund, Mfund, 1e-14);

   omega.adjoint();
   omega = -omega;
  


}


// use weight[i] = J[i] * norm(grad_r[i]) to get same results as new_omegamatrix();
// can use reduced set grad_fdel

void new_omegamatrix_fft( Matrix<complex<double> > & omega,
			  const Vector<p3vector<double> > & B,
			  const Vector<double> & weight,
			  const Matrix<p3vector<complex<double> > > & grad_fdel,
			  const unsigned int Nphi, const unsigned int Ntheta,
			  const unsigned int Nnn, const unsigned int Nmm,
			  const unsigned int Nfund, const unsigned int Mfund)
{ 
  
   const unsigned int Npts = B.size();
   const unsigned int No = grad_fdel.Ncols();
  
   //  const double dphi_dtheta= (2*PI)*(2*PI)/Npts;

   Matrix<complex<double> >  cBdotgradf(No,Npts,"new_omegamatrix_fft::cBdotgradf");
  
   for (unsigned int i = 0; i<Npts ; i++) {
      const p3vector<double> vtemp = weight[i]*B[i];
      for (unsigned int k = 0; k<No ; k++) {
	 cBdotgradf(k,i) = -conj(vtemp|grad_fdel(i,k));
      }
   }

   // perform the integration over (theta,phi) space, using FFT
   omega_half_fft(cBdotgradf, omega, Nphi, Ntheta,Nnn,Nmm, Nfund, Mfund, 1e-14);

  
}


// left out 1/(2*pi) factor here also
void new_omega_from_lambda_fft( Matrix<complex<double> > & omega,
				const Vector<double> & dlambda_dtheta,
				const Vector<double> &  dlambda_dphi,
				const Matrix<complex<double> > & fsR,
				const Vector<double> & nnR, 
				const Vector<double> & mmR,
				const unsigned int Nphi, const unsigned int Ntheta,
				const unsigned int Nnn, const unsigned int Nmm,
				const unsigned int Nfund, const unsigned int Mfund,
				const double iota, const double fluxshear)
{ 
  
   const unsigned int Npts = fsR.Nrows();
   const unsigned int NFR = fsR.Ncols();
  
   const complex<double> C0 = -complex<double>(0,1)*fluxshear;


   Matrix<complex<double> >  integrand(NFR,Npts,"new_omega_from_lambda_fft::integrand");
  
   for (unsigned int i = 0; i<Npts ; i++) {
      const complex<double> termA = C0*(iota -  dlambda_dphi[i]);
      const complex<double> termB = C0*(1 +  dlambda_dtheta[i]);
      for (unsigned int k = 0; k<NFR ; k++) {
	 integrand(k,i) = conj(fsR(i,k))*(mmR[k]*termA + nnR[k]*termB);
      }
   }

   // perform the integration over (theta,phi) space, using FFT
   omega_half_fft(integrand, omega, Nphi, Ntheta,Nnn,Nmm, Nfund, Mfund, 1e-14);

  
}


// left out 1/(2*pi) factor here also
// this is new__omega_from_lambda_fft, but divided by [fluxshear*(1+dlambda_dtheta)]
void modified_omega_from_lambda_fft( Matrix<complex<double> > & omega,
				     const Vector<double> & dlambda_dtheta,
				     const Vector<double> &  dlambda_dphi,
				     const Matrix<complex<double> > & fsR,
				     const Vector<double> & nnR, 
				     const Vector<double> & mmR,
				     const unsigned int Nphi, const unsigned int Ntheta,
				     const unsigned int Nnn, const unsigned int Nmm,
				     const unsigned int Nfund, const unsigned int Mfund,
				     const double iota, const double fluxshear)
{ 
  
   const unsigned int Npts = fsR.Nrows();
   const unsigned int NFR = fsR.Ncols();
  
   const complex<double> C0 = -complex<double>(0,1);//*fluxshear;

   Matrix<complex<double> >  integrand(NFR,Npts,"new_omega_from_lambda_fft::integrand");

   dummy(fluxshear);
  
   for (unsigned int i = 0; i<Npts ; i++) {
      const complex<double> termA = (iota -  dlambda_dphi[i]);
      const complex<double> termB = (1 +  dlambda_dtheta[i]);
      const complex<double> termC = termA/termB;
      for (unsigned int k = 0; k<NFR ; k++) {
	 integrand(k,i) = conj(fsR(i,k))*(mmR[k]*termC + nnR[k])*C0;
      }
   }

   // perform the integration over (theta,phi) space, using FFT
   omega_half_fft(integrand, omega, Nphi, Ntheta,Nnn,Nmm, Nfund, Mfund, 1e-14);

  
}






//  Normalized Omega

// gives an omega matrix whose eigenvalues are i*(n+iota*m)
// without factor of 1/(2*pi)
void NormalizedOmegamatrix_fft( Matrix<complex<double> > & omegaN,
				const Vector<p3vector<double> > & B,
				const Matrix<p3vector<complex<double> > > & grad_f,
				const Vector<double>& Bp,
				const unsigned int Nphi, const unsigned int Ntheta,
				const unsigned int Nnn, const unsigned int Nmm,
				const unsigned int Nfund, const unsigned int Mfund)
{ 
  
   const unsigned int Npts = B.size();
   const unsigned int NFgrad = grad_f.Ncols();
  
   //  const double dphi_dtheta= (2*PI)*(2*PI)/Npts;

   Matrix<complex<double> >  cBdotgradf(NFgrad,Npts,"omegamatrixW::cBdotgradf");
  
   for (unsigned int i = 0; i<Npts ; i++) {
      const p3vector<double> vtemp = B[i] / Bp[i];
      for (unsigned int k = 0; k<NFgrad ; k++) {
	 cBdotgradf(k,i) = -conj(vtemp|grad_f(i,k));
      }
   }

   // perform the integration over (theta,phi) space, using FFT
   omega_half_fft(cBdotgradf, omegaN, Nphi, Ntheta,Nnn,Nmm, Nfund, Mfund, 1e-14);

  
  
}



void NormalizedOmegamatrix_fft_fast( Matrix<complex<double> > & omegaN,
				     const Vector<double> & Bt,
				     const Vector<double>& Bp,
				     const Matrix<complex<double> > & fsR,
				     const Vector<double>& nnR,  const Vector<double>& mmR,
				     const unsigned int Nphi, const unsigned int Ntheta,
				     const unsigned int Nnn, const unsigned int Nmm,
				     const unsigned int Nfund, const unsigned int Mfund)
{ 
  
   const unsigned int Npts = fsR.Nrows();
   const unsigned int NFR = fsR.Ncols();
  
   //  const double dphi_dtheta= (2*PI)*(2*PI)/Npts;


  
   Matrix<complex<double> >  cBdotgradf(NFR,Npts,"omegamatrixW::cBdotgradf");
  
   for (unsigned int i = 0; i<Npts ; i++) {
      const double temp1 =   Bt[i] / Bp[i];
      for (unsigned int k = 0; k<NFR ; k++) {
	 const double temp2 = nnR[k] + mmR[k]*temp1;
	 cBdotgradf(k,i) = std::complex<double>(0,1)*conj(fsR(i,k)) * temp2;
      }
   }

   // perform the integration over (theta,phi) space, using FFT
   omega_half_fft(cBdotgradf, omegaN, Nphi, Ntheta,Nnn,Nmm, Nfund, Mfund, 1e-14);


}





void OmegaN_PF_from_lambda_fft( Matrix<complex<double> > & PF,
			    const  Vector<double> & lambda,
			    const Vector<double> & nn, 
			    const Vector<double> & mm,
			    const Vector<double>& thetas, 
			    const Vector<double>& phis,
			    const unsigned int Nphi, const unsigned int Ntheta,
		            const unsigned int Nnn, const unsigned int Nmm,
			    const unsigned int Nfund, const unsigned int Mfund)
{ 
  
   const unsigned int Npts = thetas.size();
   const unsigned int NF = mm.size();
  
//   const double dphi_dtheta= (2*PI)*(2*PI)/Npts;
//   const double coef = 1.0/(2.0*PI)*dphi_dtheta;

   Matrix<complex<double> >  FuncF(Npts,NF,"OmegaN_PF_from_lambda::FuncF");
  
   for (unsigned int k = 0; k<NF ; k++) {
      Vector<double> phase(Npts,"phase");
      phase= nn[k]*phis + mm[k]*(thetas+lambda);
      for (unsigned int i = 0; i<Npts ; i++) {
	 FuncF(i,k) = std::complex<double>(cos(phase[i]),sin(phase[i])); 
      }
   }

omega_half_fft2(FuncF, PF, Nphi, Ntheta, Nnn, Nmm, Nfund, Mfund,1e-10);
  
}
