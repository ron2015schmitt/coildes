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

#include "omegamatrix.hpp"



void omegamatrix( Matrix<complex<double> > & omega,
		  const LAvector<p3vector<double> > & B,
		  const Matrix<p3vector<complex<double> > > & grad_f,
		  const Matrix<complex<double> > & f_del,
		  const LAvector<p3vector<double> > dx_dtheta,
		  const LAvector<p3vector<double> > dx_dphi){ 
  
   const unsigned int Npts = B.size();
   const unsigned int NF = grad_f.Ncols();
   const unsigned int No = f_del.Ncols();
   omega.resize(NF,No);
  
   const double dphi_dtheta = (2*PI)*(2*PI)/double(Npts);
  
   Matrix<complex<double> >  Bdotgradf(Npts,NF,"omegamatrix::Bdotgradf");
  
   LAvector<double> da(Npts,"omegamatrix::da");
   for (unsigned int i = 0; i<Npts ; i++) 
      da[i] = norm(cross(dx_dtheta[i],dx_dphi[i]))*dphi_dtheta;

   for (unsigned int k = 0; k<NF ; k++) {
      for (unsigned int i = 0; i<Npts ; i++) {
	 Bdotgradf(i,k) = (B[i]|grad_f(i,k))*da[i];
	 //      disp(Bdotgradf(i,k));disp(B[i]);dispcr(grad_f(i,k));
      }
   }


   // perform the integration

   omega = -(adj(Bdotgradf)|f_del);
  
}



// different implementation, but gives same result.

void omegamatrixB( Matrix<complex<double> > & omega,
		   const LAvector<p3vector<double> > & B,
		   const Matrix<p3vector<complex<double> > > & grad_f,
		   const Matrix<complex<double> > & f_del,
		   const LAvector<p3vector<double> > dx_dr,
		   const LAvector<p3vector<double> > dx_dtheta,
		   const LAvector<p3vector<double> > dx_dphi,
		   const LAvector<p3vector<double> > grad_r)
{ 

   const unsigned int Npts = B.size();
   const unsigned int NF = grad_f.Ncols();
   const unsigned int No = f_del.Ncols();
   const double dphi_dtheta= (2*PI)*(2*PI)/Npts;
   omega.resize(NF,No);
  
   LAvector<double> J(Npts,"J");
   for (unsigned int i =0; i<Npts; i++)
      J[i] = dot(dx_dr[i], cross(dx_dtheta[i],dx_dphi[i]));

   for(unsigned int k=0; k<NF; k++){  
      for(unsigned int l=0; l<No; l++){
	 std::complex<double> temp = 0;
	 for(unsigned int i=0; i<Npts; i++){
	    const double da = J[i]*norm(grad_r[i])*dphi_dtheta;
	    temp = temp  - conj(B[i]|grad_f(i,k)) * f_del(i,l) * da;
	 }
	 omega(k,l) = temp;
      }
   }
}



void new_omegamatrix( Matrix<complex<double> > & omega,
		      const LAvector<p3vector<double> > & B,
		      const LAvector<double> & weight,
		      const Matrix<p3vector<complex<double> > > & grad_fdel,
		      const Matrix<complex<double> > & f)
{ 
  
   const unsigned int Npts = B.size();
   const unsigned int No = grad_fdel.Ncols();
   const unsigned int NF = f.Ncols();
   omega.resize(NF,No);
  
   const double dphi_dtheta = (2*PI)*(2*PI)/double(Npts);
  
   Matrix<complex<double> >  Bdotgradf(Npts,No,"new_omegamatrix::Bdotgradf");

   LAvector<double> weight2(Npts,"new_omegamatrix::weight2");
  
   for (unsigned int i = 0; i<Npts ; i++) 
      weight2[i] = weight[i] * dphi_dtheta;

   for (unsigned int k = 0; k<No ; k++) {
      for (unsigned int i = 0; i<Npts ; i++) {
	 Bdotgradf(i,k) = (B[i]|grad_fdel(i,k))*weight2[i];
      }
   }

   // perform the integration
   omega = (adj(f)|Bdotgradf);
  
}



void other_omegamatrix( Matrix<complex<double> > & omega,
			const LAvector<p3vector<double> > & B,
			const Matrix<p3vector<complex<double> > > & grad_fdel,
			const Matrix<complex<double> > & f,
			const LAvector<p3vector<double> > dx_dtheta,
			const LAvector<p3vector<double> > dx_dphi)
{ 
  
   const unsigned int Npts = B.size();
   const unsigned int No = grad_fdel.Ncols();
   const unsigned int NF = f.Ncols();
   omega.resize(No,NF);
  
   const double dphi_dtheta = (2*PI)*(2*PI)/double(Npts);
  
   Matrix<complex<double> >  Bdotgradf(Npts,NF,"omegamatrix::Bdotgradf");
  
   LAvector<double> da(Npts,"omegamatrix::da");
   for (unsigned int i = 0; i<Npts ; i++) 
      da[i] = norm(cross(dx_dtheta[i],dx_dphi[i]))*dphi_dtheta;

   for (unsigned int k = 0; k<NF ; k++) {
      for (unsigned int i = 0; i<Npts ; i++) {
	 Bdotgradf(i,k) = (B[i]|grad_fdel(i,k))*da[i];
	 //      disp(Bdotgradf(i,k));disp(B[i]);dispcr(grad_f(i,k));
      }
   }


   // perform the integration

   omega = (adj(f)|Bdotgradf);
  
}









//  Normalized Omega
void NormalizedOmegamatrix( Matrix<complex<double> > & omegaN,
			    const LAvector<p3vector<double> > & B,
			    const Matrix<complex<double> > & fs,
			    const Matrix<p3vector<complex<double> > > & grad_f,
			    const LAvector<double> Bphi)
{ 
  
   const unsigned int Npts = B.size();
   const unsigned int NF = fs.Ncols();
   const unsigned int NFgrad = grad_f.Ncols();
   omegaN.resize(NF,NFgrad);
  
   const double dphi_dtheta= (2*PI)*(2*PI)/Npts;

   Matrix<complex<double> >  Bdotgradf(Npts,NFgrad,"omegamatrixN::Bdotgradf");
  
   for (unsigned int i = 0; i<Npts ; i++) {
      const double temp = dphi_dtheta / Bphi[i];
      const p3vector<double> vtemp = B[i] * temp;
      for (unsigned int k = 0; k<NFgrad; k++) {
	 Bdotgradf(i,k) = (vtemp|grad_f(i,k));
      }
   }

   // perform the integration over (theta,phi) space
   omegaN = (adj(fs)|Bdotgradf);
  
}





//  Normalized Omega
// using one loop.  this is actually slower.  COOLL matrix multiply is optimized,whereas
// this is not.
void NormalizedOmegamatrixB( Matrix<complex<double> > & omegaN,
			     const LAvector<p3vector<double> > & B,
			     const Matrix<complex<double> > & fs,
			     const Matrix<p3vector<complex<double> > > & grad_f,
			     const LAvector<p3vector<double> > grad_phi)
{ 
  
   const unsigned int Npts = B.size();
   const unsigned int NF = fs.Ncols();
   const unsigned int NFgrad = grad_f.Ncols();
   omegaN.resize(NF,NFgrad);
  
   const double dphi_dtheta= (2*PI)*(2*PI)/Npts;

   omegaN = 0;
   for (unsigned int r = 0; r<Npts ; r++) {
      for (unsigned int i = 0; i<Npts ; i++) {
	 const double temp = dphi_dtheta / (dot(B[i],grad_phi[i]));
	 const p3vector<double> vtemp = B[i] * temp;
	 for (unsigned int k = 0; k<NF ; k++) {
	    omegaN(r,k) += conj(fs(r,i)) * (vtemp|grad_f(i,k));
	 }
      }
   }
 
}



void OmegaN_PF_from_lambda( Matrix<complex<double> > & PF,
			    const  LAvector<double> & lambda,
			    const Matrix<complex<double> > & fs,
			    const LAvector<double> & nn, 
			    const LAvector<double> & mm,
			    const LAvector<double>& thetas, 
			    const LAvector<double>& phis
			    )
{ 
  
   const unsigned int Npts = fs.Nrows();
   const unsigned int NF = fs.Ncols();
   PF.resize(NF,NF);
  
   const double dphi_dtheta= (2*PI)*(2*PI)/Npts;
   const double coef = 1.0/(2.0*PI)*dphi_dtheta;

   Matrix<complex<double> >  FuncF(Npts,NF,"OmegaN_PF_from_lambda::FuncF");
  
   for (unsigned int k = 0; k<NF ; k++) {
      LAvector<double> phase(Npts,"phase");
      phase= nn[k]*phis + mm[k]*(thetas+lambda);
      for (unsigned int i = 0; i<Npts ; i++) {
	 FuncF(i,k) = std::complex<double>(coef*cos(phase[i]),coef*sin(phase[i])); 
      }
   }

   // perform the integration over (theta,phi) space
   PF = (adj(fs)|FuncF);
  
}


void OmegaN_PFinv_from_lambda( Matrix<complex<double> > & PFinv,
			       const  LAvector<double> & lambda,
			       const  LAvector<double> & dlambda_dtheta,
			       const Matrix<complex<double> > & fs,
			       const LAvector<double> & nn, 
			       const LAvector<double> & mm,
			       const LAvector<double>& thetas, 
			       const LAvector<double>& phis
			       )
{ 
  
   const unsigned int Npts = fs.Nrows();
   const unsigned int NF = fs.Ncols();
   PFinv.resize(NF,NF);
  
   const double dphi_dtheta= (2*PI)*(2*PI)/Npts;
   const double coef = 1.0/(2.0*PI)*dphi_dtheta;

   Matrix<complex<double> >  FuncF(NF,Npts,"OmegaN_PFinv_from_lambda::FuncF");
  
   for (unsigned int k = 0; k<NF ; k++) {
      LAvector<double> phase(Npts,"phase");
      phase= -nn[k]*phis - mm[k]*(thetas+lambda);
      for (unsigned int i = 0; i<Npts ; i++) {
	 FuncF(k,i) = std::complex<double>(coef*cos(phase[i]),coef*sin(phase[i]))*(1+dlambda_dtheta[i]); 
      }

   }

   // perform the integration over (theta,phi) space
   PFinv = (FuncF|fs);
  
}





void divfree_omegamatrix( Matrix<complex<double> > & omega,
			  const LAvector<complex<double> > & JBtF,
			  const LAvector<complex<double> > & JBpF,
			  const LAvector<double> & nn,
			  const LAvector<double> & mm)
{ 
  
   const unsigned int NF = nn.size();
   const unsigned int NFR = NF-1;
   omega.resize(NFR,NFR);
   const std::complex<double> i =  std::complex<double>(0,1);
  
   const LAvector<unsigned int> tmp = findtrue((nn==0)&&(mm==0));
   const unsigned int ind00 = tmp[0];
  
  
   for (unsigned int ll = 0; ll<NF ; ll++) {
      for (unsigned int kk = 0; kk<NF ; kk++) {
	 const double deln = nn[kk] - nn[ll];
	 const double delm = mm[kk] - mm[ll];
	 const LAvector<unsigned int> inds = findtrue((nn==deln)&&(mm==delm));
	 const unsigned int ind = inds[0];

	 const unsigned int r = (kk<=ind00) ? kk : (kk-1);
	 const unsigned int c = (ll<=ind00) ? ll : (ll-1);


	 if ((kk==ind00)||(ll==ind00)) {
	    ; // omega omits the 00 fourier terms for rows and columns
	 } else if (c>r) {
	    omega(r,c) = -conj(omega(c,r)); //anti-hermitian filling of upper diagonal
	 } else if (kk==ll) { //(deln==0) && (delm==0)
	    omega(r,c) = i/(2*PI)*(nn[kk]*JBpF[ind] + mm[kk]*JBtF[ind]);
	 } else if (deln==0) {
	    omega(r,c) = i/(2*PI)*nn[kk]*JBpF[ind];
	 } else if (delm==0) {
	    omega(r,c) = i/(2*PI)*mm[kk]*JBtF[ind];
	 } else { //(deln!=0) && (delm!=0)
	    omega(r,c) = i/(2*PI)*0.5*(mm[kk]*nn[ll]-mm[ll]*nn[kk])*(JBpF[ind]/delm - JBtF[ind]/deln);
	 }
      }
   }
  

}





void lambda_omegamatrix( Matrix<complex<double> > & omega,
			 const LAvector<complex<double> > & lambdaF,
			 const LAvector<double> & nn,
			 const LAvector<double> & mm)
{ 
  
   const unsigned int NF = nn.size();
   const unsigned int NFR = NF-1;
   omega.resize(NFR,NFR);
   Matrix<complex<double> >  temp(NF,NF,"lambda_omegamatrix::temp");

   const std::complex<double> i =  std::complex<double>(0,1);
  
   const LAvector<unsigned int> tmp = findtrue((nn==0)&&(mm==0));
   const unsigned int ind00 = tmp[0];


   const std::complex<double> dpsi_p = lambdaF[ind00].real()/(2*PI);
   const std::complex<double> dpsi_t  = lambdaF[ind00].imag()/(2*PI);
  
  
   for (unsigned int c = 0; c<NF ; c++) {
      for (unsigned int r = 0; r<NF ; r++) {
	 const double deln = nn[r] - nn[c];
	 const double delm = mm[r] - mm[c];
	 const LAvector<unsigned int> inds = findtrue((nn==deln)&&(mm==delm));
	 const unsigned int ind = inds[0];

	 if ((r==ind00)||(c==ind00)) {
	    temp(r,c) = 0;
	 } else if (c>r) {
	    temp(r,c) = -conj(temp(c,r)); //anti-hermitian filling of upper diagonal
	 } else if (r==c) { 
	    temp(r,c) = i/(2*PI)*(dpsi_t*nn[r] + dpsi_p*mm[r]);
	 } else { 
	    temp(r,c) = dpsi_t*(1/(2*PI*2*PI))*(mm[c]*nn[r]-mm[r]*nn[c])*lambdaF[ind];
	 }
      }
   }
  

   for (unsigned int c = 0; c<NF ; c++) {
      for (unsigned int r = 0; r<NF ; r++) {
	 unsigned int r2 = (r<=ind00) ? r : (r-1);
	 unsigned int c2 = (c<=ind00) ? c : (c-1);
	 omega(r2,c2) = temp(r,c);
      }
   }

}












void lambda_omegamatrix( Matrix<complex<double> > & omega,
			 const LAvector<complex<double> > & lambdaF,
			 const LAvector<double> & nnR,
			 const LAvector<double> & mmR,
			 const double fluxshear,
			 const double iota)
{ 
  
   const unsigned int NFR = nnR.size();
   omega.resize(NFR,NFR);

   const std::complex<double> i =  std::complex<double>(0,1);
   const double C = 1/(2.0*PI*2.0*PI)*fluxshear;
  

   for (unsigned int c = 0; c<NFR ; c++) {
      for (unsigned int r = 0; r<NFR ; r++) {

	 if (c>r) {
	    omega(r,c) = -conj(omega(c,r)); //anti-hermitian filling of upper diagonal
	 } else if (r==c) { 
	    omega(r,c) = i*(C*(nnR[r] + iota*mmR[r]));
	 } else { 
    	   const double deln = nnR[r] - nnR[c];
	   const double delm = mmR[r] - mmR[c];
       	   const unsigned int ind = find1sttrue((nnR==deln)&&(mmR==delm));
	   if (ind >= NFR) 
              omega(r,c) = 0; 
	   else 
	     omega(r,c) = C*(mmR[c]*nnR[r]-mmR[r]*nnR[c])*lambdaF[ind];
	 }
      }
   }
}




