/************************************************************************* 
 * 
 *   File Name    :  surface.cpp
 *   Platform     : gcc
 *   Author       :  Ron Schmitt
 *   Date         :  16 May 2003
 * 
 *
 *   SYNOPSIS     
 *     Generates plasma surface and surface normals from Fourier coefficients;
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

#include "matricks.hpp"

void anglevectors(Vector<double>& thetas, Vector<double>& phis, 
		  const unsigned int Ntheta, const unsigned int  Nphi)
{
   double dtheta = 2*M_PI/((double)Ntheta);
   double dphi = 2*M_PI/((double)Nphi);
  
   unsigned int i=0;
   for (unsigned int tt = 0; tt<Ntheta ; tt++){
      for (unsigned int pp = 0; pp<Nphi; pp++,i++) {
	 thetas[i] =  double(tt) * dtheta;
	 phis[i] = double(pp) * dphi;
      }
   }
  
}





void modevectors(unsigned int& NF,  Vector<double>& nn,  Vector<double>& mm,
		 const unsigned int n_max, const unsigned int m_max,
		 const bool mode00) 
{
   NF = (2*n_max+1)*(2*m_max+1);

   if (!mode00)
      NF = NF - 1;

   nn.resize(NF);
   mm.resize(NF);

   unsigned int k =0;
   for( int n=-n_max; n <= int(n_max); n++){
      for( int m = -m_max; m <= int(m_max); m++){
	 if ( !( (n==0) && (m==0) && !mode00 ) ) {
	    nn[k] = double(n);
	    mm[k] = double(m);
	    k++;
	 }
      }
   }
} 


void modevectors(unsigned int& NF,  Vector<double>& nn,  Vector<double>& mm,
		 unsigned int& n_max, unsigned int& m_max,
		 const unsigned int N_harmonics, const unsigned int M_harmonics,
		 const bool mode00) 
{

   unsigned int n_max_new = 0;
   unsigned int m_max_new = 0;


   for( int n=0; n <= int(n_max); n += N_harmonics){
      n_max_new = n;
   }

   if (n_max != n_max_new ) {
      n_max = n_max_new;
      std::cerr << "Max toroidal mode changed to multiple of harmonic (Nharm="<<N_harmonics<<"): Nmax ="<<n_max<<endl;
   }

   for( int m=0; m <= int(m_max); m += M_harmonics){
      m_max_new = m;
   }

   if (m_max != m_max_new ) {
      m_max = m_max_new;
      std::cerr << "Max poloidal mode changed to multiple of harmonic (Mharm="<<M_harmonics<<"): Mmax ="<<m_max<<endl;
   }

   NF = (2*n_max/N_harmonics+1)*(2*m_max/M_harmonics+1);

   if (!mode00)
      NF = NF - 1;

   nn.resize(NF);
   mm.resize(NF);


   unsigned int k =0;
   for( int n=-n_max; n <= int(n_max); n += N_harmonics){
      for( int m = -m_max; m <= int(m_max); m += M_harmonics){
	 if ( !( (n==0) && (m==0) && !mode00 ) ) {
	    nn[k] = double(n);
	    mm[k] = double(m);
	    k++;
	 }
      }
   }
} 



void modevectorsR(unsigned int& NF,  Vector<double>& nn,  Vector<double>& mm,
		  const unsigned int n_max, const unsigned int m_max) 
{
   NF = ((2*n_max+1)*(2*m_max+1))-1;

   nn.resize(NF);
   mm.resize(NF);

   unsigned int k =0;
   for( int n=-n_max; n <= int(n_max); n++){
      for( int m = -m_max; m <= int(m_max); m++){
	 if (!((n==0) && (m==0))) {
	    nn[k] = double(n);
	    mm[k] = double(m);
	    k++;
	 }
      }
   }
} 



// B >= A!!
// i.e. B must have at least all the modes that are in A!!!!

void make_mode_map(unsigned int& NF_A, const unsigned int& Nnn_A, const unsigned int& Nmm_A,
	      const unsigned int Ndelta_A, const unsigned int Mdelta_A,
	      const bool mode00_A,
	      unsigned int& NF_B, const unsigned int& Nnn_B, const unsigned int& Nmm_B,
	      const unsigned int Ndelta_B, const unsigned int Mdelta_B,
	      const bool mode00_B,
	      Vector<unsigned int>& Fmap) 
{

   unsigned int NnnA = Nnn_A;
   unsigned int NmmA = Nmm_A;
   unsigned int NnnB = Nnn_B;
   unsigned int NmmB = Nmm_B;
   


   Vector<double> nn_A("nn_A");
   Vector<double> mm_A("mm_A");
   modevectors(NF_A,nn_A,mm_A,NnnA,NmmA,Ndelta_A,Mdelta_A,mode00_A);
   Vector<double> nn_B("nn_B");
   Vector<double> mm_B("mm_B");
   modevectors(NF_B,nn_B,mm_B,NnnB,NmmB,Ndelta_B,Mdelta_B,mode00_B);

   //   disp(NF_A);disp(NnnA);disp(NmmA); disp(Ndelta_A);dispcr(Mdelta_A);
   //   disp(NF_B);disp(NnnB);disp(NmmB); disp(Ndelta_B);dispcr(Mdelta_B);

   Fmap.resize(NF_A);
   for (unsigned int kA=0; kA<NF_A; kA++) {
      double m_A = mm_A[kA];
      double n_A = nn_A[kA];
      int count = 0;
      for (unsigned int kB=0; kB<NF_B; kB++) {
	 if ( (m_A == mm_B[kB]) && (n_A == nn_B[kB]) ) {
	    Fmap[kA] = kB;
	    count++;
	 }
      }
      if (count !=1) {
	 printcr ("******BUG BUG BUG: bug occured during mapping.");
	 disp(count); disp(kA); cr();
	 disp(m_A);disp(n_A);cr();
	 disp(NF_A);disp(NF_B);
      }
   }


   //   disp(NF_A);disp(NnnA);disp(NmmA); disp(Ndelta_A);dispcr(Mdelta_A);
   //   disp(NF_B);disp(NnnB);disp(NmmB); disp(Ndelta_B);dispcr(Mdelta_B);

   //    dispcr(nn_A);
   //    dispcr(mm_A);
   //    dispcr(nn_B);
   //    dispcr(mm_B);
   //    disp(NF_A);disp(NF_B);
   //    dispcr(Fmap);
  
};



// generate set of complex *orthonormal* Fourier functions

void fseries(const Vector<double>& nn,  const Vector<double>& mm,
	     const Vector<double>& thetas, const Vector<double>& phis,
	     Matrix<complex<double> >& f)
{

   const unsigned int NF = mm.size();
 const unsigned int Npts = thetas.size();

   //fortho is the starting orthonormal basis of NF funstions
   f.resize(Npts,NF);

 
   for (unsigned int k = 0; k<NF ; k++) {
      Vector<double> phase(Npts,"phase");
      phase= nn[k]*phis + mm[k]*thetas;
      f.col(k) = vcomplex(cos(phase),sin(phase)); 
   }

   //normalize
   const double coef = 1.0/(2.0*PI);
   f=f*coef;
}


void remove_mode00(const Vector<complex<double> >& vF, Vector<complex<double> >& vFR) {

   const unsigned int NF = vF.size();
   vFR.resize(NF-1);
   const unsigned int mode00 = NF/2;

 
   for (unsigned int k = 0; k<NF ; k++) {
      vFR[(k<=mode00) ?  k : (k-1)] = vF[k];
   
   }
}



void remove_mode00(const Matrix<complex<double> >& fs, Matrix<complex<double> >& fsR) {

   const unsigned int Npts = fs.Nrows();
   const unsigned int NF = fs.Ncols();
   fsR.resize(Npts,NF-1);
   const unsigned int mode00 = NF/2;
 
   for (unsigned int i = 0; i<Npts ; i++) 
      for (unsigned int k = 0; k<NF ; k++) 
	 fsR( i, (k<=mode00) ?  k : (k-1) ) = fs(i,k);

}




double calcAspectRatio(const FourierSurface& fsurface) 
{


   const Vector<double>& nn = fsurface.nn();
   const Vector<double>& mm = fsurface.mm();
   const Vector<double>& RF = fsurface.RF();
   const Vector<double>& ZF = fsurface.ZF();
   const unsigned int Nmodes = mm.size();

   Vector<bool> mask(Nmodes,"mask");
   mask = (nn==0.0) && (mm==0.0);
   Vector<double> result(Nmodes,"result");
   result=RF[mask];

   const double Ro = result[0];


   mask = (nn==0.0) && (mm==1.0);
   result=RF[mask];
   const double ra = result[0];
   result=ZF[mask];
   const double rb = result[0];
   const double r0 = sqrt(ra*ra + rb*rb)/sqrt(2.);

   const double aspectRatio = Ro/r0;

   dispcr(Ro);
   dispcr(ra);
   dispcr(rb);
   dispcr(r0);
   dispcr(aspectRatio);


   return aspectRatio;

}




// **should add check for nyquist criterion (i.e. that Ntheta and Nphi are AT LEAST
// twice largest nn,mm and then not use larger nn,mm modes)


void expandsurfaceandbases (
 Vector<p3vector<double> >&X, Vector<p3vector<double> >&dA_dthetadphi,
 Vector<p3vector<double> >& dx_dr, Vector<p3vector<double> >& dx_dtheta, Vector<p3vector<double> >& dx_dphi, 
 Vector<p3vector<double> >& grad_r, Vector<p3vector<double> >& grad_theta, Vector<p3vector<double> >& grad_phi, 
 const FourierSurface& fsurface,
 const Vector<double>& thetas, const Vector<double>& phis 
 )
{

  
   const Vector<double>& nn = fsurface.nn();
   const Vector<double>& mm = fsurface.mm();
   const Vector<double>& RF = fsurface.RF();
   const Vector<double>& ZF = fsurface.ZF();
   const unsigned int Nmodes = mm.size();
   const unsigned int Npts = phis.size();
//   const double dphi_by_dtheta = 2.0*PI* 2.0*PI/double(Npts);

   Vector<double> angle(Nmodes,"angle");
   Vector<double> sinkern(Nmodes,"sinkern");
   Vector<double> coskern(Nmodes,"coskern"); 

 


   for (unsigned int ii = 0; ii<Npts;ii++){
      const double theta=thetas[ii];
      const double phi=phis[ii];
      
      // create sin,cos kernals for fourier transform
      angle = phi*nn+theta*mm;
      sinkern = sin(angle);
      coskern = cos(angle);
    
      const double R = (RF|coskern);
      X[ii].x() = R * cos(phi);
      X[ii].y() = R * sin(phi);
      X[ii].z() = (ZF|sinkern);
      

      // Now find the tangent vectors

      
      const double dRdTheta = -(RF|mm*sinkern);
      const double dRdPhi = -(RF|nn*sinkern);
      const double dZdTheta = (ZF|mm*coskern);
      const double dZdPhi = (ZF|nn*coskern);

      const p3vector<double> dxdtheta( dRdTheta*cos(phi), dRdTheta*sin(phi), dZdTheta );
      const p3vector<double> dxdphi( dRdPhi*cos(phi)- R*sin(phi), dRdPhi*sin(phi) + R*cos(phi), dZdPhi );
    
      const p3vector<double> dAdthetadphi = cross(dxdtheta, dxdphi);
      const p3vector<double> dxdr = dAdthetadphi/norm(dAdthetadphi);
    
      const double J = dot(dxdr,dAdthetadphi);
    
      const p3vector<double> gradr = cross(dxdtheta,dxdphi)/J;
      const p3vector<double> gradtheta = cross(dxdphi,dxdr)/J;
      const p3vector<double> gradphi = cross(dxdr,dxdtheta)/J;
    
      dA_dthetadphi[ii] = dAdthetadphi;
    
      dx_dr[ii] = dxdr;
      dx_dtheta[ii] = dxdtheta;
      dx_dphi[ii] = dxdphi;
    
      grad_r[ii] = gradr;
      grad_theta[ii] = gradtheta;
      grad_phi[ii] = gradphi;

      
   }
}


//*********************************************
//this function does not work roperly.  has not been fully debugged
//**********************************************

void perturbsurface(
		    Vector<p3vector<double> >&X, Vector<p3vector<double> >&dA_dthetadphi,
		    Vector<p3vector<double> >& dx_dr, Vector<p3vector<double> >& dx_dtheta, Vector<p3vector<double> >& dx_dphi, 
		    Vector<p3vector<double> >& grad_r, Vector<p3vector<double> >& grad_theta, Vector<p3vector<double> >& grad_phi, 
		    const unsigned int n, const unsigned int m,
		    const double delRF, const double delZF,
		    const Vector<double>& thetas, const Vector<double>& phis 
		    )
{
  
   const unsigned int Npts = phis.size();
//   const double dphi_by_dtheta = 2.0*PI* 2.0*PI/double(Npts);

   for (unsigned int ii = 0; ii<Npts;ii++){
      const double theta=thetas[ii];
      const double phi=phis[ii];
      
      // create sin,cos kernals for fourier transform
      double angle = phi*n+theta*m;
      double sinval = sin(angle);
      double cosval = cos(angle);
    
      const double delR = (delRF*cosval);
      const double delX = delR * cos(phi);
      const double delY = delR * sin(phi);
      X[ii].x() += delX; 
      X[ii].y() += delY; 
      X[ii].z() += (delZF*sinval);
      
      // Now find the tangent vectors
      
      const double dRdTheta = -(delRF*m*sinval);
      const double dRdPhi = -(delRF*n*sinval);
      const double dZdTheta = (delZF*m*cosval);
      const double dZdPhi = (delZF*n*cosval);

      const p3vector<double> delx_dtheta( dRdTheta*cos(phi), dRdTheta*sin(phi), dZdTheta );
      const p3vector<double> delx_dphi( dRdPhi*cos(phi)- delY, dRdPhi*sin(phi) + delX, dZdPhi );
    
      dx_dtheta[ii] = dx_dtheta[ii] + delx_dtheta;
      dx_dphi[ii] = dx_dphi[ii] + delx_dphi;
    
      dA_dthetadphi[ii] = cross(dx_dtheta[ii], dx_dphi[ii]);
      dx_dr[ii] = dA_dthetadphi[ii]/norm(dA_dthetadphi[ii]);
    
      const double J = dot(dx_dr[ii],dA_dthetadphi[ii]);
    
      grad_r[ii] = cross(dx_dtheta[ii],dx_dphi[ii])/J;
      grad_theta[ii] = cross(dx_dphi[ii],dx_dr[ii])/J;
      grad_phi[ii] = cross(dx_dr[ii],dx_dtheta[ii])/J;
    
      
   }
}
  






void expandsurface(
		   Vector<p3vector<double> >&X, Vector<p3vector<double> >&dA_dthetadphi,
		   const FourierSurface& fsurface,
		   const Vector<double>& thetas, const Vector<double>& phis 
		   ) {
   const unsigned int Npts = phis.size();

   Vector<p3vector<double> > dx_dr(Npts, "plasmasurface::dx_dr");
   Vector<p3vector<double> > dx_dtheta(Npts,"plasmasurface::dx_dtheta");
   Vector<p3vector<double> > dx_dphi(Npts,"plasmasurface::dx_dphi");
   Vector<p3vector<double> > grad_r(Npts,"plasmasurface::grad_r");
   Vector<p3vector<double> > grad_theta(Npts,"plasmasurface::grad_theta");
   Vector<p3vector<double> > grad_phi(Npts,"plasmasurface::grad_phi");

   expandsurfaceandbases(X,dA_dthetadphi,dx_dr,dx_dtheta,dx_dphi,grad_r,grad_theta,grad_phi,fsurface,thetas,phis);
}






void transformsurface(
		      const Vector<p3vector<double> >&X,
		      FourierSurface& fsurface,
		      const Matrix<complex<double> >& fs,
		      const Vector<double>& nn, const Vector<double>& mm
		      ) {

   const unsigned int Npts = X.size();
   const unsigned int NF = fs.Ncols();

   // find R and Z
  
   Vector<double> R(Npts, "transformsurface::R");
   Vector<double> Z(Npts, "transformsurface::Z");

   for (unsigned int i = 0; i<Npts;i++) {
      R[i] = sqrt( sqr(X[i].x()) + sqr(X[i].y()) );
      Z[i] = X[i].z();
   }


   Vector<complex<double> > RF(NF,"transformsurface::RF");
   transformfunction(R,RF,fs);

   Vector<complex<double> > ZF(NF,"transformsurface::ZF");
   transformfunction(Z,ZF,fs);

   Vector<double> nRZ("nRZ");
   Vector<double> mRZ("mRZ");
   Vector<double> RZCOS("RZCOS");
   Vector<double> RZSIN("RZSIN");

   vector<double> nnTEMP;
   vector<double> mmTEMP;
   vector<double> cosTEMP;
   vector<double> sinTEMP;

   convert_exp2sincos(nRZ,mRZ,RZCOS,RZSIN,nn,mm,RF);
   for(unsigned int i = 0; i<nRZ.size(); i++) {
      nnTEMP.push_back(nRZ[i]);
      mmTEMP.push_back(mRZ[i]);
      cosTEMP.push_back(RZCOS[i]);
      sinTEMP.push_back(0);
   }

   mRZ.resize(0);
   nRZ.resize(0);
   RZCOS.resize(0);
   RZSIN.resize(0);
   convert_exp2sincos(nRZ,mRZ,RZCOS,RZSIN,nn,mm,ZF);
   for(unsigned int i = 0; i<nRZ.size(); i++) {
      nnTEMP.push_back(nRZ[i]);
      mmTEMP.push_back(mRZ[i]);
      cosTEMP.push_back(0);
      sinTEMP.push_back(RZSIN[i]);
   }
    
   mRZ.resize() = vcast<double>(mmTEMP);
   nRZ.resize() = vcast<double>(nnTEMP);
   RZCOS.resize() = vcast<double>(cosTEMP);
   RZSIN.resize() = vcast<double>(sinTEMP);

   Vector<complex<double> > RZF(NF,"RZF");
   Vector<double> nntemp = nn;
   Vector<double> mmtemp = mm;
   convert_sincos2exp(nRZ,mRZ,RZCOS,RZSIN,nntemp,mmtemp,RZF,false);
   convert_exp2sincos(nRZ,mRZ,RZCOS,RZSIN,nn,mm,RZF);

   fsurface.nn() = nRZ;
   fsurface.mm() = mRZ;
   fsurface.RF() = RZCOS;
   fsurface.ZF() = RZSIN;


}






void expandfunction(
		    Vector<complex<double> >&func,
		    const Vector<complex<double> > &funcF,
		    const Matrix<complex<double> >& f_ortho
		    ) { 
  
//   const unsigned int NF = funcF.size();
//   const unsigned int Npts = func.size();
  

   //   for (unsigned int i = 0; i<Npts;i++) {
   //     func[i] = 0;
   //     for (unsigned int k = 0; k<NF ; k++) { 
   //       func[i] += funcF[k] * f_ortho(i,k);
   //     }
   //   }

   func = (f_ortho|funcF);

}


void expandfunction(
		    Vector<double>&funcREAL,
		    const Vector<complex<double> > &funcF,
		    const Matrix<complex<double> >& f_ortho
		    ) { 
   const unsigned int Npts = funcREAL.size();
   Vector<complex<double> > funcCMPLX(Npts,"expandsurface::funcCMPLX");
   expandfunction(funcCMPLX,funcF,f_ortho);
   funcREAL = real(funcCMPLX);
}


// generate coef's for orthonormal expansion of given function
void transformfunction(
		       const Vector<double>&func,
		       Vector<complex<double> > &funcF,
		       const Matrix<complex<double> >& f_ortho
		       ) {
  
   const unsigned int NF = f_ortho.Ncols();
   const unsigned int Npts = func.size();
   funcF.resize(NF);

   //   for (unsigned int k = 0; k<NF ; k++) { 
   //     funcF[k] = 0;
   //     for (unsigned int i = 0; i<Npts;i++) {
   //       funcF[k] +=  conj(f_ortho(i,k)) * func[i];
   //     }
   //   }
  
   funcF = (adj(f_ortho)|func);

   // this comes about when approximating integral usign a summation
   // i.e. dtheta = 2*pi/Ntheta and dphi = 2*pi/Nphi
   // Npts = Ntheta*Nphi
   const double coef = sqr(2*PI) / double(Npts);
   funcF = funcF *coef;
}





// insert coef's for all conjugate modes (as needed)

void completetheseries(const Vector<double>& nn, const Vector<double>& mm,
		       Vector<complex<double> >& xF)
{
   //   int nmax = int(max(abs(nn)));
   //   int mmax = int(max(abs(mm)));
   //   unsigned int mt = 2*mmax +1;
   //   for( int n = -nmax; n <= 0; n++){
   //     for( int m = -mmax; m <= mmax; m++){
   //       unsigned int k = (n+nmax)*mt + (m+mmax);
   //       unsigned int knegate = (-n+nmax)*mt + (-m+mmax);

   unsigned int NF=nn.size();
   Vector<unsigned int> knegate(NF,"knegate");
   for( unsigned int k = 0; k<NF; k++){
      for( unsigned int j = 0; j<NF; j++){
	 if ( (mm[j] == - mm[k]) && (nn[j] == - nn[k]) ) {
	    knegate[k] = j;
	 }
      }
   }

   for( unsigned int k = 0; k<NF; k++){
      if (nn[k] != - nn[knegate[k]]) {
	 std::cerr <<"BUG ENCOUNTERED at line "<<__LINE__<<" in file "<<__FILE__<<std::endl;
      }
      if (mm[k] != - mm[knegate[k]]) {
	 std::cerr <<"BUG ENCOUNTERED at line "<<__LINE__<<" in file "<<__FILE__<<std::endl;
      }
      
      if (xF[k] == std::complex<double>(0.0,0.0)) {
	 xF[k] = conj(xF[knegate[k]]);
      } else if (xF[knegate[k]] == std::complex<double>(0.0,0.0)) {
	 xF[knegate[k]] = conj(xF[k]);
      }
   }

}



// If mode == CondenseMode_zeros_only, this function removes all zero coef's
// If mode == CondenseMode_zeros_conjugs, this function removes all zero coefs and conjugate coef's

void condensetheseries(const Vector<double>& nn, const Vector<double>& mm, 
		       const Vector<complex<double> >& xF, 
		       Matrix<double>& Aout,
		       const CondenseMode mode, const double cutoff) {
  
   const unsigned int NF = xF.size();

   std::vector<double> nnew;
   std::vector<double> mnew;
   std::vector<std::complex<double> > xnew;

   Vector<unsigned int> knegate(NF,"knegate");
   for( unsigned int k = 0; k<NF; k++){
      for( unsigned int j = 0; j<NF; j++){
	 if ( (mm[j] == - mm[k]) && (nn[j] == - nn[k]) ) {
	    knegate[k] = j;
	 }
      }
   }


   for (unsigned int k = 0; k<NF; k++){
      const int n = int(nn[k]);
      const int m = int(mm[k]);

      if ((abs(xF[k].real()) > cutoff) || (abs(xF[k].imag()) > cutoff)) {
	 if (mode==CondenseMode_zeros_only) {
	    nnew.push_back(nn[k]);
	    mnew.push_back(mm[k]);
	    xnew.push_back(xF[k]);
	 } else if (mode==CondenseMode_zeros_conjugates) {
	    if ((n==0) && (m==0)) {
	       nnew.push_back(nn[k]);
	       mnew.push_back(mm[k]);
	       xnew.push_back(xF[k]);
	    }
	    if ((n==0) && (m>=1)){
	       nnew.push_back(nn[k]);
	       mnew.push_back(mm[k]);
	       xnew.push_back(0.5*(xF[k]+conj(xF[knegate[k]])));
	    }
	    if (n>=1) {
	       nnew.push_back(nn[k]);
	       mnew.push_back(mm[k]);
	       xnew.push_back(0.5*(xF[k]+conj(xF[knegate[k]])));
	    }
	 }
      }
   }

   const unsigned int NF2 = xnew.size();

   Aout.resize(NF2,4);
   for (unsigned int k = 0; k<NF2; k++){
      Aout(k,0) = nnew[k];
      Aout(k,1) = mnew[k];
      Aout(k,2) = xnew[k].real();
      Aout(k,3) = xnew[k].imag();
   }

}





// convert 
// from condensed sin/cos series:
//      g(theta,phi) = sum{ Amn*cos(m*theta + n*phi) + Bmn*sin(m*theta +n*phi) }n=0,1,...; m=0,+/-1,+/-2,...
// to complex orthonormal exponential series
//      g(theta,phi) = sum{ Gmn * (1/2pi)*exp(i*m*theta + i*n*phi) }n=0,+/-1,+/-2...; m=0,+/-1,+/-2,...

int convert_sincos2exp(const Vector<double>& ncosin, const Vector<double>& mcosin,
		       const Vector<double>& xCOS, const Vector<double>& xSIN,
		       Vector<double>& nn, Vector<double>& mm,
		       Vector<complex<double> >& xF, const bool giveAllWarnings)
{

   const unsigned int Nmodes = ncosin.size();
   int nmax;
   int mmax;
   unsigned int NF;
   bool warned_already = false;

   if ( nn.size()==0) {
      // create mode vectors if not given
      std::cerr <<"BUG ENCOUNTERED at line "<<__LINE__<<" in file "<<__FILE__<<std::endl;
      nmax = Matricks::round2int(max(abs(ncosin)));
      mmax = Matricks::round2int(max(abs(mcosin)));
      modevectors(NF,nn,mm,nmax,mmax);
   } else {
      // otherwise work with what we're given
      nmax = Matricks::round2int(max(abs(nn)));
      mmax = Matricks::round2int(max(abs(mm)));
      NF = nn.size();
   }


   xF.resize(NF);
   xF = 0;

   for(unsigned int j = 0; j<Nmodes; j++) {
      int n=Matricks::round2int(ncosin[j]);
      int m=Matricks::round2int(mcosin[j]);
      double Acos = xCOS[j];
      double Asin = xSIN[j];

      bool found = false;

      for(unsigned int k = 0; k<NF; k++) {
	 int nn_k = Matricks::round2int(nn[k]); 
	 int mm_k = Matricks::round2int(mm[k]); 
	 if ( (nn_k==n) && (mm_k==m) ){
	    if ( (n==0) && (m==0) ) {
	       xF[k] += complex<double>(Acos,Asin);
	       found=true;
	    } else {
	       xF[k] += complex<double>(Acos/2.0,-Asin/2.0);
	       found=true;
	    }
	 }
      }
      if (!found) {
	 if (giveAllWarnings || !warned_already) {
	    //	  dispcr(giveAllWarnings);
	    //      dispcr(warned_already);
	    std::cerr<<"WARNING: mode (n="<<n<<",m="<<m<<") was ignored"<<std::endl;
	    std::cerr<<"       during conversion from sin/cos to complex exponential series"<<std::endl;
	    if (!giveAllWarnings)
	       std::cerr<<"       FURTHER WARNINGS WILL BE SUPPRESSED"<<std::endl;
	    warned_already = true;
	 }
      }
   }
  
   
   // insert coef's for all conjugate modes (as needed)
   completetheseries(nn,mm,xF);

   const double coef = (2 * PI);
   xF = xF *coef;
   return 0;
}




// convert 
// from complex orthonormal exponential series
//      g(theta,phi) = sum{ Gmn * (1/2pi)*exp(i*m*theta + i*n*phi) }n=0,+/-1,+/-2...; m=0,+/-1,+/-2,...
// to  condensed sin/cos series:
//      g(theta,phi) = sum{ Amn*cos(m*theta + n*phi) + Bmn*sin(m*theta +n*phi) }n=0,1,...; m=0,+/-1,+/-2,...



int convert_exp2sincos(Vector<double>& ncosin, Vector<double>& mcosin,
		       Vector<double>& xCOS,  Vector<double>& xSIN,
		       const Vector<double>& nn, const Vector<double>& mm,
		       const Vector<complex<double> >& xF)
{

   Matrix<double> Atemp("Atemp");
   condensetheseries(nn,mm,xF,Atemp,CondenseMode_zeros_conjugates);
  
   const unsigned int Nmodes=Atemp.Nrows();
   ncosin.resize(Nmodes);
   mcosin.resize(Nmodes);
   xCOS.resize(Nmodes);
   xSIN.resize(Nmodes);

   for(uint i = 0; i<Nmodes; i++) {
      ncosin[i] = Atemp(i,0);
      mcosin[i] = Atemp(i,1);
      int n = int(ncosin[i]);
      int m = int(mcosin[i]);
      if ((n==0) && (m==0)) {
	 xCOS[i] = Atemp(i,2);
	 xSIN[i] = Atemp(i,3);
      } else {
	 xCOS[i] = 2.0*Atemp(i,2);
	 xSIN[i] = -2.0*Atemp(i,3);
      }
   }
   const double coef = 1/(2 * PI);
   xCOS = xCOS *coef;
   xSIN = xSIN *coef;
  
   return 0;
}







