/************************************************************************* 
 * 
 *   File Name    :  calcgamma.cpp
 *   Platform     :  gcc compiler (v3.2.2)
 *   Author       :  Ron Schmitt
 *   Date         :  16 May 2003
 * 
 *
 *   SYNOPSIS     
 *     Generates gamma matrix (shape parameter displacement)
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

#include "gammamatrix.hpp"



void gammamatrix( Matrix<complex<double> > & gamma,
		  const Matrix<complex<double> >& fs,
		  const FourierSurface& fsurface, 
		  const LAvector<p3vector<double> > A,
		  const LAvector<double>& thetas, const LAvector<double>& phis,
		  const unsigned int Ntheta, const unsigned int Nphi) {


  const LAvector<double>& nnsurf = fsurface.nn();
  const LAvector<double>& mmsurf = fsurface.mm();
  const LAvector<double>& RF = fsurface.RF();
  const LAvector<double>& ZF = fsurface.ZF();
  const unsigned int Npmodes = mmsurf.size();
  const unsigned int Npts = A.size();
  const unsigned int NF = fs.Ncols();

  const double dphi = 2*PI/double(Nphi);
  const double dtheta = 2*PI/double(Ntheta);

  LAvector<double> nr(Npts,"nr");
  LAvector<double> nz(Npts,"nz");
  
  for (unsigned int i = 0; i<Npts ; i++) {
    double Amag= sqrt(sqr(A[i].x())+sqr(A[i].y())+sqr(A[i].z()));
    double nx = A[i].x()/Amag;
    double ny = A[i].y()/Amag;
    nr[i]=nx*cos(phis[i]) + ny*sin(phis[i]);
    nz[i]= A[i].z()/Amag;
  }


  unsigned int Nz = 0;
  unsigned int Nr = 0;

  for (unsigned int f = 0; f<Npmodes ; f++) {
    if (RF[f]!=0) 
      Nr++;
    if (ZF[f]!=0) 
      Nz++;
  }
  const unsigned int Ns=Nr+Nz;

  LAvector<double> m_R(Nr,"m_R");
  LAvector<double> n_R(Nr,"n_R");
  LAvector<double> m_Z(Nz,"m_Z");
  LAvector<double> n_Z(Nz,"n_Z");
  unsigned int jr=0;
  unsigned int jz=0;
  for (unsigned int f = 0; f<Npmodes ; f++) {
    if (RF[f]!=0) {
      m_R[jr]= mmsurf[f];
      n_R[jr]= nnsurf[f];
      jr++;
    }
    if (ZF[f]!=0) {
      m_Z[jz]= mmsurf[f];
      n_Z[jz]= nnsurf[f];
      jz++;
    }
  }


  Matrix<double> gp(Npts,Ns,"gp");

  for (unsigned int i = 0; i<Nr ; i++) {
    double m=m_R[i];
    double n=n_R[i];
    gp.col(i) = nr*cos(n*phis+m*thetas);
  }
  for (unsigned int i = 0; i<Nz ; i++) {
    double m=m_Z[i];
    double n=n_Z[i];
    gp.col(i+Nr) = nz*sin(n*phis+m*thetas);
  }


  // perform the integration
  gamma.resize(NF,Ns);
  gamma = (adj(fs)|gp)*dphi*dtheta;
  
}




void fullgammamatrix( Matrix<complex<double> > & gamma,
		      const LAvector<double>& nn, const LAvector<double>& mm, 
		      const Matrix<complex<double> >& fsIN,
		      const Matrix<complex<double> >& fsOUT,
		      const LAvector<p3vector<double> > A,
		      const LAvector<double>& thetas, const LAvector<double>& phis,
		      const unsigned int Ntheta, const unsigned int Nphi) 
{

  const unsigned int Npts = A.size();
  const unsigned int NFout = fsOUT.Ncols();
  const unsigned int NFin = fsIN.Ncols();

  const double dphi = 2*PI/double(Nphi);
  const double dtheta = 2*PI/double(Ntheta);

  LAvector<double> nr(Npts,"nr");
  LAvector<double> nz(Npts,"nz");
  
  for (unsigned int i = 0; i<Npts ; i++) {
    double Amag= sqrt(sqr(A[i].x())+sqr(A[i].y())+sqr(A[i].z()));
    double nx = A[i].x()/Amag;
    double ny = A[i].y()/Amag;
    nr[i]=nx*cos(phis[i]) + ny*sin(phis[i]);
    nz[i]= A[i].z()/Amag;
  }


  Matrix<double> gp(Npts,2*NFin,"gp");

  for (unsigned int i = 0; i<NFin ; i++) {
    double m=mm[i];
    double n=nn[i];
    gp.col(i) = nr*cos(n*phis+m*thetas);
  }
  for (unsigned int i = 0; i<NFin ; i++) {
    double m=mm[i];
    double n=nn[i];
    gp.col(i+NFin) = nz*sin(n*phis+m*thetas);
  }


  // perform the integration
  gamma.resize(NFout,2*NFin);
  gamma = (adj(fsOUT)|gp)*dphi*dtheta;

}
