/************************************************************************* 
 * 
 *   File Name    :  createsurface.cpp
 *   Platform     :  Red Hat LINUX 
 *   Author       :  Ron Schmitt
 *   Date         :  16 May 2003
 * 
 *
 *   SYNOPSIS     
 *     Generates a surface that encloses the given surface
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
#include "createsurface.hpp"


// fprime=df/dtheta using a third order derivative formula
// f is a (theta, phi) formatted matrix

inline double dfdtheta(const Matrix<double>& f, const unsigned int t, const unsigned int p, const double dtheta, const unsigned int Ntheta) {

  double  fprime = ( f( (t+3)%Ntheta, p) - 9* f( (t+2)%Ntheta, p) + 45* f( (t+1)%Ntheta, p) 
		     - 45*f( (t-1+Ntheta)%Ntheta, p) + 9*f( (t-2+Ntheta)%Ntheta, p) - f( (t-3+Ntheta)%Ntheta, p) ) / (60*dtheta); 
  return fprime;

}




// fp=df/dphi using a third order derivative formula
// f is a (theta, phi) formatted matrix

inline double dfdphi(const Matrix<double>& f, const unsigned int t, const unsigned int p, const double dphi, const unsigned int Nphi) {

  double fprime = ( f( t, (p+3)%Nphi ) - 9* f( t, (p+2)%Nphi ) + 45* f( t, (p+1)%Nphi ) 
		    - 45*f( t, (p-1+Nphi)%Nphi ) + 9*f( t, (p-2+Nphi)%Nphi ) - f( t, (p-3+Nphi)%Nphi) ) / (60*dphi); 
  return fprime;

}





void createsurface(const Vector<p3vector<double> >& XX, const Vector<p3vector<double> >& AA,
		   Vector<p3vector<double> >& XXnew,  Vector<p3vector<double> >& AAnew,
		   const unsigned int Ntheta, const unsigned int Nphi, 
		   const double dtheta, const double dphi,
		   const double scale)
{
  dispcr(scale);
  Matrix<double> Xnew(Ntheta,Nphi,"Xnew");
  Matrix<double> Ynew(Ntheta,Nphi,"Ynew");
  Matrix<double> Znew(Ntheta,Nphi,"Znew");
  unsigned int i = 0;
  for (unsigned int tt = 0; tt<Ntheta; tt++){
    for (unsigned int pp = 0; pp<Nphi; pp++){
      double ax = AA[i].x();
      double ay = AA[i].y();
      double az = AA[i].z();
      double  temp = scale/sqrt(ax*ax + ay*ay + az*az);
      Xnew(tt,pp) = XX[i].x() + ax * temp;
      Ynew(tt,pp) = XX[i].y() + ay * temp;
      Znew(tt,pp) = XX[i].z() + az * temp;
      i++;
    }
  }

 
  Matrix<double> AXnew(Ntheta,Nphi,"AXnew");
  Matrix<double> AYnew(Ntheta,Nphi,"AYnew");
  Matrix<double> AZnew(Ntheta,Nphi,"AZnew");
  
  for (unsigned int tt = 0; tt<Ntheta; tt++){
    for (unsigned int pp = 0; pp<Nphi; pp++){
      p3vector<double> dxdtheta, dxdphi, Atemp;
      
      dxdtheta.x() = dfdtheta(Xnew, tt, pp, dtheta, Ntheta);
      dxdtheta.y() = dfdtheta(Ynew, tt, pp, dtheta, Ntheta);
      dxdtheta.z() = dfdtheta(Znew, tt, pp, dtheta, Ntheta);
      
      dxdphi.x() = dfdphi(Xnew, tt, pp, dphi, Nphi);
      dxdphi.y() = dfdphi(Ynew, tt, pp, dphi, Nphi);
      dxdphi.z() = dfdphi(Znew, tt, pp, dphi, Nphi);
      
      Atemp = cross( dxdphi, dxdtheta )*dphi*dtheta;
      
      AXnew(tt,pp) = Atemp.x();
      AYnew(tt,pp) = Atemp.y();
      AZnew(tt,pp) = Atemp.z();

    }
  }


  i = 0;
  for (unsigned int tt = 0; tt<Ntheta; tt++){
    for (unsigned int pp = 0; pp<Nphi; pp++){
      XXnew[i].x() = Xnew(tt,pp);
      XXnew[i].y() = Ynew(tt,pp);
      XXnew[i].z() = Znew(tt,pp);
      AAnew[i].x() = AXnew(tt,pp);
      AAnew[i].y() = AYnew(tt,pp);
      AAnew[i].z() = AZnew(tt,pp);
      i++;

    }
  }
  

}  
