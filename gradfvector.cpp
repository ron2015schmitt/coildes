/************************************************************************* 
 * 
*   File Name    :  
 *   Platform     :  gcc compiler (v3.2.2)
 *   Author       : 
 *   Date         : 
 * 
 *
 *   SYNOPSIS     
 *     Generates magentic field for a straight wire
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

#include "gradfvector.hpp"



void gradfvector( const Vector<double>& thetas, const Vector<double>& phis,
		  const Vector<double>& mm, const Vector<double>& nn, 
		  const Vector<p3vector<double> >& grad_theta, 
		  const Vector<p3vector<double> >& grad_phi,
		  Matrix<p3vector<complex<double> > >& grad_f )
{
  
  const unsigned int NF = mm.size();
  const unsigned int Npts = phis.size();

  grad_f.resize(Npts,NF);

  for (unsigned int k=0; k<NF; k++ ) {
    double m = mm[k];
    double n = nn[k];
    Vector<double> angle = m*thetas + n*phis;
 
      
    // create sin,cos kernals for fourier transform
    // this coef makes the series orthonormal
    const double coef = 1/(2*PI);
    Vector<double> sinkern = sin(angle)*coef;
    Vector<double> coskern = cos(angle)*coef;

    Vector<double> df_dthetaR = - m*sinkern;
    Vector<double> df_dthetaI = m*coskern;
    Vector<double> df_dphiR = - n*sinkern;
    Vector<double> df_dphiI = n*coskern;

    

    for (unsigned int i=0; i<Npts; i++ ) {
      p3vector<double> gradfR =  df_dthetaR[i] * grad_theta[i] + df_dphiR[i] * grad_phi[i];
      p3vector<double> gradfI =  df_dthetaI[i] * grad_theta[i] + df_dphiI[i] * grad_phi[i];
      grad_f(i,k).x() = complex<double>(gradfR.x(),gradfI.x());
      grad_f(i,k).y() = complex<double>(gradfR.y(),gradfI.y());
      grad_f(i,k).z() = complex<double>(gradfR.z(),gradfI.z());


    }
  }
  
}



void gradlambda( const Vector<double>& thetas, const Vector<double>& phis,
		 const Vector<double>& mm, const Vector<double>& nn, 
		 const Vector<complex<double> > lambdaF,
		 Vector<double>& dlambda_dtheta, 
		 Vector<double>& dlambda_dphi)
{
  
  const unsigned int NF = mm.size();
  const unsigned int Npts = phis.size();


  for (unsigned int k=0; k<NF; k++ ) {
    double m = mm[k];
    double n = nn[k];
    Vector<double> angle = m*thetas + n*phis;
 
      
    // create sin,cos kernals for fourier transform
    // this coef makes the series orthonormal
    const double coef = 1/(2*PI);
    Vector<double> sinkern = sin(angle)*coef;
    Vector<double> coskern = cos(angle)*coef;

    Vector<double> df_dthetaR = - m*sinkern;
    Vector<double> df_dthetaI = m*coskern;
    Vector<double> df_dphiR = - n*sinkern;
    Vector<double> df_dphiI = n*coskern;
    Vector<complex<double> >df_dtheta(Npts);
    Vector<complex<double> >df_dphi(Npts);
    

    for (unsigned int i=0; i<Npts; i++ ) {
       df_dtheta[i] = complex<double>(df_dthetaR[i],df_dthetaI[i]);
        complex<double> temp_theta =  lambdaF[k]*df_dtheta[i];
       dlambda_dtheta[i] += temp_theta.real();

       df_dphi[i] = complex<double>(df_dphiR[i],df_dphiI[i]);
       complex<double> temp_phi =  lambdaF[k]*df_dphi[i];
       dlambda_dphi[i] += temp_phi.real();
    }
  }
  
}






