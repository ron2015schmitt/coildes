/************************************************************************* 
 * 
 *   File Name    :  inductancematrix.cpp
 *   Platform     :  Red Hat LINUX 
 *   Author       :  Ron Schmitt
 *   Date         :  16 May 2003
 * 
 *
 *   SYNOPSIS     
 *     Generates mutual inductance matrix
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













void inductancematrix(const Vector<p3vector<double> >& X1, const Vector<p3vector<double> >& n1,
		      const Vector<p3vector<double> >& X2, const Vector<p3vector<double> >& n2, 
		      Matrix<double>& M) {
  
  const unsigned int Npts1 = X1.size();
  const unsigned int Npts2 = X2.size();

 
  for(unsigned int i=0; i<Npts1; i++) {
    for(unsigned int j=0; j<Npts2; j++) {

      const p3vector<double> rij(X1[i]- X2[j]);
      const p3vector<double> n1i(n1[i]);
      const p3vector<double> n2j(n2[j]);

// this is optimized so that we only have to perform one multiplicative inverse, as 
// a floatiing point divide typically takes ~50 times longer than a multiply
      const double r = norm(rij);
      const double invr = 1/r;
      const double invr3 = invr*invr*invr;
      const double invr5 = invr3*invr*invr;

     
      M(i,j) = mu0div4pi * ( 3*dot(n1i,rij)*dot(rij,n2j)*invr5 - dot(n1i,n2j)*invr3 );

    }
  }
}



