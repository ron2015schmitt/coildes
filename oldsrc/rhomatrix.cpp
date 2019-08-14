/************************************************************************* 
 * 
 *   File Name    :  
 *   Platform     :  
 *   Author       :  Ron Schmitt
 *   Date         :  
 * 
 *
 *   SYNOPSIS     
 *
 **************************************************************************/





#include "rhomatrix.hpp"




 void rhomatrix(const double alpha_theta, const double alpha_phi, 
		const Vector<double>& mm, const Vector<double>& nn, 
		Vector<double>& rho_sqrt)
{

  
   const unsigned int NF = mm.size();
   // rho matrices should be of size [NFxNF]

   // intialize entire matrix to zero
   rho_sqrt.resize(NF);

  for (unsigned int f = 0; f<NF ; f++) {
     double mode = abs(alpha_theta*mm[f]) + abs(alpha_phi*nn[f]); 
     // divide by two in exponent performs the sqrt function
     double temp = exp(mode/2);
     rho_sqrt[f] = temp;
  }
  
}




 void rhomatrix(const double alpha_theta, const double alpha_phi, 
		const Vector<double>& mm, const Vector<double>& nn, 
		Vector<double>& rho_sqrt, Vector<double>& rho_sqrt_inv)
{

  
   const unsigned int NF = mm.size();
   // rho matrices should be of size [NFxNF]

   // intialize entire matrix to zero
   rho_sqrt.resize(NF);
   rho_sqrt_inv.resize(NF);

  for (unsigned int f = 0; f<NF ; f++) {
     double mode = abs(alpha_theta*mm[f]) + abs(alpha_phi*nn[f]); 
     // divide by two in exponent performs the sqrt function
     double temp = exp(mode/2);
     rho_sqrt[f] = temp;
     rho_sqrt_inv[f] = 1/temp;
  }
  
}




void rhomatrix2(const double alpha_theta, const double alpha_phi, 
		const Vector<double>& mm, const Vector<double>& nn, 
		Vector<double>& rho_sqrt, Vector<double>& rho_sqrt_inv)
{

  
   const unsigned int NF = mm.size();
   // rho matrices should be of size [NFxNF]

   // intialize entire matrix to zero
   rho_sqrt.resize(NF);
   rho_sqrt_inv.resize(NF);

  for (unsigned int f = 0; f<NF ; f++) {
     double mode = sqr(alpha_theta*mm[f]) + sqr(alpha_phi*nn[f]); 
     // divide by two in exponent performs the sqrt function
     double temp = exp(mode/2);
     rho_sqrt[f] = temp;
     rho_sqrt_inv[f] = 1/temp;
  }
  
}




