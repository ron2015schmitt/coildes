/************************************************************************* 
 * 
 *   File Name    :  
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  16 May 2003
 * 
 *
 *   SYNOPSIS     
 *
 *
 **************************************************************************/



// Standard C libraries
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>

// Standard C++ libraries

#include <iostream>
#include <complex>

using namespace std;



// coil libraries
#include "coils.hpp"
#include "coilio.hpp"
#include "coils_cmdline.hpp"
#include "surface.hpp"


const double NEGLECT =  1e-12;




// Main Function for code

int main (int argc, char *argv[])
{

  disable_all_options();
  enable_option(opt_ff);
  enable_option(opt_pff);
  enable_option(opt_Nphi);
  enable_option(opt_Ntheta);
  enable_option(opt_Nnn);
  enable_option(opt_Nmm);

 

  // parse command line input
 
  if (!parse_cmd(argc, argv))
     return 1;


  // variables for measuring times
  struct tms tbuff;
  clock_t ckstart;

  // display Matricks mode
  cout << endl;
  display_execution_mode();
  cout << endl;

  // Create angle grid
  const unsigned int Npts = Ntheta*Nphi;
  Vector<double> thetas(Npts,"thetas");
  Vector<double> phis(Npts,"phis");
  anglevectors(thetas, phis, Ntheta, Nphi);

  // coefficient C is the integration coef for the fourier transform
  // C = dtheta*dphi
  //   = (2*pi/Ntheta)*(2*pi/Nphi)
  //   = (2*pi*2*pi/Npts)
  const double C = (2*PI*2*PI/double(Npts));
  const double Csqr = C*C;


  // Create Fourier Mode vectors
  Vector<double> nn("nn");
  Vector<double> mm("mm");
  unsigned int NF;
  modevectors(NF,nn,mm,Nnn,Nmm);


  // Create ortho normal series


  cout << endl;
  cout<<"$ Generate orthonormal series matrix ("<<Npts<<" x "<<NF<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  Matrix<complex<double> > fs(Npts,NF,"fs");
  fseries(nn,mm,thetas,phis,fs);


  STOPTIME(tbuff,ckstart);

  

 // load the plasma flux fourier coef's
  Vector<complex<double> > FluxF(NF,"FluxF");
  cout <<endl<< "$ Loading original Plasma Flux sin/cos fourier coefficients from " << flux_filename << endl;
  if (load_coefs(flux_filename,CoefFileFormat_sincos,nn,mm,FluxF))
    return 3;
  Vector<double> Flux(Npts, "Flux");
  expandfunction(Flux,FluxF,fs);
 

 
  Vector<complex<double> > Flux2F(NF,"Flux2F");
  cout <<endl<< "$ Loading perturbed Plasma Flux sin/cos fourier coefficients from " << pert_flux_filename << endl;
  if (load_coefs(pert_flux_filename,CoefFileFormat_sincos,nn,mm,Flux2F))
    return 3;
  Vector<double> Flux2(Npts, "Flux2");
  expandfunction(Flux2,Flux2F,fs);




  Vector<double> delFlux(Npts, "delFlux");
  for (unsigned int i =0; i<Npts; i++)
    delFlux[i] = (Flux2[i]-Flux[i]);
  Vector<complex<double> > delFluxF(NF,"delFluxF");
  transformfunction(delFlux,delFluxF,fs);

  STOPTIME(tbuff,ckstart);

  

  // ALL DONE, NOW SAVE TO FILES


  ostringstream strm;
  string fname;

  massage(delFluxF,1e-8);
  fname = "delFluxF.out";
  save_coefs(fname,CoefFileFormat_sincos,nn,mm,delFluxF);

 

 
  return 0;
} // main()





