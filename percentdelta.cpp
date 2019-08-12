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
  enable_option(opt_pf);
  enable_option(opt_pf2);
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



  // Create Fourier Mode vectors
  Vector<double> nn("nn");
  Vector<double> mm("mm");
  unsigned int NF;
  modevectors(NF,nn,mm,Nnn,Nmm);


  // load the plasma surface fourier coef's

  cout << "$ Loading PLASMA SURFACE fourier coefficients from " << plasma_filename << endl;
  Vector<complex<double> > surf1F(NF,"surf1F");
  if (load_coefs(plasma_filename,CoefFileFormat_sincos,nn,mm,surf1F))
    return 1;


  cout << "$ Loading PLASMA SURFACE fourier coefficients from " << plasma2_filename << endl;
  Vector<complex<double> > surf2F(NF,"surf2F");
  if (load_coefs(plasma2_filename,CoefFileFormat_sincos,nn,mm,surf2F))
    return 1;

  //  NEED TO PERFORM THIS IN SIN/COS SPACE because
  // WE OUTPUT IN SIN COS SPACE

  Vector<complex<double> > diffF(NF,"diffF");
  Vector<double> diffFR(NF,"diffFR");
  Vector<double> diffFI(NF,"diffFI");
  diffF = surf2F-surf1F;
  diffFR = real(diffF);
  diffFI = imag(diffF);

  const double TINY = 1e-300;
  Vector<double> surf1FR(NF,"surf1FR");
  Vector<double> surf1FI(NF,"surf1FI");
  surf1FR = real(surf1F);
  surf1FR = surf1FR + TINY;
  surf1FI = imag(surf1F);
  surf1FI = surf1FI + TINY;

  Vector<double> deltaFR(NF,"deltaFR");
  Vector<double> deltaFI(NF,"deltaFI");
  deltaFR = diffFR/surf1FR*100;
  deltaFI = diffFI/surf1FI*100;

  Vector<complex<double> > deltaF(NF,"deltaF");
  deltaF = vcomplex(deltaFR, deltaFI);

  cout<<endl;
  cout << "Saving % difference coefs to percentdiff.out"<<endl;
  save_coefs("percentdiff.out",CoefFileFormat_sincos,nn,mm,deltaF);
  p3vectorformat::textformat(text_nobraces);

 
  return 0;
} // main()





