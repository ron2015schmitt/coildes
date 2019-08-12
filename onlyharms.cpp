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
  enable_option(opt_Nnn);
  enable_option(opt_Nmm);
  enable_option(opt_Nharm);
  enable_option(opt_Mharm);

 

  // parse command line input
 
  if (!parse_cmd(argc, argv))
     return 1;


  // variables for measuring times
  struct tms tbuff;
  clock_t ckstart;

  // display COOLL mode
  cout << endl;
  display_execution_mode();
  cout << endl;


  dispcr(Nnn);
  dispcr(Nmm);
  dispcr(Nharm);
  dispcr(Mharm);

  // Create Fourier Mode vectors
  LAvector<double> nn("nn");
  LAvector<double> mm("mm");
  unsigned int NF;
  bool mode00 = true;
  if ( (Nharm >1) ||(Mharm>1) )
     modevectors(NF,nn,mm,Nnn,Nmm,Nharm,Mharm,mode00);
  else
     modevectors(NF,nn,mm,Nnn,Nmm,mode00);


  dispcr(Nnn);
  dispcr(Nmm);
  dispcr(nn.size());
  dispcr(mm.size());
  dispcr(NF);

  // load the plasma surface fourier coef's

  cout << "$ Loading PLASMA SURFACE fourier coefficients from " << plasma_filename << endl;
  LAvector<complex<double> > surf1F(NF,"surf1F");
  bool giveWarning = false;
  if (load_coefs(plasma_filename,CoefFileFormat_sincos,nn,mm,surf1F,giveWarning))
    return 1;



  cout<<endl;
  cout << "Saving difference coefs to onlyharms.out"<<endl;
  save_coefs("onlyharms.out",CoefFileFormat_sincos,nn,mm,surf1F,false);
  p3vectorformat::textformat(text_nobraces);

 
  return 0;
} // main()





