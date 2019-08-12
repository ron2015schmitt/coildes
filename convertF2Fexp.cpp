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

  // display Matricks mode
  cout << endl;
  display_execution_mode();
  cout << endl;


  // Create angle grid


   ostringstream strmtemp;
   strmtemp <<Nnn;
   string Nnn_str(strmtemp.str());
   strmtemp.str("");
   strmtemp <<Nmm;
   string Nmm_str(strmtemp.str());
   strmtemp.str("");


   string nnmm_str(".Nn="+Nnn_str+".Nm="+Nmm_str);
   string fname;




  // Create Fourier Mode vectors
  Vector<double> nn("nn");
  Vector<double> mm("mm");
  unsigned int NF;
  bool mode00 = true;
  modevectors(NF,nn,mm,Nnn,Nmm,Nharm,Mharm,mode00);

  
  Vector<complex<double> > FuncF(NF,"FuncF");

  cout <<endl<< "$ Loading Function sin/cos fourier coefficients from " << flux_filename << endl;
  if (load_coefs(flux_filename,CoefFileFormat_sincos,nn,mm,FuncF))
    return 3;


  FuncF.perline(1);
  FuncF.textformat(text_nobraces);
  save(FuncF,"FuncF.out");



  return 0;
} // main()





