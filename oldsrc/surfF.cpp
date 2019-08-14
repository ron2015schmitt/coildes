/************************************************************************* 
 * 
 *   File Name    : 
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  16 May 2003
 * 
 *
 *   SYNOPSIS     
 *     This is the main source file for the coil code.
 *
 **************************************************************************/



// Standard C libraries
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <ctime>

// Standard C++ libraries

#include <iostream>
#include <complex>

using namespace std;



// coil libraries
#include "coils.hpp"
#include "coilio.hpp"
#include "coils_cmdline.hpp"
#include "surface.hpp"
#include "createsurface.hpp"




// Main Function for code

int main (int argc, char *argv[])
{

  disable_all_options();
  enable_option(opt_pf);
  enable_option(opt_Nphi);
  enable_option(opt_Ntheta);
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
  const unsigned int Npts = Ntheta*Nphi;
  Vector<double> thetas(Npts,"thetas");
  Vector<double> phis(Npts,"phis");
  anglevectors(thetas, phis, Ntheta, Nphi);

   // Create Fourier Mode vectors
   Vector<double> nn("nn");
   Vector<double> mm("mm");
   unsigned int NF;
   bool mode00 = true;
   modevectors(NF,nn,mm,Nnn,Nmm,Nharm,Mharm,mode00);


  // load the surface data

  

  // lay plasma surface onto grid 

  cout << endl;
  cout <<"$ Loading surface data onto "<<Ntheta<<" x "<<Nphi<<" (theta by phi) grid"<<endl;

  Vector<p3vector<double> > X(Npts, "X");



  p3vectorformat::textformat(text_nobraces);
  X.textformat(text_nobraces);
  load(X,plasma_filename);



  // Create fourier series


  cout << endl;
  cout<<"$ Generate fourier series matrix ("<<Npts<<" x "<<NF<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  Matrix<complex<double> > fs(Npts,NF,"fs");
  fseries(nn,mm,thetas,phis,fs);
  STOPTIME(tbuff,ckstart);

  cout << endl;
  printcr("Find fourier coef's of surface RF and ZF");

  STARTTIME(tbuff,ckstart);
  FourierSurface surffourier;
  transformsurface(X,surffourier,fs,nn,mm);

  STOPTIME(tbuff,ckstart);



  string plasma_rootname;
  string fn_ext;
  disect_filename(plasma_filename,plasma_rootname,fn_ext);



  // ALL DONE, NOW SAVE TO FILES
  ostringstream strm;
  strm.str("");
  strm << plasma_rootname << "_fourier" <<".out";
  print("Saving surface Fourier coefs to ");
  printcr(strm.str());

  save_fourier_surface(strm.str(),surffourier);

  return 0;
} // main()
