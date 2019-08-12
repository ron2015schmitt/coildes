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
  enable_option(opt_ff2);
  enable_option(opt_Nphi);
  enable_option(opt_Ntheta);
  enable_option(opt_Nnn);
  enable_option(opt_Nmm);
  enable_option(opt_Nharm);
  enable_option(opt_Mharm);

  // parse command line input
 
  if (!parse_cmd(argc, argv))
     return 1;

  string fext = "out";
  //  string fext = plasma_extname;
 //  ios_base::fmtflags flags = ios_base::right | ios_base::scientific;
  string ftemp;
  p3vectorformat::textformat(text_nobraces);

  // variables for measuring times
  struct tms tbuff;
  clock_t ckstart;

  // display COOLL mode
  cout << endl;
  display_execution_mode();
  cout << endl;


  // Create angle grid
  const unsigned int Npts = Ntheta*Nphi;
  LAvector<double> thetas(Npts,"thetas");
  LAvector<double> phis(Npts,"phis");
  anglevectors(thetas, phis, Ntheta, Nphi);


  // Create Fourier Mode vectors
  LAvector<double> nn("nn");
  LAvector<double> mm("mm");
  unsigned int NF;
  bool mode00 = true;
  if ( (Nharm >1) ||(Mharm>1) )
     modevectors(NF,nn,mm,Nnn,Nmm,Nharm,Mharm,mode00);
  else
     modevectors(NF,nn,mm,Nnn,Nmm,mode00);

  cout << endl;
  cout<<"$ Generate orthonormal series matrix ("<<Npts<<" x "<<NF<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  Matrix<complex<double> > fs(Npts,NF,"fs");
  fseries(nn,mm,thetas,phis,fs);

  STOPTIME(tbuff,ckstart);

  
  LAvector<complex<double> > Func1F(NF,"Func1F");

  cout <<endl<< "$ Loading numerator Function sin/cos fourier coefficients from " << flux_filename << endl;
  if (load_coefs(flux_filename,CoefFileFormat_sincos,nn,mm,Func1F))
    return 3;


  LAvector<double> func1(Npts,"func");
  expandfunction(func1,Func1F,fs);


  LAvector<complex<double> > Func2F(NF,"Func2F");

  cout <<endl<< "$ Loading denomenator Function sin/cos fourier coefficients from " << flux2_filename << endl;
  if (load_coefs(flux2_filename,CoefFileFormat_sincos,nn,mm,Func2F))
    return 3;

  LAvector<double> func2(Npts,"func");
  expandfunction(func2,Func2F,fs);


  func1 = func1 / func2;


  transformfunction(func1,Func1F,fs);
  //  massage(func1F,1e-10);

  LAvector<unsigned int> temp = findtrue((nn==0)&&(mm==0));
print("index[n=0,m=0] = ");dispcr(temp);

  LAvector <double> datavec(NF,"datavec");
  datavec.perline(1);
  datavec.textformat(text_nobraces);
  string fname("");

  datavec = real(Func1F);
  fname = "divfuncs.R.out";
  dispcr(fname);
  save(datavec,fname);

  datavec = imag(Func1F);
  fname = "divfuncs.I.out";
  dispcr(fname);
  save(datavec,fname);

  ostringstream strm;
  strm.str("divfuncs.out");
  //  strm << flux_filename << "_" << Ntheta<<"x"<<Nphi<<".out";
  print("Saving function data to ");
  printcr(strm.str());
  save_coefs(strm.str(),CoefFileFormat_sincos,nn,mm,Func1F);


  return 0;
} // main()





