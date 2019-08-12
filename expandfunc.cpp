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

  // display COOLL mode
  cout << endl;
  display_execution_mode();
  cout << endl;


  // Create angle grid
  const unsigned int Npts = Ntheta*Nphi;
  LAvector<double> thetas(Npts,"thetas");
  LAvector<double> phis(Npts,"phis");
  anglevectors(thetas, phis, Ntheta, Nphi);


   ostringstream strmtemp;
   strmtemp <<Nnn;
   string Nnn_str(strmtemp.str());
   strmtemp.str("");
   strmtemp <<Nmm;
   string Nmm_str(strmtemp.str());
   strmtemp.str("");

   strmtemp <<Ntheta;
   string Ntheta_str(strmtemp.str());
   strmtemp.str("");
   strmtemp <<Nphi;
   string Nphi_str(strmtemp.str());
   strmtemp.str("");

   string nnmm_str(".Nn="+Nnn_str+".Nm="+Nmm_str);
   string NtNp_str(".Ntheta="+Ntheta_str+".Nphi="+Nphi_str);
   string fname;
   // save angle vectors
   thetas.perline(1);
   thetas.textformat(text_nobraces);
   phis.perline(1);
   phis.textformat(text_nobraces);
   fname =  "thetas"+NtNp_str +".out";
   save(thetas,fname);
   fname =  "phis"+NtNp_str +".out";
   save(phis,fname);




  // Create Fourier Mode vectors
  LAvector<double> nn("nn");
  LAvector<double> mm("mm");
  unsigned int NF;
  bool mode00 = true;
  modevectors(NF,nn,mm,Nnn,Nmm,Nharm,Mharm,mode00);

  cout << endl;
  cout<<"$ Generate orthonormal series matrix ("<<Npts<<" x "<<NF<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  Matrix<complex<double> > fs(Npts,NF,"fs");
  fseries(nn,mm,thetas,phis,fs);

  STOPTIME(tbuff,ckstart);

  
  LAvector<complex<double> > FuncF(NF,"FuncF");

  cout <<endl<< "$ Loading Function sin/cos fourier coefficients from " << flux_filename << endl;
  if (load_coefs(flux_filename,CoefFileFormat_sincos,nn,mm,FuncF))
    return 3;



   LAvector <double> datavec("datavec");
   datavec.perline(1);
   datavec.textformat(text_nobraces);


  LAvector<double> func(Npts,"func");
  expandfunction(func,FuncF,fs);


  // ALL DONE, NOW SAVE TO FILES
  p3vectorformat::textformat(text_nobraces);
  func.perline(1);
  func.textformat(text_nobraces);

  ostringstream strm;

  strm.str("");
  strm << flux_filename << "_" << Ntheta<<"x"<<Nphi<<".out";
  print("Saving function data to ");
  printcr(strm.str());
  save(func,strm.str());


  strm.str("");
  strm << flux_filename << ".exp.R.out";
  print("Saving Fourier exp function data to ");
  printcr(strm.str());
  datavec = real(FuncF);
  save(datavec,strm.str());

  strm.str("");
  strm << flux_filename << ".exp.I.out";
  print("Saving Fourier exp function data to ");
  printcr(strm.str());
  datavec = imag(FuncF);
  save(datavec,strm.str());



  return 0;
} // main()





