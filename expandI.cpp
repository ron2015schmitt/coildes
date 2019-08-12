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
#include "coilfft.hpp"

const double NEGLECT =  1e-12;




// Main Function for code

int main (int argc, char *argv[])
{
  disable_all_options();
  enable_option(opt_if);
  enable_option(opt_if2);
  enable_option(opt_Nphi);
  enable_option(opt_Ntheta);
  enable_option(opt_Nnn);
  enable_option(opt_Nmm);
  enable_option(opt_Itoroidal);
  enable_option(opt_Ipoloidal);

  enable_option(opt_Nharm);
  enable_option(opt_Mharm);
 
  // parse command line input
 
  if (!parse_cmd(argc, argv))
     return 1;


  // variables for measuring times
//  struct tms tbuff;
//  clock_t ckstart;

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
  LAvector<double> nnR("nn");
  LAvector<double> mmR("mm");
  unsigned int NFR;
  bool mode00 = false;
  if ( (Nharm >1) ||(Mharm>1) )
     modevectors(NFR,nnR,mmR,Nnn,Nmm,Nharm,Mharm,mode00);
  else
     modevectors(NFR,nnR,mmR,Nnn,Nmm,mode00);


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



  cout <<endl<< "$ Loading COIL CURRENT fourier coefficients from " << current_filename << endl;
  LAvector<complex<double> > IF(NFR,"IF");
  if (load_coefs(current_filename,CoefFileFormat_sincos,nnR,mmR,IF))
    return 3;  


  if (!current_filename2.empty()){
    cout <<endl<< "$ Loading 2nd set of COIL CURRENT fourier coefficients from " << current_filename2 << endl;
    LAvector<complex<double> > IF2(NFR,"IF2");
    if (load_coefs(current_filename2,CoefFileFormat_sincos,nnR,mmR,IF2))
      return 4;  
    IF=IF+IF2;
    save_coefs("IF1plus2.out",CoefFileFormat_sincos,nnR,mmR,IF);
  }


  LAvector<complex<double> > IF_LHC(NFR,"IF_LHC");
  if (load_coefs(current_filename,CoefFileFormat_sincos_RHC2LHC,nnR,mmR,IF_LHC))
    return 4;  

  if (!current_filename2.empty()){
    LAvector<complex<double> > IF2_LHC(NFR,"IF2_LHC");
    if (load_coefs(current_filename2,CoefFileFormat_sincos_RHC2LHC,nnR,mmR,IF2_LHC))
      return 5;  
    IF_LHC=IF_LHC+IF2_LHC;
  }

//  cout << endl;
//  cout<<"$ Generate orthonormal series matrix ("<<Npts<<" x "<<NFR<<")"<<endl;

//  STARTTIME(tbuff,ckstart);

//  Matrix<complex<double> > fs(Npts,NFR,"fs");
//  fseries(nnR,mmR,thetas,phis,fs);

//  STOPTIME(tbuff,ckstart);
 



  // expand current onto theta,phi grid 

  LAvector<double> kappa0(Npts,"kappa0");
  LAvector<double> kappa(Npts,"kappa");

//  expandfunction(kappa0,IF,fs);
    mode00 = false;
    ifft2d(kappa0,IF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,1/(2*PI),mode00);


  LAvector<double> kappa_pol(Npts,"kappa_pol");
  LAvector<double> kappa_tor(Npts,"kappa_tor");

double kappa_RMS =0;
  for (unsigned int i = 0; i<Npts; i++) {
     kappa_pol[i] =  Ipoloidal/(2*PI)*phis[i];
     kappa_tor[i] =  -Itoroidal/(2*PI)*thetas[i];
     kappa[i] =  kappa0[i] + kappa_pol[i] + kappa_tor[i];
     kappa_RMS += kappa0[i]*kappa0[i];
  }
  kappa_RMS = sqrt(kappa_RMS/double(Npts));
  dispcr(kappa_RMS);


  // ALL DONE, NOW SAVE TO FILES
  p3vectorformat::textformat(text_nobraces);
  kappa.perline(1);
  kappa.textformat(text_nobraces);
  kappa0.perline(1);
  kappa0.textformat(text_nobraces);

  string fn_root;
  string fn_ext;
  disect_filename(current_filename,fn_root,fn_ext);
  
  ostringstream strm;

  strm.str("");
  strm << fn_root << ".RMS"<<"."<<fn_ext;
  print("Saving RMS data to ");
  printcr(strm.str());
  save(kappa_RMS,strm.str());

  strm.str("");
  strm << fn_root << ".grid=" << Ntheta<<"x"<<Nphi<<"."<<fn_ext;
  print("Saving function data to ");
  printcr(strm.str());
  save(kappa,strm.str());


  strm.str("");
  strm << fn_root << ".0.grid=" << Ntheta<<"x"<<Nphi<<"."<<fn_ext;
  print("Saving function data to ");
  printcr(strm.str());
  save(kappa0,strm.str());




  // REPEAT FOR REVERSE COORDINATE HANDEDNESS

//  expandfunction(kappa0,IF_LHC,fs);
    mode00 = false;
    ifft2d(kappa0,IF_LHC,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,1/(2*PI),mode00);

  for (unsigned int i = 0; i<Npts; i++) {
     kappa_pol[i] =  Ipoloidal/(2*PI)*phis[i];
     kappa_tor[i] =  -Itoroidal/(2*PI)*thetas[i];
     kappa[i] =  kappa0[i] + kappa_pol[i] + kappa_tor[i];
  }



  strm.str("");
  strm << fn_root << ".grid=" << Ntheta<<"x"<<Nphi<<".handsw."<<fn_ext;
  print("Saving function data to ");
  printcr(strm.str());
  save(kappa,strm.str());

  strm.str("");
  strm << fn_root << ".0.grid=" << Ntheta<<"x"<<Nphi<<".handsw."<<fn_ext;
  print("Saving function data to ");
  printcr(strm.str());
  save(kappa0,strm.str());


  return 0;
} // main()





