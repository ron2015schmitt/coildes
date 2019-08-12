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




unsigned int mod(int j, unsigned int N) {

   unsigned int k;

   if (j>=int(N))
      k=j % N;
   else if (j<0) {
      int m = j%N;
      if (m<0)
	 m = - m;
      dispcr(m);
      k = m;
   } else
      k = j;

   return k;
}



// Main Function for code

int main (int argc, char *argv[])
{
  disable_all_options();
  enable_option(opt_if);
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
  LAvector<double> nnR("nn");
  LAvector<double> mmR("mm");
  unsigned int NFR;
  bool mode00 = false;
  if ( (Nharm >1) ||(Mharm>1) )
     modevectors(NFR,nnR,mmR,Nnn,Nmm,Nharm,Mharm,mode00);
  else
     modevectors(NFR,nnR,mmR,Nnn,Nmm,mode00);


  cout <<endl<< "$ Loading COIL CURRENT fourier coefficients from " << current_filename << endl;
  LAvector<complex<double> > IF(NFR,"IF");
  if (load_coefs(current_filename,CoefFileFormat_sincos,nnR,mmR,IF))
    return 3;  
  LAvector<double> IFreal(NFR,"IFreal");
  LAvector<double> IFimag(NFR,"IFimag");
  IFreal = real(IF);
  IFimag = imag(IF);
  // print coef's
  //  printfouriercoefs(nnR,mmR,IFreal,IFimag,10,18);


  cout << endl;
  cout<<"$ Generate orthonormal series matrix ("<<Npts<<" x "<<NFR<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  Matrix<complex<double> > fs(Npts,NFR,"fs");
  fseries(nnR,mmR,thetas,phis,fs);

  STOPTIME(tbuff,ckstart);
 



  // expand current onto theta,phi grid 

  LAvector<double> kappa(Npts,"kappa");
  expandfunction(kappa,IF,fs);

  LAvector<double> kappa_pol(Npts,"kappa_pol");
  LAvector<double> kappa_tor(Npts,"kappa_tor");

  for (unsigned int i = 0; i<Npts; i++) {
     kappa_pol[i] =  Ipoloidal/(2*PI)*phis[i];
     kappa_tor[i] =  -Itoroidal/(2*PI)*thetas[i];
     kappa[i] =  kappa[i] + kappa_pol[i] + kappa_tor[i];
  }


  // cutting algortihm
  unsigned int cutfactor = 8;
  unsigned int Ncuts = Nphi/cutfactor;
  LAvector<double> coil(Npts,"coil");


   double dtheta = 2*M_PI/((double)Ntheta);
  double dphi = 2*M_PI/((double)Nphi);

  vector<double> T;
  vector<double> P;
  LAvector<int> dtt(8,'dtt');
  LAvector<int> dpp(8,'dpp');
    "{-1,-1,-1, 0, 0, 1, 1, 1}" >> dtt;
    "{-1, 0, 1,-1, 1,-1, 0, 1" >> dpp;
  LAvector<int> pp_next(8,'pp_next');
  LAvector<int> tt_next(8,'tt_next');
  dispcr(dtt);
  dispcr(dpp);
  for (unsigned int n = 0; n<0; n++) {
     unsigned int tt= 0;
     unsigned int pp = n*cutfactor;
     unsigned int N = 0;
     bool done =true;
     while (!done) {
	T.push_back(dtheta*tt);
	P.push_back(dphi*pp);
	//	pp_next = mod(pp+dpp,Nphi);
	//	tt_next = mod(tt+dtt,Ntheta);

     }
  }

  unsigned int k1 ;
  int k2 = 13;
  k1 = mod(k2,5);
  disp(k2);dispcr(k1);

  k2=-3;
  k1 = mod(k2,5);
  disp(k2);dispcr(k1);


  k2=-3;
  k1 = mod(k2,5);
  disp(k2);dispcr(k1);


  k2=4;
  k1 = mod(k2,5);
  disp(k2);dispcr(k1);


  int k3 ;
  k3 = k2 % 5;
  disp(k2);dispcr(k3);


  // ALL DONE, NOW SAVE TO FILES

  string fn_root;
  string fn_ext;
  disect_filename(current_filename,fn_root,fn_ext);
  
  ostringstream strm;

  strm.str("");
  strm << fn_root << ".grid=" << Ntheta<<"x"<<Nphi<<"."<<fn_ext;
  print("Saving function data to ");
  printcr(strm.str());
  //  save(kappa,strm.str());

  return 0;
} // main()





