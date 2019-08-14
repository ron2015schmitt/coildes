/************************************************************************* 
 * 
 *   File Name    :  current2flux_fft.cpp
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  16 May 2003
 * 
 *
 *   SYNOPSIS     
 *     This file finds the plasma surface flux for a given coil current.
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
#include "inductancematrix.hpp"
#include "coilfft.hpp"

const double NEGLECT =  1e-12;




// Main Function for code

int main (int argc, char *argv[])
{

  disable_all_options();
  enable_option(opt_pf);
  enable_option(opt_cf);
  enable_option(opt_if);
  enable_option(opt_Nphi);
  enable_option(opt_Ntheta);
  enable_option(opt_Nnn);
  enable_option(opt_Nmm);
   enable_option(opt_Nharm);
   enable_option(opt_Mharm);

   enable_option(opt_Mf);


   // parse command line input
 
   if (!parse_cmd(argc, argv))
      return 1;

   printcr("----------------------------------------------------");
   for(int i=0; i<argc; i++){
      print(argv[i]);
      print(" ");
   }
   cr();
   printcr("----------------------------------------------------");


   // 
   //  ios_base::fmtflags flags = ios_base::right | ios_base::scientific;
   string fext = "txt";
   string ftemp;
   p3vectorformat::textformat(text_nobraces);

   Vector <double> datavec("datavec");
   datavec.perline(1);
   datavec.textformat(text_nobraces);
   Matrix <double> data("data");
   data.perline(1);
   data.textformat(text_nobraces);

   string fname;

   // variables for measuring times
   struct tms tbuff;
   clock_t ckstart;

   // display Matricks mode
   cout << endl;
   display_execution_mode();
   cout << endl;

   ostringstream strmtemp;

   string methstr(".f2c");
   string plasma_name;
   {
      unsigned int j = plasma_filename.rfind("_");
      plasma_name = plasma_filename.substr(0,j);
      plasma_name = "."+plasma_name;
   }


   strmtemp <<Ntheta;
   string Ntheta_str(strmtemp.str());
   strmtemp.str("");
   Ntheta_str = ".Ntheta="+Ntheta_str;

   strmtemp <<Nphi;
   string Nphi_str(strmtemp.str());
   strmtemp.str("");
   Nphi_str = ".Nphi="+Nphi_str;

   strmtemp <<Nnn;
   string Nnn_str(strmtemp.str());
   strmtemp.str("");
   Nnn_str = ".Nn="+Nnn_str;

   strmtemp <<Nmm;
   string Nmm_str(strmtemp.str());
   strmtemp.str("");
   Nmm_str = ".Nm="+Nmm_str;


   dispcr(plasma_name);
   dispcr(methstr);
   dispcr(Nphi_str+Ntheta_str);
   dispcr(Nnn_str+Nmm_str);



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


   // these exclude the n=0,m=0 case
   Vector<double> nnR("nnR");
   Vector<double> mmR("mmR");
   unsigned int NFR;
   bool mode00 = false;
   modevectors(NFR,nnR,mmR,Nnn,Nmm,Nharm,Mharm,mode00);


   // load the plasma surface fourier coef's

   cout << "$ Loading PLASMA SURFACE fourier coefficients from " << plasma_filename << endl;
   FourierSurface plasmafourier;
   if (load_fourier_surface(plasma_filename,plasmafourier)){
      printcr("Above ERROR occurred in "+myname+".");
      return 1;
   }

   // load the coil surface fourier coef's

   cout << "$ Loading COIL SURFACE fourier coefficients from " << coil_filename << endl;
   FourierSurface coilfourier;
   if (load_fourier_surface(coil_filename,coilfourier)){
      printcr("Above ERROR occurred in "+myname+".");
      return 2;
   }


  // load the coil CURRENT fourier coef's
  // at some point, add code so that user can select type of coef's from command line

  cout <<endl<< "$ Loading COIL CURRENT fourier coefficients from " << current_filename << endl;
  Vector<complex<double> > IF(NFR,"IF");
  if (load_coefs(current_filename,CoefFileFormat_sincos,nnR,mmR,IF)){
      printcr("Above ERROR occurred in "+myname+".");
      return 3;
   }

 
 
  // lay plasma surface onto grid 
  
  Vector<p3vector<double> > X(Npts, "X");
  Vector<p3vector<double> > dA_dtdp(Npts, "dA_dtdp");

  cout << endl;
  cout <<"$ Mapping plasma surface fourier coefficients to "<<Ntheta<<" x "<<Nphi<<" (theta by phi) grid"<<endl;

  STARTTIME(tbuff,ckstart);

  expandsurface(X,dA_dtdp,plasmafourier,thetas,phis);

  STOPTIME(tbuff,ckstart);

 
  

  // lay coil surface onto grid 
  
  Vector<p3vector<double> > Xcoil(Npts, "Xcoil");
  Vector<p3vector<double> > dA_dtdp_coil(Npts, "dA_dtdp_coil");

  cout << endl;
  cout <<"$ Mapping coil surface fourier coefficients to "<<Ntheta<<" x "<<Nphi<<" (theta by phi) grid"<<endl;

  STARTTIME(tbuff,ckstart);

  expandsurface(Xcoil,dA_dtdp_coil,coilfourier,thetas,phis);

  STOPTIME(tbuff,ckstart);



  // print grid spacing info
  {
    double del_plasma = gridspacing(X);
    cout <<endl<< "  estimate of plasma grid poloidal spacing = " << del_plasma <<endl;
    double del_coil = gridspacing(Xcoil);
    cout << "  estimate of coil grid poloidal spacing = " << del_coil <<endl;
    double del_c2p = coil2plasmaspacing(X,Xcoil);
    cout << "  estimate of spacing between coil and plasma grid = " << del_c2p <<endl<<endl;
    double gridfactor = 1.0;
    if ( (del_c2p < del_plasma*gridfactor) || (del_c2p < del_coil*gridfactor) ) {
      cout << "**WARNING: spacing between plasma and coils is less than (" << 1 << " * the grid spacing)"<< endl;
      cout <<endl;
    }
  }





   Matrix<complex<double> > MF(NFR,NFR,"MF"); 
   if ( M_filename.empty() ) {

      // Create mutual inductance matrix

      cout << endl;
      cout<<"$ Generating inductance matrix ("<<Npts<<" x "<<Npts<<")"<<endl;

      STARTTIME(tbuff,ckstart);

      Matrix<double> M(Npts, Npts, "M");
      inductancematrix(X,dA_dtdp,Xcoil,dA_dtdp_coil,M);
   
      STOPTIME(tbuff,ckstart);

      // fft  M matrix
      cout << endl;
      cout<<"$ FFT of inductance matrix ("<<Npts<<" x "<<NFR<<")"<<endl;
      
      STARTTIME(tbuff,ckstart);
  
      fft_of_M(M, MF, Nphi, Ntheta,Nnn,Nmm,  Nharm,Mharm,1e-12);
      M.clear();
      STOPTIME(tbuff,ckstart);

      fname = "MF"+plasma_name+Nnn_str+Nmm_str+".R.out";
      data.resize() = real(MF);
      data.perline(1);
      data.textformat(text_nobraces);
      dispcr(fname);
      save(data,fname);
      fname = "MF"+plasma_name+Nnn_str+Nmm_str+".I.out";
      data = imag(MF);
      dispcr(fname);
      save(data,fname);
      data.clear();

   } else {
      cout << endl;
      cout<<"$ Load M (inductance) matrix("<<NFR<<"x"<<NFR<<")"<<endl;
      
      STARTTIME(tbuff,ckstart);
      
      data.resize(NFR,NFR);
      data.perline(1);
      data.textformat(text_nobraces);
      fname = M_filename +Nnn_str+Nmm_str+ ".R.out";
      dispcr(fname);
      load(data,fname);
      Matrix <double> data2("data2");
      data2.resize(NFR,NFR);
      data2.perline(1);
      data2.textformat(text_nobraces);
      fname =  M_filename +Nnn_str+Nmm_str+ ".I.out";
      dispcr(fname);
      load(data2,fname);
      MF = mcomplex(data,data2);
      data.clear();
      data2.clear();

      STOPTIME(tbuff,ckstart);

   }






  //*************************************

  //////////////////////////////////////////////////////

  double maxMFR = max(abs(real(MF)));
  dispcr(maxMFR);
  double maxMFI = max(abs(imag(MF)));
  dispcr(maxMFI);
 
  /////////////////////////////////////////////////////



  // calculate flux
  
  cout << endl;
  cout<<"$ Calculate plasma flux ("<<NFR<<" x 1)"<<endl;

  STARTTIME(tbuff,ckstart);

  Vector<complex<double> > FluxF(NFR,"FluxF");
  FluxF = (MF|IF);

  STOPTIME(tbuff,ckstart);


  //////////////////////////////////////////////////////

  double maxFluxR = max(abs(real(FluxF)));
  dispcr(maxFluxR);
  double maxFluxI = max(abs(imag(FluxF)));
  dispcr(maxFluxI);
 
 
  /////////////////////////////////////////////////////







  // ALL DONE, NOW SAVE TO FILES

  massage(FluxF,1e-10);

  save_coefs("FluxF.current2flux_fft.out",CoefFileFormat_sincos,nnR,mmR,FluxF);

  // save in same format as David L., for comparison
  //  const double scale = 1.0/mu0;
  //  save_coefs("fluxcoefs.exp.scaled.out",CoefFileFormat_complexexp,nn,mm,FluxF*scale);


  return 0;
} // main()





