/************************************************************************* 
 * 
 *   File Name    :  coilbwd.cpp
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  16 May 2003
 * 
 *
 *   SYNOPSIS     
 *     This file finds the coil current from given plasma surface flux 
 *    using smple Merkel method.
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
   bool mode00 = true;
   if ( (Nharm >1) ||(Mharm>1) )
      modevectors(NF,nn,mm,Nnn,Nmm,Nharm,Mharm,mode00);
   else
      modevectors(NF,nn,mm,Nnn,Nmm,mode00);



   // these exclude the n=0,m=0 case
   Vector<double> nnR("nnR");
   Vector<double> mmR("mmR");
   unsigned int NFR;
   mode00 = false;
   if ( (Nharm >1) ||(Mharm>1) )
      modevectors(NFR,nnR,mmR,Nnn,Nmm,Nharm,Mharm,mode00);
   else
      modevectors(NFR,nnR,mmR,Nnn,Nmm,mode00);


  // load the plasma surface fourier coef's

  cout << "$ Loading PLASMA SURFACE fourier coefficients from " << plasma_filename << endl;
  FourierSurface plasmafourier;
  if (load_fourier_surface(plasma_filename,plasmafourier))
    return 1;
  plasmafourier.RF().name("p.RF");
  plasmafourier.ZF().name("p.ZF");

  // print coef's
  printfouriercoefs(plasmafourier.nn(),plasmafourier.mm(),plasmafourier.RF(),plasmafourier.ZF(),10,18);
  


  // load the coil surface fourier coef's

  cout << "$ Loading COIL SURFACE fourier coefficients from " << coil_filename << endl;
  FourierSurface coilfourier;
  if (load_fourier_surface(coil_filename,coilfourier))
    return 2;

  coilfourier.RF().name("c.RF");
  coilfourier.ZF().name("c.ZF");

  // print coef's
  printfouriercoefs(coilfourier.nn(),coilfourier.mm(),coilfourier.RF(),coilfourier.ZF(),10,18);
    

  // load the plasma flux fourier coef's
  // at some point, add code so that user can select type of coef's from command line

  Vector<complex<double> > BnF(NF,"BnF");
  cout <<endl<< "$ Loading Plasma Bnormal sin/cos fourier coefficients from " << flux_filename << endl;
  if (load_coefs(flux_filename,CoefFileFormat_sincos,nn,mm,BnF))
    return 3;


 
   // lay plasma surface onto grid 
  
  Vector<p3vector<double> > X(Npts, "X");
  Vector<p3vector<double> > dAdtdp(Npts, "dAdtdp");
  Vector<p3vector<double> > dx_dr(Npts, "dx_dr");
  Vector<p3vector<double> > dx_dtheta(Npts,"dx_dtheta");
  Vector<p3vector<double> > dx_dphi(Npts,"dx_dphi");
  Vector<p3vector<double> > grad_r(Npts,"grad_r");
  Vector<p3vector<double> > grad_theta(Npts,"grad_theta");
  Vector<p3vector<double> > grad_phi(Npts,"grad_phi");
 
  cout << endl;
  cout <<"$ Mapping plasma surface fourier coefficients to "<<Ntheta<<" x "<<Nphi<<" (theta by phi) grid"<<endl;

  STARTTIME(tbuff,ckstart);

  expandsurfaceandbases(X,dAdtdp,dx_dr,dx_dtheta,dx_dphi,grad_r,grad_theta,grad_phi,plasmafourier,thetas,phis);

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




  // Create mutual inductance matrix

  cout << endl;
  cout<<"$ Generating inductance matrix ("<<Npts<<" x "<<Npts<<")"<<endl;

  STARTTIME(tbuff,ckstart);


  Matrix<double> M(Npts, Npts, "M");
  inductancematrix(X,grad_r,Xcoil,dA_dtdp_coil,M);

  STOPTIME(tbuff,ckstart);





  // Fourier transform mutual inductance matrix

  cout << endl;
  cout<<"$ FFT of  inductance matrix ("<<NFR<<" x "<<NFR<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  Matrix<complex<double> > MF(NFR,NFR,"MF"); 

  // note: this clobbers M and releases its memory

  fft_of_M(M, MF, Nphi, Ntheta,Nnn,Nmm,  Nharm,Mharm,1e-12);

  STOPTIME(tbuff,ckstart);




  //*************************************

  //////////////////////////////////////////////////////

  double maxMFR = max(abs(real(MF)));
  dispcr(maxMFR);
  double maxMFI = max(abs(imag(MF)));
  dispcr(maxMFI);
 
  /////////////////////////////////////////////////////



  // calculate SVD of MF
  cout << endl;
  cout<<"$ SVD of fourier inductance matrix ("<<NFR<<" x "<<NFR<<")"<<endl;

  STARTTIME(tbuff,ckstart);

   
  Matrix<complex<double> > U(NFR,NFR,"U");
  Matrix<complex<double> > V(NFR,NFR,"V");
  Vector<double> S(NFR);
  matricks_lapack::svd(MF,U,S,V);

  Vector<double> Snorm = S/S[0];
                                                                                
  STOPTIME(tbuff,ckstart);
   

 cout<<"Largest SV="<<S[0]<<endl;
 cout<<"Smallest SV="<<S[NFR-1]<<endl;
 double condition_number = S[0]/S[NFR-1];
 dispcr(condition_number);



  // calculate inverse of MF
  cout << endl;
  cout<<"$ Pseudo Inverse of MF("<<NFR<<" x "<<NFR<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  //   Matrix<complex<double> > MFinv(NFR,NFR,"MFinv");
  //   Matrix<double> Sinv(NFR,NFR,"Sinv");
  //   Sinv=0;
  //   for(unsigned int k=0; k<(NFR);k++)
  //     Sinv(k,k) = 1/S[k];
  //   MFinv = (V|Sinv|adj(U));


  double mincond = 1e-7;
  Matrix<complex<double> > SU(NFR,NFR,"SU");
  for(unsigned int r=0; r<NFR; r++) {
     double temp;
     if (Snorm[r] >= mincond)
	temp = 1/S[r];
     else
	temp = 0;
    for(unsigned int c=0; c<NFR; c++) {
      SU(r,c) = conj(U(c,r)) * temp;
    }
  }
  U.resize(0,0);
  S.resize(0);
  Matrix<complex<double> > MFinv(NFR,NFR,"MFinv");
  MFinv = (V|SU);
  
  STOPTIME(tbuff,ckstart);


  // calculate current
  
  cout << endl;
  cout<<"$ Calculate coil current ("<<NFR<<" x 1)"<<endl;

  STARTTIME(tbuff,ckstart);

  Vector<complex<double> > IF(NFR,"IF");
  IF = (MFinv|BnF);

  STOPTIME(tbuff,ckstart);


  massage(IF,1e-10);

  // ALL DONE, NOW SAVE TO FILES

  save_coefs("IF.Bnormal2current.out",CoefFileFormat_sincos,nnR,mmR,IF);

  return 0;
} // main()





