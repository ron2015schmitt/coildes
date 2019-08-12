/************************************************************************* 
 * 
 *   File Name    :  coilfwd.cpp
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  16 May 2003
 * 
 *
 *   SYNOPSIS     
 *     This file finds the plasma surface Bnormal for a given coil current.
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
#include "rhomatrix.hpp"


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

  // coefficient C is the integration coef for the fourier transform
  // C = dtheta*dphi
  //   = (2*pi/Ntheta)*(2*pi/Nphi)
  //   = (2*pi*2*pi/Npts)
  const double C = (2*PI*2*PI/double(Npts));
  const double Csqr = C*C;


  // Create Fourier Mode vectors
  LAvector<double> nn("nn");
  LAvector<double> mm("mm");
  unsigned int NF;
  modevectors(NF,nn,mm,Nnn,Nmm);


  // these exclude the n=0,m=0 case
  LAvector<double> nnR("nnR");
  LAvector<double> mmR("mmR");
  unsigned int NFR;
  modevectorsR(NFR,nnR,mmR,Nnn,Nmm);

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
    


  // load the coil CURRENT fourier coef's
  // at some point, add code so that user can select type of coef's from command line

  cout <<endl<< "$ Loading COIL CURRENT fourier coefficients from " << current_filename << endl;
  LAvector<complex<double> > IF(NFR,"IF");
  if (load_coefs(current_filename,CoefFileFormat_sincos,nnR,mmR,IF))
    return 3;  
  LAvector<double> IFreal(NFR,"IFreal");
  LAvector<double> IFimag(NFR,"IFimag");
  IFreal = real(IF);
  IFimag = imag(IF);

  // print coef's
  printfouriercoefs(nnR,mmR,IFreal,IFimag,10,18);

 
 
  // lay plasma surface onto grid 
  
  LAvector<p3vector<double> > X(Npts, "X");
  LAvector<p3vector<double> > dA_dtdp(Npts, "dA_dtdp");

  cout << endl;
  cout <<"$ Mapping plasma surface fourier coefficients to "<<Ntheta<<" x "<<Nphi<<" (theta by phi) grid"<<endl;

  STARTTIME(tbuff,ckstart);

  expandsurface(X,dA_dtdp,plasmafourier,thetas,phis);

  STOPTIME(tbuff,ckstart);

 
  

  // lay coil surface onto grid 
  
  LAvector<p3vector<double> > Xcoil(Npts, "Xcoil");
  LAvector<p3vector<double> > dA_dtdp_coil(Npts, "dA_dtdp_coil");

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





  // Create ortho normal series


  cout << endl;
  cout<<"$ Generate orthonormal series matrix ("<<Npts<<" x "<<NF<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  Matrix<complex<double> > fs(Npts,NF,"fs");
  fseries(nn,mm,thetas,phis,fs);


  STOPTIME(tbuff,ckstart);
  cout << endl;
  cout<<"$ Generate reduced orthonormal series matrix ("<<Npts<<" x "<<NFR<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  Matrix<complex<double> > fsR(Npts,NFR,"fsR");
  fseries(nnR,mmR,thetas,phis,fsR);

  STOPTIME(tbuff,ckstart);


  // Create mutual inductance matrix

  cout << endl;
  cout<<"$ Generating inductance matrix ("<<Npts<<" x "<<Npts<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  LAvector<p3vector<double> > n(Npts, "n");
  for (unsigned int j =0; j<Npts; j++)
    n[j] = dA_dtdp[j] / norm(dA_dtdp[j]);

  Matrix<double> M(Npts, Npts, "M");
  inductancematrix(X,n,Xcoil,dA_dtdp_coil,M);

  STOPTIME(tbuff,ckstart);


  // Fourier transform mutual inductance matrix, exp method

  cout << endl;
  cout<<"$ Fourier transform inductance matrix ("<<NF<<" x "<<NFR<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  Matrix<complex<double> > MF(NF,NFR,"MF"); 
  MF = (adj(fs)|M|fsR)*Csqr;

  STOPTIME(tbuff,ckstart);

  // clear M and its release its memory

  M.resize(0,0);

  // calculate Bnormal
  
  cout << endl;
  cout<<"$ Calculate plasma B normal ("<<NF<<" x 1)"<<endl;

  STARTTIME(tbuff,ckstart);

  LAvector<complex<double> > BnF(NF,"BnF");
  BnF = (MF|IF);



  STOPTIME(tbuff,ckstart);

  // ALL DONE, NOW SAVE TO FILES

  massage(BnF,1e-10);

  save_coefs("BnF.current2Bnormal.out",CoefFileFormat_sincos,nn,mm,BnF);

  // save in same format as David L., for comparison
  //  const double scale = 1.0/mu0;
  //  save_coefs("fluxcoefs.exp.scaled.out",CoefFileFormat_complexexp,nn,mm,FluxF*scale);


  return 0;
} // main()





