/************************************************************************* 
 * 
 *   File Name    :  coilbwd2.cpp
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  
 * 
 *
 *   SYNOPSIS     
 *     This file finds the coil current from given plasma surface flux 
 *    using simple SDV inversion of M
 *
 * This version uses the required flux as the input.
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
  enable_option(opt_ff);
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
    

  // load the plasma flux fourier coef's
  // at some point, add code so that user can select type of coef's from command line

  LAvector<complex<double> > FluxF(NF,"FluxF");
  cout <<endl<< "$ Loading Plasma Flux sin/cos fourier coefficients from " << flux_filename << endl;
  if (load_coefs(flux_filename,CoefFileFormat_sincos,nn,mm,FluxF))
    return 3;


 
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
  inductancematrix(X,dA_dtdp,Xcoil,dA_dtdp_coil,M);

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

  // calculate SVD of MF
  cout << endl;
  cout<<"$ SVD of fourier inductance matrix ("<<NF<<" x "<<NFR<<")"<<endl;

  STARTTIME(tbuff,ckstart);

   
  Matrix<complex<double> > U(NF,NF,"U");
  Matrix<complex<double> > V(NFR,NFR,"V");
  LAvector<double> S(NFR);
  cooll_lapack::svd(MF,U,S,V);
                                                                                
  STOPTIME(tbuff,ckstart);
   


  // calculate inverse of MF
  cout << endl;
  cout<<"$ Pseudo Inverse of MF("<<NFR<<" x "<<NF<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  Matrix<complex<double> > MFinv(NFR,NF,"MFinv");
  Matrix<double> Sinv(NFR,NF,"Sinv");
  Sinv=0;
  for(unsigned int k=0; k<(NFR);k++)
    Sinv(k,k) = 1/S[k];
  MFinv = (V|Sinv|adj(U));
                                                                                
  STOPTIME(tbuff,ckstart);



  // calculate current
  
  cout << endl;
  cout<<"$ Calculate coil current ("<<NFR<<" x 1)"<<endl;

  STARTTIME(tbuff,ckstart);

  LAvector<complex<double> > IF(NFR,"IF");
  IF = (MFinv|FluxF);

  STOPTIME(tbuff,ckstart);


  massage(IF,1e-10);

  // ALL DONE, NOW SAVE TO FILES

  save_coefs("IF.flux2current.out",CoefFileFormat_sincos,nnR,mmR,IF);

  return 0;
} // main()





