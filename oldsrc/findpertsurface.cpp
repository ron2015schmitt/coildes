/************************************************************************* 
 * 
 *   File Name    :  findpertsurface.cpp
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  
 * 
 *
 *   SYNOPSIS     
 *
 *  This takes the differnece between the necessary fourier flux coef's 
 * and the acutal flux coef's and finds the location of the perturbed magnetic
 * surface.
 *    
 *
 * VERSION NOTES:

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
#include "omegamatrix.hpp"
#include "gradfvector.hpp"


// This is the file that defines the B field configuration


// Main Function for code

int main (int argc, char *argv[])
{

  disable_all_options();
  enable_option(opt_pf);
  enable_option(opt_ff);
  enable_option(opt_Nphi);
  enable_option(opt_Ntheta);
  enable_option(opt_Nnn);
  enable_option(opt_Nmm);
  enable_option(opt_Btf);
  enable_option(opt_Bpf);
   enable_option(opt_Nharm);
   enable_option(opt_Mharm);



  // parse command line input
 
  if (!parse_cmd(argc, argv))
     return 1;


  // 
  //  ios_base::fmtflags flags = ios_base::right | ios_base::scientific;
  string fext = "txt";
  string ftemp;
  p3vectorformat::textformat(text_nobraces);


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




  nn.perline(1);
  nn.textformat(text_nobraces);
  save(nn,"nn.out");
  mm.perline(1);
  mm.textformat(text_nobraces);
  save(mm,"mm.out");

  nnR.perline(1);
  nnR.textformat(text_nobraces);
  save(nnR,"nnR.out");
  mmR.perline(1);
  mmR.textformat(text_nobraces);
  save(mmR,"mmR.out");

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

  // load the plasma surface fourier coef's

  cout << "$ Loading PLASMA SURFACE fourier coefficients from " << plasma_filename << endl;
  FourierSurface plasmafourier;
  if (load_fourier_surface(plasma_filename,plasmafourier))
    return 1;
  plasmafourier.RF().name("p.RF");
  plasmafourier.ZF().name("p.ZF");

  // print coef's
  //  printfouriercoefs(plasmafourier.nn(),plasmafourier.mm(),plasmafourier.RF(),plasmafourier.ZF(),10,18);
  



  // load the plasma flux fourier coef's
  Vector<complex<double> > delFluxF(NF,"delFluxF");
  cout <<endl<< "$ Loading  perturbation plasma Flux sin/cos fourier coefficients from " << flux_filename << endl;
  if (load_coefs(flux_filename,CoefFileFormat_sincos,nn,mm,delFluxF))
    return 3;


 
  // lay plasma surface onto grid 
  
  Vector<p3vector<double> > X(Npts, "X");
  Vector<p3vector<double> > dA_dtdp(Npts, "dA_dtdp");

  Vector<p3vector<double> > dx_dr(Npts, "dx_dr");
  Vector<p3vector<double> > dx_dtheta(Npts,"dx_dtheta");
  Vector<p3vector<double> > dx_dphi(Npts,"dx_dphi");
  Vector<p3vector<double> > grad_r(Npts,"grad_r");
  Vector<p3vector<double> > grad_theta(Npts,"grad_theta");
  Vector<p3vector<double> > grad_phi(Npts,"grad_phi");

    

  cout << endl;
  cout <<"$ Mapping plasma surface fourier coefficients to "<<Ntheta<<" x "<<Nphi<<" (theta by phi) grid"<<endl;

  STARTTIME(tbuff,ckstart);

  expandsurfaceandbases(X,dA_dtdp,dx_dr,dx_dtheta,dx_dphi,grad_r,grad_theta,grad_phi,plasmafourier,thetas,phis);
  Vector<double> J(Npts, "J");
  for (unsigned int j =0; j<Npts; j++) {
     J[j] = dot(dx_dr[j],cross(dx_dtheta[j], dx_dphi[j]));
  }


  STOPTIME(tbuff,ckstart);


  // load B field
  cout << endl;
  cout<<"$ Load tangent B field "<<endl;

  STARTTIME(tbuff,ckstart);

  Vector<p3vector<double> > B(Npts, "B");

  cout <<endl<< "$ Loading BTOTAL_theta fourier coefficients from " << Bt_filename << endl;
  cout <<endl<< "$ Loading BTOTAL_phi fourier coefficients from " << Bp_filename << endl;

  Vector<complex<double> > BtF(NF,"BtF");
  if (load_coefs( Bt_filename,CoefFileFormat_sincos,nn,mm,BtF))
    return 5;
  Vector<complex<double> > BpF(NF,"BpF");
  if (load_coefs( Bp_filename,CoefFileFormat_sincos,nn,mm,BpF))
    return 6;

  Vector<double> Bt(Npts, "Bt");
  expandfunction(Bt,BtF,fs);
  Vector<double>Bp(Npts, "Bp");
  expandfunction(Bp,BpF,fs);
  
  for (unsigned int j =0; j<Npts; j++)
    B[j] =  Bt[j] * dx_dtheta[j] + Bp[j] * dx_dphi[j];



  STOPTIME(tbuff,ckstart);

     // Calculate Omega matrix

     cout << endl;
     cout<<"$ Calculate Omega matrix ("<<NF<<"x"<<NFR<<")"<<endl;
     
     STARTTIME(tbuff,ckstart);
     
     Matrix<complex<double> > Omega(NF,NFR,"Omega");
     Matrix<p3vector<complex<double> > >  grad_fsR(Npts,NFR,"grad_fsR");
     gradfvector(thetas,phis,mmR,nnR,grad_theta,grad_phi,grad_fsR);
     
     Vector<double> w(Npts, "w");
     for (unsigned int i = 0; i<Npts; i++) 
	w[i] = J[i];
     
     new_omegamatrix(Omega,B, w, grad_fsR,fs);

     STOPTIME(tbuff,ckstart);

   // Invert Omega Matrix
  
   cout << endl;
   cout<<"$ Invert Omega Matrix. Result is ("<<NFR<<"x"<<NF<<")"<<endl;
   Matrix<complex<double> > OmegaInv(NFR,NF,"OmegaInv");
     
//    STARTTIME(tbuff,ckstart);
//    matricks_lapack::inv(Omega,OmegaInv); 
//    STOPTIME(tbuff,ckstart);


   STARTTIME(tbuff,ckstart);

   Matrix<complex<double> > U(NF,NF,"U");
   Matrix<complex<double> > V(NFR,NFR,"V");
   Vector<double> S(NFR,"S");
  
   matricks_lapack::svd(Omega,U,S,V);


   Matrix<double> Sinv(NFR,NF,"Sinv");
   Sinv=0.0;
   dispcr(S);
   for(unsigned int k=0; k<NFR;k++) {
      Sinv(k,k) = 1/S[k];
   }
   OmegaInv = (V|Sinv|adj(U));

   V.clear();
   U.clear();
               
   STOPTIME(tbuff,ckstart);



  // calculate perturbations
  
  cout << endl;
  cout<<"$ Calculate perturbations to plasma surface"<<endl;

  STARTTIME(tbuff,ckstart);


  Vector<complex<double> > deltaF(NFR,"deltaF");
  deltaF = (OmegaInv|(delFluxF));

  Vector<double> delta(Npts,"delta");
  expandfunction(delta,deltaF,fsR);

  
  STOPTIME(tbuff,ckstart);


  // calculate new surface coords
  
  cout << endl;
  cout<<"$ Calculate new plasma surface"<<endl;

  STARTTIME(tbuff,ckstart);

  Vector<p3vector<double> > ksi(Npts,"ksi");
  Vector<p3vector<double> > Xnew(Npts,"Xnew");
   for(unsigned int i=0; i<Npts;i++) {
    ksi[i] = delta[i] * dx_dr[i];
    Xnew[i] = X[i] + ksi[i];
  }
  

  // ALL DONE, NOW SAVE TO FILES


  ostringstream strm;
  string fname;


  massage(deltaF,1e-8);
  fname = "deltaF.out";
  save_coefs(fname,CoefFileFormat_sincos,nnR,mmR,deltaF);

  FourierSurface newfourier;
  transformsurface(Xnew,newfourier,fs,nn,mm);
  massage(newfourier,1e-8);
  fname = "perturbedsurface.out";
  save_fourier_surface(fname,newfourier);


  return 0;
} // main()





