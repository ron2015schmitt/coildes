/************************************************************************* 
 * 
 *   File Name    :  findpertsurf_N.cpp
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

 * This version uses Normalized Omega matrix.
 * shoudl give same resutls as findpertsurf_N
 * should work regardless of what method was used to perturb surface

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
#include "omegamatrix_fft.hpp"
#include "gradfvector.hpp"
#include "coilfft.hpp"


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
   enable_option(opt_Nharm);
   enable_option(opt_Mharm);
   enable_option(opt_Bpf);
  enable_option(opt_Btf);

   // parse command line input
 
   if (!parse_cmd(argc, argv))
      return 1;

   // 
   //  ios_base::fmtflags flags = ios_base::right | ios_base::scientific;
   string fext = "txt";
   string ftemp;
   p3vectorformat::textformat(text_nobraces);

   LAvector <double> datavec("datavec");
   datavec.perline(1);
   datavec.textformat(text_nobraces);
   Matrix <double> data("data");
   data.perline(1);
   data.textformat(text_nobraces);

   string fname;

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



   ostringstream strmtmp1;
   strmtmp1 <<Nnn;
   string Nnn_str(strmtmp1.str());
   ostringstream strmtmp2;
   strmtmp2 <<Nmm;
   string Nmm_str(strmtmp2.str());

   // Create Fourier Mode vectors
   LAvector<double> nn("nn");
   LAvector<double> mm("mm");
   unsigned int NF;
   bool mode00 = true;
   modevectors(NF,nn,mm,Nnn,Nmm,Nharm,Mharm,mode00);
   


   // these exclude the n=0,m=0 case
   LAvector<double> nnR("nnR");
   LAvector<double> mmR("mmR");
   unsigned int NFR;
   mode00 = false;
   modevectors(NFR,nnR,mmR,Nnn,Nmm,Nharm,Mharm,mode00);

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
   LAvector<complex<double> > delFluxF(NFR,"delFluxF");
   cout <<endl<< "$ Loading perturbed Plasma Flux sin/cos fourier coefficients from " << flux_filename << endl;
   if (load_coefs(flux_filename,CoefFileFormat_sincos,nnR,mmR,delFluxF,false))
      return 3;
 

 



 
   // lay plasma surface onto grid 
  
   LAvector<p3vector<double> > X(Npts, "X");
   LAvector<p3vector<double> > dA_dtdp(Npts, "dA_dtdp");

   LAvector<p3vector<double> > dx_dr(Npts, "dx_dr");
   LAvector<p3vector<double> > dx_dtheta(Npts,"dx_dtheta");
   LAvector<p3vector<double> > dx_dphi(Npts,"dx_dphi");
   LAvector<p3vector<double> > grad_r(Npts,"grad_r");
   LAvector<p3vector<double> > grad_theta(Npts,"grad_theta");
   LAvector<p3vector<double> > grad_phi(Npts,"grad_phi");

    

   cout << endl;
   cout <<"$ Mapping plasma surface fourier coefficients to "<<Ntheta<<" x "<<Nphi<<" (theta by phi) grid"<<endl;

   STARTTIME(tbuff,ckstart);

   expandsurfaceandbases(X,dA_dtdp,dx_dr,dx_dtheta,dx_dphi,grad_r,grad_theta,grad_phi,plasmafourier,thetas,phis);

   LAvector<double> J(Npts, "J");
   for (unsigned int j =0; j<Npts; j++) {
      J[j] = dot(dx_dr[j],cross(dx_dtheta[j], dx_dphi[j]));
   }

   STOPTIME(tbuff,ckstart);




  // load B field
  cout << endl;
  cout<<"$ Load tangent B field "<<endl;

  STARTTIME(tbuff,ckstart);

  LAvector<p3vector<double> > B(Npts, "B");

  cout <<endl<< "$ Loading BTOTAL_theta fourier coefficients from " << Bt_filename << endl;
  cout <<endl<< "$ Loading BTOTAL_phi fourier coefficients from " << Bp_filename << endl;

  LAvector<complex<double> > BtF(NF,"BtF");
  if (load_coefs( Bt_filename,CoefFileFormat_sincos,nn,mm,BtF))
    return 5;
  LAvector<complex<double> > BpF(NF,"BpF");
  if (load_coefs( Bp_filename,CoefFileFormat_sincos,nn,mm,BpF))
    return 6;

  LAvector<double> Bt(Npts, "Bt");
  expandfunction(Bt,BtF,fs);
  LAvector<double>Bp(Npts, "Bp");
  expandfunction(Bp,BpF,fs);
  
  for (unsigned int j =0; j<Npts; j++)
    B[j] =  Bt[j] * dx_dtheta[j] + Bp[j] * dx_dphi[j];



  STOPTIME(tbuff,ckstart);




   cout << endl;
   cout<<"$ Calculate Omega matrix ("<<NFR<<"x"<<NFR<<")"<<endl;

   STARTTIME(tbuff,ckstart);

   //   NormalizedOmegamatrix_fft(OmegaN, Bt,Bp, fsR,nnR,mmR, Nphi,Ntheta, Nnn,Nmm,Nharm,Mharm);


     Matrix<complex<double> > Omega(NFR,NFR,"Omega");
     Matrix<p3vector<complex<double> > >  grad_fsR(Npts,NFR,"grad_fsR");
     gradfvector(thetas,phis,mmR,nnR,grad_theta,grad_phi,grad_fsR);
     
     NormalizedOmegamatrix_fft_fast(Omega, Bt,Bp, fsR,nnR,mmR, Nphi,Ntheta, Nnn,Nmm,Nharm,Mharm);

   STOPTIME(tbuff,ckstart);



   // Invert Omega Matrix
  
   cout << endl;
   cout<<"$ Invert Omega Matrix ("<<NFR<<"x"<<NFR<<")"<<endl;
     
   STARTTIME(tbuff,ckstart);

   Matrix<complex<double> > OmegaInv(NFR,NFR,"OmegaInv");
   cooll_lapack::inv(Omega,OmegaInv);
    
   STOPTIME(tbuff,ckstart);




   cout << endl;
   cout<<"$ Calculate weighted flux"<<endl;

   LAvector<double> delFlux(Npts, "delFlux");
   expandfunction(delFlux,delFluxF,fsR);

   for(unsigned int i=0; i<Npts;i++) {
      // NOrmalize the flux
      delFlux[i] = delFlux[i]/J[i]/Bp[i];
   }

   transformfunction(delFlux,delFluxF,fsR);


   STOPTIME(tbuff,ckstart);

   // calculate perturbations
  
   cout << endl;
   cout<<"$ Calculate perturbations to plasma surface"<<endl;

   STARTTIME(tbuff,ckstart);


   LAvector<complex<double> > deltaF(NFR,"deltaF");
   deltaF = (OmegaInv|(delFluxF));
   LAvector<double> delta(Npts,"delta");
   expandfunction(delta,deltaF,fsR);


  
   STOPTIME(tbuff,ckstart);


   // calculate new surface coords
  
   cout << endl;
   cout<<"$ Calculate new plasma surface"<<endl;

   STARTTIME(tbuff,ckstart);

   LAvector<p3vector<double> > ksi(Npts,"ksi");
   LAvector<p3vector<double> > Xnew(Npts,"Xnew");
   for(unsigned int i=0; i<Npts;i++) {
      ksi[i] = delta[i] * dx_dr[i];
      Xnew[i] = X[i] + ksi[i];
   }
   STOPTIME(tbuff,ckstart);



  
 

   ostringstream strm;


   massage(deltaF,1e-8);
   fname = "deltaF_N.out";
   save_coefs(fname,CoefFileFormat_sincos,nnR,mmR,deltaF);


   FourierSurface newfourier;
   transformsurface(Xnew,newfourier,fs,nn,mm);
   //  massage(newfourier,1e-8);
   fname = "perturbedsurface_N.out";
   save_fourier_surface(fname,newfourier);


   return 0;
} // main()





