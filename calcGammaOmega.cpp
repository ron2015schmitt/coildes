/************************************************************************* 
 * 
 *   File Name    :  
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  
 * 
 *
 *   SYNOPSIS     
 *     Stellarator COIL Design 
 *     This file finds the coil current from given plasma surface flux 
 *    
 *
 * VERSION NOTES:
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
#include "gammamatrix.hpp"
#include "omegamatrix.hpp"
#include "gradfvector.hpp"
#include "coilfft.hpp"

// This is the file that defines the B field configuration



// Main Function for code

int main (int argc, char *argv[])
{

   disable_all_options();
   enable_option(opt_pf);
   enable_option(opt_Nphi);
   enable_option(opt_Ntheta);
   enable_option(opt_Nnn);
   enable_option(opt_Nmm);
   enable_option(opt_Btf);
   enable_option(opt_Bpf);
   enable_option(opt_Nharm);
   enable_option(opt_Mharm);
   enable_option(opt_omegaf);
   enable_option(opt_gammaf);

   // parse command line input
 
   if (!parse_cmd(argc, argv))
      return 1;


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


   // Create angle grid
   const unsigned int Npts = Ntheta*Nphi;
   Vector<double> thetas(Npts,"thetas");
   Vector<double> phis(Npts,"phis");
   anglevectors(thetas, phis, Ntheta, Nphi);



  
   ostringstream strmtmp1;
   strmtmp1 <<Nnn;
   string Nnn_str(strmtmp1.str());
   ostringstream strmtmp2;
   strmtmp2 <<Nmm;
   string Nmm_str(strmtmp2.str());


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


   // coefficient C is the integration coef for the fourier transform
   // C = dtheta*dphi
   //   = (2*pi/Ntheta)*(2*pi/Nphi)
   //   = (2*pi*2*pi/Npts)
   const double C = (2*PI*2*PI/double(Npts));
   const double Csqr = C*C;

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


   STOPTIME(tbuff,ckstart);


 
   // Create Gamma Matrix
  
   Matrix<complex<double> > Gamma("Gamma");
  
   // if gamma filename given, then load gamma (SVD form) from files
   if (!gamma_filename.empty()) {
      cout << endl;
      cout<<"$ Load Gamma matrix ("<<NFR<<"x"<<(2*NF)<<")";
      cout << " from "<< gamma_filename <<".R.out";
      cout << " and "<< gamma_filename <<".R.out";
      cout<<endl;

      STARTTIME(tbuff,ckstart);

      data.perline(1);
      data.textformat(text_nobraces);
      data.resize(NFR,2*NF);
      fname = gamma_filename + ".R.out";
      load(data,fname);

      Matrix <double> data2("data2");
      data2.perline(1);
      data2.textformat(text_nobraces);
      data2.resize(NFR,2*NF);
      fname = gamma_filename + ".I.out";
      load(data2,fname);
      Gamma = mcomplex(data,data2);
      data.clear();
      data2.clear();

      STOPTIME(tbuff,ckstart);

      ///////////////////////////////////////////////////////
      // ********DEBUG*************************
      // save Gamma 
      printcr("$ Saving Gamma Matrix");
      data.resize() = real(Gamma);
      fname = "gamma.Nn="+Nnn_str+".Nm="+Nmm_str+".R.debug.out";
      dispcr(fname);
      save(data,fname);
      data.resize() = imag(Gamma);
      fname = "gamma.Nn="+Nnn_str+".Nm="+Nmm_str+".I.debug.out";
      dispcr(fname);
      save(data,fname);
      //////////////////////////////////////////////////////


   } else {
      // Calculate Gamma matrix

      cout << endl;
      cout<<"$ Generate Gamma matrix "<<endl;

      STARTTIME(tbuff,ckstart);

      fullgammamatrix(Gamma,nn,mm,fs,fsR,dA_dtdp,thetas,phis,Ntheta,Nphi);

      const unsigned int Ns = Gamma.Ncols();
      dispcr(Gamma.Nrows());
      dispcr(Gamma.Ncols());
      dispcr(NFR);
      dispcr(Ns);
      STOPTIME(tbuff,ckstart);
  

 
      // save Gamma 
      printcr("$ Saving Gamma Matrix");
      data.resize() = real(Gamma);
      fname = "gamma.Nn="+Nnn_str+".Nm="+Nmm_str+".R.out";
      dispcr(fname);
      save(data,fname);
      data.resize() = imag(Gamma);
      fname = "gamma.Nn="+Nnn_str+".Nm="+Nmm_str+".I.out";
      dispcr(fname);
      save(data,fname);


   }





   Matrix<complex<double> > Omega(NFR,NFR,"Omega");

   // if omega filename given, then load omega (SVD form) from files
   if (!omega_filename.empty()) {
      cout << endl;
      cout<<"$ Load Omega matrix ("<<NFR<<"x"<<NFR<<")";
      cout << " from "<< omega_filename <<".R.out";
      cout << " and "<< omega_filename <<".R.out";
      cout<<endl;

      STARTTIME(tbuff,ckstart);

      data.perline(1);
      data.textformat(text_nobraces);
      data.resize(NFR,NFR);
      fname = omega_filename + ".R.out";
      load(data,fname);

      Matrix <double> data2("data2");
      data2.perline(1);
      data2.textformat(text_nobraces);
      data2.resize(NFR,NFR);
      fname = omega_filename + ".I.out";
      load(data2,fname);
      Omega = mcomplex(data,data2);
      data.clear();
      data2.clear();

      STOPTIME(tbuff,ckstart);

      ///////////////////////////////////////////////////////
      // ********DEBUG*************************
      // save Omega 
      printcr("$ Saving Omega Matrix");
      data.resize() = real(Omega);
      fname = "omega.Nn="+Nnn_str+".Nm="+Nmm_str+".R.debug.out";
      dispcr(fname);
      save(data,fname);
      data.resize() = imag(Omega);
      fname = "omega.Nn="+Nnn_str+".Nm="+Nmm_str+".I.debug.out";
      dispcr(fname);
      save(data,fname);
      //////////////////////////////////////////////////////


   } else {
      // Calculate Omega matrix

      // first we need to load B field
      cout << endl;
      cout<<"$ Load tangent B field "<<endl;

      STARTTIME(tbuff,ckstart);
     
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
      Vector<double> Bp(Npts, "Bp");
      expandfunction(Bp,BpF,fs);
     
      Vector<p3vector<double> > B(Npts, "B");
      for (unsigned int j =0; j<Npts; j++)
	 B[j] =  Bt[j] * dx_dtheta[j] + Bp[j] * dx_dphi[j];

      STOPTIME(tbuff,ckstart);

     
      cout << endl;
      cout<<"$ Calculate Omega matrix ("<<NFR<<"x"<<NFR<<")"<<endl;
     
      STARTTIME(tbuff,ckstart);
     
      Matrix<p3vector<complex<double> > >  grad_fsR(Npts,NFR,"grad_fsR");
      gradfvector(thetas,phis,mmR,nnR,grad_theta,grad_phi,grad_fsR);
     
      //  omegamatrix(Omega,B,grad_fsR,f_delta,dx_dr,dx_dtheta,dx_dphi,grad_r,grad_theta,grad_phi);
      omegamatrix(Omega,B,grad_fsR,fsR,dx_dr,dx_dtheta,dx_dphi,grad_r,grad_theta,grad_phi);
      

      STOPTIME(tbuff,ckstart);
   }

 
   // save Omega 
   printcr("$ Saving Omega Matrix");
   data.resize() = real(Omega);
   fname = "omega.Nn="+Nnn_str+".Nm="+Nmm_str+".R.out";
   dispcr(fname);
   save(data,fname);
   data.resize() = imag(Omega);
   fname = "omega.Nn="+Nnn_str+".Nm="+Nmm_str+".I.out";
   dispcr(fname);
   save(data,fname);

   // clear some variables

   fsR.clear();
   fs.clear();


   // note that Gamma seems to always be real and Omega seems to always be imaginary.
   // using this symmetry would greatly reduce computation

   cout << endl;
   cout<<"$ Calculate (Omega|Gamma) matrix ("<<NFR<<"x"<<(2*NF)<<")"<<endl;
     
   STARTTIME(tbuff,ckstart);

   Matrix<complex<double> > OG(NFR,(2*NF),"OG");

   OG=(Omega|Gamma);

   STOPTIME(tbuff,ckstart);

 
   // save OmegaGamma 
   printcr("$ Saving OmegaGamma Matrix");
   data.resize() = real(OG);
   fname = "omegagamma.Nn="+Nnn_str+".Nm="+Nmm_str+".R.out";
   dispcr(fname);
   save(data,fname);
   data.resize() = imag(OG);
   fname = "omegagamma.Nn="+Nnn_str+".Nm="+Nmm_str+".I.out";
   dispcr(fname);
   save(data,fname);


   return 0;
} // main()





