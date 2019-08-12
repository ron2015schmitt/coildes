/************************************************************************* 
 * 
 *   File Name    :  residualmodes_XVII
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
 *  This version was based on scoild_fft_XVII.cpp on 2006 oct 18
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
   enable_option(opt_ff);
   enable_option(opt_Nphi);
   enable_option(opt_Ntheta);
   enable_option(opt_Nnn);
   enable_option(opt_Nmm);
   enable_option(opt_Btf);
   enable_option(opt_Bpf);
   enable_option(opt_Nharm);
   enable_option(opt_Mharm);
   enable_option(opt_omegaf);
   enable_option(opt_lambdaf);
   //    enable_option(opt_Npert);
   enable_option(opt_Nprin);
   //    enable_option(opt_Ncomb);

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


   if (Nprin>NFR) {
      cout << "Nprin="<<Nprin<<" is larger than the number of Fourier modes ("<<NFR<<")."<<endl;
      cout << "Nprin has been changed to equal the number of Fourier modes (Nprin="<<NFR<<")."<<endl;
      Nprin = NFR;
   }
   strmtemp  <<Nprin;
   string Nprin_str(strmtemp.str());
   strmtemp.str("");
   dispcr(Nprin);

   // Create angle grid
   const unsigned int Npts = Ntheta*Nphi;
   LAvector<double> thetas(Npts,"thetas");
   LAvector<double> phis(Npts,"phis");
   anglevectors(thetas, phis, Ntheta, Nphi);

   const double dtheta = 2*M_PI/((double)Ntheta);
   const double dphi = 2*M_PI/((double)Nphi);
   const double dphi_by_dtheta = 2.0*PI * 2.0*PI/double(Npts);



  
   ostringstream strmtemp;
   strmtemp <<Nnn;
   string Nnn_str(strmtemp.str());
   strmtemp.str("");
   strmtemp <<Nmm;
   string Nmm_str(strmtemp.str());
   strmtemp.str("");


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


   // coefficient C is the integration coef for the fourier transform
   // C = dtheta*dphi
   //   = (2*pi/Ntheta)*(2*pi/Nphi)
   //   = (2*pi*2*pi/Npts)
   const double C = (2*PI*2*PI/double(Npts));
   const double Csqr = C*C;





   LAvector<complex<double> > FluxF(NFR,"FluxF");

   cout <<endl<< "$ Loading Plasma Flux sin/cos fourier coefficients from " << flux_filename << endl;
   if (load_coefs(flux_filename,CoefFileFormat_sincos,nnR,mmR,FluxF,false))
      return 3;


 
   Matrix<complex<double> > P(NFR,NFR,"V");
   LAvector<complex<double> >  W(NFR,"W");

   // if omega filename given, then load omega (EVD form) from files
   if (!omega_filename.empty()) {
      cout << endl;
      cout<<"$ Load Omega matrix ("<<NFR<<"x"<<NFR<<")"<<endl;

      STARTTIME(tbuff,ckstart);

      // UR seems to always be effectively zero.  do we need it?

      // PRINT OUT FILENAMES!!!
      data.perline(1);
      data.textformat(text_nobraces);
      data.resize(NFR,NFR);
      fname = omega_filename + ".PR.out";
      load(data,fname);

      Matrix <double> data2("data2");
      data2.perline(1);
      data2.textformat(text_nobraces);
      data2.resize(NFR,NFR);
      fname = omega_filename + ".PI.out";
      load(data2,fname);
      P = mcomplex(data,data2);
      data.clear();
      data2.clear();

      datavec.perline(1);
      datavec.textformat(text_nobraces);
      datavec.resize(NFR);
      fname = omega_filename + ".WR.out";
      load(datavec,fname);

      LAvector <double> datavec2("datavec2");
      datavec2.perline(1);
      datavec2.textformat(text_nobraces);
      datavec2.resize(NFR);
      fname = omega_filename + ".WI.out";
      load(datavec2,fname);
      W = vcomplex(datavec,datavec2);
      datavec.clear();
      datavec2.clear();

      STOPTIME(tbuff,ckstart);

   } else {
      // Calculate Omega matrix

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

      // print coef's
      //  printfouriercoefs(plasmafourier.nn(),plasmafourier.mm(),plasmafourier.RF(),plasmafourier.ZF(),10,18);
  
      // load B field
      cout << endl;
      cout<<"$ Load tangent B field (total equilibrium B field)"<<endl;

      STARTTIME(tbuff,ckstart);

      LAvector<p3vector<double> > B(Npts, "B");

      cout <<endl<< "$ Loading BTOTAL_theta fourier coefficients from " << Bt_filename << endl;

      LAvector<complex<double> > BtF(NF,"BtF");
      if (load_coefs( Bt_filename,CoefFileFormat_sincos,nn,mm,BtF,false))
	 return 5;


      cout <<endl<< "$ Loading BTOTAL_phi fourier coefficients from " << Bp_filename << endl;
      LAvector<complex<double> > BpF(NF,"BpF");
      if (load_coefs( Bp_filename,CoefFileFormat_sincos,nn,mm,BpF,false))
	 return 6;


      LAvector<complex<double> > lambdaF(NF,"lambdaF");
      cout <<endl<< "$ Loading lambda sin/cos fourier coefficients from " << lambda_filename << endl;
      if (load_coefs(lambda_filename,CoefFileFormat_sincos,nn,mm,lambdaF,false))
	 return 3;

      LAvector<double> Bt(Npts, "Bt");
      expandfunction(Bt,BtF,fs);
      LAvector<double> Bp(Npts, "Bp");
      expandfunction(Bp,BpF,fs);

      for (unsigned int j =0; j<Npts; j++)
	 B[j] =  Bt[j] * dx_dtheta[j] + Bp[j] * dx_dphi[j];

      STOPTIME(tbuff,ckstart);



      cr();
      printcr("Find contravariant * J vector components.");

      STARTTIME(tbuff,ckstart);
      LAvector<double> JBt(Npts,"JBt");
      LAvector<double> JBp(Npts,"JBp");
 
      for (unsigned int j =0; j<Npts; j++) {
	 JBt[j] = J[j] * Bt[j];
	 JBp[j] = J[j] * Bp[j];
      }
      STOPTIME(tbuff,ckstart);


      LAvector<complex<double> > JBtF(NF,"JBtF");
      transformfunction(JBt,JBtF,fs);
      LAvector<complex<double> > JBpF(NF,"JBpF");
      transformfunction(JBp,JBpF,fs);

      cout << endl;
      cout<<"$ Calculate Omega matrix ("<<NFR<<"x"<<NFR<<")"<<endl;
     
      STARTTIME(tbuff,ckstart);
     
      const LAvector<unsigned int> tmp = findtrue((nn==0)&&(mm==0));
      const unsigned int ind00 = tmp[0];

      dispcr(JBtF[ind00].real());
      dispcr(JBpF[ind00].real());

      lambdaF[ind00] = 2*PI*complex<double>(JBtF[ind00].real(),JBpF[ind00].real());

      Matrix<complex<double> > Omega(NFR,NFR,"Omega");
      lambda_omegamatrix(Omega,lambdaF,nn,mm);

      STOPTIME(tbuff,ckstart);


      //SAVE OMEGA MATRIX
      printcr("$ Saving Omega Matrix");
      data.resize() = real(Omega);
      fname = "omega.XVII.Nn="+Nnn_str+".Nm="+Nmm_str+".R.out";
      dispcr(fname);
      save(data,fname);
      data.resize() = imag(Omega);
      fname = "omega.XVII.Nn="+Nnn_str+".Nm="+Nmm_str+".I.out";
      dispcr(fname);
      save(data,fname);




      // EV decomposition of Omega Matrix
  
      cout << endl;
      cout<<"$ EV decomposition of Omega matrix ("<<NFR<<"x"<<NFR<<")"<<endl;
     
      STARTTIME(tbuff,ckstart);


      cooll_lapack::eigwincsort(Omega,P,W);
      Omega.clear();
     
      STOPTIME(tbuff,ckstart);


      // save  Omega in eigen form
      printcr("$ Saving  Omega Matrix in eigen form");
      data.resize() = real(P);
      fname = "scoild_XVII.Nn="+Nnn_str+".Nm="+Nmm_str+".PR.out";
      dispcr(fname);
      save(data,fname);
      data.resize() = imag(P);
      fname = "scoild_XVII.Nn="+Nnn_str+".Nm="+Nmm_str+".PI.out";
      dispcr(fname);
      save(data,fname);

      datavec.resize() = real(W);
      datavec.perline(1);
      fname = "scoild_XVII.Nn="+Nnn_str+".Nm="+Nmm_str+".WR.out";
      dispcr(fname);
      save(datavec,fname);
      datavec.resize() = imag(W);
      datavec.perline(1);
      fname = "scoild_XVII.Nn="+Nnn_str+".Nm="+Nmm_str+".WI.out";
      dispcr(fname);
      save(datavec,fname);

   }


   // Invert Omega Eigenvector Matrix
  
   cout << endl;
   cout<<"$ Invert Omega Eigenvector Matrix ("<<NFR<<"x"<<NFR<<")"<<endl;
     
   STARTTIME(tbuff,ckstart);

   Matrix<complex<double> > Pinv(NFR,NFR,"Pinv");
   cooll_lapack::inv(P,Pinv);
     
   STOPTIME(tbuff,ckstart);




   // "prin" modes

   // Project FluxF into perturbation eigenspace
  
   cout << endl;
   cout<<"$ Project FluxF into perturbation eigenspace"<<endl;
     
   STARTTIME(tbuff,ckstart);

   LAvector<complex<double> > FluxMF(NFR,"FluxMF");
   FluxMF = (Pinv|FluxF);

   STOPTIME(tbuff,ckstart);

   LAvector<unsigned int> ii_prin(NFR,"ii_prin");
   LAvector<double> magn_prin(NFR,"magn_prin");
   LAvector<double> magn_percent_prin(NFR,"magn_percent_prin");
   LAvector<double> magn_cumpercent_prin(NFR,"magn_cumpercent_prin");

   magn_prin = abs(FluxMF);
   double ss_prin = COOLL::sum(magn_prin*magn_prin);
   ii_prin = sortrevwind(magn_prin);

   LAvector<complex<double> > residualFluxMF(NFR,"residualFluxMF");

   for (unsigned int j =0; j<Nprin; j++) 
      residualFluxMF[ii_prin[j]] = 0;

   for (unsigned int j =Nprin; j<NFR; j++) 
      residualFluxMF[ii_prin[j]] = FluxMF[ii_prin[j]];


   /////////////////////////////////
   FluxMF.perline(1);
   FluxMF.textformat(text_nobraces);
   save(FluxMF,"FluxMF.out");

   residualFluxMF.perline(1);
   residualFluxMF.textformat(text_nobraces);
   save(residualFluxMF,"residualFluxMF.out");
   /////////////////////////////////   			   

   magn_prin = abs(residualFluxMF);
   ss_prin = COOLL::sum(magn_prin*magn_prin);
   ii_prin = sortrevwind(magn_prin);


   magn_percent_prin = (magn_prin*magn_prin)/ss_prin;

   for (unsigned int j =0; j<NFR; j++) {
      if (j==0)
	 magn_cumpercent_prin[j] = magn_percent_prin[j];
      else
	 magn_cumpercent_prin[j] = magn_cumpercent_prin[j-1] + magn_percent_prin[j];
   }
   

   ii_prin.perline(1);
   ii_prin.textformat(text_nobraces);
   save(ii_prin,"residualFlux_prin_modenum_pertsubspace.out");

   magn_prin.perline(1);
   magn_prin.textformat(text_nobraces);
   save(magn_prin,"residualFlux_prin_mag_pertsubspace.out");


   magn_percent_prin.perline(1);
   magn_percent_prin.textformat(text_nobraces);
   save(magn_percent_prin,"residualFlux_prin_magpercent_pertsubspace.out");


   magn_cumpercent_prin.perline(1);
   magn_cumpercent_prin.textformat(text_nobraces);
   save(magn_cumpercent_prin,"residualFlux_prin_magcumper_pertsubspace.out");



   // "comb" modes


   // Multiply residualFluxF in magnetic coords by (inverted) perturbation eigenvalues
  
   cout << endl;
   cout<<"$ Multiply residualFluxF in magnetic coords by perturbation eigenvalues"<<endl;
   cout<<"$ This gives predicted surface pertubation"<<endl;

   LAvector<unsigned int> ii_comb(NFR,"ii_comb");
   LAvector<double> magn_comb(NFR,"magn_comb");
   LAvector<double> magn_percent_comb(NFR,"magn_percent_comb");
   LAvector<double> magn_cumpercent_comb(NFR,"magn_cumpercent_comb");

     
   STARTTIME(tbuff,ckstart);

   LAvector<complex<double> > residualFluxMpertF(NFR,"residualFluxMpertF");
   residualFluxMpertF = residualFluxMF / W;

   STOPTIME(tbuff,ckstart);


   magn_comb = abs(residualFluxMpertF);
   double ss_comb = COOLL::sum(magn_comb*magn_comb);
   ii_comb = sortrevwind(magn_comb);
   magn_percent_comb = (magn_comb*magn_comb)/ss_comb;

   for (unsigned int j =0; j<NFR; j++) {
      if (j==0)
	 magn_cumpercent_comb[j] = magn_percent_comb[j];
      else
	 magn_cumpercent_comb[j] = magn_cumpercent_comb[j-1] + magn_percent_comb[j];
   }

   ii_comb.perline(1);
   ii_comb.textformat(text_nobraces);
   save(ii_comb,"residualFlux_comb_modenum_pertsubspace.out");

   magn_comb.perline(1);
   magn_comb.textformat(text_nobraces);
   save(magn_comb,"residualFlux_comb_mag_pertsubspace.out");


   magn_percent_comb.perline(1);
   magn_percent_comb.textformat(text_nobraces);
   save(magn_percent_comb,"residualFlux_comb_magpercent_pertsubspace.out");


   magn_cumpercent_comb.perline(1);
   magn_cumpercent_comb.textformat(text_nobraces);
   save(magn_cumpercent_comb,"residualFlux_comb_magcumper_pertsubspace.out");



   



   return 0;
} // main()





