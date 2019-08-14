/************************************************************************* 
 * 
 *   File Name    :  
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  16 May 2003
 * 
 *
 *   SYNOPSIS     
 *     Calculate flux perturbation on surface exactly using numerical techniques, 
 *     capturing nonlinear behavior.
 *
 **************************************************************************/



// Standard C libraries
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <ctime>

// Standard C++ libraries

#include <iostream>
#include <complex>

using namespace std;



// coil libraries
#include "coils.hpp"
#include "coilio.hpp"
#include "coils_cmdline.hpp"
#include "surface.hpp"

// generic file for generating the coil+plasma bfield
#include "bfield_coils.hpp"



// Main Function for code

int main (int argc, char *argv[])
{

   // display Matricks mode
   cout << endl;
   display_execution_mode();
   cout << endl;

   disable_all_options();
   enable_option(opt_pf);
   enable_option(opt_cf);
   enable_option(opt_if);
   enable_option(opt_Nphi);
   enable_option(opt_Ntheta);
   enable_option(opt_Nnn);
   enable_option(opt_Nmm);
   enable_option(opt_Itoroidal);
   enable_option(opt_Ipoloidal);
   enable_option(opt_Nharm);
   enable_option(opt_Mharm);
   enable_option(opt_pertvalZ);
   enable_option(opt_pertvalR);



   // parse command line input
 
   if (!parse_cmd(argc, argv))
      return 1;

   string fext = "out";
   //  string fext = plasma_extname;
   //  ios_base::fmtflags flags = ios_base::right | ios_base::scientific;
   string ftemp;
   p3vectorformat::textformat(text_nobraces);

   string fname;
   Matrix <double> dataM("dataM");
   dataM.perline(1);
   dataM.textformat(text_nobraces);
   Vector <double> dataV("dataV");
   dataV.perline(1);
   dataV.textformat(text_nobraces);

   ostringstream strmtmp1;
   strmtmp1 <<Nnn;
   string Nnn_str(strmtmp1.str());
   ostringstream strmtmp2;
   strmtmp2 <<Nmm;
   string Nmm_str(strmtmp2.str());

   ostringstream strmtmp3;
   strmtmp3 << perturbation_valueR;
   string pertR_str(strmtmp3.str());

   ostringstream strmtmp4;
   strmtmp4 << perturbation_valueZ;
   string pertZ_str(strmtmp4.str());

   string tail_str = "Nn="+Nnn_str+".Nm="+Nmm_str+".pR="+pertR_str+".pZ="+pertZ_str;

   // variables for measuring times
   struct tms tbuff;
   clock_t ckstart;




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

   // load coil surface and coil currents. run setup for biot-savrt calcs

   setupcoils(coil_filename, current_filename, Itoroidal, Ipoloidal, 
	      Nnn, Nmm, Nphi, Ntheta, Nharm, Mharm);
   dispcr(Itoroidal);
   dispcr(Ipoloidal);

   // load the plasma surface fourier coef's

   cout << "$ Loading PLASMA SURFACE fourier coefficients from " << plasma_filename << endl;
   FourierSurface plasmafourier0;
   if (load_fourier_surface(plasma_filename,plasmafourier0))
      return 1;
   plasmafourier0.RF().name("p0.RF");
   plasmafourier0.ZF().name("p0.ZF");

 
   // lay plasma surface onto grid 
  
   Vector<p3vector<double> > X0(Npts, "X0");
   Vector<p3vector<double> > dA_dtdp0(Npts, "dAdtdp0");
   Vector<p3vector<double> > dx_dr0(Npts, "dx_dr0");
   Vector<p3vector<double> > dx_dtheta0(Npts,"dx_dtheta0");
   Vector<p3vector<double> > dx_dphi0(Npts,"dx_dphi0");
   Vector<p3vector<double> > grad_r0(Npts,"grad_r0");
   Vector<p3vector<double> > grad_theta0(Npts,"grad_theta0");
   Vector<p3vector<double> > grad_phi0(Npts,"grad_phi0");
 


   cout << endl;
   cout <<"$ Mapping plasma surface fourier coefficients to "<<Ntheta<<" x "<<Nphi<<" (theta by phi) grid"<<endl;

   STARTTIME(tbuff,ckstart);

   expandsurfaceandbases(X0,dA_dtdp0,dx_dr0,dx_dtheta0,dx_dphi0,grad_r0,grad_theta0,grad_phi0,plasmafourier0,thetas,phis);

   STOPTIME(tbuff,ckstart);


   ///////////////////////////////////////////////////////////
   //////////////////DEBUGGING//////////////////////////////
   // save surface to files
   //    fname = "X0.out";
   //    dispcr(fname);
   //    X0.perline(1);
   //    X0.textformat(text_nobraces);
   //    save(X0,fname);

   //    fname = "dA_dtdp0.out";
   //    dispcr(fname);
   //    dA_dtdp0.perline(1);
   //    dA_dtdp0.textformat(text_nobraces);
   //    save(dA_dtdp0,fname);

   //    fname = "dx_dr0.out";
   //    dispcr(fname);
   //    dx_dr0.perline(1);
   //    dx_dr0.textformat(text_nobraces);
   //    save(dx_dr0,fname);

   //    fname = "dx_dtheta0.out";
   //    dispcr(fname);
   //    dx_dtheta0.perline(1);
   //    dx_dtheta0.textformat(text_nobraces);
   //    save(dx_dtheta0,fname);

   //    fname = "dx_dphi0.out";
   //    dispcr(fname);
   //    dx_dphi0.perline(1);
   //    dx_dphi0.textformat(text_nobraces);
   //    save(dx_dphi0,fname);

   //    fname = "grad_r0.out";
   //    dispcr(fname);
   //    grad_r0.perline(1);
   //    grad_r0.textformat(text_nobraces);
   //    save(grad_r0,fname);

   //    fname = "grad_theta0.out";
   //    dispcr(fname);
   //    grad_theta0.perline(1);
   //    grad_theta0.textformat(text_nobraces);
   //    save(grad_theta0,fname);

   //    fname = "grad_phi0.out";
   //    dispcr(fname);
   //    grad_phi0.perline(1);
   //    grad_phi0.textformat(text_nobraces);
   //    save(grad_phi0,fname);
   ///////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////



   // calculate magnetic field on plasma surface
    
   cout << endl;
   cout<<"$ Calc B field of unperturbed field"<<endl;
 
   STARTTIME(tbuff,ckstart);
     
   Vector<p3vector<double> > Bcoils0(Npts, "Bcoils0");
   Vector<p3vector<double> > Btotal0(Npts, "Btotal0");
   unsigned int count0 = 0;
   unsigned int checkup0 = static_cast<unsigned int>(Npts*0.01);
   for (unsigned int j =0; j<Npts; j++) {
      if (++count0 == checkup0) {
	 print(Matricks::round(double(j)/Npts*100));cout <<" %"<<endl;
	 count0 =0;
      }
      bTotalandbCoils(X0[j],Btotal0[j], Bcoils0[j]);
   }

   STOPTIME(tbuff,ckstart);


   Vector<p3vector<double> > v_normal0(Npts, "v_normal0");
   for (unsigned int j =0; j<Npts; j++)
      v_normal0[j] = dA_dtdp0[j]/norm(dA_dtdp0[j]);

   Vector<double> Flux0(Npts, "Flux0");
 
   for (unsigned int j =0; j<Npts; j++) {
      Flux0[j] = dot(Btotal0[j], dA_dtdp0[j]);
   }
   Vector<complex<double> > FluxF0(NFR,"FluxF0");
   transformfunction(Flux0,FluxF0,fsR);


   printcr("$ Saving unperturbed flux");
   dataV.resize() = real(FluxF0);
   fname = "FluxF0."+tail_str+".R.out";
   dispcr(fname);
   save(dataV,fname);
   dataV.resize() = imag(FluxF0);
   fname = "FluxF0."+tail_str+".I.out";
   dispcr(fname);
   save(dataV,fname);



   double Flux_RMS0 = 0 ;
   for (unsigned int j =0; j<Npts; j++) {
      double temp;
      Flux_RMS0 += Flux0[j]*Flux0[j];
   }
   Flux_RMS0 = sqrt(Flux_RMS0/Npts);
   dispcr(Flux_RMS0);
   double FluxF_RMS0 = 0 ;
   for (unsigned int k =0; k<NFR; k++) {
      double temp;
      temp = abs(FluxF0[k]);
      FluxF_RMS0 += temp*temp;
   }
   FluxF_RMS0 = sqrt(FluxF_RMS0)/(2*PI);
   dispcr(FluxF_RMS0);cr();

   printcr("$ Saving unperturbed RMS flux");
   dataV.resize(2);
   dataV[0] =  Flux_RMS0;
   dataV[1] =  FluxF_RMS0;
   fname = "flux.rms."+tail_str+".out";
   dispcr(fname);
   save(dataV,fname);


   Matrix<complex<double> > OG(NFR,(2*NF),"OG");
   Vector<double> dFlux_RMS(2*NF, "dFlux_RMS");
   Vector<double> dFluxF_RMS(2*NF, "dFluxF_RMS");


   Vector<double> kvec(2*NF, "kvec");
   Vector<double> mvec(2*NF, "mvec");
   Vector<double> nvec(2*NF, "nvec");

   cout << endl;
   cout<<"$ Now perturb each fourier mode..."<<endl;
   dispcr(perturbation_valueR);
   dispcr(perturbation_valueZ);


   // Loop for each perturbation mode

   for (unsigned int RZ =0; RZ<=1; RZ++) {
      for (unsigned int k =0; k<NF; k++) {

	 double npert = nn[k];

	 double mpert = mm[k];


	 const unsigned int mode_index= k+RZ*NF;

	 nvec[mode_index] = npert;
	 mvec[mode_index] = mpert;

	 Vector<p3vector<double> > X(Npts, "X");
	 Vector<p3vector<double> > dA_dtdp(Npts, "dA_dtdp");
	 Vector<p3vector<double> > dx_dr(Npts, "dx_dr");
	 Vector<p3vector<double> > dx_dtheta(Npts,"dx_dtheta");
	 Vector<p3vector<double> > dx_dphi(Npts,"dx_dphi");
	 Vector<p3vector<double> > grad_r(Npts,"grad_r");
	 Vector<p3vector<double> > grad_theta(Npts,"grad_theta");
	 Vector<p3vector<double> > grad_phi(Npts,"grad_phi");

	 const unsigned int NFS = plasmafourier0.mm().size();
	 FourierSurface plasmafourier;
	 plasmafourier.RF().name("p.RF");
	 plasmafourier.ZF().name("p.ZF");

	 cout << endl;
	 if (RZ==0) {
	    cout<<"$ Calc surface perturbations for cosine mode (RF) n="<<npert<<", m="<<mpert<<endl;
	 }else {
	    cout<<"$ Calc surface perturbations for sine mode (ZF) n="<<npert<<", m="<<mpert<<endl;
	 }

	 STARTTIME(tbuff,ckstart);

	 bool found = false;
	 for (unsigned int ff = 0; ff<NFS; ff++) {
	    if ( (mpert ==  plasmafourier0.mm()[ff]) && (npert ==  plasmafourier0.nn()[ff]) ) 
	       found = true;
	 }

	 if (!found) {
	    plasmafourier.mm().resize(NFS+1);
	    plasmafourier.nn().resize(NFS+1);
	    plasmafourier.RF().resize(NFS+1);
	    plasmafourier.ZF().resize(NFS+1);
	    plasmafourier.mm()[NFS] = mpert;
	    plasmafourier.nn()[NFS] = npert;
	    if (RZ==0) 
	       plasmafourier.RF()[NFS] = perturbation_valueR;
	    else
	       plasmafourier.ZF()[NFS] = perturbation_valueZ;
	 } else {
	    plasmafourier.mm().resize(NFS);
	    plasmafourier.nn().resize(NFS);
	    plasmafourier.RF().resize(NFS);
	    plasmafourier.ZF().resize(NFS);
	 }

	 for (unsigned int ff = 0; ff<NFS; ff++) {
	    plasmafourier.mm()[ff] =  plasmafourier0.mm()[ff];
	    plasmafourier.nn()[ff] =  plasmafourier0.nn()[ff];
	    plasmafourier.RF()[ff] =  plasmafourier0.RF()[ff];
	    plasmafourier.ZF()[ff] =  plasmafourier0.ZF()[ff];
	    if ( (mpert ==  plasmafourier0.mm()[ff]) && (npert ==  plasmafourier0.nn()[ff]) ) {
	       if (RZ==0) 
		  plasmafourier.RF()[ff] += perturbation_valueR;
	       else
		  plasmafourier.ZF()[ff] += perturbation_valueZ;
	    }
	 }

	 printfouriercoefs( plasmafourier.nn(), plasmafourier.mm(), plasmafourier.RF(), plasmafourier.ZF(), 5);

	 expandsurfaceandbases(X,dA_dtdp,dx_dr,dx_dtheta,dx_dphi,grad_r,grad_theta,grad_phi,plasmafourier,thetas,phis);
	 
	 STOPTIME(tbuff,ckstart);
    

	 // calculate magnetic field on plasma surface
    
	 cout << endl;
	 if (RZ==0) {
	    cout<<"$ Calc B field for cosine mode (RF) n="<<npert<<", m="<<mpert<<endl;
	 }else {
	    cout<<"$ Calc B field for sine mode (ZF) n="<<npert<<", m="<<mpert<<endl;
	 }
 
	 STARTTIME(tbuff,ckstart);
     
	 Vector<p3vector<double> > Bcoils(Npts, "Bcoils");
	 Vector<p3vector<double> > Btotal(Npts, "Btotal");
	 unsigned int count = 0;
	 unsigned int checkup = static_cast<unsigned int>(Npts*0.01);
	 for (unsigned int j =0; j<Npts; j++) {
// 	    if (++count == checkup) {
// 	       print(Matricks::round(double(j)/Npts*100));cout <<" %"<<endl;
// 	       count =0;
// 	    }
	    bTotalandbCoils(X[j],Btotal[j], Bcoils[j]);
	 }

	 STOPTIME(tbuff,ckstart);


	 Vector<p3vector<double> > v_normal(Npts, "v_normal");
	 for (unsigned int j =0; j<Npts; j++)
	    v_normal[j] = dA_dtdp[j]/norm(dA_dtdp[j]);

	 Vector<double> Flux(Npts, "Flux");
	 Vector<double> dFlux(Npts, "dFlux");
 
	 // calculate flux, then delta flux
	 for (unsigned int j =0; j<Npts; j++) {
	    Flux[j] = dot(Btotal[j], dA_dtdp[j]);
	    dFlux[j] = Flux[j] - Flux0[j];
	 }
  


	 // Create fourier series for flux
	 Vector<complex<double> > dFluxF(NFR,"dFluxF");
	 transformfunction(dFlux,dFluxF,fsR);

	 for (unsigned int r =0; r<NFR; r++) 
	    OG(r,mode_index) = dFluxF[r];

	 //************CALC RMS VALUES TOO
	 double temp = 0 ;
	 for (unsigned int j =0; j<Npts; j++) {
	    temp += dFlux[j]*dFlux[j];
	 }
	 temp = sqrt(temp/Npts);
	 dFlux_RMS[mode_index] = temp;

	 dispcr(dFlux_RMS[mode_index]);

	 temp = 0;
	 for (unsigned int k =0; k<NFR; k++) {
	    double fluxsqr = abs(dFluxF[k]);
	    temp += fluxsqr*fluxsqr;
	 }
	 temp = sqrt(temp)/(2*PI);
	 dFluxF_RMS[mode_index] = temp;

	 dispcr(dFluxF_RMS[mode_index]);

      } // end of mode for loop
   } // end of RZ for loop





   // save OmegaGamma 
   printcr("$ Saving ExactOmegaGamma Matrix");
   dataM.resize() = real(OG);
   fname = "exact.omegagamma."+tail_str+".R.out";
   dispcr(fname);
   save(dataM,fname);
   dataM.resize() = imag(OG);
   fname = "exact.omegagamma."+tail_str+".I.out";
   dispcr(fname);
   save(dataM,fname);

   //save RMS values
   printcr("$ Saving RMS flux perturbations");
   dataV.resize() = dFlux_RMS;
   dataV.perline(1);
   fname = "dflux.rms."+tail_str+".out";
   dispcr(fname);
   save(dataV,fname);
   dataV.resize() = dFluxF_RMS;
   dataV.perline(1);
   fname = "dFluxF.rms."+tail_str+".out";
   dispcr(fname);
   save(dataV,fname);


   dataV = nvec;
   dataV.perline(1);
   fname = "nvec."+tail_str+".out";
   dispcr(fname);
   save(dataV,fname);

   dataV = mvec;
   dataV.perline(1);
   fname = "mvec."+tail_str+".out";
   dispcr(fname);
   save(dataV,fname);
   
   
   return 0;
} // main()





