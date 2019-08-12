/************************************************************************* 
 * 
 *   File Name    :  nescoil.cpp
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
#include "inductancematrix.hpp"
#include "omegamatrix_fft.hpp"
#include "coilfft.hpp"

// This is the file that defines the B field configuration



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

   enable_option(opt_NphiI);
   enable_option(opt_NthetaI);
   enable_option(opt_NmmIF);
   enable_option(opt_NnnIF);

   // set this flag to compute and store extensive data
   enable_option(opt_extradata);


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

   ostringstream strmtemp;

   string methstr(".nescoil");
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


   strmtemp <<NnnIF;
   string NnnIF_str(strmtemp.str());
   strmtemp.str("");
   NnnIF_str = ".NnIF="+NnnIF_str;

   strmtemp <<NmmIF;
   string NmmIF_str(strmtemp.str());
   strmtemp.str("");
   NmmIF_str = ".NmIF="+NmmIF_str;

   dispcr(plasma_name);
   dispcr(methstr);
   dispcr(Nphi_str+Ntheta_str);
   dispcr(Nnn_str+Nmm_str);
   dispcr(NnnIF_str+NmmIF_str);

   dispcr(condnum);



   // Create angle grid
   const unsigned int Npts = Ntheta*Nphi;
   LAvector<double> thetas(Npts,"thetas");
   LAvector<double> phis(Npts,"phis");
   anglevectors(thetas, phis, Ntheta, Nphi);

   const unsigned int NptsI = NthetaI*NphiI;
   LAvector<double> thetasI(NptsI,"thetasI");
   LAvector<double> phisI(NptsI,"phisI");
   anglevectors(thetasI, phisI, NthetaI, NphiI);

   // coefficient C is the integration coef for the fourier transform
   // C = dtheta*dphi
   //   = (2*pi/Ntheta)*(2*pi/Nphi)
   //   = (2*pi*2*pi/Npts)
   //   const double C = (2*PI*2*PI/double(Npts));
   //   const double Csqr = C*C;


   // these exclude the n=0,m=0 case
   LAvector<double> nnR("nnR");
   LAvector<double> mmR("mmR");
   unsigned int NFR;
   bool mode00 = false;
   modevectors(NFR,nnR,mmR,Nnn,Nmm,Nharm,Mharm,mode00);


   // these exclude the n=0,m=0 case
   LAvector<double> nnIF("nnIF");
   LAvector<double> mmIF("mmIF");
   unsigned int NIF;
   mode00 = false;
   modevectors(NIF,nnIF,mmIF,NnnIF,NmmIF,Nharm,Mharm,mode00);
   dispcr(NIF);


   // save mode vectors

   fname = "nnR"+Nnn_str+Nmm_str+".out";
   nnR.perline(1);
   nnR.textformat(text_nobraces);
   dispcr(fname);
   save(nnR,fname);
   fname = "mmR"+Nnn_str+Nmm_str+".out";
   mmR.perline(1);
   mmR.textformat(text_nobraces);
   dispcr(fname);
   save(mmR,fname);

   strmtemp.str("");
   strmtemp <<NnnIF;
   fname = "nnR.Nn="+strmtemp.str();
   strmtemp.str("");
   strmtemp <<NmmIF;
   fname = fname +".Nm="+strmtemp.str()+".out";
   nnIF.perline(1);
   nnIF.textformat(text_nobraces);
   dispcr(fname);
   save(nnIF,fname);
   strmtemp.str("");
   strmtemp <<NnnIF;
   fname = "mmR.Nn="+strmtemp.str();
   strmtemp.str("");
   strmtemp <<NmmIF;
   fname = fname +".Nm="+strmtemp.str()+".out";
   mmIF.perline(1);
   mmIF.textformat(text_nobraces);
   dispcr(fname);
   save(mmIF,fname);



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



   LAvector<complex<double> > FluxF(NFR,"FluxF");

   cout <<endl<< "$ Loading Plasma Flux sin/cos fourier coefficients from " << flux_filename << endl;
   if (load_coefs(flux_filename,CoefFileFormat_sincos,nnR,mmR,FluxF,false)) {
      printcr("Above ERROR occurred in "+myname+".");
      return 3;
   }

   LAvector<double> Flux(Npts, "Flux");
   mode00=false;
   ifft2d(Flux,FluxF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,1/(2*PI),mode00);

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

   LAvector<p3vector<double> > rootn(Npts, "rootn");
   LAvector<double> J(Npts, "J");
   for (unsigned int j =0; j<Npts; j++) {
      J[j] = dot(dx_dr[j],cross(dx_dtheta[j], dx_dphi[j]));
      double temp = 1/sqrt(norm(dA_dtdp[j]));
      rootn[j] = dA_dtdp[j]*temp;
      Flux[j] = Flux[j]*temp;
   }


   STOPTIME(tbuff,ckstart);


 
  
   // lay coil surface onto grid 
  
   LAvector<p3vector<double> > Xcoil(NptsI, "Xcoil");
   LAvector<p3vector<double> > dA_dtdp_coil(NptsI, "dA_dtdp_coil");

   LAvector<p3vector<double> > dx_dr_coil(NptsI, "dx_dr_coil");
   LAvector<p3vector<double> > dx_dtheta_coil(NptsI,"dx_dtheta_coil");
   LAvector<p3vector<double> > dx_dphi_coil(NptsI,"dx_dphi_coil");
   LAvector<p3vector<double> > grad_r_coil(NptsI,"grad_r_coil");
   LAvector<p3vector<double> > grad_theta_coil(NptsI,"grad_theta_coil");
   LAvector<p3vector<double> > grad_phi_coil(NptsI,"grad_phi_coil");

   cout << endl;
   cout <<"$ Mapping coil surface fourier coefficients to "<<NthetaI<<" x "<<NphiI<<" (theta by phi) grid"<<endl;

   STARTTIME(tbuff,ckstart);

   expandsurfaceandbases(Xcoil,dA_dtdp_coil,
			 dx_dr_coil,dx_dtheta_coil,dx_dphi_coil,
			 grad_r_coil,grad_theta_coil,grad_phi_coil,
			 coilfourier,thetasI,phisI);


 



   // Create mutual inductance matrix

   cout << endl;
   cout<<"$ Generating inductance matrix ("<<Npts<<" x "<<NptsI<<")"<<endl;

   STARTTIME(tbuff,ckstart);

   Matrix<double> M(Npts, NptsI, "M");
   inductancematrix(X,rootn,Xcoil,dA_dtdp_coil,M);
   
   STOPTIME(tbuff,ckstart);


  
   // fft of coil side of M matrix

   Matrix<complex<double> > FG(Npts,NIF,"FG"); 

   cout << endl;
   cout<<"$ FFT coil side of inductance matrix -- Fourier-Green Functions ("<<Npts<<" x "<<NIF<<")"<<endl;
   
   STARTTIME(tbuff,ckstart);
      
   half_fft_of_M(M, FG, NphiI, NthetaI, NnnIF, NmmIF, Nharm,Mharm,1e-14);

   STOPTIME(tbuff,ckstart);

   M.clear();

   if (extradata_flag){
      fname="FG"+plasma_name+Nphi_str+Ntheta_str+NnnIF_str+NmmIF_str+".R.out";
      data.resize() = real(FG);
      data.perline(1);
      data.textformat(text_nobraces);
      dispcr(fname);
      save(data,fname);
      fname="FG"+plasma_name+Nphi_str+Ntheta_str+NnnIF_str+NmmIF_str+".I.out";
      data = imag(FG);
      dispcr(fname);
      save(data,fname);
   }

   cout << endl;
   cout<<"$ Project inductance matrix onto Fourier-Green functions("<<NIF<<" x "<<NIF<<")"<<endl;

   STARTTIME(tbuff,ckstart);

   Matrix<complex<double> > MF(NIF,NIF,"MF"); 
   MF = (adj(FG)|FG);
  
   STOPTIME(tbuff,ckstart);

   if (extradata_flag){

      fname="FGadjFG"+plasma_name+Nphi_str+Ntheta_str+NnnIF_str+NmmIF_str+".R.out";
      data.resize() = real(MF);
      data.perline(1);
      data.textformat(text_nobraces);
      dispcr(fname);
      save(data,fname);
      fname="FGadjFG"+plasma_name+Nphi_str+Ntheta_str+NnnIF_str+NmmIF_str+".I.out";
      data = imag(MF);
      dispcr(fname);
      save(data,fname);
   }

   dispcr(MF.Nrows());
   dispcr(MF.Ncols());



   cout << endl;
   cout<<"$ Project Flux onto Fourier-Green functions"<<endl;
     
   STARTTIME(tbuff,ckstart);

   LAvector<complex<double> > FluxFG(NIF,"FluxFG");
   FluxFG = (adj(FG)|Flux);

   STOPTIME(tbuff,ckstart);

   massage(FluxFG,1e-8);
   fname="FluxFG"+plasma_name+NnnIF_str+NmmIF_str+".out";
   dispcr(fname);
   save_coefs(fname,CoefFileFormat_sincos,nnIF,mmIF,FluxFG);

   FG.clear();





   // calculate inverse of MF
   cout << endl;
   cout<<"$ invert MF ("<<NIF<<" x "<<NIF<<")"<<endl;

   STARTTIME(tbuff,ckstart);

   Matrix<complex<double> > MFinv(NIF,NIF,"MFinv");
   cooll_lapack::inv(MF,MFinv);
   MF.clear();

   STOPTIME(tbuff,ckstart);



   // calculate current
  
   cout << endl;
   cout<<"$ Calculate coil current ("<<NIF<<" x 1)"<<endl;

   STARTTIME(tbuff,ckstart);
   LAvector<complex<double> > IF(NIF,"IF");
   IF = (MFinv|FluxFG);

   STOPTIME(tbuff,ckstart);

   MFinv.clear();

   // ALL DONE, NOW SAVE TO FILES


   massage(IF,1e-8);
   fname = "IF"+plasma_name+methstr+NnnIF_str+NmmIF_str+".out";
   dispcr(fname);
   save_coefs(fname,CoefFileFormat_sincos,nnIF,mmIF,IF);

   return 0;
} // main()





