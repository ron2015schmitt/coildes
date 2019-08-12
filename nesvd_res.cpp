/************************************************************************* 
 * 
 *   File Name    : nesvd.cpp
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  
 * 
 *
 *   SYNOPSIS     
 *     This file finds the coil current from given plasma surface flux 
 *    using simple SDV inversion of M
 *
 *
 * VERSION NOTES:
 * adds rho matrix filtering
 *
 * based on  nesvd.cpp of 2007-11-02
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
   enable_option(opt_Nharm);
   enable_option(opt_Mharm);

   enable_option(opt_NphiI);
   enable_option(opt_NthetaI);
   enable_option(opt_NmmIF);
   enable_option(opt_NnnIF);

   enable_option(opt_condnum);
   enable_option(opt_No);

  enable_option(opt_alpha_theta);  
  enable_option(opt_alpha_phi);


//   enable_option(opt_Mf);

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

   string methstr(".nesvd_res");
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


   strmtemp <<alpha_theta;
   string alpha_theta_str(strmtemp.str());
   strmtemp.str("");
   alpha_theta_str = ".at="+alpha_theta_str;
   strmtemp <<alpha_phi;
   string alpha_phi_str(strmtemp.str());
   strmtemp.str("");
   alpha_phi_str = ".ap="+alpha_phi_str;

   string alpha_str =  alpha_theta_str+ alpha_phi_str;

   dispcr(plasma_name);
   dispcr(methstr);
   dispcr(Nphi_str+Ntheta_str);
   dispcr(Nnn_str+Nmm_str);
   dispcr(NnnIF_str+NmmIF_str);

   dispcr(condnum);
   dispcr(No);
   bool using_No =false;

   if (No>0)
      using_No =true;



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
   const double C = (2*PI*2*PI/double(Npts));
   const double Csqr = C*C;


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

   fname = "nnR"+NnnIF_str+NmmIF_str+".out";
   nnIF.perline(1);
   nnIF.textformat(text_nobraces);
   save(nnIF,fname);
   fname = "mmR"+NnnIF_str+NmmIF_str+".out";
   mmIF.perline(1);
   mmIF.textformat(text_nobraces);
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
      inductancematrix(X,dA_dtdp,Xcoil,dA_dtdp_coil,M);
   
      STOPTIME(tbuff,ckstart);




      // fft  M matrix
      cout << endl;
      cout<<"$ FFT of inductance matrix ("<<NFR<<" x "<<NIF<<")"<<endl;
      STARTTIME(tbuff,ckstart);
      Matrix<complex<double> > MF(NFR,NIF,"MF"); 
      fft_of_M(M, MF, Nphi,Ntheta,Nnn,Nmm,  NphiI,NthetaI,NnnIF,NmmIF, Nharm,Mharm,1e-14);
      M.clear();
      STOPTIME(tbuff,ckstart);



   //*************************************

   //////////////////////////////////////////////////////

   double maxMFR = max(abs(real(MF)));
   dispcr(maxMFR);
   double maxMFI = max(abs(imag(MF)));
   dispcr(maxMFI);
 
   dispcr(MF.Nrows());
   dispcr(MF.Ncols());
 
   /////////////////////////////////////////////////////



   



   ///////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////

 // Create rho matrix

  cout << endl;
  cout<<"$ Generating rho matrices"<<endl;

  STARTTIME(tbuff,ckstart);

  LAvector<double> rhorootC(NIF, "rhoroot");
  rhomatrix(alpha_theta, alpha_phi, mmIF,  nnIF, rhorootC);

  LAvector<double> rhorootP(NFR, "rhorootP");
  rhomatrix(alpha_theta, alpha_phi, mmR,  nnR, rhorootP);

  STOPTIME(tbuff,ckstart);



  // multiply  flux by rhomatrix^(-1/2) 
//  cout << endl;
//  cout<<"$ multiply flux by inv(root(rho)) matrix "<<endl;

//  STARTTIME(tbuff,ckstart);

//  for (unsigned int k =0; k<NFR; k++) 
//     FluxF[k] = (1/rhorootP[k])*FluxF[k];
  
//  STOPTIME(tbuff,ckstart);


  // multiply  flux side of M matrix by rhomatrix^(-1/2)  
//  cout << endl;
 // cout<<"$ multiply plasma side of inductance matrix by inv(root(rho)) matrix "<<endl;

//  STARTTIME(tbuff,ckstart);

//  for (unsigned int r =0; r<NFR; r++) {
//     double temp = (1/rhorootP[r]);
//    for (unsigned int c =0; c<NIF; c++) {
//      MF(r,c) = temp*MF(r,c);
//    }
//  }

//  STOPTIME(tbuff,ckstart);

  // project coil side of M matrix onto inverse root rho matrix
  cout << endl;
  cout<<"$ multiply coil side of inductance matrix by inv(root(rho)) matrix "<<endl;

  STARTTIME(tbuff,ckstart);

  for (unsigned int c =0; c<NIF; c++) {
     double temp = (1/rhorootC[c]);
     for (unsigned int r =0; r<NFR; r++) {
	MF(r,c) = MF(r,c)*temp;
     }
  }
  STOPTIME(tbuff,ckstart);

   ///////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////




   STOPTIME(tbuff,ckstart);

   // calculate SVD of MF
   cout << endl;
   cout<<"$ SVD of fourier inductance matrix ("<<NFR<<" x "<<NIF<<")"<<endl;

   STARTTIME(tbuff,ckstart);

   
   Matrix<complex<double> > U(NFR,NFR,"U");
   Matrix<complex<double> > V(NIF,NIF,"V");
   const unsigned int Nmin = min(NFR,NIF);
   disp(Nmin);
   LAvector<double> S(Nmin);
   cooll_lapack::svd(MF,U,S,V);

   LAvector<double> Snorm = S/S[0];
                                                                                
   STOPTIME(tbuff,ckstart);
   


   cout<<"Largest SV="<<S[0]<<endl;
   cout<<"Smallest SV="<<S[Nmin-1]<<endl;
   double condition_number_of_matrix = S[0]/S[Nmin-1];
   dispcr(condition_number_of_matrix);



   // calculate inverse of MF
   cout << endl;
   cout<<"$ Pseudo Inverse of MF("<<NFR<<" x "<<NIF<<")"<<endl;

   STARTTIME(tbuff,ckstart);

   unsigned int Nmodes = 0;

   LAvector<double> Sinv(Nmin);

   for(unsigned int k=0; k<Nmin; k++) {
      if (Snorm[k] >= condnum) {
	 Sinv[k] = 1/S[k];
	 Nmodes++;
      } else {
	 Sinv[k] = 0;
      }
   }
   dispcr(Nmodes);
   if ( (No !=0) && (No<Nmodes)) {
      for(unsigned int k=No; k<Nmodes; k++) {
	 Sinv[k] = 0;
      }
      Nmodes = No;
      condnum = S[No-1]/S[0];
      dispcr(condnum);
   }
   
   Matrix<complex<double> > MFinv(NIF,NFR,"MFinv");

disp(U(20,13));
dispcr(U(13,20));

   // in place adjoint
   U.adjoint();

disp(U(20,13));
dispcr(U(13,20));

dispcr(V.Nrows());
dispcr(V.Ncols());
dispcr(U.Nrows());
dispcr(U.Ncols());
dispcr(MFinv.Nrows());
dispcr(MFinv.Ncols());
dispcr(S.size());


   // compute MFinv = (V|S|adj(U));
   for(unsigned int r=0; r<NIF; r++) {
     for(unsigned int c=0; c<NFR; c++) {
      complex<double> temp = 0;
      for(unsigned int k=0; k<Nmodes; k++) {
         temp  += V(r,k) * U(k,c) * Sinv[k];
      }
      MFinv(r,c) = temp;
     }
   }
  
   U.clear();
   V.clear();
   STOPTIME(tbuff,ckstart);


   strmtemp <<Nmodes;
   string Nmodes_str(strmtemp.str());
   strmtemp.str("");
   Nmodes_str = ".Nmodes="+Nmodes_str;
   dispcr(Nmodes_str);

   strmtemp <<condnum;
   string condnum_str(strmtemp.str());
   strmtemp.str("");
   condnum_str = ".cond="+condnum_str;
   dispcr(condnum_str);





   ///////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////

  // project invMF matrix onto rho matrix
  cout << endl;
  cout<<"$ multiply inv(rho) by  inverse inductance matrix ("<<NFR<<" x "<<N2<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  for (unsigned int r =0; r<NIF; r++) {
     double temp = (1/rhorootC[r]);
     for (unsigned int c =0; c<NFR; c++) {
	MFinv(r,c) = MFinv(r,c)*temp;
     }
  }

   STOPTIME(tbuff,ckstart);
   ///////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////


   // calculate current
  
   cout << endl;
   cout<<"$ Calculate coil current ("<<NIF<<" x 1)"<<endl;

   STARTTIME(tbuff,ckstart);

   LAvector<complex<double> > IF(NIF,"IF");
   IF = (MFinv|FluxF);

   STOPTIME(tbuff,ckstart);



   massage(IF,1e-10);


   if (using_No) {
      fname = "IF"+plasma_name+methstr+NnnIF_str+NmmIF_str+Nmodes_str+alpha_str+".out";
      save_coefs(fname,CoefFileFormat_sincos,nnIF,mmIF,IF);
   }else {
      fname = "IF"+plasma_name+methstr+NnnIF_str+NmmIF_str+condnum_str+alpha_str+".out";
      save_coefs(fname,CoefFileFormat_sincos,nnIF,mmIF,IF);
   }

   return 0;
} // main()





