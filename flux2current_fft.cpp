/************************************************************************* 
 * 
 *   File Name    :  flux2current_fft.cpp
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
 * -No No 
 *     limit dimension of SVD pseudo inverse to No 
 * -cond condnum
 *     treat any singular values that are below condnum*(largest SV) as zero
 *     when performing the pseudo inverse.
 *  
 * When No and condnum are both specified, which ever gives the smaller number of
 * modes is used.
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

   enable_option(opt_condnum);
   enable_option(opt_No);
   enable_option(opt_Mf);

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

   string methstr(".f2c");
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


   dispcr(plasma_name);
   dispcr(methstr);
   dispcr(Nphi_str+Ntheta_str);
   dispcr(Nnn_str+Nmm_str);

   dispcr(condnum);
   dispcr(No);

   bool using_No =false;



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


   // these exclude the n=0,m=0 case
   LAvector<double> nnR("nnR");
   LAvector<double> mmR("mmR");
   unsigned int NFR;
   bool mode00 = false;
   modevectors(NFR,nnR,mmR,Nnn,Nmm,Nharm,Mharm,mode00);


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
  
   LAvector<p3vector<double> > Xcoil(Npts, "Xcoil");
   LAvector<p3vector<double> > dA_dtdp_coil(Npts, "dA_dtdp_coil");

   LAvector<p3vector<double> > dx_dr_coil(Npts, "dx_dr_coil");
   LAvector<p3vector<double> > dx_dtheta_coil(Npts,"dx_dtheta_coil");
   LAvector<p3vector<double> > dx_dphi_coil(Npts,"dx_dphi_coil");
   LAvector<p3vector<double> > grad_r_coil(Npts,"grad_r_coil");
   LAvector<p3vector<double> > grad_theta_coil(Npts,"grad_theta_coil");
   LAvector<p3vector<double> > grad_phi_coil(Npts,"grad_phi_coil");

   cout << endl;
   cout <<"$ Mapping coil surface fourier coefficients to "<<Ntheta<<" x "<<Nphi<<" (theta by phi) grid"<<endl;

   STARTTIME(tbuff,ckstart);

   expandsurfaceandbases(Xcoil,dA_dtdp_coil,
			 dx_dr_coil,dx_dtheta_coil,dx_dphi_coil,
			 grad_r_coil,grad_theta_coil,grad_phi_coil,
			 coilfourier,thetas,phis);

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


   Matrix<complex<double> > MF(NFR,NFR,"MF"); 
   if ( M_filename.empty() ) {

      // Create mutual inductance matrix

      cout << endl;
      cout<<"$ Generating inductance matrix ("<<Npts<<" x "<<Npts<<")"<<endl;

      STARTTIME(tbuff,ckstart);

      Matrix<double> M(Npts, Npts, "M");
      inductancematrix(X,dA_dtdp,Xcoil,dA_dtdp_coil,M);
   
      STOPTIME(tbuff,ckstart);

      // fft  M matrix
      cout << endl;
      cout<<"$ FFT of inductance matrix ("<<Npts<<" x "<<NFR<<")"<<endl;
      
      STARTTIME(tbuff,ckstart);
  
      fft_of_M(M, MF, Nphi, Ntheta,Nnn,Nmm,  Nharm,Mharm,1e-12);
      M.clear();
      STOPTIME(tbuff,ckstart);

      fname = "MF"+plasma_name+Nnn_str+Nmm_str+".R.out";
      data.resize() = real(MF);
      data.perline(1);
      data.textformat(text_nobraces);
      dispcr(fname);
      save(data,fname);
      fname = "MF"+plasma_name+Nnn_str+Nmm_str+".I.out";
      data = imag(MF);
      dispcr(fname);
      save(data,fname);
      data.clear();

   } else {
      cout << endl;
      cout<<"$ Load M (inductance) matrix("<<NFR<<"x"<<NFR<<")"<<endl;
      
      STARTTIME(tbuff,ckstart);
      
      data.resize(NFR,NFR);
      data.perline(1);
      data.textformat(text_nobraces);
      fname = M_filename +Nnn_str+Nmm_str+ ".R.out";
      dispcr(fname);
      load(data,fname);
      Matrix <double> data2("data2");
      data2.resize(NFR,NFR);
      data2.perline(1);
      data2.textformat(text_nobraces);
      fname =  M_filename +Nnn_str+Nmm_str+ ".I.out";
      dispcr(fname);
      load(data2,fname);
      MF = mcomplex(data,data2);
      data.clear();
      data2.clear();

      STOPTIME(tbuff,ckstart);

   }

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
   LAvector<double> S(NFR);
   cooll_lapack::svd(MF,U,S,V);

   LAvector<double> Snorm = S/S[0];
                                                                                
   STOPTIME(tbuff,ckstart);
   

   S.perline(1);
   S.textformat(text_nobraces);
   fname = "S"+plasma_name+methstr+Nnn_str+Nmm_str+".out";
   save(S,fname);

   cout<<"Largest SV="<<S[0]<<endl;
   cout<<"Smallest SV="<<S[NFR-1]<<endl;
   double condition_number_of_matrix = S[0]/S[NFR-1];
   dispcr(condition_number_of_matrix);



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


   unsigned int Nmodes = 0;
   Matrix<complex<double> > SU(NFR,NFR,"SU");
   for(unsigned int r=0; r<NFR; r++) {
      double temp;
      if (Snorm[r] >= condnum) {
	 temp = 1/S[r];
	 Nmodes++;
      } else {
	 temp = 0;
      }
      for(unsigned int c=0; c<NFR; c++) {
	 SU(r,c) = conj(U(c,r)) * temp;
      }
   }
   dispcr(Nmodes);
   if ( (No !=0) && (No<Nmodes)) {
      for(unsigned int r=No; r<Nmodes; r++) {
	 for(unsigned int c=0; c<NFR; c++) {
	    SU(r,c) = 0;
	 }
      }
      Nmodes = No;
      condnum = S[No-1]/S[0];
      dispcr(condnum);
      using_No =true;
   }

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


   
   U.resize(0,0);
   S.resize(0);
   Matrix<complex<double> > MFinv(NFR,NFR,"MFinv");
   MFinv = (V|SU);
  
   V.clear();
   STOPTIME(tbuff,ckstart);



   // calculate current
  
   cout << endl;
   cout<<"$ Calculate coil current ("<<NFR<<" x 1)"<<endl;

   STARTTIME(tbuff,ckstart);

   LAvector<complex<double> > IF(NFR,"IF");
   IF = (MFinv|FluxF);

   STOPTIME(tbuff,ckstart);

   fname="Nmodes"+plasma_name+methstr+Nnn_str+Nmm_str+condnum_str+".out";
   save(Nmodes,fname);

   fname="condnum"+plasma_name+methstr+Nnn_str+Nmm_str+Nmodes_str+".out";
   save(condnum,fname);


   massage(IF,1e-10);

   if (using_No) {
      fname = "IF"+plasma_name+methstr+Nnn_str+Nmm_str+Nmodes_str+".out";
      save_coefs(fname,CoefFileFormat_sincos,nnR,mmR,IF);
   }else {
      fname = "IF"+plasma_name+methstr+Nnn_str+Nmm_str+condnum_str+".out";
      save_coefs(fname,CoefFileFormat_sincos,nnR,mmR,IF);
   }

   return 0;
} // main()





