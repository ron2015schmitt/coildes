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
 *
 * This version allows different size fourier space for current.
 *
 * based on  flux2current_fft.cpp of 2007-10-18
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

   enable_option(opt_NphiI);
   enable_option(opt_NthetaI);
   enable_option(opt_NmmIF);
   enable_option(opt_NnnIF);

   enable_option(opt_condnum);
   enable_option(opt_No);


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

   ostringstream strmtemp;

   string methstr(".nesvd");
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
   dispcr(No);
   bool using_No =false;

   if (No>0)
      using_No =true;



   // Create angle grid
   const unsigned int Npts = Ntheta*Nphi;
   Vector<double> thetas(Npts,"thetas");
   Vector<double> phis(Npts,"phis");
   anglevectors(thetas, phis, Ntheta, Nphi);

   const unsigned int NptsI = NthetaI*NphiI;
   Vector<double> thetasI(NptsI,"thetasI");
   Vector<double> phisI(NptsI,"phisI");
   anglevectors(thetasI, phisI, NthetaI, NphiI);

   // coefficient C is the integration coef for the fourier transform
   // C = dtheta*dphi
   //   = (2*pi/Ntheta)*(2*pi/Nphi)
   //   = (2*pi*2*pi/Npts)
//   const double C = (2*PI*2*PI/double(Npts));
//   const double Csqr = C*C;


   // these exclude the n=0,m=0 case
   Vector<double> nnR("nnR");
   Vector<double> mmR("mmR");
   unsigned int NFR;
   bool mode00 = false;
   modevectors(NFR,nnR,mmR,Nnn,Nmm,Nharm,Mharm,mode00);


   // these exclude the n=0,m=0 case
   Vector<double> nnIF("nnIF");
   Vector<double> mmIF("mmIF");
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



   Vector<complex<double> > FluxF(NFR,"FluxF");

   cout <<endl<< "$ Loading Plasma Flux sin/cos fourier coefficients from " << flux_filename << endl;
   if (load_coefs(flux_filename,CoefFileFormat_sincos,nnR,mmR,FluxF,false)) {
      printcr("Above ERROR occurred in "+myname+".");
      return 3;
   }




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


 
  
   // lay coil surface onto grid 
  
   Vector<p3vector<double> > Xcoil(NptsI, "Xcoil");
   Vector<p3vector<double> > dA_dtdp_coil(NptsI, "dA_dtdp_coil");

   Vector<p3vector<double> > dx_dr_coil(NptsI, "dx_dr_coil");
   Vector<p3vector<double> > dx_dtheta_coil(NptsI,"dx_dtheta_coil");
   Vector<p3vector<double> > dx_dphi_coil(NptsI,"dx_dphi_coil");
   Vector<p3vector<double> > grad_r_coil(NptsI,"grad_r_coil");
   Vector<p3vector<double> > grad_theta_coil(NptsI,"grad_theta_coil");
   Vector<p3vector<double> > grad_phi_coil(NptsI,"grad_phi_coil");

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



   
   // Compute magnitude of FGs

   cout << endl;
   cout<<"$ Compute magnitude of Fourier-Green functions"<<endl;

   STARTTIME(tbuff,ckstart);

   Vector<double> FGmag(NIF,"FGmag");
   for(unsigned int k=1; k<NIF; k++) {
       double sumsqrs = 0;
       for(unsigned int L=1; L<NFR; L++) {
	   complex<double> temp = MF(L,k);
	   sumsqrs += temp.real()*temp.real() + temp.imag()*temp.imag();
       }
       FGmag[k] = sqrt(sumsqrs);
   }
   
   STOPTIME(tbuff,ckstart);

   FGmag.perline(1);
   FGmag.textformat(text_nobraces);
   fname="fgmag"+plasma_name+NnnIF_str+NmmIF_str+".out";
   save(FGmag,fname);



   // Compute inner product of FGs with largest order mode

//    unsigned int ind = maxind(FGmag);
//    Vector <int> tesv = vcast<int>(nnIF);

//    ind = find1sttrue( (nnIF==24) && (mmIF==7) );
//    disp(nnIF[ind]);dispcr(mmIF[ind]);

//    cout << endl;
//    cout<<"$ Compute inner product of Fourier-Green functions with largest mode"<<endl;

//    STARTTIME(tbuff,ckstart);

//    Vector<double> FGip(NIF,"FGmip");
//    Vector<complex<double> > FGbig(NIF,"FGbig");
//    for(unsigned int L=1; L<NFR; L++) {
//        FGbig[L] = conj(MF(L,ind));
//    }

//    for(unsigned int k=1; k<NIF; k++) {
//        complex<double> temp(0,0);
//        for(unsigned int L=1; L<NFR; L++) {
// 	   temp += FGbig[L]*MF(L,k);
//        }
//        FGip[k] = sqrt(temp.real()*temp.real() + temp.imag()*temp.imag());
//    }
   
//    STOPTIME(tbuff,ckstart);

//    FGip.perline(1);
//    FGip.textformat(text_nobraces);
//    fname="fgip"+plasma_name+NnnIF_str+NmmIF_str+".out";
//    save(FGip,fname);





   // calculate SVD of MF
   cout << endl;
   cout<<"$ SVD of fourier inductance matrix ("<<NFR<<" x "<<NIF<<")"<<endl;

   STARTTIME(tbuff,ckstart);

   
   Matrix<complex<double> > U(NFR,NFR,"U");
   Matrix<complex<double> > V(NIF,NIF,"V");
   const unsigned int Nmin = min(NFR,NIF);
   disp(Nmin);
   Vector<double> S(Nmin);
   matricks_lapack::svd(MF,U,S,V);

   Vector<double> Snorm = S/S[0];
                                                                                
   STOPTIME(tbuff,ckstart);
   

   S.perline(1);
   S.textformat(text_nobraces);
   fname = "S"+plasma_name+methstr+NnnIF_str+NmmIF_str+".out";
   save(S,fname);

   cout<<"Largest SV="<<S[0]<<endl;
   cout<<"Smallest SV="<<S[Nmin-1]<<endl;
   double condition_number_of_matrix = S[0]/S[Nmin-1];
   dispcr(condition_number_of_matrix);



   // calculate inverse of MF
   cout << endl;
   cout<<"$ Pseudo Inverse of MF("<<NFR<<" x "<<NIF<<")"<<endl;

   STARTTIME(tbuff,ckstart);

   unsigned int Nmodes = 0;

   Vector<double> Sinv(Nmin);

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
      condnum = S[Nmodes-1]/S[0];
      dispcr(condnum);
   }
   if (condnum==0) {
      condnum = S[Nmodes-1]/S[0];
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


   // compute MFinv = (V|S|adj(U));vb
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


   data.resize() = real(MFinv);
   data.perline(1);
   data.textformat(text_nobraces);
fname="MFinv"+plasma_name+NnnIF_str+NmmIF_str+".R.out";
 save(data,fname);
   data.resize() = imag(MFinv);
   data.perline(1);
   data.textformat(text_nobraces);
fname="MFinv"+plasma_name+NnnIF_str+NmmIF_str+".I.out";
 save(data,fname);




   
   // Compute magnitude of invFGs

   cout << endl;
   cout<<"$ Compute magnitude of Fourier-Green functions"<<endl;

   STARTTIME(tbuff,ckstart);

   Vector<double> FGinvmag(NIF,"FGinvmag");
   for(unsigned int k=1; k<NIF; k++) {
       double sumsqrs = 0;
       for(unsigned int L=1; L<NFR; L++) {
	   complex<double> temp = MFinv(k,L);
	   sumsqrs += temp.real()*temp.real() + temp.imag()*temp.imag();
       }
       FGinvmag[k] = sqrt(sumsqrs);
   }
   
   STOPTIME(tbuff,ckstart);

   FGinvmag.perline(1);
   FGinvmag.textformat(text_nobraces);
   fname="fgmaginv"+plasma_name+NnnIF_str+NmmIF_str+".out";
   save(FGinvmag,fname);



   // calculate current
  
   cout << endl;
   cout<<"$ Calculate coil current ("<<NIF<<" x 1)"<<endl;

   STARTTIME(tbuff,ckstart);

   Vector<complex<double> > IF(NIF,"IF");
   IF = (MFinv|FluxF);

   STOPTIME(tbuff,ckstart);

   fname="Nmodes"+plasma_name+methstr+NnnIF_str+NmmIF_str+condnum_str+".out";
   save(Nmodes,fname);

   fname="condnum"+plasma_name+methstr+NnnIF_str+NmmIF_str+Nmodes_str+".out";
   save(condnum,fname);


   massage(IF,1e-10);

   if (using_No) {
      fname = "IF"+plasma_name+methstr+NnnIF_str+NmmIF_str+Nmodes_str+".out";
      save_coefs(fname,CoefFileFormat_sincos,nnIF,mmIF,IF);
   }else {
      fname = "IF"+plasma_name+methstr+NnnIF_str+NmmIF_str+condnum_str+".out";
      save_coefs(fname,CoefFileFormat_sincos,nnIF,mmIF,IF);
   }

   return 0;
} // main()





