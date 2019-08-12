/************************************************************************* 
 * 
 *   File Name    :  
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  
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
#include "coilfft.hpp"
#include "coils_cmdline.hpp"
#include "surface.hpp"
#include "omegamatrix_fft.hpp"


// Main Function for code

int main (int argc, char *argv[])
{

   disable_all_options();
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



 
   Matrix<complex<double> > P(NFR,NFR,"V");
   LAvector<complex<double> >  W(NFR,"W");


   cout << endl;
   cout<<"$ Generate reduced orthonormal series matrix ("<<Npts<<" x "<<NFR<<")"<<endl;

   STARTTIME(tbuff,ckstart);

   Matrix<complex<double> > fsR(Npts,NFR,"fsR");
   fseries(nnR,mmR,thetas,phis,fsR);
      
   STOPTIME(tbuff,ckstart);


  
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

   LAvector<double> Bt(Npts, "Bt");
   ifft2d(Bt,BtF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,1/(2*PI));
   //  expandfunction(Bt,BtF,fs);
   LAvector<double> Bp(Npts, "Bp");
   ifft2d(Bp,BpF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,1/(2*PI));
   //  expandfunction(Bp,BpF,fs);

   STOPTIME(tbuff,ckstart);



   cout << endl;
   cout<<"$ Calculate OmegaN matrix ("<<NFR<<"x"<<NFR<<")"<<endl;
     
   STARTTIME(tbuff,ckstart);
     
   Matrix<complex<double> > OmegaN(NFR,NFR,"OmegaN");
   NormalizedOmegamatrix_fft_fast(OmegaN, Bt,Bp, fsR,nnR,mmR, Nphi,Ntheta, Nnn,Nmm,Nharm,Mharm);

   fsR.clear();
   STOPTIME(tbuff,ckstart);


   //SAVE OMEGAN MATRIX
   printcr("$ Saving OmegaN Matrix");
   data.resize() = real(OmegaN);
   fname = "omegaN.Nn="+Nnn_str+".Nm="+Nmm_str+".R.out";
   dispcr(fname);
   save(data,fname);
   data.resize() = imag(OmegaN);
   fname = "omegaN.Nn="+Nnn_str+".Nm="+Nmm_str+".I.out";
   dispcr(fname);
   save(data,fname);



   // EV decomposition of OmegaN Matrix
  
   cout << endl;
   cout<<"$ EV decomposition of OmegaN matrix ("<<NFR<<"x"<<NFR<<")"<<endl;
     
   STARTTIME(tbuff,ckstart);


   cooll_lapack::eigwincsort(OmegaN,P,W);
   OmegaN.clear();
     
   STOPTIME(tbuff,ckstart);


   // save  omega in eigen form
   printcr("$ Saving  OmegaN Matrix in eigen form");
   data.resize() = real(P);
   fname = "omegaN.Nn="+Nnn_str+".Nm="+Nmm_str+".PR.out";
   dispcr(fname);
   save(data,fname);
   data.resize() = imag(P);
   fname = "omegaN.Nn="+Nnn_str+".Nm="+Nmm_str+".PI.out";
   dispcr(fname);
   save(data,fname);

   datavec.resize() = real(W);
   datavec.perline(1);
   fname = "omegaN.Nn="+Nnn_str+".Nm="+Nmm_str+".WR.out";
   dispcr(fname);
   save(datavec,fname);
   datavec.resize() = imag(W);
   datavec.perline(1);
   fname = "omegaN.Nn="+Nnn_str+".Nm="+Nmm_str+".WI.out";
   dispcr(fname);
   save(datavec,fname);
   






   return 0;
} // main()





