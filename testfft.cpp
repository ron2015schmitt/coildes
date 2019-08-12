/************************************************************************* 
 * 
 *   File Name    : 
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  
 * 
 *
 *   SYNOPSIS     
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
#include "coilfft.hpp"

// This is the file that defines the B field configuration



// Main Function for code

int main (int argc, char *argv[])
{

   const double SMALL = 5e-14;

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
   modevectors(NF,nn,mm,Nnn,Nmm,Nharm,Mharm,mode00);


   // these exclude the n=0,m=0 case
   Vector<double> nnR("nnR");
   Vector<double> mmR("mmR");
   unsigned int NFR;
   mode00 = false;
   modevectors(NFR,nnR,mmR,Nnn,Nmm,Nharm,Mharm,mode00);

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


   cout << endl;
   cout<<"$ Generate reduced orthonormal series matrix ("<<Npts<<" x "<<NFR<<")"<<endl;

   STARTTIME(tbuff,ckstart);

   Matrix<complex<double> > fsR2(Npts,NFR,"fsR");
   remove_mode00(fs,fsR2);

   STOPTIME(tbuff,ckstart);

   {
      Matrix<double> del(Npts,NFR,"del");
      del = abs(fsR - fsR2);

      double maxval = max(abs(fsR));
      dispcr(maxval);
      double maxerr = max(del);
      maxerr = (maxerr>=SMALL) ? maxerr : 0; 
      dispcr(maxerr);
   }

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

   Vector<complex<double> > BtF(NF,"BtF");
   if (load_coefs( Bt_filename,CoefFileFormat_sincos,nn,mm,BtF,false))
      return 5;

   cout <<endl<< "$ Loading BTOTAL_phi fourier coefficients from " << Bp_filename << endl;
   Vector<complex<double> > BpF(NF,"BpF");
   if (load_coefs( Bp_filename,CoefFileFormat_sincos,nn,mm,BpF,false))
      return 6;

   Vector<double> Bt(Npts, "Bt");
   expandfunction(Bt,BtF,fs);
   Vector<double> Bp(Npts, "Bp");
   expandfunction(Bp,BpF,fs);


   STOPTIME(tbuff,ckstart);


   //  Fourier Transform of complex vector

   Vector<complex<double> > ron(Npts,"ron");
   ron = vcomplex<double>(Bt,Bp);
   datavec.resize() = real(ron);
   datavec.perline(1);
   save(datavec,"ron1.R.out");
   datavec = imag(ron);
   save(datavec,"ron1.I.out");


   Vector<complex<double> > ronF1(NF,"ronF1");

   fft2d(ron,ronF1,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,2*PI);


   datavec.resize() = real(ronF1);
   datavec.perline(1);
   save(datavec,"ronF1.R.out");
   datavec = imag(ronF1);
   save(datavec,"ronF1.I.out");


   Vector<complex<double> > ronF2(NF,"ronF2");
   ronF2= (adj(fs)|ron);

   const double coeff = sqr(2*PI) / double(Npts);
   ronF2 = ronF2 *coeff;

   datavec.resize() = real(ronF2);
   datavec.perline(1);
   save(datavec,"ronF2.R.out");
   datavec = imag(ronF2);
   save(datavec,"ronF2.I.out");


   {
      Vector<double> del(NF,"del");
      del = abs(ronF2-ronF1);

      double maxval = max(abs(ronF1));
      dispcr(maxval);
      double maxerr = max(del);
      maxerr = (maxerr>=SMALL) ? maxerr : 0; 
      dispcr(maxerr);
   }




   //  Fourier Transform of real vector

   Vector<double> yat(Npts,"yat");
   yat = Bt;
   Vector<complex<double> > yatF1(NF,"yatF1");

   fft2d(yat,yatF1,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,2*PI);


   datavec.resize() = real(yatF1);
   datavec.perline(1);
   save(datavec,"yatF1.R.out");
   datavec = imag(yatF1);
   save(datavec,"yatF1.I.out");


   Vector<complex<double> > yatF2(NF,"yatF2");
   yatF2= (adj(fs)|yat);

   yatF2 = yatF2 *coeff;

   datavec.resize() = real(yatF2);
   datavec.perline(1);
   save(datavec,"yatF2.R.out");
   datavec = imag(yatF2);
   save(datavec,"yatF2.I.out");

   {
      Vector<double> del(NF,"del");
      del = abs(yatF2-yatF1);

      double maxval = max(abs(yatF1));
      dispcr(maxval);
      double maxerr = max(del);
      maxerr = (maxerr>=SMALL) ? maxerr : 0; 
      dispcr(maxerr);
   }



   ////////////////////////////////////////////
   // DEBUGGING
   ////////////////////////////////////////////


   //   Vector<complex<double> > ronF3(NF,"ronF3");

   //   testfunc(ron,ronF3,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,2*PI);


   //   datavec.resize() = real(ronF3);
   //   datavec.perline(1);
   //   save(datavec,"ronF3.R.out");
   //   datavec = imag(ronF3);
   //   save(datavec,"ronF3.I.out");

   /////////////////////////////////////////////

   //  Inverse Fourier Transform to a complex vector

   Vector<complex<double> > ron2(Npts,"ron2");
   ifft2d(ron2,ronF1,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,1.0/(2*PI));


   datavec.resize() = real(ron2);
   datavec.perline(1);
   save(datavec,"ron2.R.out");
   datavec = imag(ron2);
   save(datavec,"ron2.I.out");


   {
      Vector<double> del(Npts,"del");
      del = abs(ron-ron2);

      double maxval = max(abs(ron));
      dispcr(maxval);
      double maxerr = max(del);
      maxerr = (maxerr>=SMALL) ? maxerr : 0; 
      dispcr(maxerr);
   }



   //  Inverse Fourier Transform to a REAL vector

   Vector<double> yat2(Npts,"yat2");
   ifft2d(yat2,yatF1,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,1.0/(2*PI));

   {
      Vector<double> del(Npts,"del");
      del = abs(yat-yat2);

      double maxval = max(abs(yat));
      dispcr(maxval);
      double maxerr = max(del);
      maxerr = (maxerr>=SMALL) ? maxerr : 0; 
      dispcr(maxerr);
   }



   //  Fourier Transform of complex vector (no mode00)

   Vector<complex<double> > ronFR1(NFR,"ronFR1");

   fft2d(ron,ronFR1,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,2*PI,false);

   Vector<complex<double> > ronFR2(NFR,"ronFR2");
   ronFR2= (adj(fsR)|ron);

   ronFR2 = ronFR2 *coeff;

   {
      Vector<double> del(NFR,"del");
      del = abs(ronFR2-ronFR1);

      double maxval = max(abs(ronFR1));
      dispcr(maxval);
      double maxerr = max(del);
      maxerr = (maxerr>=SMALL) ? maxerr : 0; 
      dispcr(maxerr);
   }

   

   //  Fourier Transform of real vector (no mode00)

   Vector<complex<double> > yatFR1(NF,"yatFR1");

   fft2d(yat,yatFR1,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,2*PI,false);

   Vector<complex<double> > yatFR2(NF,"yatFR2");
   yatFR2= (adj(fsR)|yat);
   yatFR2 = yatFR2 *coeff;

   {
      Vector<double> del(NF,"del");
      del = abs(yatFR2-yatFR1);

      double maxval = max(abs(yatFR1));
      dispcr(maxval);
      double maxerr = max(del);
      maxerr = (maxerr>=SMALL) ? maxerr : 0; 
      dispcr(maxerr);
   }



   //  Inverse Fourier Transform (no mode00) to a complex vector

   Vector<complex<double> > ron3(Npts,"ron3");
   ifft2d(ron3,ronFR2,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,1.0/(2*PI),false);

   ron3 = ron3 + mean(ron);

   datavec.resize() = real(ron3);
   datavec.perline(1);
   save(datavec,"ron3.R.out");
   datavec = imag(ron3);
   save(datavec,"ron3.I.out");


   {
      Vector<double> del(Npts,"del");
      del = abs(ron-ron3);

      double maxval = max(abs(ron));
      dispcr(maxval);
      double maxerr = max(del);
      maxerr = (maxerr>=SMALL) ? maxerr : 0; 
      dispcr(maxerr);
   }


   //  Inverse Fourier Transform (no mode00) to a REAL vector

   Vector<double> yat3(Npts,"yat3");
   ifft2d(yat3,yatFR1,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,1.0/(2*PI),false);
   yat3 = yat3 + mean(yat);

   {
      Vector<double> del(Npts,"del");
      del = abs(yat-yat3);

      double maxval = max(abs(yat));
      dispcr(maxval);
      double maxerr = max(del);
      maxerr = (maxerr>=SMALL) ? maxerr : 0; 
      dispcr(maxerr);
   }



   return 0;
} // main()





