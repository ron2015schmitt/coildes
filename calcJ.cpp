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

  disable_all_options();
  enable_option(opt_pf);
  enable_option(opt_Nphi);
  enable_option(opt_Ntheta);
  enable_option(opt_Nnn);
  enable_option(opt_Nmm);
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




  // coefficient C is the integration coef for the fourier transform
  // C = dtheta*dphi
  //   = (2*pi/Ntheta)*(2*pi/Nphi)
  //   = (2*pi*2*pi/Npts)
  const double C = (2*PI*2*PI/double(Npts));
  const double Csqr = C*C;

  // Create ortho normal series

//   cout << endl;
//   cout<<"$ Generate orthonormal series matrix ("<<Npts<<" x "<<NF<<")"<<endl;

//   STARTTIME(tbuff,ckstart);

//   Matrix<complex<double> > fs(Npts,NF,"fs");
//   fseries(nn,mm,thetas,phis,fs);

//   STOPTIME(tbuff,ckstart);


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


  LAvector<double> Xx(Npts, "Xx");
  LAvector<double> Xy(Npts, "Xy");
  LAvector<double> Xz(Npts, "Xz");
  LAvector<double> Xr(Npts, "Xr");

  LAvector<double> dx_dtheta_x(Npts, "dx_dtheta_x");
  LAvector<double> dx_dtheta_y(Npts, "dx_dtheta_y");
  LAvector<double> dx_dtheta_z(Npts, "dx_dtheta_z");
  LAvector<double> dx_dtheta_r(Npts, "dx_dtheta_r");
  LAvector<double> dx_dtheta_p(Npts, "dx_dtheta_p");


  LAvector<double> dx_dphi_x(Npts, "dx_dphi_x");
  LAvector<double> dx_dphi_y(Npts, "dx_dphi_y");
  LAvector<double> dx_dphi_z(Npts, "dx_dphi_z");
  LAvector<double> dx_dphi_r(Npts, "dx_dphi_r");
  LAvector<double> dx_dphi_p(Npts, "dx_dphi_p");


  LAvector<double> J(Npts, "J");
  LAvector<double> Jsqr(Npts, "Jsqr");
  for (unsigned int j =0; j<Npts; j++) {
     const double phi=phis[j];
     Xx[j] = X[j].x();
     Xy[j] = X[j].y();
     Xz[j] = X[j].z();
     Xr[j] = sqrt( sqr(Xx[j]) + sqr(Xy[j]) );

     dx_dtheta_x[j] = dx_dtheta[j].x();
     dx_dtheta_y[j] = dx_dtheta[j].y();
     dx_dtheta_z[j] = dx_dtheta[j].z();
     dx_dtheta_r[j] = dx_dtheta_x[j]*cos(phi) + dx_dtheta_y[j]*sin(phi);
     dx_dtheta_p[j] = -dx_dtheta_x[j]*sin(phi) + dx_dtheta_y[j]*cos(phi);

     dx_dphi_x[j] = dx_dphi[j].x();
     dx_dphi_y[j] = dx_dphi[j].y();
     dx_dphi_z[j] = dx_dphi[j].z();
     dx_dphi_r[j] = dx_dphi_x[j]*cos(phi) + dx_dphi_y[j]*sin(phi);
     dx_dphi_p[j] = -dx_dphi_x[j]*sin(phi) + dx_dphi_y[j]*cos(phi);

     J[j] = dot(dx_dr[j],cross(dx_dtheta[j], dx_dphi[j]));
     Jsqr[j] = J[j]*J[j];
  }


  STOPTIME(tbuff,ckstart);

  dispcr(stddev(Xx));
  dispcr(stddev(Xy));
  dispcr(stddev(Xr));
  dispcr(stddev(Xz));
  dispcr(stddev(J));
  dispcr(stddev(Jsqr));



   LAvector<complex<double> > XxF(NF,"XxF");
   LAvector<complex<double> > XyF(NF,"XyF");
   LAvector<complex<double> > XzF(NF,"XzF");
   LAvector<complex<double> > XrF(NF,"XrF");
   LAvector<complex<double> > XpF(NF,"XpF");


   LAvector<complex<double> > dx_dtheta_xF(NF,"dx_dtheta_xF");
   LAvector<complex<double> > dx_dtheta_yF(NF,"dx_dtheta_yF");
   LAvector<complex<double> > dx_dtheta_zF(NF,"dx_dtheta_zF");
   LAvector<complex<double> > dx_dtheta_rF(NF,"dx_dtheta_rF");
   LAvector<complex<double> > dx_dtheta_pF(NF,"dx_dtheta_pF");

   LAvector<complex<double> > dx_dphi_xF(NF,"dx_dphi_xF");
   LAvector<complex<double> > dx_dphi_yF(NF,"dx_dphi_yF");
   LAvector<complex<double> > dx_dphi_zF(NF,"dx_dphi_zF");
   LAvector<complex<double> > dx_dphi_rF(NF,"dx_dphi_rF");
   LAvector<complex<double> > dx_dphi_pF(NF,"dx_dphi_pF");


   LAvector<complex<double> > JF(NF,"JF");
   LAvector<complex<double> > JsqrF(NF,"JsqrF");

   fft2d(Xx,XxF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,2*PI);
   fft2d(Xy,XyF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,2*PI);
   fft2d(Xz,XzF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,2*PI);
   fft2d(Xr,XrF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,2*PI);

   fft2d( dx_dtheta_x, dx_dtheta_xF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,2*PI);
   fft2d( dx_dtheta_y, dx_dtheta_yF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,2*PI);
   fft2d( dx_dtheta_z, dx_dtheta_zF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,2*PI);
   fft2d( dx_dtheta_r, dx_dtheta_rF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,2*PI);
   fft2d( dx_dtheta_p, dx_dtheta_pF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,2*PI);

   fft2d( dx_dphi_x, dx_dphi_xF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,2*PI);
   fft2d( dx_dphi_y, dx_dphi_yF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,2*PI);
   fft2d( dx_dphi_z, dx_dphi_zF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,2*PI);
   fft2d( dx_dphi_r, dx_dphi_rF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,2*PI);
   fft2d( dx_dphi_p, dx_dphi_pF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,2*PI);

   fft2d(J,JF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,2*PI);
   fft2d(Jsqr,JsqrF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,2*PI);

  
   dispcr(norm(XxF)/(2*PI));
   dispcr(norm(XyF)/(2*PI));
   dispcr(norm(XrF)/(2*PI));
   dispcr(norm(XzF)/(2*PI));

   dispcr(norm( dx_dtheta_xF)/(2*PI));
   dispcr(norm( dx_dtheta_yF)/(2*PI));
   dispcr(norm( dx_dtheta_rF)/(2*PI));
   dispcr(norm( dx_dtheta_pF)/(2*PI));
   dispcr(norm( dx_dtheta_zF)/(2*PI));

   dispcr(norm( dx_dphi_xF)/(2*PI));
   dispcr(norm( dx_dphi_yF)/(2*PI));
   dispcr(norm( dx_dphi_rF)/(2*PI));
   dispcr(norm( dx_dphi_pF)/(2*PI));
   dispcr(norm( dx_dphi_zF)/(2*PI));

   dispcr(norm(JF)/(2*PI));
   dispcr(norm(JsqrF)/(2*PI));

   save_coefs("XxF.out",CoefFileFormat_sincos,nn,mm,XxF);
   save_coefs("XyF.out",CoefFileFormat_sincos,nn,mm,XyF);
   save_coefs("XzF.out",CoefFileFormat_sincos,nn,mm,XzF);
   save_coefs("XrF.out",CoefFileFormat_sincos,nn,mm,XrF);

   save_coefs("dx_dtheta_xF.out",CoefFileFormat_sincos,nn,mm, dx_dtheta_xF);
   save_coefs("dx_dtheta_yF.out",CoefFileFormat_sincos,nn,mm, dx_dtheta_yF);
   save_coefs("dx_dtheta_zF.out",CoefFileFormat_sincos,nn,mm, dx_dtheta_zF);
   save_coefs("dx_dtheta_rF.out",CoefFileFormat_sincos,nn,mm, dx_dtheta_rF);
   save_coefs("dx_dtheta_pF.out",CoefFileFormat_sincos,nn,mm, dx_dtheta_pF);

   save_coefs("dx_dphi_xF.out",CoefFileFormat_sincos,nn,mm, dx_dphi_xF);
   save_coefs("dx_dphi_yF.out",CoefFileFormat_sincos,nn,mm, dx_dphi_yF);
   save_coefs("dx_dphi_zF.out",CoefFileFormat_sincos,nn,mm, dx_dphi_zF);
   save_coefs("dx_dphi_rF.out",CoefFileFormat_sincos,nn,mm, dx_dphi_rF);
   save_coefs("dx_dphi_pF.out",CoefFileFormat_sincos,nn,mm, dx_dphi_pF);

   save_coefs("JF.out",CoefFileFormat_sincos,nn,mm,JF);
   save_coefs("JsqrF.out",CoefFileFormat_sincos,nn,mm,JsqrF);



  return 0;
} // main()





