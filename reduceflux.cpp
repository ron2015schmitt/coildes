/************************************************************************* 
 * 
 *   File Name    :  
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  16 May 2003
 * 
 *
 *   SYNOPSIS     
 *     Stellarator COIL Design 
 *     This file finds the coil current from given plasma surface flux 
 *    
 *
 * VERSION NOTES:
 * latest version as of Jan 17, 2005
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
  enable_option(opt_ff2);
  enable_option(opt_Nphi);
  enable_option(opt_Ntheta);
  enable_option(opt_Nnn);
  enable_option(opt_Nmm);
  enable_option(opt_No);
  enable_option(opt_Btf);
  enable_option(opt_Bpf);

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


  // Create Fourier Mode vectors
  Vector<double> nn("nn");
  Vector<double> mm("mm");
  unsigned int NF;
  modevectors(NF,nn,mm,Nnn,Nmm);

  // these exclude the n=0,m=0 case
  Vector<double> nnR("nnR");
  Vector<double> mmR("mmR");
  unsigned int NFR;
  modevectorsR(NFR,nnR,mmR,Nnn,Nmm);


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
  



  Vector<complex<double> > FluxF(NFR,"FluxF");

  cout <<endl<< "$ Loading Plasma Flux sin/cos fourier coefficients from " << flux_filename << endl;
  if (load_coefs(flux_filename,CoefFileFormat_sincos,nnR,mmR,FluxF))
    return 3;

  Vector<complex<double> > FluxF2(NFR,"FluxF2");

  cout <<endl<< "$ Loading second ste of Flux sin/cos fourier coefficients from " << flux2_filename << endl;
  if (load_coefs(flux2_filename,CoefFileFormat_sincos,nnR,mmR,FluxF2))
    return 4;

 
 
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

  cout << endl;
  cout<<"$ Generate Gamma matrix "<<endl;

  STARTTIME(tbuff,ckstart);

  Matrix<complex<double> > Gamma("Gamma");
  gammamatrix(Gamma,fsR,plasmafourier,dA_dtdp,thetas,phis,Ntheta,Nphi);

  const unsigned int Ns = Gamma.Ncols();
  dispcr(Gamma.Nrows());
  dispcr(Ns);
  STOPTIME(tbuff,ckstart);

  if (NFR<=Ns) {
    printcr("*****ERROR: NFR is too small.");
    return(4);
  }

  const unsigned int Nmin1 = min(NFR,Ns);
  if (No==0)
    No = Nmin1;

  ostringstream strm;
  strm <<No;
  string No_str(strm.str());
  dispcr(No);

  if (No>Nmin1) {
    printcr("*****ERROR: No is too large.");
    return(4);
  }

  // SVD of Gamma Matrix
  cout << endl;
  cout<<"$ SVD of Gamma matrix ("<<NFR<<"x"<<Ns<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  Matrix<complex<double> > U(NFR,NFR,"U");
  Matrix<complex<double> > V(Ns,Ns,"V");
  Vector<double> S(Nmin1,"S");
  matricks_lapack::svd(Gamma,U,S,V);

  STOPTIME(tbuff,ckstart);

  //************************************************
//   data.resize() = real(U);
//   fname = "gamma.scoild.No="+No_str+".UR.out";
//   save(data,fname);
//   data.resize() = imag(U);
//   fname = "gamma.scoild.No="+No_str+".UI.out";
//   save(data,fname);

//   fname = "gamma.scoild.No="+No_str+".S.out";
//   save(S,fname);

//   data.resize() = real(V);
//   fname = "gamma.scoild.No="+No_str+".VR.out";
//   save(data,fname);
//   data.resize() = imag(V);
//   fname = "gamma.scoild.No="+No_str+".VI.out";
//   save(data,fname);
  //************************************************



  S.perline(1);
  S.textformat(text_nobraces);

  printcr("gamma SVs: ");
  dispcr(S);
                                                                                
  // REDUCE SVD to No significant sv's

  cout << endl;
  print("$ Reduce to No=");
  print(No);
  printcr(" significant vlaues and find new basis ");

  STARTTIME(tbuff,ckstart);

  Matrix<complex<double> > Ur(NFR,No,"Ur");
  for (unsigned int c=0; c<No; c++)
    Ur.col(c)=U.col(c);  

  Matrix<complex<double> > f_delta(Npts,No,"f_delta");
  f_delta = (fsR|Ur);
                                                                                
  STOPTIME(tbuff,ckstart);


  ////////// debuging
//   Matrix<complex<double> > Ao(No,No,"C");
//   Ao = (adj(f_delta)|f_delta);
//   for (unsigned int c=0; c<No; c++)
//     disp(Ao(c,c));
//   double maxoffval = 0.0;
//   for (unsigned int r=0; r<No; r++) {
//     for (unsigned int c=0; c<No; c++) {
//       if (r!=c) {
// 	maxoffval = max(maxoffval,abs(Ao(r,c)));
//       }
//     }
//   }
//   dispcr(maxoffval);
  //////////

  // load or calculate B field
  cout << endl;
  cout<<"$ Load tangent B field "<<endl;

  STARTTIME(tbuff,ckstart);

  Vector<p3vector<double> > B(Npts, "B");

  //  for (unsigned int j =0; j<Npts; j++)
  //    bTotal(X[j], B[j]);

//   Vector<complex<double> > BrF(NF,"BrF");
//   if (load_coefs("BrTotal.sincos.out",CoefFileFormat_sincos,nn,mm,BrF))
//     return 4;
//   Vector<complex<double> > Br(Npts, "Br");
//   expandfunction(Br,BrF,fs);


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

  for (unsigned int j =0; j<Npts; j++)
    B[j] =  Bt[j] * dx_dtheta[j] + Bp[j] * dx_dphi[j];

  STOPTIME(tbuff,ckstart);


  // Calculate Omega matrix

  cout << endl;
  cout<<"$ Calculate Omega matrix ("<<NFR<<"x"<<No<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  Matrix<p3vector<complex<double> > >  grad_fsR(Npts,NFR,"grad_fsR");
  gradfvector(thetas,phis,mmR,nnR,grad_theta,grad_phi,grad_fsR);
  Matrix<complex<double> > Omega(NFR,No,"Omega");
  omegamatrix(Omega,B,grad_fsR,f_delta,dx_dr,dx_dtheta,dx_dphi,grad_r,grad_theta,grad_phi);

  STOPTIME(tbuff,ckstart);

 
  // SVD of Omega Matrix

  cout << endl;
  cout<<"$ SVD of Omega matrix ("<<NFR<<"x"<<No<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  U.resize(NFR,NFR);
  V.resize(No,No);
  S.resize(min(NFR,No));
  matricks_lapack::svd(Omega,U,S,V);

  STOPTIME(tbuff,ckstart);

  //************************************************
//   data.resize() = real(U);
//   fname = "omega.scoild.No="+No_str+".UR.out";
//   save(data,fname);
//   data.resize() = imag(U);
//   fname = "omega.scoild.No="+No_str+".UI.out";
//   save(data,fname);

//   fname = "omega.scoild.No="+No_str+".S.out";
//   save(S,fname);

//   data.resize() = real(V);
//   fname = "omega.scoild.No="+No_str+".VR.out";
//   save(data,fname);
//   data.resize() = imag(V);
//   fname = "omega.scoild.No="+No_str+".VI.out";
//   save(data,fname);
  //************************************************


                                                                                


 // Calculate reduced basis of No functions in "Gamma space"


  cout << endl;
  cout<<"$ Calculate optimal basis for flux "<<endl;

  STARTTIME(tbuff,ckstart);

  Ur.resize(NFR,No);
  for (unsigned int c=0; c<No; c++)
    Ur.col(c)=U.col(No-1-c);  
  //    Ur.col(c)=U.col(c);  

  Matrix<complex<double> > f_phix(Npts,No,"f_phix");
  f_phix = (fsR|Ur);

  STOPTIME(tbuff,ckstart);

  // project flux into optimal space
  cout << endl;
  cout<<"$ Project flux on basis "<<endl;

  STARTTIME(tbuff,ckstart);

  Vector<complex<double> > FluxFP(No,"FluxFP");
  FluxFP = (adj(Ur)|FluxF);

  Vector<complex<double> > FluxFP2(No,"FluxFP2");
  FluxFP2 = (adj(Ur)|FluxF2);

  STOPTIME(tbuff,ckstart);

 
//   massage(Ur,1e-10);
//   Ur.perline(No);
//   dispcr(Ur);
//   cr();
//   cr();

//   Matrix<complex<double> > UU(NFR,NFR,"Ur|adj(Ur)");
//   UU= (Ur|adj(Ur));

//   massage(UU,1e-10);
//   FluxFP.perline(NFR);
//   dispcr(UU);

  cr();
  cr();
  massage(FluxFP,1e-10);
  FluxFP.perline(1);
  dispcr(FluxFP);

  cr();
  cr();
  massage(FluxFP2,1e-10);
  FluxFP2.perline(1);
  dispcr(FluxFP2);


  Vector<complex<double> > FluxFPdelta(No,"FluxFPdelta");
  FluxFPdelta = FluxFP2 - FluxFP;
 
  cr();
  cr();
  FluxFPdelta.perline(1);
  dispcr(FluxFPdelta);
   

   datavec.resize() = real(FluxFP);
   fname = "FluxF.projected.RE.No="+No_str+".out";
   save(datavec,fname);
   datavec.resize() = imag(FluxFP);
   fname = "FluxF.projected.IM.No="+No_str+".out";
   save(datavec,fname);

   datavec.resize() = real(FluxFP2);
   fname = "FluxF2.projected.RE.No="+No_str+".out";
   save(datavec,fname);
   datavec.resize() = imag(FluxFP2);
   fname = "FluxF2.projected.IM.No="+No_str+".out";
   save(datavec,fname);

   datavec.resize() = real(FluxFPdelta);
   fname = "FluxF.diff.projected.RE.No="+No_str+".out";
   save(datavec,fname);
   datavec.resize() = imag(FluxFPdelta);
   fname = "FluxF.diff.projected.IM.No="+No_str+".out";
   save(datavec,fname);



  return 0;
} // main()





