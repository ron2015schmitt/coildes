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
#include "omegasubspace.hpp"

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


  // Create Fourier Mode vectors
  LAvector<double> nn("nn");
  LAvector<double> mm("mm");
  unsigned int NF;
  modevectors(NF,nn,mm,Nnn,Nmm);

  // these exclude the n=0,m=0 case
  LAvector<double> nnR("nnR");
  LAvector<double> mmR("mmR");
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
  


  // load the coil surface fourier coef's

  cout << "$ Loading COIL SURFACE fourier coefficients from " << coil_filename << endl;
  FourierSurface coilfourier;
  if (load_fourier_surface(coil_filename,coilfourier))
    return 2;

  coilfourier.RF().name("c.RF");
  coilfourier.ZF().name("c.ZF");

  // print coef's
  //  printfouriercoefs(coilfourier.nn(),coilfourier.mm(),coilfourier.RF(),coilfourier.ZF(),10,18);
    


  LAvector<complex<double> > FluxF(NFR,"FluxF");

  cout <<endl<< "$ Loading Plasma Flux sin/cos fourier coefficients from " << flux_filename << endl;
  if (load_coefs(flux_filename,CoefFileFormat_sincos,nnR,mmR,FluxF))
    return 3;


 
 
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


  STOPTIME(tbuff,ckstart);


 
  

  // lay coil surface onto grid 
  
  LAvector<p3vector<double> > Xcoil(Npts, "Xcoil");
  LAvector<p3vector<double> > dA_dtdp_coil(Npts, "dA_dtdp_coil");

  cout << endl;
  cout <<"$ Mapping coil surface fourier coefficients to "<<Ntheta<<" x "<<Nphi<<" (theta by phi) grid"<<endl;

  STARTTIME(tbuff,ckstart);

  expandsurface(Xcoil,dA_dtdp_coil,coilfourier,thetas,phis);

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
  LAvector<double> S(Nmin1,"S");
  cooll_lapack::svd(Gamma,U,S,V);

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

  LAvector<p3vector<double> > B(Npts, "B");

  //  for (unsigned int j =0; j<Npts; j++)
  //    bTotal(X[j], B[j]);

//   LAvector<complex<double> > BrF(NF,"BrF");
//   if (load_coefs("BrTotal.sincos.out",CoefFileFormat_sincos,nn,mm,BrF))
//     return 4;
//   LAvector<complex<double> > Br(Npts, "Br");
//   expandfunction(Br,BrF,fs);


  cout <<endl<< "$ Loading BTOTAL_theta fourier coefficients from " << Bt_filename << endl;
  cout <<endl<< "$ Loading BTOTAL_phi fourier coefficients from " << Bp_filename << endl;

  LAvector<complex<double> > BtF(NF,"BtF");
  if (load_coefs( Bt_filename,CoefFileFormat_sincos,nn,mm,BtF))
    return 5;
  LAvector<complex<double> > BpF(NF,"BpF");
  if (load_coefs( Bp_filename,CoefFileFormat_sincos,nn,mm,BpF))
    return 6;

  LAvector<double> Bt(Npts, "Bt");
  expandfunction(Bt,BtF,fs);
  LAvector<double> Bp(Npts, "Bp");
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
  cooll_lapack::svd(Omega,U,S,V);

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

  LAvector<complex<double> > FluxFP(No,"FluxFP");
  FluxFP = (adj(Ur)|FluxF);

  STOPTIME(tbuff,ckstart);


  // Create mutual inductance matrix

  cout << endl;
  cout<<"$ Generating inductance matrix ("<<Npts<<" x "<<Npts<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  Matrix<double> M(Npts, Npts, "M");
  inductancematrix(X,dA_dtdp,Xcoil,dA_dtdp_coil,M);

  STOPTIME(tbuff,ckstart);



  // project M matrix into optimal space (from configuration space)
  cout << endl;
  cout<<"$ Project inductance matrix ("<<No<<" x "<<NFR<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  Matrix<complex<double> > MF(No,NFR,"MF"); 
  MF = (adj(f_phix)|M|fsR)*Csqr;

  STOPTIME(tbuff,ckstart);

  // clear M and its release its memory

  M.resize(0,0);



  // calculate SVD of MF
  cout << endl;
  cout<<"$ SVD of fourier inductance matrix ("<<No<<" x "<<NFR<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  
  U.resize(No,No);
  V.resize(NFR,NFR);
  S.resize(min(NFR,No));
  cooll_lapack::svd(MF,U,S,V);
               
  STOPTIME(tbuff,ckstart);

  //************************************************
//   data.resize() = real(U);
//   fname = "MF.scoild.No="+No_str+".UR.out";
//   save(data,fname);
//   data.resize() = imag(U);
//   fname = "MF.scoild.No="+No_str+".UI.out";
//   save(data,fname);

//   fname = "MF.scoild.No="+No_str+".S.out";
//   save(S,fname);

//   data.resize() = real(V);
//   fname = "MF.scoild.No="+No_str+".VR.out";
//   save(data,fname);
//   data.resize() = imag(V);
//   fname = "MF.scoild.No="+No_str+".VI.out";
//   save(data,fname);
  //************************************************






  // calculate inverse of MF
  cout << endl;
  cout<<"$ Pseudo Inverse of MF =("<<NFR<<" x "<<No<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  Matrix<complex<double> > MFinv(NFR,No,"MFinv");
  Matrix<double> Sinv(NFR,No,"Sinv");
  Sinv=0.0;

  const unsigned int Nmin = S.size();
  for(unsigned int k=0; k<Nmin;k++) {
    Sinv(k,k) = 1/S[k];
  }
  MFinv = (V|Sinv|adj(U));

  STOPTIME(tbuff,ckstart);



  // calculate current
  
  cout << endl;
  cout<<"$ Calculate coil current ("<<NFR<<" x 1)"<<endl;

  STARTTIME(tbuff,ckstart);

  LAvector<complex<double> > IF(NFR,"IF");
  IF = (MFinv|FluxFP);

  STOPTIME(tbuff,ckstart);


  // ALL DONE, NOW SAVE TO FILES


  massage(IF,1e-8);


  fname = "IF.scoild.No="+No_str+".out";
  save_coefs(fname,CoefFileFormat_sincos,nnR,mmR,IF);

  nn.perline(1);
  nn.textformat(text_nobraces);
  mm.perline(1);
  mm.textformat(text_nobraces);
  save(nn,"nn.out");
  save(mm,"mm.out");


  nnR.perline(1);
  nnR.textformat(text_nobraces);
  mmR.perline(1);
  mmR.textformat(text_nobraces);
  save(mmR,"mmR.out");
  save(nnR,"nnR.out");


  // auxilliary files


   

//   data.resize() = real(fs);
//   fname = "fsR.scoild.No="+No_str+".exp.out";
//   save(data,fname);

//   data.resize() = imag(fs);
//   fname = "fsI.scoild.No="+No_str+".exp.out";
//   save(data,fname);


//   data.resize() = real(fsR2);
//   fname = "fsR2R.scoild.No="+No_str+".exp.out";
//   save(data,fname);

//   data.resize() = imag(fsR2);
//   fname = "fsR2I.scoild.No="+No_str+".exp.out";
//   save(data,fname);



  return 0;
} // main()





