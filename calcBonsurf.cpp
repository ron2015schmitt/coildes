/************************************************************************* 
 * 
 *   File Name    :  
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  16 May 2003
 * 
 *
 *   SYNOPSIS     
 *     Calculate B field on teh given plasma surface
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

#include "bfield_ext.hpp"



 




// Main Function for code

int main (int argc, char *argv[])
{

  disable_all_options();
  enable_option(opt_pf);
  enable_option(opt_Nphi);
  enable_option(opt_Ntheta);
  enable_option(opt_Nnn);
  enable_option(opt_Nmm);

  // parse command line input
 
   if (!parse_cmd(argc, argv))
     return 1;

  string fext = "out";
  //  string fext = plasma_extname;
 //  ios_base::fmtflags flags = ios_base::right | ios_base::scientific;
  string ftemp;
  p3vectorformat::textformat(text_nobraces);

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


  // load the plasma surface fourier coef's

  cout << "$ Loading PLASMA SURFACE fourier coefficients from " << plasma_filename << endl;
  FourierSurface plasmafourier;
  if (load_fourier_surface(plasma_filename,plasmafourier))
    return 1;
  plasmafourier.RF().name("p.RF");
  plasmafourier.ZF().name("p.ZF");

  // print coef's
  printfouriercoefs(plasmafourier.nn(),plasmafourier.mm(),plasmafourier.RF(),plasmafourier.ZF(),10,18);
  
 
  // lay plasma surface onto grid 
  
  Vector<p3vector<double> > X(Npts, "X");
  Vector<p3vector<double> > dAdtdp(Npts, "dAdtdp");
  Vector<p3vector<double> > dx_dr(Npts, "dx_dr");
  Vector<p3vector<double> > dx_dtheta(Npts,"dx_dtheta");
  Vector<p3vector<double> > dx_dphi(Npts,"dx_dphi");
  Vector<p3vector<double> > grad_r(Npts,"grad_r");
  Vector<p3vector<double> > grad_theta(Npts,"grad_theta");
  Vector<p3vector<double> > grad_phi(Npts,"grad_phi");
 
  cout << endl;
  cout <<"$ Mapping plasma surface fourier coefficients to "<<Ntheta<<" x "<<Nphi<<" (theta by phi) grid"<<endl;

  STARTTIME(tbuff,ckstart);

  expandsurfaceandbases(X,dAdtdp,dx_dr,dx_dtheta,dx_dphi,grad_r,grad_theta,grad_phi,plasmafourier,thetas,phis);

  STOPTIME(tbuff,ckstart);


  // calculate magnetic field on plasma surface


  cout << endl;
  cout<<"$ Calc B field "<<endl;

  STARTTIME(tbuff,ckstart);

  Vector<p3vector<double> > Btotal(Npts, "Btotal");
  for (unsigned int j =0; j<Npts; j++)
    bTotal(X[j], Btotal[j]);
 
  Vector<p3vector<double> > Bplasma(Npts, "Bplasma");
  for (unsigned int j =0; j<Npts; j++)
    bplasma(X[j], Bplasma[j]);


  Vector<p3vector<double> > Bext(Npts, "Bext");
  for (unsigned int j =0; j<Npts; j++)
    bext(X[j], Bext[j]);

  STOPTIME(tbuff,ckstart);

  // calculate magnetic field flux

  cout << endl;
  cout<<"$ Calc total Bnormal "<<endl;

  STARTTIME(tbuff,ckstart);

  Vector<p3vector<double> > n(Npts, "n");
  for (unsigned int j =0; j<Npts; j++)
    n[j] = dAdtdp[j] / norm(dAdtdp[j]);

  Vector<double> Bn(Npts, "Bn");
  for (unsigned int j =0; j<Npts; j++)
    Bn[j] = dot(Btotal[j], n[j]);


  Vector<double> Bn_plasma(Npts, "Bn_plasma");
  for (unsigned int j =0; j<Npts; j++)
    Bn_plasma[j] = dot(Bplasma[j], n[j]);

  Vector<double> Flux_plasma(Npts, "Flux_plasma");
  for (unsigned int j =0; j<Npts; j++)
    Flux_plasma[j] = dot(Bplasma[j], dAdtdp[j]);

  Vector<double> Bn_ext(Npts, "Bn_ext");
  for (unsigned int j =0; j<Npts; j++)
    Bn_ext[j] = dot(Bext[j], n[j]);

  STOPTIME(tbuff,ckstart);

  

  // Create fourier series


  cout << endl;
  cout<<"$ Generate fourier series matrix ("<<Npts<<" x "<<NF<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  Matrix<complex<double> > fs(Npts,NF,"fs");
  fseries(nn,mm,thetas,phis,fs);
  STOPTIME(tbuff,ckstart);




  cout << endl;
  printcr("Find fourier coef's of magnetic field");
  STARTTIME(tbuff,ckstart);

  Vector<double> Br(Npts, "Br");
  for (unsigned int j =0; j<Npts; j++)
    Br[j] =  dot ( Btotal[j], grad_r[j] );
  Vector<complex<double> > BrF(NF,"BrF");
  transformfunction(Br,BrF,fs);
 
  Vector<double> Bt(Npts, "Bt");
  for (unsigned int j =0; j<Npts; j++)
    Bt[j] =  dot ( Btotal[j], grad_theta[j] );
  Vector<complex<double> > BtF(NF,"BtF");
  transformfunction(Bt,BtF,fs);

 Vector<double> Bp(Npts, "Bp");
  for (unsigned int j =0; j<Npts; j++)
    Bp[j] =  dot ( Btotal[j], grad_phi[j] );
  Vector<complex<double> > BpF(NF,"BpF");
  transformfunction(Bp,BpF,fs);


  Vector<complex<double> > BnF(NF,"BnF");
  transformfunction(Bn,BnF,fs);

  Vector<complex<double> > Bn_plasmaF(NF,"Bn_plasmaF");
  transformfunction(Bn_plasma,Bn_plasmaF,fs);

  Vector<complex<double> > Flux_plasmaF(NF,"Flux_plasmaF");
  transformfunction(Flux_plasma,Flux_plasmaF,fs);

  Vector<complex<double> > Bn_extF(NF,"Bn_extF");
  transformfunction(Bn_ext,Bn_extF,fs);


  STOPTIME(tbuff,ckstart);

  massage(BrF,1e-10);
  massage(BtF,1e-10);
  massage(BpF,1e-10);
  massage(BnF,1e-10);
  massage(Bn_plasmaF,1e-10);
  massage(Flux_plasmaF,1e-10);
  massage(Bn_extF,1e-10);
  save_coefs("BrTOTALF.out",CoefFileFormat_sincos,nn,mm,BrF);
  save_coefs("BtTOTALF.out",CoefFileFormat_sincos,nn,mm,BtF);
  save_coefs("BpTOTALF.out",CoefFileFormat_sincos,nn,mm,BpF);
  save_coefs("BnTOTALF.out",CoefFileFormat_sincos,nn,mm,BnF);
  save_coefs("BnF_negplasmaF.out",CoefFileFormat_sincos,nn,mm,-Bn_plasmaF);
  save_coefs("FluxF_negplasmaF.out",CoefFileFormat_sincos,nn,mm,-Flux_plasmaF);
  save_coefs("BnF_extF.out",CoefFileFormat_sincos,nn,mm,Bn_extF);


  Bn_plasma.textformat(text_nobraces);
  Bn_plasma.perline(1);
  Bn_plasma =- Bn_plasma;
  save(Bn_plasma,"Flux_negplasma.out");
  Bn_plasma =- Bn_plasma;

  Vector<double> Bx(Npts, "Bx");
  Vector<double> By(Npts, "By");
  Vector<double> Bz(Npts, "Bz");

  for (unsigned int j =0; j<Npts; j++) {
     Bx[j] = Btotal[j].x();
     By[j] = Btotal[j].y();
     Bz[j] = Btotal[j].z();
  }

  save(Bx,"BTOTALx.out");
  save(By,"BTOTALy.out");
  save(Bz,"BTOTALz.out");



  for (unsigned int j =0; j<Npts; j++) {
     p3vector<double> B = Br[j] * dx_dr[j];
     Bx[j] = B.x();
     By[j] = B.y();
     Bz[j] = B.z();
  }

  save(Bx,"Brx.out");
  save(By,"Bry.out");
  save(Bz,"Brz.out");


  for (unsigned int j =0; j<Npts; j++) {
     p3vector<double> B = Bt[j] * dx_dtheta[j];
     Bx[j] = B.x();
     By[j] = B.y();
     Bz[j] = B.z();
  }

  save(Bx,"Btx.out");
  save(By,"Bty.out");
  save(Bz,"Btz.out");

  for (unsigned int j =0; j<Npts; j++) {
     p3vector<double> B = Bp[j] * dx_dphi[j];
     Bx[j] = B.x();
     By[j] = B.y();
     Bz[j] = B.z();
  }

  save(Bx,"Bpx.out");
  save(By,"Bpy.out");
  save(Bz,"Bpz.out");


  // calculate RMS values and display
  double BrF_RMS = 0 ;
  double BtF_RMS = 0 ;
  double BpF_RMS = 0 ;
  double Bn_plasmaF_RMS = 0 ;
  double Flux_plasmaF_RMS = 0 ;
  for (unsigned int k =0; k<NF; k++) {
    double temp;
    
    temp = abs(BrF[k]);
    BrF_RMS += temp*temp;

    temp = abs(BtF[k]);
    BtF_RMS += temp*temp;

    temp = abs(BpF[k]);
    BpF_RMS += temp*temp;


    temp = abs(Bn_plasmaF[k]);
    Bn_plasmaF_RMS += temp*temp;

    temp = abs(Flux_plasmaF[k]);
    Flux_plasmaF_RMS += temp*temp;

  }
  BrF_RMS = sqrt(BrF_RMS)/(2*PI);
  BtF_RMS = sqrt(BtF_RMS)/(2*PI);
  BpF_RMS = sqrt(BpF_RMS)/(2*PI);
  Bn_plasmaF_RMS = sqrt(Bn_plasmaF_RMS)/(2*PI);
  Flux_plasmaF_RMS = sqrt(Flux_plasmaF_RMS)/(2*PI);

  double B_RMS = 0 ;
  double Br_RMS = 0 ;
  double Bt_RMS = 0 ;
  double Bp_RMS = 0 ;
  double Bn_plasma_RMS = 0 ;
  double Flux_plasma_RMS = 0 ;
  for (unsigned int j =0; j<Npts; j++) {
    double temp;

    temp = dot(Btotal[j],Btotal[j]);
    B_RMS += temp;

    temp = abs(Br[j]);
    Br_RMS += temp*temp;

    temp = abs(Bt[j]);
    Bt_RMS += temp*temp;

    temp = abs(Bp[j]);
    Bp_RMS += temp*temp;

    Bn_plasma_RMS += Bn_plasma[j]*Bn_plasma[j];
    Flux_plasma_RMS += Flux_plasma[j]*Flux_plasma[j];
  }
  B_RMS = sqrt(B_RMS/Npts);
  Br_RMS = sqrt(Br_RMS/Npts);
  Bt_RMS = sqrt(Bt_RMS/Npts);
  Bp_RMS = sqrt(Bp_RMS/Npts);
  Bn_plasma_RMS = sqrt(Bn_plasma_RMS/Npts);
  Flux_plasma_RMS = sqrt(Flux_plasma_RMS/Npts);


  dispcr(B_RMS);cr();
  dispcr(Br_RMS);
  dispcr(BrF_RMS);cr();
  dispcr(Bt_RMS);
  dispcr(BtF_RMS);cr();
  dispcr(Bp_RMS);
  dispcr(BpF_RMS);cr();
  dispcr(Bn_plasma_RMS);
  dispcr(Bn_plasmaF_RMS);cr();
  dispcr(Flux_plasma_RMS);
  dispcr(Flux_plasmaF_RMS);cr();



  save(B_RMS,"BTOTAL.RMS.out");
  save(Br_RMS,"BrTOTAL.RMS.out");
  
  return 0;
} // main()





