/************************************************************************* 
 * 
 *   File Name    :  coilfwd2.cpp
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  16 May 2003
 * 
 *
 *   SYNOPSIS     
 *     This file finds the plasma surface flux for a given coil current.
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
#include "coilfft.hpp"

const double NEGLECT =  1e-12;




// Main Function for code

int main (int argc, char *argv[])
{

  disable_all_options();
  enable_option(opt_pf);
  enable_option(opt_cf);
  enable_option(opt_if);
  enable_option(opt_Nphi);
  enable_option(opt_Ntheta);
  enable_option(opt_Nnn);
  enable_option(opt_Nmm);
   enable_option(opt_Nharm);
   enable_option(opt_Mharm);
  enable_option(opt_Itoroidal);
  enable_option(opt_Ipoloidal);

  // parse command line input
 
  if (!parse_cmd(argc, argv))
     return 1;


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
    bool mode00 = true;
    if ( (Nharm >1) ||(Mharm>1) )
       modevectors(NF,nn,mm,Nnn,Nmm,Nharm,Mharm,mode00);
    else
       modevectors(NF,nn,mm,Nnn,Nmm,mode00);



   // these exclude the n=0,m=0 case
   LAvector<double> nnR("nnR");
   LAvector<double> mmR("mmR");
   unsigned int NFR;
   bool modeR00 = false;
   if ( (Nharm >1) ||(Mharm>1) )
      modevectors(NFR,nnR,mmR,Nnn,Nmm,Nharm,Mharm,modeR00);
   else
      modevectors(NFR,nnR,mmR,Nnn,Nmm,modeR00);


  // coefficient C is the integration coef for the fourier transform
  // C = dtheta*dphi
  //   = (2*pi/Ntheta)*(2*pi/Nphi)
  //   = (2*pi*2*pi/Npts)
  const double C = (2*PI*2*PI/double(Npts));
  const double Csqr = C*C;



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
    


  // load the coil CURRENT fourier coef's
  // at some point, add code so that user can select type of coef's from command line

  cout <<endl<< "$ Loading COIL CURRENT fourier coefficients from " << current_filename << endl;
  LAvector<complex<double> > IF(NFR,"IF");
  if (load_coefs(current_filename,CoefFileFormat_sincos,nnR,mmR,IF))
    return 3;  
  LAvector<double> IFreal(NFR,"IFreal");
  LAvector<double> IFimag(NFR,"IFimag");
  IFreal = real(IF);
  IFimag = imag(IF);

  // print coef's
  //  printfouriercoefs(nnR,mmR,IFreal,IFimag,10,18);

 
 
  // lay plasma surface onto grid 
  
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








  ////////////////////////////////////////

  // Create theta mutual inductance matrix

  cout << endl;
  cout<<"$ Generating theta inductance matrix ("<<Npts<<" x "<<Npts<<")"<<endl;

  STARTTIME(tbuff,ckstart);
  LAvector<p3vector<double> > Jgradtheta(Npts,"Jgradtheta");
  for (unsigned int j =0; j<Npts; j++) {
    Jgradtheta[j] = J[j] * grad_theta[j];
  }
  Matrix<double> Mtheta(Npts, Npts, "Mtheta");
  inductancematrix(X,Jgradtheta,Xcoil,dA_dtdp_coil,Mtheta);
  STOPTIME(tbuff,ckstart);


  // fft of coil side of Mtheta matrix
  cout << endl;
  cout<<"$ FFT of theta inductance matrix ("<<Npts<<" x "<<NFR<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  Matrix<complex<double> > MFtheta1(Npts,NFR,"MFtheta1"); 

  half_fft_of_M(Mtheta, MFtheta1, Nphi, Ntheta, Nnn, Nmm, Nharm,Mharm,1e-14);

  STOPTIME(tbuff,ckstart);


  // find FFT (only n=0,m=0) of plasma side of Mtheta matrix
  cout << endl;
  cout<<"$ find FFT (only n=0,m=0) of plasma side of Mtheta matrix ("<<1<<" x "<<NFR<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  LAvector<complex<double> > MFtheta00(NFR,"MFtheta00"); 
  for (unsigned int k =0; k<NFR; k++) {
    double coef = 1/double(Npts);
     complex<double> temp = 0;
     for (unsigned int j =0; j<Npts; j++) {
	temp = temp + MFtheta1(j,k);
     }
     MFtheta00[k] = coef*temp;
  }
  

  STOPTIME(tbuff,ckstart);




   // Create phi mutual inductance matrix

  cout << endl;
  cout<<"$ Generating phi inductance matrix ("<<Npts<<" x "<<Npts<<")"<<endl;

  STARTTIME(tbuff,ckstart);
  LAvector<p3vector<double> > Jgradphi(Npts,"Jgradphi");
  for (unsigned int j =0; j<Npts; j++) {
    Jgradphi[j] = J[j] * grad_phi[j];
  }
  Matrix<double> Mphi(Npts, Npts, "Mphi");
  inductancematrix(X,Jgradphi,Xcoil,dA_dtdp_coil,Mphi);
  STOPTIME(tbuff,ckstart);


  // fft of coil side of Mphi matrix
  cout << endl;
  cout<<"$ FFT of phi inductance matrix ("<<Npts<<" x "<<NFR<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  Matrix<complex<double> > MFphi1(Npts,NFR,"MFphi1"); 

  half_fft_of_M(Mphi, MFphi1, Nphi, Ntheta, Nnn, Nmm, Nharm,Mharm,1e-14);

  STOPTIME(tbuff,ckstart);


  // find FFT (only n=0,m=0) of plasma side of Mphi matrix
  cout << endl;
  cout<<"$ find FFT (only n=0,m=0) of plasma side of Mphi matrix ("<<1<<" x "<<NFR<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  LAvector<complex<double> > MFphi00(NFR,"MFphi00"); 
  for (unsigned int k =0; k<NFR; k++) {
     double coef = 1/double(Npts);
     complex<double> temp = 0;
     for (unsigned int j =0; j<Npts; j++) {
	temp = temp + MFphi1(j,k);
     }
     MFphi00[k] = temp*coef;
  }

  STOPTIME(tbuff,ckstart);


  // calculate B components due to IF (without net currents)
  
  cout << endl;
  cout<<"$ calculate B components due to IF (without net currents)"<<endl;

  STARTTIME(tbuff,ckstart);

  double JB1t00;
  double JB1p00;
  complex<double> temp;
  temp = (MFtheta00|IF);
  dispcr(temp);
  JB1t00 = temp.real();

  temp = (MFphi00|IF);
  dispcr(temp);
  JB1p00 = temp.real();

  dispcr(JB1t00);
  dispcr(JB1p00);

  STOPTIME(tbuff,ckstart);



  // calculate Bfield from Ipol and Itor
  
  cout << endl;
  cout<<"$ B field due to Ipoloidal and Itoroidal (net currents)"<<endl;

  STARTTIME(tbuff,ckstart);

  // jcoil is actually j*Jacobian
  LAvector<p3vector<double> > jcoil(Npts,"jcoil");
  
  double j_theta =  (Ipoloidal)/(2*PI);
  double j_phi =  (Itoroidal)/(2*PI);
  for (unsigned int i = 0; i<Npts; i++) {
     jcoil[i] = j_theta * dx_dtheta_coil[i] + j_phi * dx_dphi_coil[i];
  }

  
  //  save(jcoil,string("iotafromIF.jcoil.out"));


  LAvector<p3vector<double> > B2(Npts,"B2");
  const double dphi_by_dtheta = 2.0*PI * 2.0*PI/double(Npts);
  const double CBS= mu0div4pi * dphi_by_dtheta;
  for (unsigned int jj = 0; jj<Npts; jj++) {
     p3vector<double> Btemp(0,0,0);
     for (unsigned int ii = 0; ii<Npts; ii++) {
	const p3vector<double> R = X[jj] - Xcoil[ii];
	const double r3 = pow(norm(R),3.0);
	Btemp = Btemp + cross(jcoil[ii],R) / r3;
     }
     B2[jj] = CBS*Btemp;
  }

  STOPTIME(tbuff,ckstart);


  cr();
  printcr("Find contravariant * J vector components for net current B");

  STARTTIME(tbuff,ckstart);
  LAvector<double > JB2t(Npts,"JB2t");
  LAvector<double > JB2p(Npts,"JB2p"); 
  for (unsigned int j = 0; j<Npts; j++) {
     JB2t[j] = J[j] * dot(B2[j],grad_theta[j]);
     JB2p[j] = J[j] * dot(B2[j],grad_phi[j]);
  }
  STOPTIME(tbuff,ckstart);


 
  // find FFT (only n=0,m=0) of J*Bt and JBp
  cout << endl;
  cout<<"$ find FFT (only n=0,m=0) of J*Bt"<<endl;

  STARTTIME(tbuff,ckstart);
  double JB2t00 = 0 ;
  double JB2p00 = 0 ;
  for (unsigned int j =0; j<Npts; j++) {
    JB2t00 = JB2t00 + JB2t[j];
    JB2p00 = JB2p00 + JB2p[j];
  }
  JB2t00 = JB2t00/double(Npts);
  JB2p00 = JB2p00/double(Npts);
  STOPTIME(tbuff,ckstart);


  dispcr(JB2t00);
  dispcr(JB2p00);





  // calculate total B and iota
  
  cout << endl;
  cout<<"$ Calculate iota"<<endl;

  STARTTIME(tbuff,ckstart);

  double JBTotalt00= JB1t00 + JB2t00;
  double JBTotalp00= JB1p00 + JB2p00;

  dispcr(JBTotalt00);
  dispcr(JBTotalp00);

  double iota = JBTotalt00/JBTotalp00;
  dispcr(iota);

  STOPTIME(tbuff,ckstart);


  ////////////////////////////////////////////////////////////





 

  return 0;
} // main()





