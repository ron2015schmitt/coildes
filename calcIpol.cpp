/************************************************************************* 
 * 
 *   File Name    :  
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  16 May 2003
 * 
 *
 *   SYNOPSIS     
 *
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


const double NEGLECT =  1e-12;




// Main Function for code

int main (int argc, char *argv[])
{
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
  bool mode00 = true;
  if ( (Nharm >1) ||(Mharm>1) )
     modevectors(NF,nn,mm,Nnn,Nmm,Nharm,Mharm,mode00);
  else
     modevectors(NF,nn,mm,Nnn,Nmm,mode00);



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


  STOPTIME(tbuff,ckstart);



  // load or calculate B field
  cout << endl;
  cout<<"$ Load tangent B field "<<endl;

  STARTTIME(tbuff,ckstart);

  LAvector<p3vector<double> > B(Npts, "B");


  cout <<endl<< "$ Loading BTOTAL_theta fourier coefficients from " << Bt_filename << endl;
  cout <<endl<< "$ Loading BTOTAL_phi fourier coefficients from " << Bp_filename << endl;

  LAvector<complex<double> > BtF(NF,"BtF");
  if (load_coefs( Bt_filename,CoefFileFormat_sincos,nn,mm,BtF))
    return 5;
  LAvector<complex<double> > BpF(NF,"BpF");
  if (load_coefs( Bp_filename,CoefFileFormat_sincos,nn,mm,BpF))
    return 6;

  // contravariant components
  LAvector<double> Bt(Npts, "Bt");
  expandfunction(Bt,BtF,fs);
  LAvector<double> Bp(Npts, "Bp");
  expandfunction(Bp,BpF,fs);

  // covariant components
  LAvector<double> B_t(Npts, "B_t");
  LAvector<double> B_p(Npts, "B_p");

  for (unsigned int j =0; j<Npts; j++) {
    B[j] =  Bt[j] * dx_dtheta[j] + Bp[j] * dx_dphi[j];
    B_t[j] = dot( B[j], dx_dtheta[j] );
    B_p[j] = dot( B[j], dx_dphi[j] );
  }



  double dtheta = 2*M_PI/((double)Ntheta);
  double dphi = 2*M_PI/((double)Nphi);

  // calculate Ntheta estimates of Ipol
  LAvector<double> Ipol(Ntheta, "Ipol");

  for (unsigned int tt = 0; tt<Ntheta ; tt++){
     double temp =0;
    for (unsigned int pp = 0; pp<Nphi; pp++) {
       unsigned int ind = tt*Nphi + pp;
       temp += B_p[ind];
    }
     Ipol[tt] = temp*dphi/mu0;
     //     dispcr(Ipol[tt]);
  }

  cout.precision(10);

  double Ipol_avg = sum(Ipol)/double(Ntheta);
  dispcr(Ipol_avg);

  double Ipol_stddev = norm(Ipol-Ipol_avg)/sqrt(double(Ntheta));
  dispcr(Ipol_stddev);

  double Ipol_rms_err = Ipol_stddev/Ipol_avg;
  dispcr(Ipol_rms_err);

  double Ipol_min = min(Ipol);
  dispcr(Ipol_min);

  double Ipol_max = max(Ipol);
  dispcr(Ipol_max);






  // calculate Ntheta estimates of Itor
  LAvector<double> Itor(Nphi, "Itor");

  for (unsigned int pp = 0; pp<Nphi; pp++) {
     double temp =0;
     for (unsigned int tt = 0; tt<Ntheta ; tt++){
       unsigned int ind = tt*Nphi + pp;
       temp += B_t[ind];// * norm(dx_dphi[ind]);
    }
     Itor[pp] = temp*dtheta/mu0;
  }

  double Itor_avg = sum(Itor)/double(Nphi);
  dispcr(Itor_avg);

  double Itor_stddev = norm(Itor-Itor_avg)/sqrt(double(Nphi));
  dispcr(Itor_stddev);

  double Itor_rms_err = Itor_stddev/Itor_avg;
  dispcr(Itor_rms_err);

  double Itor_min = min(Itor);
  dispcr(Itor_min);

  double Itor_max = max(Itor);
  dispcr(Itor_max);

  double Itor_div_Ipol = Itor_avg/Ipol_avg;
  dispcr(Itor_div_Ipol);

  return 0;
} // main()





