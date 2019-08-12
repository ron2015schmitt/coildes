/************************************************************************* 
 * 
 *   File Name    : 
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  16 May 2003
 * 
 *
 *   SYNOPSIS     
 *     This is the main source file for the coil code.
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
#include "createsurface.hpp"




// Main Function for code

int main (int argc, char *argv[])
{

  disable_all_options();
  enable_option(opt_pf);
  enable_option(opt_ff);
  enable_option(opt_Nphi);
  enable_option(opt_Ntheta);
  enable_option(opt_Nnn);
  enable_option(opt_Nmm);
  enable_option(opt_Nharm);
  enable_option(opt_Mharm);
  enable_option(opt_zero);
  // parse command line input
 
   if (!parse_cmd(argc, argv))
     return 1;


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
   bool mode00 = true;
   modevectors(NF,nn,mm,Nnn,Nmm,Nharm,Mharm,mode00);


  // load the surface data
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


  // lay function onto grid 

  cout << endl;
  cout <<"$ Loading surface data onto "<<Ntheta<<" x "<<Nphi<<" (theta by phi) grid"<<endl;

  //load vector field, B, in cartesian coordinates
  Vector<p3vector<double> > B(Npts, "B");

  p3vectorformat::textformat(text_nobraces);
  B.textformat(text_nobraces);
  load(B,flux_filename);

  //find contravariant components of B
  Vector<double> Br(Npts, "Br");
  Vector<double> Bt(Npts, "Bt");
  Vector<double> Bp(Npts, "Bp");

  printcr("Find contravariant vector components.");

  STARTTIME(tbuff,ckstart);
  for (unsigned int j =0; j<Npts; j++) {
     const double J = dot(dx_dr[j],cross(dx_dtheta[j], dx_dphi[j]));
     Br[j] = dot(B[j],grad_r[j]);
     Bt[j] = dot(B[j],grad_theta[j]);
     Bp[j] = dot(B[j],grad_phi[j]);
  }
  STOPTIME(tbuff,ckstart);

  // Create fourier series
  cout << endl;
  cout<<"$ Generate fourier series matrix ("<<Npts<<" x "<<NF<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  Matrix<complex<double> > fs(Npts,NF,"fs");
  fseries(nn,mm,thetas,phis,fs);
  STOPTIME(tbuff,ckstart);

  cout << endl;
  printcr("Find fourier coef's for each contravariant component");

  STARTTIME(tbuff,ckstart);
  Vector<complex<double> > BrF(NF,"BrF");
  transformfunction(Br,BrF,fs);
  Vector<complex<double> > BtF(NF,"BtF");
  transformfunction(Bt,BtF,fs);
  Vector<complex<double> > BpF(NF,"BpF");
  transformfunction(Bp,BpF,fs);
  STOPTIME(tbuff,ckstart);



  string flux_rootname;
  string fn_ext;
  disect_filename(flux_filename,flux_rootname,fn_ext);


  printcr("Saving Fourier coefs");


  //  remove small modes that are (presumably) due to numerical error

  double maxR = max(max(abs(real(BrF))),max(abs(imag(BrF))));
  double maxT = max(max(abs(real(BtF))),max(abs(imag(BtF))));
  double maxP = max(max(abs(real(BpF))),max(abs(imag(BpF))));

  double largestmode;
  largestmode = max(maxR,maxT);
  largestmode = max(largestmode,maxP);
  dispcr(maxR);
  dispcr(maxT);
  dispcr(maxP);
  dispcr(largestmode);

  const double zeroval = relative_zero*largestmode;

  massage_absolute(BrF,zeroval);
  massage_absolute(BtF,zeroval);
  massage_absolute(BpF,zeroval);


  // ALL DONE, NOW SAVE TO FILES
  ostringstream strm_r;
  strm_r << flux_rootname << "_sup_r_F" <<".out";
  printcr(strm_r.str());
  save_coefs(strm_r.str(),CoefFileFormat_sincos,nn,mm,BrF);

  ostringstream strm_theta;
  strm_theta << flux_rootname << "_sup_theta_F" <<".out";
  printcr(strm_theta.str());
  save_coefs(strm_theta.str(),CoefFileFormat_sincos,nn,mm,BtF);

  ostringstream strm_phi;
  strm_phi << flux_rootname << "_sup_phi_F" <<".out";
  printcr(strm_phi.str());
  save_coefs(strm_phi.str(),CoefFileFormat_sincos,nn,mm,BpF);

ostringstream strm_r2;
  strm_r2 << flux_rootname << "_sup_r_Fexp" <<".out";
  printcr(strm_r2.str());
  Vector<double> BrF_Re(NF,"BrF_Re");
  BrF_Re.perline(1);
  BrF_Re.textformat(text_nobraces);
  BrF_Re =real(BrF);
  save(BrF_Re,strm_r2.str());

  ostringstream strm_theta2;
  strm_theta2 << flux_rootname << "_sup_theta_Fexp" <<".out";
  printcr(strm_theta2.str());
  Vector<double> BtF_Re(NF,"BtF_Re");
  BtF_Re.perline(1);
  BtF_Re.textformat(text_nobraces);
  BtF_Re =real(BtF);
  save(BtF_Re,strm_theta2.str());

  ostringstream strm_phi2;
  strm_phi2 << flux_rootname << "_sup_phi_Fexp" <<".out";
  printcr(strm_phi2.str());
  Vector<double> BpF_Re(NF,"BpF_Re");
  BpF_Re.perline(1);
  BpF_Re.textformat(text_nobraces);
  BpF_Re =real(BpF);
  save(BpF_Re,strm_phi2.str());

  // calculate RMS values and display
  double BrF_RMS = 0 ;
  double BtF_RMS = 0 ;
  double BpF_RMS = 0 ;
  for (unsigned int k =0; k<NF; k++) {
    double temp;
    
    temp = abs(BrF[k]);
    BrF_RMS += temp*temp;

    temp = abs(BtF[k]);
    BtF_RMS += temp*temp;

    temp = abs(BpF[k]);
    BpF_RMS += temp*temp;
  }

  BrF_RMS = sqrt(BrF_RMS)/(2*PI);
  BtF_RMS = sqrt(BtF_RMS)/(2*PI);
  BpF_RMS = sqrt(BpF_RMS)/(2*PI);

  double Br_RMS = 0 ;
  double Bt_RMS = 0 ;
  double Bp_RMS = 0 ;
  double Bmag_RMS = 0 ;
  double Bmag_RMS2 = 0 ;
  for (unsigned int j =0; j<Npts; j++) {
    double temp;

    temp = abs(Br[j]);
    Br_RMS += temp*temp;

    temp = abs(Bt[j]);
    Bt_RMS += temp*temp;

    temp = abs(Bp[j]);
    Bp_RMS += temp*temp;

    Bmag_RMS += dot(B[j],B[j]);

    Bmag_RMS2 += abs(Br[j])*abs(Br[j]) * dot(dx_dr[j],dx_dr[j]);
    Bmag_RMS2 += abs(Br[j])*abs(Bt[j]) * dot(dx_dr[j],dx_dtheta[j]);
    Bmag_RMS2 += abs(Br[j])*abs(Bp[j]) * dot(dx_dr[j],dx_dphi[j]);

    Bmag_RMS2 += abs(Bt[j])*abs(Br[j]) * dot(dx_dtheta[j],dx_dr[j]);
    Bmag_RMS2 += abs(Bt[j])*abs(Bt[j]) * dot(dx_dtheta[j],dx_dtheta[j]);
    Bmag_RMS2 += abs(Bt[j])*abs(Bp[j]) * dot(dx_dtheta[j],dx_dphi[j]);

    Bmag_RMS2 += abs(Bp[j])*abs(Br[j]) * dot(dx_dphi[j],dx_dr[j]);
    Bmag_RMS2 += abs(Bp[j])*abs(Bt[j]) * dot(dx_dphi[j],dx_dtheta[j]);
    Bmag_RMS2 += abs(Bp[j])*abs(Bp[j]) * dot(dx_dphi[j],dx_dphi[j]);

  }
  Br_RMS = sqrt(Br_RMS/Npts);
  Bt_RMS = sqrt(Bt_RMS/Npts);
  Bp_RMS = sqrt(Bp_RMS/Npts);
  Bmag_RMS = sqrt(Bmag_RMS/Npts);
  Bmag_RMS2 = sqrt(Bmag_RMS2/Npts);

  dispcr(Br_RMS);
  dispcr(BrF_RMS);cr();
  dispcr(Bt_RMS);
  dispcr(BtF_RMS);cr();
  dispcr(Bp_RMS);
  dispcr(BpF_RMS);cr();
  dispcr(Bmag_RMS);cr();
  dispcr(Bmag_RMS2);cr();

  return 0;
} // main()
