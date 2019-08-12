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

  LAvector<p3vector<double> > X(Npts, "X");
  LAvector<p3vector<double> > dAdtdp(Npts, "dAdtdp");
  LAvector<p3vector<double> > dx_dr(Npts, "dx_dr");
  LAvector<p3vector<double> > dx_dtheta(Npts,"dx_dtheta");
  LAvector<p3vector<double> > dx_dphi(Npts,"dx_dphi");
  LAvector<p3vector<double> > grad_r(Npts,"grad_r");
  LAvector<p3vector<double> > grad_theta(Npts,"grad_theta");
  LAvector<p3vector<double> > grad_phi(Npts,"grad_phi");

  cout << endl;
  cout <<"$ Mapping plasma surface fourier coefficients to "<<Ntheta<<" x "<<Nphi<<" (theta by phi) grid"<<endl;

  STARTTIME(tbuff,ckstart);

  expandsurfaceandbases(X,dAdtdp,dx_dr,dx_dtheta,dx_dphi,grad_r,grad_theta,grad_phi,plasmafourier,thetas,phis);

  STOPTIME(tbuff,ckstart);


  // lay function onto grid 

  cout << endl;
  cout <<"$ Loading surface data onto "<<Ntheta<<" x "<<Nphi<<" (theta by phi) grid"<<endl;

  //load vector field, B, in cartesian coordinates
  LAvector<p3vector<double> > B(Npts, "B");

  p3vectorformat::textformat(text_nobraces);
  B.textformat(text_nobraces);
  load(B,flux_filename);

  //find contravariant components of B
  LAvector<double> JBr(Npts, "JBr");
  LAvector<double> JBt(Npts, "JBt");
  LAvector<double> JBp(Npts, "JBp");
  LAvector<double> J(Npts, "J");

  printcr("Find contravariant*J vector components.");

  STARTTIME(tbuff,ckstart);
  for (unsigned int j =0; j<Npts; j++) {
     J[j] = dot(dx_dr[j],cross(dx_dtheta[j], dx_dphi[j]));
     JBr[j] = J[j]*dot(B[j],grad_r[j]);
     JBt[j] = J[j]*dot(B[j],grad_theta[j]);
     JBp[j] = J[j]*dot(B[j],grad_phi[j]);
  }
  STOPTIME(tbuff,ckstart);


  //as a test, recreate B
  LAvector<p3vector<double> > B2(Npts, "B2");
  for (unsigned int j =0; j<Npts; j++) {
     B2[j] = JBr[j]*dx_dr[j] + JBt[j]*dx_dtheta[j] + JBp[j]*dx_dphi[j];
     B2[j] = B2[j]/J[j];
  }

  LAvector<p3vector<double> > Bdiff(Npts, "Bdiff");
  Bdiff = B2 - B;

  LAvector<double> Bmag(Npts, "Bmag");
  LAvector<double> B2mag(Npts, "B2mag");
  LAvector<double> Bdiffmag(Npts, "Bmag");
  for (unsigned int j =0; j<Npts; j++) {
     Bmag[j] = norm(B[j]);
     B2mag[j] = norm(B2[j]);
     Bdiffmag[j] = norm(Bdiff[j]);
  }
  dispcr(max(Bmag));
  dispcr(max(B2mag));
  dispcr(max(Bdiffmag));
  dispcr(mean(Bmag));
  dispcr(mean(B2mag));
  dispcr(mean(Bdiffmag));

  


  // Create fourier series
  cout << endl;
  cout<<"$ Generate fourier series matrix ("<<Npts<<" x "<<NF<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  Matrix<complex<double> > fs(Npts,NF,"fs");
  fseries(nn,mm,thetas,phis,fs);
  STOPTIME(tbuff,ckstart);

  cout << endl;
  printcr("Find fourier coef's for each contravariant*J component");

  STARTTIME(tbuff,ckstart);
  LAvector<complex<double> > JBrF(NF,"JBrF");
  transformfunction(JBr,JBrF,fs);
  LAvector<complex<double> > JBtF(NF,"JBtF");
  transformfunction(JBt,JBtF,fs);
  LAvector<complex<double> > JBpF(NF,"JBpF");
  transformfunction(JBp,JBpF,fs);
  STOPTIME(tbuff,ckstart);





  string flux_rootname;
  string fn_ext;
  disect_filename(flux_filename,flux_rootname,fn_ext);


  printcr("Saving Fourier coefs");


  //  remove small modes that are (presumably) due to numerical error

  double maxR = max(max(abs(real(JBrF))),max(abs(imag(JBrF))));
  double maxT = max(max(abs(real(JBtF))),max(abs(imag(JBtF))));
  double maxP = max(max(abs(real(JBpF))),max(abs(imag(JBpF))));

  double largestmode;
  largestmode = max(maxR,maxT);
  largestmode = max(largestmode,maxP);
  dispcr(maxR);
  dispcr(maxT);
  dispcr(maxP);
  dispcr(largestmode);

  const double zeroval = relative_zero*largestmode;

  massage_absolute(JBrF,zeroval);
  massage_absolute(JBtF,zeroval);
  massage_absolute(JBpF,zeroval);


  // ALL DONE, NOW SAVE TO FILES
  ostringstream strm_r;
  strm_r << flux_rootname << "_sup_rJ_F" <<".out";
  printcr(strm_r.str());
  save_coefs(strm_r.str(),CoefFileFormat_sincos,nn,mm,JBrF);

  ostringstream strm_theta;
  strm_theta << flux_rootname << "_sup_thetaJ_F" <<".out";
  printcr(strm_theta.str());
  save_coefs(strm_theta.str(),CoefFileFormat_sincos,nn,mm,JBtF);

  ostringstream strm_phi;
  strm_phi << flux_rootname << "_sup_phiJ_F" <<".out";
  printcr(strm_phi.str());
  save_coefs(strm_phi.str(),CoefFileFormat_sincos,nn,mm,JBpF);


  ostringstream strm_r2;
  strm_r2 << flux_rootname << "_sup_rJ_Fexp" <<".out";
  printcr(strm_r2.str());
  LAvector<double> JBrF_Re(NF,"JBrF_Re");
  JBrF_Re.perline(1);
  JBrF_Re.textformat(text_nobraces);
  JBrF_Re =real(JBrF);
  save(JBrF_Re,strm_r2.str());

  ostringstream strm_theta2;
  strm_theta2 << flux_rootname << "_sup_thetaJ_Fexp" <<".out";
  printcr(strm_theta2.str());
  LAvector<double> JBtF_Re(NF,"JBtF_Re");
  JBtF_Re.perline(1);
  JBtF_Re.textformat(text_nobraces);
  JBtF_Re =real(JBtF);
  save(JBtF_Re,strm_theta2.str());

  ostringstream strm_phi2;
  strm_phi2 << flux_rootname << "_sup_phiJ_Fexp" <<".out";
  printcr(strm_phi2.str());
  LAvector<double> JBpF_Re(NF,"JBpF_Re");
  JBpF_Re.perline(1);
  JBpF_Re.textformat(text_nobraces);
  JBpF_Re =real(JBpF);
  save(JBpF_Re,strm_phi2.str());





  LAvector<complex<double> > JF(NF,"JF");
  transformfunction(J,JF,fs);
  massage(JF,1e-10);
  ostringstream strm_JF;
  strm_JF << flux_rootname << "_JF" <<".out";
  printcr(strm_JF.str());
  save_coefs(strm_JF.str(),CoefFileFormat_sincos,nn,mm,JF);



  LAvector<complex<double> > JdivBF(NF,"JdivBF");
  for (unsigned int k =0; k<NF; k++)
     JdivBF[k] = complex<double>(0.0,1.0)*(mm[k]*JBtF[k] + nn[k]*JBpF[k]);
  massage(JdivBF,1e-10);
  ostringstream strm_JdivBF;
  strm_JdivBF << flux_rootname << "_JdivBF" <<".out";
  printcr(strm_JdivBF.str());
  save_coefs(strm_JdivBF.str(),CoefFileFormat_sincos,nn,mm,JdivBF);

  //  for (unsigned int k =0; k<NF; k++) {
  // if (abs(nn(k)) <=4.00000001){
  //disp(mm[k]); disp(nn[k]); disp(mm[k]*JBtF[k]); disp(nn[k]*JJBpF[k]);cr();
  // }
     //     disp(mm[k]); disp(nn[k]); disp(JF[k]); disp(JF[k]);cr();
  //}

  ///////////////////////////////////////////////////////////
//   JdivBF.perline(1);
//   JdivBF.textformat(text_nobraces); 
//   save(JdivBF,"JdivBF.out");
  ///////////////////////////////////////////////////////////

  LAvector<double> JdivB(Npts, "JdivB");
  expandfunction(JdivB,JdivBF,fs);

  LAvector<double> divB(Npts, "divB");
  for (unsigned int j =0; j<Npts; j++)
    divB[j] =  JdivB[j]/J[j];
  LAvector<complex<double> > divBF(NF,"divBF");
  transformfunction(divB,divBF,fs);
  massage(divBF,1e-10);
  ostringstream strm_divBF;
  strm_divBF << flux_rootname << "_divBF" <<".out";
  printcr(strm_divBF.str());
  save_coefs(strm_divBF.str(),CoefFileFormat_sincos,nn,mm,divBF);









  // calculate RMS values and display
  double JBrF_RMS = 0 ;
  double JBtF_RMS = 0 ;
  double JBpF_RMS = 0 ;
  for (unsigned int k =0; k<NF; k++) {
    double temp;
    
    temp = abs(JBrF[k]);
    JBrF_RMS += temp*temp;

    temp = abs(JBtF[k]);
    JBtF_RMS += temp*temp;

    temp = abs(JBpF[k]);
    JBpF_RMS += temp*temp;
  }

  JBrF_RMS = sqrt(JBrF_RMS)/(2*PI);
  JBtF_RMS = sqrt(JBtF_RMS)/(2*PI);
  JBpF_RMS = sqrt(JBpF_RMS)/(2*PI);

  double JBr_RMS = 0 ;
  double JBt_RMS = 0 ;
  double JBp_RMS = 0 ;
  double Bmag_RMS = 0 ;
  double Bmag_RMS2 = 0 ;
  for (unsigned int j =0; j<Npts; j++) {
    double temp;

    temp = abs(JBr[j]);
    JBr_RMS += temp*temp;

    temp = abs(JBt[j]);
    JBt_RMS += temp*temp;

    temp = abs(JBp[j]);
    JBp_RMS += temp*temp;

    Bmag_RMS += dot(B[j],B[j]);

    Bmag_RMS2 += abs(JBr[j])*abs(JBr[j]) * dot(dx_dr[j],dx_dr[j]);
    Bmag_RMS2 += abs(JBr[j])*abs(JBt[j]) * dot(dx_dr[j],dx_dtheta[j]);
    Bmag_RMS2 += abs(JBr[j])*abs(JBp[j]) * dot(dx_dr[j],dx_dphi[j]);

    Bmag_RMS2 += abs(JBt[j])*abs(JBr[j]) * dot(dx_dtheta[j],dx_dr[j]);
    Bmag_RMS2 += abs(JBt[j])*abs(JBt[j]) * dot(dx_dtheta[j],dx_dtheta[j]);
    Bmag_RMS2 += abs(JBt[j])*abs(JBp[j]) * dot(dx_dtheta[j],dx_dphi[j]);

    Bmag_RMS2 += abs(JBp[j])*abs(JBr[j]) * dot(dx_dphi[j],dx_dr[j]);
    Bmag_RMS2 += abs(JBp[j])*abs(JBt[j]) * dot(dx_dphi[j],dx_dtheta[j]);
    Bmag_RMS2 += abs(JBp[j])*abs(JBp[j]) * dot(dx_dphi[j],dx_dphi[j]);

  }
  JBr_RMS = sqrt(JBr_RMS/Npts);
  JBt_RMS = sqrt(JBt_RMS/Npts);
  JBp_RMS = sqrt(JBp_RMS/Npts);
  Bmag_RMS = sqrt(Bmag_RMS/Npts);
  Bmag_RMS2 = sqrt(Bmag_RMS2/Npts);

  dispcr(JBr_RMS);
  dispcr(JBrF_RMS);cr();
  dispcr(JBt_RMS);
  dispcr(JBtF_RMS);cr();
  dispcr(JBp_RMS);
  dispcr(JBpF_RMS);cr();
  dispcr(Bmag_RMS);cr();
  dispcr(Bmag_RMS2);cr();

  return 0;
} // main()
