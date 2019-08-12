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
  printcr("(contravariant * J vector components = 2 form components).");

  STARTTIME(tbuff,ckstart);

  cout <<endl<< "$ Loading BTOTAL_theta fourier coefficients from " << Bt_filename << endl;

  Vector<complex<double> > JBtF(NF,"JBtF");
  if (load_coefs( Bt_filename,CoefFileFormat_sincos,nn,mm,JBtF,false))
    return 5;

  cout <<endl<< "$ Loading BTOTAL_phi fourier coefficients from " << Bp_filename << endl;
  Vector<complex<double> > JBpF(NF,"JBpF");
  if (load_coefs( Bp_filename,CoefFileFormat_sincos,nn,mm,JBpF,false))
    return 6;


  const Vector<unsigned int> tmp = findtrue((nn==0)&&(mm==0));
  const unsigned int ind00 = tmp[0];

  Vector<complex<double> > lambdaF(NF,"lambdaF");

  for(unsigned int k = 0; k<NF; k++) {
     const std::complex<double> i =  std::complex<double>(0,1);
     if ((mm[k]==0) && (nn[k]==0)) {
	// lambda00 is arbitrary constant
	// use this to save JBt00 and JBp00 so that they can
	// be used if needed.  recall that iota=JBt00/JBp00
	// note that in the file as cos/sin constants, the values will be
	// JBtF[ind00]/(2*pi) JBpF[ind00]/(2*pi)
	lambdaF[k] = 2*PI*complex<double>(JBtF[ind00].real(),JBpF[ind00].real());
     } else if (mm[k]==0) {
	lambdaF[k] = 2*PI*i*JBtF[k]/(nn[k]*JBpF[ind00].real());
     } else if (nn[k]==0) {
	lambdaF[k] = -2*PI*i*JBpF[k]/(mm[k]*JBpF[ind00].real());
     } else {
	const std::complex<double> tmp1 =  2*PI*i*JBtF[k]/(nn[k]*JBpF[ind00].real());
	const std::complex<double> tmp2 = -2*PI*i*JBpF[k]/(mm[k]*JBpF[ind00].real());
	lambdaF[k] = 0.5*(tmp1+tmp2);
     }
  }
   
  massage(lambdaF,1e-10);
  save_coefs("lambda.out",CoefFileFormat_sincos,nn,mm,lambdaF);


  //calcuate the divergence
  Vector<complex<double> > JdivBF(NF,"JdivBF");
  for(unsigned int k = 0; k<NF; k++) {
     JdivBF[k] = complex<double>(0,1)*(mm[k]*JBtF[k] + nn[k]*JBpF[k]);
  }

  Vector<double> JdivB(Npts,"JdivB");
  ifft2d(JdivB,JdivBF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,2*PI);

  Vector<double> divB(Npts,"divB");
  divB = JdivB / J;

  Vector<complex<double> > divBF(NF,"divBF");
  fft2d(divB,divBF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,2*PI);

  massage(divBF,1e-10);
  save_coefs("divBF.b2lambda.out",CoefFileFormat_sincos,nn,mm,divBF);




  double JdivBF_RMS = 0 ;
  double divBF_RMS = 0 ;
  for (unsigned int k =0; k<NF; k++) {
    double temp;
    temp = abs(JdivBF[k]);
    JdivBF_RMS += temp*temp;

    temp = abs(divBF[k]);
    divBF_RMS += temp*temp;
  }
  JdivBF_RMS = sqrt(JdivBF_RMS)/(2*PI);
  divBF_RMS = sqrt(divBF_RMS)/(2*PI);
  double JdivB_RMS = 0 ;
  double divB_RMS = 0 ;

  for (unsigned int j =0; j<Npts; j++) {
    double temp;
    temp = abs(JdivB[j]);
    JdivB_RMS += temp*temp;

    temp = abs(divB[j]);
    divB_RMS += temp*temp;

  }
  JdivB_RMS = sqrt(JdivB_RMS/Npts);
  divB_RMS = sqrt(divB_RMS/Npts);
  dispcr(JdivB_RMS);
  dispcr(JdivBF_RMS);cr();
  dispcr(divB_RMS);
  dispcr(divBF_RMS);cr();

  complex<double> iota = JBtF[ind00]/JBpF[ind00];
  dispcr(iota);


  return 0;
} // main()





