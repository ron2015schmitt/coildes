/************************************************************************* 
 * 
 *   File Name    :  
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  16 May 2003
 * 
 *
 *   SYNOPSIS     
 *     Calculate B field due to teh coils, on teh given plasma surface
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

// generic file for generating the coil+plasma bfield
#include "bfield_coils.hpp"



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
  enable_option(opt_Itoroidal);
  enable_option(opt_Ipoloidal);
  enable_option(opt_Nharm);
  enable_option(opt_Mharm);


   enable_option(opt_NphiI);
   enable_option(opt_NthetaI);

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

   const unsigned int NptsI = NthetaI*NphiI;
   Vector<double> thetasI(NptsI,"thetasI");
   Vector<double> phisI(NptsI,"phisI");
   anglevectors(thetasI, phisI, NthetaI, NphiI);



  // Create Fourier Mode vectors
  Vector<double> nn("nn");
  Vector<double> mm("mm");
  unsigned int NF;
  bool mode00 = true;
  if ( (Nharm >1) ||(Mharm>1) )
     modevectors(NF,nn,mm,Nnn,Nmm,Nharm,Mharm,mode00);
  else
     modevectors(NF,nn,mm,Nnn,Nmm,mode00);


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


  dispcr(current_filename);


  setupcoils(coil_filename, current_filename, Itoroidal, Ipoloidal, 
	     Nnn, Nmm, NphiI, NthetaI, Nharm, Mharm);


  dispcr(Itoroidal);
  dispcr(Ipoloidal);

  // calculate magnetic field on plasma surface

    
  cout << endl;
  cout<<"$ Calc B field "<<endl;
 
  STARTTIME(tbuff,ckstart);

  Vector<p3vector<double> > Bcoils(Npts, "Bcoils");
  Vector<p3vector<double> > Btotal(Npts, "Btotal");
  unsigned int count = 0;
  unsigned int checkup = static_cast<unsigned int>(Npts*0.01);
  for (unsigned int j =0; j<Npts; j++) {
    if (++count == checkup) {
      print(Matricks::round(double(j)/Npts*100));cout <<" %"<<endl;
      count =0;
    }
    bTotalandbCoils(X[j],Btotal[j], Bcoils[j]);
  }

  STOPTIME(tbuff,ckstart);





  Vector<double> Br(Npts, "Bn");
  Vector<double> Bt(Npts, "Bt");
  Vector<double> Bp(Npts, "Bp");

  Vector<double> J(Npts, "J");

  Vector<double> JBr(Npts, "JBr");
  Vector<double> JBt(Npts, "JBt");
  Vector<double> JBp(Npts, "JBp");

  Vector<double> Flux(Npts, "Flux");
  Vector<double> Br_coils(Npts, "Br_coils");
  Vector<double> Flux_coils(Npts, "Flux_coils");
 
  for (unsigned int j =0; j<Npts; j++) {
     // calculate contravariant components
    Br[j] = dot(Btotal[j], grad_r[j]);
    Bt[j] =  dot ( Btotal[j], grad_theta[j] );
    Bp[j] =  dot ( Btotal[j], grad_phi[j] );

    J[j] = dot(dx_dr[j],cross(dx_dtheta[j], dx_dphi[j]));

    JBr[j] =  J[j]*Br[j];
    JBt[j] =  J[j]*Bt[j];
    JBp[j] =  J[j]*Bp[j];

    Flux[j] = dot(Btotal[j], dAdtdp[j]);
    Br_coils[j] = dot(Bcoils[j], grad_r[j]);
    Flux_coils[j] = dot(Bcoils[j], dAdtdp[j]);

  }
  


  // Create fourier series


  cout << endl;
  cout<<"$ Generate fourier series matrix ("<<Npts<<" x "<<NF<<")"<<endl;

  STARTTIME(tbuff,ckstart);

  Matrix<complex<double> > fs(Npts,NF,"fs");
  fseries(nn,mm,thetas,phis,fs);

  STOPTIME(tbuff,ckstart);

  cout << endl;
  cout<<"$ save data in various formats" <<endl;


  Vector<complex<double> > FluxF(NF,"FluxF");
  transformfunction(Flux,FluxF,fs);
  massage(FluxF,1e-10);
  save_coefs("FluxF.out",CoefFileFormat_sincos,nn,mm,FluxF);

  Flux.perline(1);
  Flux.textformat(text_nobraces);
  save(Flux,"Flux.out");


  Vector<complex<double> > Br_coilsF(NF,"Br_coilsF");
  transformfunction(Br_coils,Br_coilsF,fs);
  massage(Br_coilsF,1e-10);
  save_coefs("Br_coilsF.out",CoefFileFormat_sincos,nn,mm,Br_coilsF);

  Vector<complex<double> > Flux_coilsF(NF,"Flux_coilsF");
  transformfunction(Flux_coils,Flux_coilsF,fs);
  massage(Flux_coilsF,1e-10);
  save_coefs("Flux_coilsF.out",CoefFileFormat_sincos,nn,mm,Flux_coilsF);


  Vector<complex<double> > BrF(NF,"BrF");
  transformfunction(Br,BrF,fs);
  massage(BrF,1e-10);
  save_coefs("BrTOTALF.fromcoils.out",CoefFileFormat_sincos,nn,mm,BrF);

  Vector<complex<double> > BtF(NF,"BtF");
  transformfunction(Bt,BtF,fs);
  massage(BtF,1e-10);
  save_coefs("BtTOTALF.fromcoils.out",CoefFileFormat_sincos,nn,mm,BtF);

  Vector<complex<double> > BpF(NF,"BpF");
  transformfunction(Bp,BpF,fs);
  massage(BpF,1e-10);
  save_coefs("BpTOTALF.fromcoils.out",CoefFileFormat_sincos,nn,mm,BpF);



  Vector<complex<double> > JBrF(NF,"JBrF");
  transformfunction(JBr,JBrF,fs);
  massage(JBrF,1e-10);
  save_coefs("JBrTOTALF.fromcoils.out",CoefFileFormat_sincos,nn,mm,JBrF);

  Vector<complex<double> > JBtF(NF,"JBtF");
  transformfunction(JBt,JBtF,fs);
  massage(JBtF,1e-10);
  save_coefs("JBtTOTALF.fromcoils.out",CoefFileFormat_sincos,nn,mm,JBtF);

  Vector<complex<double> > JBpF(NF,"JBpF");
  transformfunction(JBp,JBpF,fs);
  massage(JBpF,1e-10);
  save_coefs("JBpTOTALF.fromcoils.out",CoefFileFormat_sincos,nn,mm,JBpF);


  Vector<complex<double> > JF(NF,"JF");
  transformfunction(J,JF,fs);
  massage(JF,1e-10);
  save_coefs("JF.out",CoefFileFormat_sincos,nn,mm,JF);



  Vector<complex<double> > JdivBF(NF,"JdivBF");
  for (unsigned int k =0; k<NF; k++)
     JdivBF[k] = complex<double>(0.0,1.0)*(mm[k]*JBtF[k] + nn[k]*JBpF[k]);
  massage(JdivBF,1e-10);
  save_coefs("JdivBTOTALF.fromcoils.out",CoefFileFormat_sincos,nn,mm,JdivBF);

  //  for (unsigned int k =0; k<NF; k++) {
  // if (abs(nn(k)) <=4.00000001){
  //disp(mm[k]); disp(nn[k]); disp(mm[k]*JBtF[k]); disp(nn[k]*JBpF[k]);cr();
  // }
     //     disp(mm[k]); disp(nn[k]); disp(JF[k]); disp(JF[k]);cr();
  //}



//   Bn_plasma.textformat(text_nobraces);
//   Bn_plasma.perline(1);
//   Bn_plasma =- Bn_plasma;
//   save(Bn_plasma,"Flux_negplasma.out");
//   Bn_plasma =- Bn_plasma;



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





  ///////////////////////////////////////////////////////////
//   JdivBF.perline(1);
//   JdivBF.textformat(text_nobraces); 
//   save(JdivBF,"JdivBF.out");
  ///////////////////////////////////////////////////////////

  Vector<double> JdivB(Npts, "JdivB");
  expandfunction(JdivB,JdivBF,fs);

  Vector<double> divB(Npts, "divB");
  for (unsigned int j =0; j<Npts; j++)
    divB[j] =  JdivB[j]/J[j];
  Vector<complex<double> > divBF(NF,"divBF");
  transformfunction(divB,divBF,fs);
  massage(divBF,1e-10);
  save_coefs("divBTOTALF.fromcoils.out",CoefFileFormat_sincos,nn,mm,divBF);




  double FluxTotal = 0;
  double absFluxTotal = 0;
  for (unsigned int j =0; j<Npts; j++) {
    FluxTotal += JBr[j];
    absFluxTotal += abs(JBr[j]);
  }
   // C = dtheta*dphi
   //   = (2*pi/Ntheta)*(2*pi/Nphi)
   //   = (2*pi*2*pi/Npts)
  const double C = (2*PI*2*PI/double(Npts));
  FluxTotal = FluxTotal * C;
  absFluxTotal = absFluxTotal * C;


  Vector<unsigned int> temp = findtrue((nn==0)&&(mm==0));
print("index[n=0,m=0] = ");dispcr(temp);
  double JBTheta00 = JBtF[temp[0]].real()/(2*PI);
  double JBPhi00 = JBpF[temp[0]].real()/(2*PI);
  double iota = JBTheta00/JBPhi00;


  dispcr(FluxTotal);
  dispcr(absFluxTotal);
dispcr(JBTheta00);
dispcr(JBPhi00);
dispcr(iota);

 

  // calculate RMS values and display
  double BrF_RMS = 0 ;
  double BtF_RMS = 0 ;
  double BpF_RMS = 0 ;
  double Br_coilsF_RMS = 0 ;
  double Flux_coilsF_RMS = 0 ;
  double JF_RMS = 0 ;
  double JdivBF_RMS = 0 ;
  double divBF_RMS = 0 ;
  for (unsigned int k =0; k<NF; k++) {
    double temp;
    
    temp = abs(BrF[k]);
    BrF_RMS += temp*temp;

    temp = abs(BtF[k]);
    BtF_RMS += temp*temp;

    temp = abs(BpF[k]);
    BpF_RMS += temp*temp;


    temp = abs(Br_coilsF[k]);
    Br_coilsF_RMS += temp*temp;

    temp = abs(Flux_coilsF[k]);
    Flux_coilsF_RMS += temp*temp;

    temp = abs(JF[k]);
    JF_RMS += temp*temp;

    temp = abs(JdivBF[k]);
    JdivBF_RMS += temp*temp;

    temp = abs(divBF[k]);
    divBF_RMS += temp*temp;


  }
  BrF_RMS = sqrt(BrF_RMS)/(2*PI);
  BtF_RMS = sqrt(BtF_RMS)/(2*PI);
  BpF_RMS = sqrt(BpF_RMS)/(2*PI);
  Br_coilsF_RMS = sqrt(Br_coilsF_RMS)/(2*PI);
  Flux_coilsF_RMS = sqrt(Flux_coilsF_RMS)/(2*PI);

  JF_RMS = sqrt(JF_RMS)/(2*PI);
  JdivBF_RMS = sqrt(JdivBF_RMS)/(2*PI);
  divBF_RMS = sqrt(divBF_RMS)/(2*PI);

  double B_RMS = 0 ;
  double Br_RMS = 0 ;
  double Bt_RMS = 0 ;
  double Bp_RMS = 0 ;
  double Br_coils_RMS = 0 ;
  double Flux_coils_RMS = 0 ;
  double J_RMS = 0 ;
  double JdivB_RMS = 0 ;
  double divB_RMS = 0 ;

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

    Br_coils_RMS += Br_coils[j]*Br_coils[j];
    Flux_coils_RMS += Flux_coils[j]*Flux_coils[j];

    temp = abs(J[j]);
    J_RMS += temp*temp;

    temp = abs(JdivB[j]);
    JdivB_RMS += temp*temp;

    temp = abs(divB[j]);
    divB_RMS += temp*temp;

  }
  B_RMS = sqrt(B_RMS/Npts);
  Br_RMS = sqrt(Br_RMS/Npts);
  Bt_RMS = sqrt(Bt_RMS/Npts);
  Bp_RMS = sqrt(Bp_RMS/Npts);
  Br_coils_RMS = sqrt(Br_coils_RMS/Npts);
  Flux_coils_RMS = sqrt(Flux_coils_RMS/Npts);

  J_RMS = sqrt(J_RMS/Npts);
  JdivB_RMS = sqrt(JdivB_RMS/Npts);
  divB_RMS = sqrt(divB_RMS/Npts);

  dispcr(B_RMS);cr();
  dispcr(Br_RMS);
  dispcr(BrF_RMS);cr();
  dispcr(Bt_RMS);
  dispcr(BtF_RMS);cr();
  dispcr(Bp_RMS);
  dispcr(BpF_RMS);cr();
  dispcr(Br_coils_RMS);
  dispcr(Br_coilsF_RMS);cr();
  dispcr(Flux_coils_RMS);
  dispcr(Flux_coilsF_RMS);cr();

  dispcr(J_RMS);
  dispcr(JF_RMS);cr();
  dispcr(JdivB_RMS);
  dispcr(JdivBF_RMS);cr();
  dispcr(divB_RMS);
  dispcr(divBF_RMS);cr();

  save(B_RMS,"BTOTAL.RMS.out");

  save(Br_RMS,"BrTOTAL.RMS.out");

  save(iota,"iota.out");

  return 0;
} // main()





