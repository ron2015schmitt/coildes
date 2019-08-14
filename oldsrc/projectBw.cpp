

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
   enable_option(opt_ff);
   enable_option(opt_Nphi);
   enable_option(opt_Ntheta);
   enable_option(opt_Nnn);
   enable_option(opt_Nmm);
   enable_option(opt_Bpf);
   enable_option(opt_Nharm);
   enable_option(opt_Mharm);
   enable_option(opt_omegaf);

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
   if ( (Nharm >1) ||(Mharm>1) )
      modevectors(NF,nn,mm,Nnn,Nmm,Nharm,Mharm,mode00);
   else
      modevectors(NF,nn,mm,Nnn,Nmm,mode00);



   // these exclude the n=0,m=0 case
   Vector<double> nnR("nnR");
   Vector<double> mmR("mmR");
   unsigned int NFR;
   mode00 = false;
   if ( (Nharm >1) ||(Mharm>1) )
      modevectors(NFR,nnR,mmR,Nnn,Nmm,Nharm,Mharm,mode00);
   else
      modevectors(NFR,nnR,mmR,Nnn,Nmm,mode00);


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



   Vector<complex<double> > BrF(NFR,"BrF");

   cout <<endl<< "$ Loading Plasma Br sin/cos fourier coefficients from " << flux_filename << endl;
   if (load_coefs(flux_filename,CoefFileFormat_sincos,nnR,mmR,BrF))
      return 3;

   Vector<double> Br(Npts, "Br");
   expandfunction(Br,BrF,fsR);

   // load or calculate B field
   cout << endl;
   cout<<"$ Load tangent B field "<<endl;

   STARTTIME(tbuff,ckstart);

   cout <<endl<< "$ Loading BTOTAL_phi fourier coefficients from " << Bp_filename << endl;

   Vector<complex<double> > BpF(NF,"BpF");
   if (load_coefs( Bp_filename,CoefFileFormat_sincos,nn,mm,BpF))
      return 6;

   Vector<double> Bp(Npts, "Bp");
   expandfunction(Bp,BpF,fs);

   STOPTIME(tbuff,ckstart);


   Vector<double> Bw(Npts, "Bw");
   Bw = Br / Bp;

   Vector<complex<double> > BwF(NFR,"BwF");
   transformfunction(Bw,BwF,fsR);

  datavec = real(BwF);
  fname = "BwF.Nn="+Nnn_str+".Nm="+Nmm_str+".R.out";
  dispcr(fname);
  save(datavec,fname);


  datavec = imag(BwF);
  fname = "BwF.Nn="+Nnn_str+".Nm="+Nmm_str+".I.out";
  dispcr(fname);
  save(datavec,fname);


   Matrix<complex<double> > P(NFR,NFR,"P");
   
   cout<<"$ Load Omega matrix ("<<NFR<<"x"<<NFR<<")"<<endl;

   STARTTIME(tbuff,ckstart);
   data.perline(1);
   data.textformat(text_nobraces);
   data.resize(NFR,NFR);
   fname = omega_filename + ".PR.out";
   load(data,fname);

   Matrix <double> data2("data2");
   data2.perline(1);
   data2.textformat(text_nobraces);
   data2.resize(NFR,NFR);
   fname = omega_filename + ".PI.out";
   load(data2,fname);
   P = mcomplex(data,data2);
   data.clear();
   data2.clear();

   STOPTIME(tbuff,ckstart);


   Matrix<complex<double> > Pinv(NFR,NFR,"Pinv");

   matricks_lapack::inv(P,Pinv);

   Vector<complex<double> > Bw_projF(NFR,"Bw_projF");
   Bw_projF = (Pinv|BwF);


  datavec = real(Bw_projF);
  fname = "projectBw.Nn="+Nnn_str+".Nm="+Nmm_str+".R.out";
  dispcr(fname);
  save(datavec,fname);


  datavec = imag(Bw_projF);
  fname = "projectBw.Nn="+Nnn_str+".Nm="+Nmm_str+".I.out";
  dispcr(fname);
  save(datavec,fname);


  Vector<unsigned int> ii;
  Vector<double> magn = abs(Bw_projF);
  ii = sortrevwind(magn);


   fname = "projectBw.isort.out";
   dispcr(fname);
   ii.perline(1);
   ii.textformat(text_nobraces);
   save(ii,fname);

  Vector<double> mmRsorted(NFR,"mmRsorted");
  Vector<double> nnRsorted(NFR,"nnRsorted");

  Bw_projF = Bw_projF[ii];

  fname = omega_filename + ".mmRsort.out";
  mmRsorted.textformat(text_nobraces);
  mmRsorted.perline(1);
  load(mmRsorted,fname);
  mmRsorted = mmRsorted[ii];

  fname = omega_filename + ".nnRsort.out";
  nnRsorted.textformat(text_nobraces);
  nnRsorted.perline(1);
  load(nnRsorted,fname);
  nnRsorted = nnRsorted[ii];

  datavec = real(Bw_projF);
  fname = "projectBw.Nn="+Nnn_str+".Nm="+Nmm_str+".Rsorted.out";
  dispcr(fname);
  save(datavec,fname);


  datavec = imag(Bw_projF);
  fname = "projectBw.Nn="+Nnn_str+".Nm="+Nmm_str+".Isorted.out";
  dispcr(fname);
  save(datavec,fname);

   fname = "projectBw.Nn="+Nnn_str+".Nm="+Nmm_str+".mmRsort.out";
   dispcr(fname);
   mmRsorted.perline(1);
   mmRsorted.textformat(text_nobraces);
   save(mmRsorted,fname);
   fname = "projectBw.Nn="+Nnn_str+".Nm="+Nmm_str+".nnRsort.out";
   dispcr(fname);
   nnRsorted.perline(1);
   nnRsorted.textformat(text_nobraces);
   save(nnRsorted,fname);

   Vector<double> EV(NFR,"EV");
   fname = "omegaNormd.Nn="+Nnn_str+".Nm="+Nmm_str+".EV.out";
   dispcr(fname);
   EV.perline(1);
   EV.textformat(text_nobraces);
   load(EV,fname);
   
   EV=1/EV[ii];
   

  datavec = EV;
  fname = "projectBw.Nn="+Nnn_str+".Nm="+Nmm_str+".EVinv.sorted.out";
  dispcr(fname);
  save(datavec,fname);


  Bw_projF = Bw_projF * EV;

  magn = abs(Bw_projF);
  ii = sortrevwind(magn);

  Bw_projF = Bw_projF[ii];
  mmRsorted = mmRsorted[ii];
  nnRsorted = nnRsorted[ii];


  datavec = real(Bw_projF);
  fname = "projectBw.Nn="+Nnn_str+".Nm="+Nmm_str+".BwEV.Rsorted.out";
  dispcr(fname);
  save(datavec,fname);
  datavec = imag(Bw_projF);
  fname = "projectBw.Nn="+Nnn_str+".Nm="+Nmm_str+".BwEV.Isorted.out";
  dispcr(fname);
  save(datavec,fname);

   fname = "projectBw.Nn="+Nnn_str+".Nm="+Nmm_str+".BwEV.mmRsort.out";
   dispcr(fname);
   mmRsorted.perline(1);
   mmRsorted.textformat(text_nobraces);
   save(mmRsorted,fname);
   fname = "projectBw.Nn="+Nnn_str+".Nm="+Nmm_str+".BwEV.nnRsort.out";
   dispcr(fname);
   nnRsorted.perline(1);
   nnRsorted.textformat(text_nobraces);
   save(nnRsorted,fname);

   return 0;
} // main()





