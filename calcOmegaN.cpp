/************************************************************************* 
 * 
 *   File Name    :  
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  
 * 
 *
 * VERSION NOTES:
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
#include "omegamatrix.hpp"
#include "gradfvector.hpp"


// Main Function for code

int main (int argc, char *argv[])
{

   disable_all_options();
   enable_option(opt_pf);
   enable_option(opt_cf);
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


   //////////////////////////////////////

 //   LAvector<complex<double> > v(1);
//    "{(5,2)}">>v;
//    dispcr(v);

//    Matrix<complex<double> > A(2,2);
//    "{{(1,2), (3,4)}, { (5,-6), (7,-8)}}" >>A;
//    //   "{{1,2},{3,4}}">>A;
//    dispcr(A);
   
//    Matrix<complex<double> > Ainv(2,2);
//    cooll_lapack::inv(A,Ainv);

//    dispcr(A);
//    dispcr(Ainv);

   //////////////////////////////////////

   // Create angle grid
   const unsigned int Npts = Ntheta*Nphi;
   LAvector<double> thetas(Npts,"thetas");
   LAvector<double> phis(Npts,"phis");
   anglevectors(thetas, phis, Ntheta, Nphi);
  
   ostringstream strmtmp1;
   strmtmp1 <<Nnn;
   string Nnn_str(strmtmp1.str());
   ostringstream strmtmp2;
   strmtmp2 <<Nmm;
   string Nmm_str(strmtmp2.str());


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




   // load B field
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

   LAvector<double> Bt(Npts, "Bt");
   expandfunction(Bt,BtF,fs);
   LAvector<double> Bp(Npts, "Bp");
   expandfunction(Bp,BpF,fs);

   for (unsigned int j =0; j<Npts; j++)
      B[j] =  Bt[j] * dx_dtheta[j] + Bp[j] * dx_dphi[j];

   STOPTIME(tbuff,ckstart);



   cr();
   printcr("Find contravariant * J vector components.");

   STARTTIME(tbuff,ckstart);
   LAvector<double> JBt(Npts,"JBt");
   LAvector<double> JBp(Npts,"JBp");
 
   for (unsigned int j =0; j<Npts; j++) {
      JBt[j] = J[j] * Bt[j];
      JBp[j] = J[j] * Bp[j];
   }
   STOPTIME(tbuff,ckstart);


 
   // find FFT (only n=0,m=0) of J*Bt
   cout << endl;
   cout<<"$ find FFT (only n=0,m=0) of J*Bt"<<endl;

   STARTTIME(tbuff,ckstart);
   double JBt00C;
   for (unsigned int j =0; j<Npts; j++) {
      JBt00C = JBt00C + JBt[j];
   }
   double JBt00 = JBt00C /double(Npts);
   STOPTIME(tbuff,ckstart);


   // find FFT (only n=0,m=0) of J*Bp
   cout << endl;
   cout<<"$ find FFT (only n=0,m=0) of J*Bp"<<endl;

   STARTTIME(tbuff,ckstart);
   double JBp00C;
   for (unsigned int j =0; j<Npts; j++) {
      JBp00C = JBp00C + JBp[j];
   }
   double JBp00 = JBp00C /double(Npts);
   STOPTIME(tbuff,ckstart);

   cout << endl;
   cout<<"$ calculate iota"<<endl;

   double iota = JBt00/JBp00;

   dispcr(JBt00);
   dispcr(JBp00);
   dispcr(iota);  

   cout << endl;
   cout<<"$ calculate expected eigen values"<<endl;
   LAvector<double> EV(NFR,"EV");
   EV = abs(nnR + iota*mmR);
   //dispcr(EV);

   cout << endl;
   cout<<"$ sort expected eigen values"<<endl;
   LAvector<unsigned int> ii;
   ii = sortwind(EV);
   LAvector<double> mmRsorted(NFR,"mmRsorted");
   LAvector<double> nnRsorted(NFR,"nnRsorted");
   mmRsorted = mmR[ii];
   nnRsorted = nnR[ii];

   //   dispcr(EV);
   //dispcr(mmRsorted);
   //dispcr(nnRsorted);


   fname = "omegaNormd.Nn="+Nnn_str+".Nm="+Nmm_str+".EV.out";
   dispcr(fname);
   EV.perline(1);
   EV.textformat(text_nobraces);
   save(EV,fname);
   fname = "omegaNormd.Nn="+Nnn_str+".Nm="+Nmm_str+".mmRsort.out";
   dispcr(fname);
   mmRsorted.perline(1);
   mmRsorted.textformat(text_nobraces);
   save(mmRsorted,fname);
   fname = "omegaNormd.Nn="+Nnn_str+".Nm="+Nmm_str+".nnRsort.out";
   dispcr(fname);
   nnRsorted.perline(1);
   nnRsorted.textformat(text_nobraces);
   save(nnRsorted,fname);




   Matrix<complex<double> > V(NFR,NFR,"V");
   LAvector<complex<double> >  W(NFR,"W");

   // Calculate Omega matrix
  
   cout << endl;
   cout<<"$ Calculate Omega matrix ("<<NFR<<"x"<<NFR<<")"<<endl;
     
   STARTTIME(tbuff,ckstart);
     
   Matrix<complex<double> > OmegaN(NFR,NFR,"OmegaN");
   Matrix<p3vector<complex<double> > >  grad_fsR(Npts,NFR,"grad_fsR");
   gradfvector(thetas,phis,mmR,nnR,grad_theta,grad_phi,grad_fsR);

   NormalizedOmegamatrix(OmegaN, B, fsR, grad_fsR, Bp);

   
   grad_fsR.clear();

   STOPTIME(tbuff,ckstart);

   //save omega
   data.resize() = real(OmegaN);
   fname = "omegaNormd.Nn="+Nnn_str+".Nm="+Nmm_str+".omegaR.out";
   dispcr(fname);
   save(data,fname);
   data.resize() = imag(OmegaN);
   fname = "omegaNormd.Nn="+Nnn_str+".Nm="+Nmm_str+".omegaI.out";
   dispcr(fname);
   save(data,fname);

 
   // Spectrum of Omega Matrix
     
   cout << endl;
   cout<<"$ eigen-decomposition of Omega matrix ("<<NFR<<"x"<<NFR<<")"<<endl;
     
   STARTTIME(tbuff,ckstart);

   cooll_lapack::eigwincsort(OmegaN,V,W);
   OmegaN.clear();
     
   STOPTIME(tbuff,ckstart);



   // save Normalized Omega in eigen form
   printcr("$ Saving Nromalized Omega Matrix in eigen form");
   data.resize() = real(V);
   fname = "omegaNormd.Nn="+Nnn_str+".Nm="+Nmm_str+".VR.out";
   dispcr(fname);
   save(data,fname);
   data.resize() = imag(V);
   fname = "omegaNormd.Nn="+Nnn_str+".Nm="+Nmm_str+".VI.out";
   dispcr(fname);
   save(data,fname);

   datavec.resize() = real(W);
   datavec.perline(1);
   fname = "omegaNormd.Nn="+Nnn_str+".Nm="+Nmm_str+".WR.out";
   dispcr(fname);
   save(datavec,fname);
   datavec.resize() = imag(W);
   datavec.perline(1);
   fname = "omegaNormd.Nn="+Nnn_str+".Nm="+Nmm_str+".WI.out";
   dispcr(fname);
   save(datavec,fname);






   return 0;
} // main()





