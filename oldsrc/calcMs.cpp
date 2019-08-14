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
#include "inductancematrix.hpp"
#include "coilfft.hpp"


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
   //   enable_option(opt_Btf);
   //   enable_option(opt_Bpf);
   enable_option(opt_Nharm);
   enable_option(opt_Mharm);
   enable_option(opt_Mf);

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


   //////////////////////////////////////

   //   Vector<complex<double> > v(1);
   //    "{(5,2)}">>v;
   //    dispcr(v);

   //    Matrix<complex<double> > A(2,2);
   //    "{{(1,2), (3,4)}, { (5,-6), (7,-8)}}" >>A;
   //    //   "{{1,2},{3,4}}">>A;
   //    dispcr(A);
   
   //    Matrix<complex<double> > Ainv(2,2);
   //    matricks_lapack::inv(A,Ainv);

   //    dispcr(A);
   //    dispcr(Ainv);

   //////////////////////////////////////

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

  ostringstream strmtmp3;
  strmtmp3 <<Ntheta;
  string Ntheta_str(strmtmp3.str());
  ostringstream strmtmp4;
  strmtmp4 <<Nphi;
  string Nphi_str(strmtmp4.str());


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


#ifdef DFTM
   cout << endl;
   cout<<"$ Generate reduced orthonormal series matrix ("<<Npts<<" x "<<NFR<<")"<<endl;

   STARTTIME(tbuff,ckstart);

   Matrix<complex<double> > fsR(Npts,NFR,"fsR");
   fseries(nnR,mmR,thetas,phis,fsR);

   STOPTIME(tbuff,ckstart);
#endif DFTM



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


 
  

   // lay coil surface onto grid 
  
   Vector<p3vector<double> > Xcoil(Npts, "Xcoil");
   Vector<p3vector<double> > dA_dtdp_coil(Npts, "dA_dtdp_coil");

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



   cout << endl;
   cout<<"$ Allocate memory for inductance matrix ("<<Npts<<" x "<<Npts<<")"<<endl;

   STARTTIME(tbuff,ckstart);
   
   Matrix<double> M(Npts, Npts, "M");

   STOPTIME(tbuff,ckstart);


  if (!M_filename.empty()) {
     cout << endl;
     cout<<"$ Load Inductance Matrix ("<<Npts<<"x"<<Npts<<")"<<endl;

     M.perline(1);
     M.textformat(text_nobraces);
     fname = M_filename + ".out";
     load(M,fname);

  } else {
     // Create mutual inductance matrix

     cout << endl;
     cout<<"$ Generating inductance matrix ("<<Npts<<" x "<<Npts<<")"<<endl;

     STARTTIME(tbuff,ckstart);

     inductancematrix(X,dA_dtdp,Xcoil,dA_dtdp_coil,M);

     STOPTIME(tbuff,ckstart);

     //SAVE Inductance MATRIX
//      printcr("$ Saving Inductance Matrix");
//      fname = "M.Ntheta="+Ntheta_str+".Nphi="+Nphi_str+".out";
//      dispcr(fname);
//      M.perline(1);
//      M.textformat(text_nobraces);
//      save(M,fname);

  }


   // Calculate M matrix with DFT

#ifdef DFTM

   // DFT  M matrix
   cout << endl;
   cout<<"$ DFT of inductance matrix ("<<NFR<<" x "<<NFR<<")"<<endl;

   STARTTIME(tbuff,ckstart);

   Matrix<complex<double> > MF1(NFR,NFR,"MF1"); 
   MF1 = (adj(fsR)|M|fsR)*Csqr;


   STOPTIME(tbuff,ckstart);
   
   //save MF1
   data.resize() = real(MF1);
   fname = "MF.fsR_M_fsR.Nn="+Nnn_str+".Nm="+Nmm_str+".R.out";
   dispcr(fname);
   save(data,fname);
   data.resize() = imag(MF1);
   fname = "MF.fsR_M_fsR.Nn="+Nnn_str+".Nm="+Nmm_str+".I.out";
   dispcr(fname);
   save(data,fname);

#endif DFTM



   // Calculate M matrix with FFT

   // fft  M matrix
   cout << endl;
   cout<<"$ FFT of inductance matrix ("<<NFR<<" x "<<NFR<<")"<<endl;

   STARTTIME(tbuff,ckstart);

   Matrix<complex<double> > MF2(NFR,NFR,"MF2"); 

   fft_of_M(M, MF2, Nphi, Ntheta,Nnn,Nmm,  Nharm,Mharm,1e-12);

   STOPTIME(tbuff,ckstart);
   
   //save MF2
   data.resize() = real(MF2);
   fname = "MF.fft_of_M.Nn="+Nnn_str+".Nm="+Nmm_str+".R.out";
   dispcr(fname);
   save(data,fname);
   data.resize() = imag(MF2);
   fname = "MF.fft_of_M.Nn="+Nnn_str+".Nm="+Nmm_str+".I.out";
   dispcr(fname);
   save(data,fname);




   return 0;
} // main()





