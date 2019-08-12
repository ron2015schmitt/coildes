/************************************************************************* 
 * 
 *   File Name    :  findpertsurface2_fft.cpp
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  
 * 
 *
 *   SYNOPSIS     
 *
 *  This takes the differnece between the necessary fourier flux coef's 
 * and the acutal flux coef's and finds the location of the perturbed magnetic
 * surface.
 *    
 *
 * VERSION NOTES:
 * Uses magnetic coordinates to calculate Omega
 *
 * created from findpertsurface_fft.cpp on 05.12.2007
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
#include "omegamatrix_fft.hpp"
#include "omegasubspace.hpp"
#include "coilfft.hpp"



// This is the file that defines the B field configuration


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

   enable_option(opt_lambdaf);
   enable_option(opt_iota);
   enable_option(opt_fluxshear);


   // parse command line input
 
   if (!parse_cmd(argc, argv))
      return 1;

   printcr("----------------------------------------------------");
   for(int i=0; i<argc; i++){
      print(argv[i]);
      print(" ");
   }
   cr();
   printcr("----------------------------------------------------");


   // 
   //  ios_base::fmtflags flags = ios_base::right | ios_base::scientific;
   string fext = "txt";
   string ftemp;
   p3vectorformat::textformat(text_nobraces);

   Vector <double> datavec("datavec");
   datavec.perline(1);
   datavec.textformat(text_nobraces);
   Matrix <double> data("data");
      Matrix <double> data2("data2");
   data2.perline(1);
   data2.textformat(text_nobraces);

   string fname;

   // variables for measuring times
   struct tms tbuff;
   clock_t ckstart;

   // display Matricks mode
   cout << endl;
   display_execution_mode();
   cout << endl;

   ostringstream strmtemp;

   string methstr(".pert2");
   string plasma_name;
   {
      unsigned int j = plasma_filename.rfind("_");
      plasma_name = plasma_filename.substr(0,j);
      plasma_name = "."+plasma_name;
   }

   strmtemp <<Ntheta;
   string Ntheta_str(strmtemp.str());
   strmtemp.str("");
   Ntheta_str = ".Ntheta="+Ntheta_str;

   strmtemp <<Nphi;
   string Nphi_str(strmtemp.str());
   strmtemp.str("");
   Nphi_str = ".Nphi="+Nphi_str;

   strmtemp <<Nnn;
   string Nnn_str(strmtemp.str());
   strmtemp.str("");
   Nnn_str = ".Nn="+Nnn_str;

   strmtemp <<Nmm;
   string Nmm_str(strmtemp.str());
   strmtemp.str("");
   Nmm_str = ".Nm="+Nmm_str;


   strmtemp <<iota_given;
   cout.precision(6);
   string iotaval_str(strmtemp.str());
   strmtemp.str("");
   cout.precision(15);
   iotaval_str = ".iota="+iotaval_str;

  

   dispcr(plasma_name);
   dispcr(methstr);
   dispcr(Nphi_str+Ntheta_str);
   dispcr(Nnn_str+Nmm_str);
   dispcr(iotaval_str);  
   dispcr(fluxshear_given);

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




   nn.perline(1);
   nn.textformat(text_nobraces);
   save(nn,"nn.out");
   mm.perline(1);
   mm.textformat(text_nobraces);
   save(mm,"mm.out");

   nnR.perline(1);
   nnR.textformat(text_nobraces);
   save(nnR,"nnR.out");
   mmR.perline(1);
   mmR.textformat(text_nobraces);
   save(mmR,"mmR.out");


   // load the plasma surface fourier coef's

   cout << "$ Loading PLASMA SURFACE fourier coefficients from " << plasma_filename << endl;
   FourierSurface plasmafourier;
   if (load_fourier_surface(plasma_filename,plasmafourier))
      return 1;
   plasmafourier.RF().name("p.RF");
   plasmafourier.ZF().name("p.ZF");

   // print coef's
   //  printfouriercoefs(plasmafourier.nn(),plasmafourier.mm(),plasmafourier.RF(),plasmafourier.ZF(),10,18);
  



   // load the plasma flux fourier coef's
   // WE CAN IGNORE THE (m=0,n=0) mode because the omega matrix is identically zero
   // for (m=0,n=0) (both in row and column index)
   Vector<complex<double> > delFluxF(NFR,"FluxF");
   cout <<endl<< "$ Loading perturbed Plasma Flux sin/cos fourier coefficients from " << flux_filename << endl;
   if (load_coefs(flux_filename,CoefFileFormat_sincos,nnR,mmR,delFluxF))
      return 3;

 
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




   double iota = iota_given;

   
   Vector<double>  W0(NFR,"W0");
   const double coef1 = fluxshear_given/ (2.0*PI);
   for (unsigned int j =0; j<NFR; j++) {
      W0[j] = (coef1)*(nnR[j] + iota*mmR[j]);
      //      dispcr(W0[j]);
   }

   // P converts from magnetic to geometric fourier
   // so rows are in (nnR,mmR) order and 
   // columns are in (nnR_maf,mmR_mag) order
   Matrix<complex<double> > P(NFR,NFR,"P");

   if ( omega_filename.empty() ) {

      Vector<complex<double> > lambdaF(NFR,"lambdaF");

      cout <<endl<< "$ Loading lambda sin/cos fourier coefficients from " << lambda_filename << endl;
      if (load_coefs(lambda_filename,CoefFileFormat_sincos,nnR,mmR,lambdaF,false)){
	 printcr("Above ERROR occurred in "+myname+".");
	 return 4;
      }

      cout << endl;
      cout<<"$ Calculate OmegaN matrix ("<<NFR<<"x"<<NFR<<")"<<endl;
      
      STARTTIME(tbuff,ckstart);

      Vector<double> lambda(Npts, "lambda");
      //      expandfunction(lambda,lambdaF,fsR);
      mode00=false;
      ifft2d(lambda,lambdaF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,1/(2*PI),mode00);

      OmegaN_PF_from_lambda_fft(P,lambda,nnR,mmR,thetas,phis,Nphi, Ntheta, Nnn, Nmm, Nharm, Mharm);          
      STOPTIME(tbuff,ckstart);STARTTIME(tbuff,ckstart);
  

   } else {

      cout << endl;
      cout<<"$ Load Omega eigenvector matrix("<<NFR<<"x"<<NFR<<")"<<endl;
      
      STARTTIME(tbuff,ckstart);
      
      data.resize(NFR,NFR);
      data.perline(1);
      data.textformat(text_nobraces);
      fname = "PF."+omega_filename +".R.out";
      dispcr(fname);
      Matricks::load(data,fname);

      data2.resize(NFR,NFR);
      data2.perline(1);
      data2.textformat(text_nobraces);
      fname = "PF."+omega_filename +".I.out";
      dispcr(fname);
      Matricks::load(data2,fname);

      P = mcomplex(data,data2);
      data.clear();
      data2.clear();

      STOPTIME(tbuff,ckstart);

   }


   // Invert Omega Eigenvector Matrix
  
   cout << endl;
   cout<<"$ Invert Omega Eigenvector Matrix ("<<NFR<<"x"<<NFR<<")"<<endl;
     
   STARTTIME(tbuff,ckstart);

   Matrix<complex<double> > Pinv(NFR,NFR,"Pinv");
   matricks_lapack::inv(P,Pinv);
     
   STOPTIME(tbuff,ckstart);



   // Calculate Omega Inv matrix

   cout << endl;
   cout<<"$ Calculate Omega Inv matrix ("<<NFR<<"x"<<NFR<<")"<<endl;
     
   STARTTIME(tbuff,ckstart);

   for(unsigned int r=0; r<NFR;r++) {
      double temp1 = 1.0/W0[r];
      complex<double> temp2 = std::complex<double>(0,-temp1);
      for(unsigned int c=0; c<NFR;c++) {
	 Pinv(r,c) = temp2*Pinv(r,c);
      }
   }
     
   Matrix<complex<double> > OmegaInv(NFR,NFR,"OmegaInv");
   OmegaInv = (P|Pinv);
   P.clear();
   Pinv.clear();

   STOPTIME(tbuff,ckstart);

   // calculate perturbations
  
   cout << endl;
   cout<<"$ Calculate perturbations to plasma surface"<<endl;

   STARTTIME(tbuff,ckstart);


   Vector<complex<double> > deltaF(NFR,"deltaF");
   deltaF = (OmegaInv|(delFluxF));

   Vector<double> delta(Npts,"delta");
   mode00=false;
   ifft2d(delta,deltaF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,1/(2*PI),mode00);
  
   STOPTIME(tbuff,ckstart);


   // calculate new surface coords
  
   cout << endl;
   cout<<"$ Calculate new plasma surface"<<endl;

   STARTTIME(tbuff,ckstart);

   Vector<p3vector<double> > ksi(Npts,"ksi");
   Vector<p3vector<double> > Xnew(Npts,"Xnew");
   for(unsigned int i=0; i<Npts;i++) {
      ksi[i] = delta[i] * dx_dr[i];
      Xnew[i] = X[i] + ksi[i];
   }
  

   // ALL DONE, NOW SAVE TO FILES

   OmegaInv.clear();


   massage(deltaF,1e-8);
   fname = "deltaF"+plasma_name+methstr+".out";
   dispcr(fname);
   save_coefs(fname,CoefFileFormat_sincos,nnR,mmR,deltaF);

   cout << endl;
   cout<<"$ Generate orthonormal series matrix ("<<Npts<<" x "<<NF<<")"<<endl;

   STARTTIME(tbuff,ckstart);

   Matrix<complex<double> > fs(Npts,NF,"fs");
   fseries(nn,mm,thetas,phis,fs);

   STOPTIME(tbuff,ckstart);


   FourierSurface newfourier;
   transformsurface(Xnew,newfourier,fs,nn,mm);
   massage(newfourier,1e-8);
   fname = "perturbedsurface"+plasma_name+methstr+".out";
   dispcr(fname);
   save_fourier_surface(fname,newfourier);

   return 0;
} // main()





