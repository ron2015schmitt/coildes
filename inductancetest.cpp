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
#include "rhomatrix.hpp"


const double NEGLECT =  1e-12;




// Main Function for code

int main (int argc, char *argv[])
{

   disable_all_options();
   enable_option(opt_pf);
   enable_option(opt_cf);
   enable_option(opt_Nphi);
   enable_option(opt_Ntheta);
   //   enable_option(opt_Itoroidal);
   enable_option(opt_Ipoloidal);

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
   data.perline(1);
   data.textformat(text_nobraces);
   Matrix <double> data2("data2");
   data2.perline(1);
   data2.textformat(text_nobraces);
   

   string fname;

   // variables for measuring times
   struct tms tbuff;
   clock_t ckstart;


   // display command name
   //   cout << string(myname)<<endl;

   // display Matricks mode
   cout << endl;
   display_execution_mode();
   cout << endl;

 
   ostringstream strmtemp;

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



   // Create angle grid
   const unsigned int Npts = Ntheta*Nphi;
   Vector<double> thetas(Npts,"thetas");
   Vector<double> phis(Npts,"phis");
   anglevectors(thetas, phis, Ntheta, Nphi);

   const double dtheta = 2*M_PI/((double)Ntheta);
   const double dphi = 2*M_PI/((double)Nphi);
   const double dphi_by_dtheta = 2.0*PI * 2.0*PI/double(Npts);
   // coefficient C is the integration coef for the fourier transform
   // C = dtheta*dphi
   //   = (2*pi/Ntheta)*(2*pi/Nphi)
   //   = (2*pi*2*pi/Npts)
   const double C = (2*PI*2*PI/double(Npts));
   const double Csqr = C*C;





   // save angle vectors

   thetas.perline(1);
   thetas.textformat(text_nobraces);
   phis.perline(1);
   phis.textformat(text_nobraces);
   fname =  "thetas"+Ntheta_str+Nphi_str +".out";
   save(thetas,fname);
   fname =  "phis"+Ntheta_str+Nphi_str +".out";
   save(phis,fname);


   
   // load the plasma surface fourier coef's

   cout << "$ Loading PLASMA SURFACE fourier coefficients from " << plasma_filename << endl;
   FourierSurface plasmafourier;
   if (load_fourier_surface(plasma_filename,plasmafourier)) {
      printcr("Above ERROR occurred in "+myname+".");
      return 1;
   }
   plasmafourier.RF().name("p.RF");
   plasmafourier.ZF().name("p.ZF");



   // load the coil surface fourier coef's

   cout << "$ Loading COIL SURFACE fourier coefficients from " << coil_filename << endl;
   FourierSurface coilfourier;
   if (load_fourier_surface(coil_filename,coilfourier)){
      printcr("Above ERROR occurred in "+myname+".");
      return 2;
   }
      
   coilfourier.RF().name("c.RF");
   coilfourier.ZF().name("c.ZF");




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

   Vector<p3vector<double> > nhat(Npts, "nhat");
   Vector<double> J(Npts, "J");

   Vector<double> dA(Npts, "dA");

   for (unsigned int p =0; p<Npts; p++) {
      const double temp = sqrt(dot(dA_dtdp[p],dA_dtdp[p]));
      dA[p] = temp*dphi_by_dtheta;
      nhat[p] = dA_dtdp[p]/temp;
      J[p] = dot(dx_dr[p],cross(dx_dtheta[p], dx_dphi[p]));
   }


   STOPTIME(tbuff,ckstart);


 
  
   // lay coil surface onto grid 
  
   Vector<p3vector<double> > Xcoil(Npts, "Xcoil");
   Vector<p3vector<double> > dA_dtdp_coil(Npts, "dA_dtdp_coil");

   Vector<p3vector<double> > dx_dr_coil(Npts, "dx_dr_coil");
   Vector<p3vector<double> > dx_dtheta_coil(Npts,"dx_dtheta_coil");
   Vector<p3vector<double> > dx_dphi_coil(Npts,"dx_dphi_coil");
   Vector<p3vector<double> > grad_r_coil(Npts,"grad_r_coil");
   Vector<p3vector<double> > grad_theta_coil(Npts,"grad_theta_coil");
   Vector<p3vector<double> > grad_phi_coil(Npts,"grad_phi_coil");

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




   
   Itoroidal = 1;


   Vector<double> j_theta(Npts,"j_theta");
   Vector<double> j_phi(Npts,"j_phi");
   j_theta = 0;
   j_phi = 0;
   
   Vector<p3vector<double> >jcoil(Npts,"jcoil"); 
   Vector<p3vector<double> >jcoil_d(Npts,"jcoil_d"); 

   for (unsigned int p = 0; p<Npts; p++) {
      jcoil[p] = Itoroidal/(2*PI) * dx_dphi_coil[p];
      jcoil_d[p] = -Ipoloidal/(2*PI) * dx_dtheta_coil[p];
   }



   Vector<p3vector<double> > B(Npts, "B");
   Vector<p3vector<double> > Bd(Npts, "Bd");

   unsigned int count = 0;
   unsigned int checkup = static_cast<unsigned int>(Npts*0.01);
   for (unsigned int p =0; p<Npts; p++) {
      if (++count == checkup) {
	 print(Matricks::round(double(p)/Npts*100));cout <<" %"<<endl;
	 count =0;
      }
      p3vector<double> Btemp = p3vector<double>(0,0,0);
      p3vector<double> Bdtemp = p3vector<double>(0,0,0);
      for (unsigned int c = 0; c<Npts; c++) {
	 const p3vector<double> R = X[p] - Xcoil[c];
	 const double r3inv = 1/pow(norm(R),3.0);
	 Btemp = Btemp + mu0div4pi * dphi_by_dtheta  * cross(jcoil[c],R) * r3inv;
	 Bdtemp = Bdtemp + mu0div4pi * dphi_by_dtheta  * cross(jcoil_d[c],R) * r3inv;
      }
      B[p] = Btemp;
      Bd[p] = Bdtemp;
   }

   Vector<double> Bn(Npts, "Bn");
   Vector<double> Bdn(Npts, "Bdn");
   for (unsigned int p =0; p<Npts; p++) {
      Bn[p] = dot(nhat[p],B[p]);
      Bdn[p] = dot(nhat[p],Bd[p]);
   }

   double Integral1 = 0;
   double Integral2 = 0;
   for (unsigned int p =0; p<Npts; p++) {
      Integral1 += dA[p]*Bn[p]*Bdn[p];
      Integral2 += dA[p]*Bn[p]*Bn[p];
   }
   dispcr(Integral1);
   dispcr(Integral2);
   double Itor_mse= Integral1/Integral2;


   dispcr(Itor_mse);

   save(Itor_mse,"Itoroidal.out");

  
   // TESTING
   // these exclude the n=0,m=0 case
   Vector<double> nnR("nnR");
   Vector<double> mmR("mmR");
   unsigned int NFR;
   mode00 = false;
   modevectors(NFR,nnR,mmR,Nnn,Nmm,Nharm,Mharm,mode00);


   Vector<complex<double> > BdnF(NFR,"BdnF");
   mode00=false;
   fft2d(Bdn,BdnF,Nphi,Ntheta,Nnnp,Nmmp,Nharm,Mharm,1e-10,1/(2*PI),mode00);

   massage(BdnF,1e-10);
   save_coefs("BdnF.out",CoefFileFormat_sincos,nnR,mmR,BdnF);


   return 0;
} // main()
