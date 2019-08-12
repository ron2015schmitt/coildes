/************************************************************************* 
 * 
 *   File Name    :  fg_mag.cpp
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  
 * 
 *
 *   SYNOPSIS     
 *    
 *
 * project fourer green functions into fourier space
 * and work with flux fourier space instead of flux physical space.
 *
 * In addition to limiting the FG functions using N3 size U and J
 * matrices (technique used in IIb and IIc),
 * we get null any modes that have magnitude less than 
 * condnum * magnitude of largest orthogonal mode, using the Dinv matrix.
 * Also make some name changes so as to allow easy merger with scoild_XXIV
 *
 * DIFFERENCES FROM fgcoils_IId
 *
 * Only calculates first N3 orthonogal FGs, to save time.
 * 
 * removed unecessary calculation of fsR
 *
 * FIXED BUG IN RANK ORDERING 
 *
 * VERSION NOTES:
 *
 * based on fgcoils_IIe.cpp 2007-10-12
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

// This is the file that defines the B field configuration



// Main Function for code

int main (int argc, char *argv[])
{

   disable_all_options();
   enable_option(opt_pf);
   enable_option(opt_cf);
   enable_option(opt_ff);
   enable_option(opt_Nphi);
   enable_option(opt_Ntheta);
   enable_option(opt_Nnn);
   enable_option(opt_Nmm);
   enable_option(opt_Nharm);
   enable_option(opt_Mharm);

   enable_option(opt_NmmIF);
   enable_option(opt_NnnIF);
   enable_option(opt_condnum);

   // the following file can be provided to speed things up on repeated runs
   enable_option(opt_Mf);


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

   LAvector <double> datavec("datavec");
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

   // display COOLL mode
   cout << endl;
   display_execution_mode();
   cout << endl;

 
   ostringstream strmtemp;

   string methstr(".fgmag");
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

   strmtemp <<NnnIF;
   string NnnIF_str(strmtemp.str());
   strmtemp.str("");
   NnnIF_str = ".NnIF="+NnnIF_str;

   strmtemp <<NmmIF;
   string NmmIF_str(strmtemp.str());
   strmtemp.str("");
   NmmIF_str = ".NmIF="+NmmIF_str;

   strmtemp <<condnum;
   string condnum_str(strmtemp.str());
   strmtemp.str("");
   condnum_str = ".cond="+condnum_str;


   dispcr(plasma_name);
   dispcr(methstr);
   dispcr(Nphi_str+Ntheta_str);
   dispcr(Nnn_str+Nmm_str);
   dispcr(NnnIF_str+NmmIF_str);
   dispcr(condnum_str);


   // Create angle grid
   const unsigned int Npts = Ntheta*Nphi;
   LAvector<double> thetas(Npts,"thetas");
   LAvector<double> phis(Npts,"phis");
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



   // these exclude the n=0,m=0 case
   LAvector<double> nnR("nnR");
   LAvector<double> mmR("mmR");
   unsigned int NFR;
   bool mode00 = false;
   modevectors(NFR,nnR,mmR,Nnn,Nmm,Nharm,Mharm,mode00);


   // save angle vectors

   thetas.perline(1);
   thetas.textformat(text_nobraces);
   phis.perline(1);
   phis.textformat(text_nobraces);
   fname =  "thetas"+Ntheta_str+Nphi_str +".out";
   save(thetas,fname);
   fname =  "phis"+Ntheta_str+Nphi_str +".out";
   save(phis,fname);

   // save mode vectors

   fname = "nnR"+Nnn_str+Nmm_str+".out";
   nnR.perline(1);
   nnR.textformat(text_nobraces);
   save(nnR,fname);
   fname = "mmR"+Nnn_str+Nmm_str+".out";
   mmR.perline(1);
   mmR.textformat(text_nobraces);
   save(mmR,fname);

   


   // load the plasma surface fourier coef's

   cout << "$ Loading PLASMA SURFACE fourier coefficients from " << plasma_filename << endl;
   FourierSurface plasmafourier;
   if (load_fourier_surface(plasma_filename,plasmafourier)) {
      printcr("Above ERROR occurred in "+myname+".");
      return 1;
   }
   plasmafourier.RF().name("p.RF");
   plasmafourier.ZF().name("p.ZF");

   // print coef's
   //  printfouriercoefs(plasmafourier.nn(),plasmafourier.mm(),plasmafourier.RF(),plasmafourier.ZF(),10,18);
  


   // load the coil surface fourier coef's

   cout << "$ Loading COIL SURFACE fourier coefficients from " << coil_filename << endl;
   FourierSurface coilfourier;
   if (load_fourier_surface(coil_filename,coilfourier)){
      printcr("Above ERROR occurred in "+myname+".");
      return 2;
   }
      
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

   LAvector<p3vector<double> > dx_dr_coil(Npts, "dx_dr_coil");
   LAvector<p3vector<double> > dx_dtheta_coil(Npts,"dx_dtheta_coil");
   LAvector<p3vector<double> > dx_dphi_coil(Npts,"dx_dphi_coil");
   LAvector<p3vector<double> > grad_r_coil(Npts,"grad_r_coil");
   LAvector<p3vector<double> > grad_theta_coil(Npts,"grad_theta_coil");
   LAvector<p3vector<double> > grad_phi_coil(Npts,"grad_phi_coil");

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




   cout << endl;
   cout<<"$ Generating inductance matrix ("<<Npts<<" x "<<Npts<<")"<<endl;

   STARTTIME(tbuff,ckstart);

   Matrix<double> M(Npts, Npts, "M");
   inductancematrix(X,dA_dtdp,Xcoil,dA_dtdp_coil,M);
      
   STOPTIME(tbuff,ckstart);


   LAvector<double> fgmag(NFR,"fgmag");

   LAvector<complex<double> > f(Npts,"f");
   LAvector<complex<double> > fg(Npts,"fg");

   STARTTIME(tbuff,ckstart);

  unsigned int count = 0;
  unsigned int checkup = static_cast<unsigned int>(NFR*0.001);

   for (unsigned int k = 0; k<NFR ; k++) {
    if (++count == checkup) {
      print(double(k)/NFR*100);cout <<" %";
       time_t tim=time(0);
      std::cout << ctime(&tim);
      std::cout.flush();
      count =0;
    }
      LAvector<double> phase(Npts,"phase");
      phase= nnR[k]*phis + mmR[k]*thetas;
      f = vcomplex(cos(phase),sin(phase)); 
      //normalize
      const double coef = 1.0/(2.0*PI);
      f=f*coef;
      fg=(M|f);
      double temp = 0;
      for (unsigned int j =0; j<Npts; j++) {
	 temp += fg[j].real()*fg[j].real() + fg[j].imag()*fg[j].imag();
      }      
      fgmag[k] = sqrt(temp);
   }

   STOPTIME(tbuff,ckstart);

   fname = "fg_mag"+plasma_name+Nnn_str+Nmm_str+".out";
   fgmag.perline(1);
   fgmag.textformat(text_nobraces);
   save(fgmag,fname);

   return 0;
} // main()





