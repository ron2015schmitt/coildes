/************************************************************************* 
 * 
 *   File Name    :  fourgreen
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  
 * 
 *
 *   SYNOPSIS     
 *    
 *
 * VERSION NOTES:
 *
 *  This version was based on scoild_fft_XXII.cpp on 2006 nov 29
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
//    enable_option(opt_Btf);
//    enable_option(opt_Bpf);
   enable_option(opt_Nharm);
   enable_option(opt_Mharm);

//    enable_option(opt_omegaf);
//    enable_option(opt_Npert);
//    enable_option(opt_Nprin);
   //   enable_option(opt_Ncomb);
   //   enable_option(opt_ic);
   //   enable_option(opt_modesf);
   //   enable_option(opt_alpha_theta);  
   //   enable_option(opt_alpha_phi);

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


   // display command name
   cout << string(myname)<<endl;

   // display Matricks mode
   cout << endl;
   display_execution_mode();
   cout << endl;


   // Create angle grid
   const unsigned int Npts = Ntheta*Nphi;
   Vector<double> thetas(Npts,"thetas");
   Vector<double> phis(Npts,"phis");
   anglevectors(thetas, phis, Ntheta, Nphi);

   const double dtheta = 2*M_PI/((double)Ntheta);
   const double dphi = 2*M_PI/((double)Nphi);
   const double dphi_by_dtheta = 2.0*PI * 2.0*PI/double(Npts);



  
   ostringstream strmtemp;
   strmtemp <<Nnn;
   string Nnn_str(strmtemp.str());
   strmtemp.str("");
   strmtemp <<Nmm;
   string Nmm_str(strmtemp.str());
   strmtemp.str("");

   string plasma_name;
   {
      unsigned int j = plasma_filename.rfind("_");
       plasma_name = plasma_filename.substr(0,j);
   }
   string fhead(plasma_name+".Nn="+Nnn_str+".Nm="+Nmm_str);
   dispcr(plasma_name);
   dispcr(fhead);

   // Create Fourier Mode vectors
   Vector<double> nn("nn");
   Vector<double> mm("mm");
   unsigned int NF;
   bool mode00 = true;
   modevectors(NF,nn,mm,Nnn,Nmm,Nharm,Mharm,mode00);



   // these exclude the n=0,m=0 case
   Vector<double> nnR("nnR");
   Vector<double> mmR("mmR");
   unsigned int NFR;
   mode00 = false;
   modevectors(NFR,nnR,mmR,Nnn,Nmm,Nharm,Mharm,mode00);


   fname = "nnR."+fhead+".out";
   nnR.perline(1);
   nnR.textformat(text_nobraces);
   save(nnR,fname);
   fname = "mmR."+fhead+".out";
   mmR.perline(1);
   mmR.textformat(text_nobraces);
   save(mmR,fname);


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
    


  Vector<complex<double> > FluxF(NFR,"FluxF");

   cout <<endl<< "$ Loading Plasma Flux sin/cos fourier coefficients from " << flux_filename << endl;
   if (load_coefs(flux_filename,CoefFileFormat_sincos,nnR,mmR,FluxF,false))
      return 3;

   {
      unsigned int j = flux_filename.rfind("out");
      string strtemp = flux_filename.substr(0,j);
      dispcr(strtemp);
      fname = strtemp + "exp.R.out";
      datavec.resize() = real(FluxF);
      datavec.perline(1);
      datavec.textformat(text_nobraces);
      save(datavec,fname);
      fname = strtemp + "exp.I.out";
      datavec.resize() = imag(FluxF);
      datavec.perline(1);
      datavec.textformat(text_nobraces);
      save(datavec,fname);

   }

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


   // Create mutual inductance matrix

   cout << endl;
   cout<<"$ Generating inductance matrix ("<<Npts<<" x "<<Npts<<")"<<endl;

   STARTTIME(tbuff,ckstart);

   Matrix<double> M(Npts, Npts, "M");
   inductancematrix(X,dA_dtdp,Xcoil,dA_dtdp_coil,M);

   STOPTIME(tbuff,ckstart);


   Matrix<complex<double> > FG(Npts,NFR,"FG"); 

   cout << endl;
   cout<<"$ FFT coil side of inductance matrix -- Fourier-Green Functions ("<<Npts<<" x "<<NFR<<")"<<endl;
   
   STARTTIME(tbuff,ckstart);
      
   half_fft_of_M(M, FG, Nphi, Ntheta, Nnn, Nmm, Nharm,Mharm,1e-14);

   STOPTIME(tbuff,ckstart);

   fname = "gmn."+fhead+".R.out";
   data.resize() = real(FG);
   data.perline(1);
   data.textformat(text_nobraces);
   save(data,fname);
   fname = "gmn."+fhead+".I.out";
   data = imag(FG);
   save(data,fname);

   Vector<double> FGmag(NFR,"FGmag");

   for (unsigned int k=0; k<NFR; k++) {
      double temp = 0;
      for (unsigned int j =0; j<Npts; j++) {
	 temp += FG(j,k).real()*FG(j,k).real() + FG(j,k).imag()*FG(j,k).imag();
      }      
      FGmag[k] = sqrt(temp);
   }

   massage(FGmag,1e-8);

   fname = "FGmag."+fhead+".out";
   FGmag.perline(1);
   FGmag.textformat(text_nobraces);
   save(FGmag,fname);

   Vector<unsigned int> FGrank(NFR,"FGrank");
  Vector<double> magtemp(NFR,"magtemp");
   magtemp = FGmag;
   FGrank = sortrevwind(magtemp);
   fname = "FGrank."+fhead+".out";
   FGrank.perline(1);
   FGrank.textformat(text_nobraces);
   save(FGrank,fname);

   M.clear();


//project flux onto normalized FG functions
   printcr("calc Flux as projected onto normalized FG functions");
   STARTTIME(tbuff,ckstart);

   Vector<double> Flux(Npts, "Flux");
   expandfunction(Flux,FluxF,fsR);

   double FluxRMS = norm(Flux);
dispcr(FluxRMS);


   Vector<complex<double> > FluxFG(NFR,"FluxFG");
   FluxFG = (Flux|FG)/FGmag;

   STOPTIME(tbuff,ckstart);

   fname = "FluxFG."+fhead+".R.out";
   datavec.resize() = real(FluxFG);
   datavec.perline(1);
   datavec.textformat(text_nobraces);
   save(datavec,fname);
   fname = "FluxFG."+fhead+".I.out";
   datavec = imag(FluxFG);
   save(datavec,fname);


   // next step is to orthogonalize (not normalize) using gram-schmidt
   // use rank order for this
   // then calc mag and rank of orthogonal FG functions

   printcr("create orthogonal FG functions using Gram-Schmidt process");
   STARTTIME(tbuff,ckstart);

   Matrix<double> Jmatrix(NFR,NFR,"Jmatrix");
   Jmatrix=0.0;

   Matrix<complex<double> > Up(NFR,NFR,"Up"); 
   Up=0.0;
 
   Matrix<complex<double> > FG_ortho(Npts,NFR,"FG_ortho"); 
   for (unsigned int kk=0; kk<NFR; kk++) {

      unsigned int ll = FGrank[kk];

      Jmatrix(ll,kk) = 1.0;
      Up(kk,kk) = 1.0;

      Vector<complex<double> > vold(Npts,"vold");
      vold = FG.col(ll);

      Vector<complex<double> > v(Npts,"v");
      v = vold;
      for (unsigned int h=0; h<kk; h++) {
	 Vector<complex<double> > v_h(Npts,"v_h");
	 v_h = FG_ortho.col(h);  
	 complex<double> ratio = (conj(v_h)|vold)/(conj(v_h)|v_h);
	 v = v - ratio * v_h;
	 for (unsigned int h2=0; h2<=h; h2++) {
	    Up(h2,kk) += -ratio*Up(h2,h);
	 }
      }
      FG_ortho.col(kk) = v;
   }
   FG.clear();

   STOPTIME(tbuff,ckstart);

   fname = "gmn_ortho."+fhead+".R.out";
   data.resize() = real(FG_ortho);
   data.perline(1);
   data.textformat(text_nobraces);
   save(data,fname);
   fname = "gmn_ortho."+fhead+".I.out";
   data = imag(FG_ortho);
   save(data,fname);


   fname = "Jmatrix."+fhead+".out";
   Jmatrix.perline(1);
   Jmatrix.textformat(text_nobraces);
   save(Jmatrix,fname);

   fname = "Up."+fhead+".R.out";
   data.resize() = real(Up);
   data.perline(1);
   data.textformat(text_nobraces);
   save(data,fname);
   fname = "Up."+fhead+".I.out";
   data = imag(Up);
   save(data,fname);


   printcr("normalized magnitude");
   STARTTIME(tbuff,ckstart);

   Vector<double> FG_orthomag(NFR,"FG_orthomag");

   for (unsigned int k=0; k<NFR; k++) {
      double temp = 0;
      for (unsigned int j =0; j<Npts; j++) {
	 temp += FG_ortho(j,k).real()*FG_ortho(j,k).real() + FG_ortho(j,k).imag()*FG_ortho(j,k).imag();
      }      
      FG_orthomag[k] = sqrt(temp);
   }
   STOPTIME(tbuff,ckstart);

   massage(FG_orthomag,1e-8);

   fname = "FG_orthomag."+fhead+".out";
   FG_orthomag.perline(1);
   FG_orthomag.textformat(text_nobraces);
   save(FG_orthomag,fname);


   printcr("normalize");
   STARTTIME(tbuff,ckstart);

   for (unsigned int kk=0; kk<NFR; kk++) {
 //     dispcr(kk);
      Vector<complex<double> > vold(Npts,"vold");
      for (unsigned int j =0; j<Npts; j++) {
        vold[j] = FG_ortho(j,kk) / FG_orthomag[kk];
      }
      for (unsigned int j =0; j<Npts; j++) {
        FG_ortho(j,kk) = vold[j];
      }
   }
   STOPTIME(tbuff,ckstart);


   fname = "gmn_orthoNormal."+fhead+".R.out";
   data.resize() = real(FG_ortho);
   data.perline(1);
   data.textformat(text_nobraces);
   save(data,fname);
   fname = "gmn_orthoNormal."+fhead+".I.out";
   data = imag(FG_ortho);
   save(data,fname);



   //project flux onto orthonormalized FG functions
   printcr("project flux onto orthofunctions");
   STARTTIME(tbuff,ckstart);

   FluxFG = (Flux|FG_ortho);

   STOPTIME(tbuff,ckstart);

   fname = "FluxFGorthonorm."+fhead+".R.out";
   datavec.resize() = real(FluxFG);
   datavec.perline(1);
   datavec.textformat(text_nobraces);
   save(datavec,fname);
   fname = "FluxFGorthonorm."+fhead+".I.out";
   datavec = imag(FluxFG);
   save(datavec,fname);

   return 0;
} // main()





