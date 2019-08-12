/************************************************************************* 
 * 
 *   File Name    :  fgcoils_final.cpp
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
 * DIFFERENCES FROM fgcoils_IIe
 * use N3 highest modes, instead of NnnIF, NmmIF
 *
 * VERSION NOTES:
 *
 * based on fgcoils_IIe.cpp of 2007-10-19
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

   enable_option(opt_N3);
   //   enable_option(opt_condnum);

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

   string methstr(".fgfinal");
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

   strmtemp <<N3;
   string N3_str(strmtemp.str());
   strmtemp.str("");
   N3_str = ".Nmodes="+N3_str;


   dispcr(plasma_name);
   dispcr(methstr);
   dispcr(Nphi_str+Ntheta_str);
   dispcr(Nnn_str+Nmm_str);
   dispcr(N3_str);

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
    


   LAvector<complex<double> > FluxF(NFR,"FluxF");

   cout <<endl<< "$ Loading Plasma Flux sin/cos fourier coefficients from " << flux_filename << endl;
   if (load_coefs(flux_filename,CoefFileFormat_sincos,nnR,mmR,FluxF,false)){
      printcr("Above ERROR occurred in "+myname+".");
      return 3;
   }

   LAvector<double> Flux(Npts, "Flux");
   mode00=false;
   ifft2d(Flux,FluxF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,1/(2*PI),mode00);


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

   LAvector<p3vector<double> > rootn(Npts, "rootn");
   LAvector<double> J(Npts, "J");
   for (unsigned int j =0; j<Npts; j++) {
      J[j] = dot(dx_dr[j],cross(dx_dtheta[j], dx_dphi[j]));
      double temp = 1/sqrt(norm(dA_dtdp[j]));
      rootn[j] = dA_dtdp[j]*temp;
      Flux[j] = Flux[j]*temp;
   }
   mode00=false;
   fft2d(Flux,FluxF,Nphi,Ntheta,Nnn,Nmm,Nharm,Mharm,1e-10,(2*PI),mode00);
   fname= "FluxF.scaled.out";
   save_coefs(fname,CoefFileFormat_sincos,nnR,mmR,FluxF);

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



   // Create mutual inductance matrix

   Matrix<complex<double> > MF(NFR,NFR,"MF"); 
   if ( M_filename.empty() ) {

      // Create mutual inductance matrix

      cout << endl;
      cout<<"$ Generating inductance matrix ("<<Npts<<" x "<<Npts<<")"<<endl;

      STARTTIME(tbuff,ckstart);

      Matrix<double> M(Npts, Npts, "M");
      inductancematrix(X,rootn,Xcoil,dA_dtdp_coil,M);
   
      STOPTIME(tbuff,ckstart);

      // fft  M matrix
      cout << endl;
      cout<<"$ FFT of inductance matrix ("<<Npts<<" x "<<NFR<<")"<<endl;
      
      STARTTIME(tbuff,ckstart);
  
      fft_of_M(M, MF, Nphi, Ntheta,Nnn,Nmm,  Nharm,Mharm,1e-12);
      M.clear();
      STOPTIME(tbuff,ckstart);

      fname = "MF"+plasma_name+Nnn_str+Nmm_str+".R.out";
      data.resize() = real(MF);
      data.perline(1);
      data.textformat(text_nobraces);
      dispcr(fname);
      save(data,fname);
      fname = "MF"+plasma_name+Nnn_str+Nmm_str+".I.out";
      data = imag(MF);
      dispcr(fname);
      save(data,fname);
      data.clear();

   } else {
      cout << endl;
      cout<<"$ Load M (inductance) matrix("<<NFR<<"x"<<NFR<<")"<<endl;
      
      STARTTIME(tbuff,ckstart);
      
      data.resize(NFR,NFR);
      data.perline(1);
      data.textformat(text_nobraces);
      fname = M_filename +Nnn_str+Nmm_str+ ".R.out";
      dispcr(fname);
      load(data,fname);
      Matrix <double> data2("data2");
      data2.resize(NFR,NFR);
      data2.perline(1);
      data2.textformat(text_nobraces);
      fname =  M_filename +Nnn_str+Nmm_str+ ".I.out";
      dispcr(fname);
      load(data2,fname);
      MF = mcomplex(data,data2);
      data.clear();
      data2.clear();

      STOPTIME(tbuff,ckstart);

   }

   //////////////////////////////////////////////////////
   // this (Nmodes, FluxFP etc) doesn't have any purpose other
   // than to prep the two codes for merger
   const unsigned int Nmodes = NFR;
   LAvector<complex<double> > FluxFP(Nmodes,"FluxFP");
   FluxFP = FluxF;

   Matrix<complex<double> > FG(Nmodes,NFR,"FG"); 
   FG = MF;
   MF.clear();

   STOPTIME(tbuff,ckstart);
   //////////////////////////////////////////

   dispcr(FG.Nrows());
   dispcr(FG.Ncols());


   LAvector<double> FGmag(NFR,"FGmag");

   for (unsigned int k=0; k<NFR; k++) {
       double temp = 0;
      for (unsigned int j =0; j<Nmodes; j++) {
	 temp += FG(j,k).real()*FG(j,k).real() + FG(j,k).imag()*FG(j,k).imag();
      }      
      FGmag[k] = sqrt(temp);
   }


   fname = "FGmag"+plasma_name+Nnn_str+Nmm_str+".out";
   FGmag.perline(1);
   FGmag.textformat(text_nobraces);
   save(FGmag,fname);

   //project flux onto normalized FG functions
   printcr("calc Flux as projected onto normalized FG functions");
   STARTTIME(tbuff,ckstart);

   LAvector<complex<double> > FluxFG(NFR,"FluxFG");
   FluxFG = (FluxFP|FG)/FGmag;

   STOPTIME(tbuff,ckstart);

   fname = "FluxFG"+plasma_name+methstr+Nnn_str+Nmm_str+N3_str+".R.out";
   datavec.resize() = real(FluxFG);
   datavec.perline(1);
   datavec.textformat(text_nobraces);
   save(datavec,fname);
   fname = "FluxFG"+plasma_name+methstr+Nnn_str+Nmm_str+N3_str+".I.out";
   datavec = imag(FluxFG);
   save(datavec,fname);


   LAvector<unsigned int> FGrank(NFR,"FGrank");
   LAvector<double> magtemp(NFR,"magtemp");
   magtemp = FGmag;
   FGrank = sortrevwind(magtemp);

   fname = "FGrank"+plasma_name+Nnn_str+Nmm_str+".out";
   FGrank.perline(1);
   FGrank.textformat(text_nobraces);
   save(FGrank,fname);


   // next step is to orthogonalize (not normalize) using gram-schmidt
   // use rank order for this
   // then calc mag and rank of orthogonal FG functions

   printcr("create orthogonal FG functions using Gram-Schmidt process");
   STARTTIME(tbuff,ckstart);

   Matrix<double> Jmatrix(NFR,N3,"Jmatrix");
   Jmatrix=0.0;

   Matrix<complex<double> > Up(N3,N3,"Up"); 
   Up=0.0;
 
   Matrix<complex<double> > FG_ortho(Nmodes,N3,"FG_ortho"); 
   for (unsigned int kk=0; kk<N3; kk++) {

      unsigned int ll = FGrank[kk];

      Jmatrix(ll,kk) = 1.0;
      Up(kk,kk) = 1.0;

      LAvector<complex<double> > vold(Nmodes,"vold");
      vold = FG.col(ll);

      LAvector<complex<double> > v(Nmodes,"v");
      v = vold;
      for (unsigned int h=0; h<kk; h++) {
	 LAvector<complex<double> > v_h(Nmodes,"v_h");
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

   fname = "FG_ortho"+plasma_name+Nnn_str+Nmm_str+".R.out";
   data.resize() = real(FG_ortho);
   data.perline(1);
   data.textformat(text_nobraces);
   save(data,fname);
   fname = "FG_ortho"+plasma_name+Nnn_str+Nmm_str+".I.out";
   data = imag(FG_ortho);
   save(data,fname);


   fname = "Jmatrix"+plasma_name+Nnn_str+Nmm_str+".out";
   Jmatrix.perline(1);
   Jmatrix.textformat(text_nobraces);
   save(Jmatrix,fname);

   fname = "Up"+plasma_name+Nnn_str+Nmm_str+".R.out";
   data.resize() = real(Up);
   data.perline(1);
   data.textformat(text_nobraces);
   save(data,fname);
   fname = "Up"+plasma_name+Nnn_str+Nmm_str+".I.out";
   data = imag(Up);
   save(data,fname);


   printcr("calculate magnitude of orthogonal FGs");
   STARTTIME(tbuff,ckstart);

   LAvector<double> FG_orthomag(N3,"FG_orthomag");

   for (unsigned int k=0; k<N3; k++) {
      double temp = 0;
      for (unsigned int j =0; j<Nmodes; j++) {
	 temp += FG_ortho(j,k).real()*FG_ortho(j,k).real() + FG_ortho(j,k).imag()*FG_ortho(j,k).imag();
      }      
      FG_orthomag[k] = sqrt(temp);
   }
   STOPTIME(tbuff,ckstart);

   fname = "FG_orthomag"+plasma_name+methstr+Nnn_str+Nmm_str+N3_str+
      ".out";
   FG_orthomag.perline(1);
   FG_orthomag.textformat(text_nobraces);
   save(FG_orthomag,fname);


   printcr("find Dinv and normalize orthogonal FGs");
   STARTTIME(tbuff,ckstart);

   LAvector<double> Dinv(N3,"Dinv");
   for (unsigned int kk=0; kk<N3; kk++) {
      double dtemp = 1/FG_orthomag[kk];
      if ( FG_orthomag[kk] == 0) 
	 dtemp = 0;
      if (kk>=Nmodes)
	 dtemp = 0;
      Dinv[kk] = dtemp;
   }

   for (unsigned int kk=0; kk<N3; kk++) {
      double dtemp = Dinv[kk];
      for (unsigned int j =0; j<Nmodes; j++) {
	 FG_ortho(j,kk) = dtemp * FG_ortho(j,kk);
      }
   }
   STOPTIME(tbuff,ckstart);

   fname = "Dinv"+plasma_name+methstr+Nnn_str+Nmm_str+N3_str+
      ".out";
   Dinv.perline(1);
   Dinv.textformat(text_nobraces);
   save(Dinv,fname);


   fname = "FG_orthoNormal"+plasma_name+Nnn_str+Nmm_str+".R.out";
   data.resize() = real(FG_ortho);
   data.perline(1);
   data.textformat(text_nobraces);
   save(data,fname);
   fname = "FG_orthoNormal"+plasma_name+Nnn_str+Nmm_str+".I.out";
   data = imag(FG_ortho);
   save(data,fname);



   //project flux onto orthonormalized FG functions
   printcr("project flux onto orthofunctions");
   STARTTIME(tbuff,ckstart);

   LAvector<complex<double> > FluxFGo(N3,"FluxFGo");
   FluxFGo = (FluxFP|FG_ortho);

   STOPTIME(tbuff,ckstart);

   fname = "FluxFGorthonorm"+plasma_name+methstr+Nnn_str+Nmm_str+N3_str+
      ".R.out";
   datavec.resize() = real(FluxFGo);
   datavec.perline(1);
   datavec.textformat(text_nobraces);
   save(datavec,fname);
   fname = "FluxFGorthonorm"+plasma_name+methstr+Nnn_str+Nmm_str+N3_str+
      ".I.out";
   datavec = imag(FluxFGo);
   save(datavec,fname);


   /////////////////////////////////////////////
   // above is same fourgreen.cpp
   // below is code to calculate Ginv, which is 
   // equivalent ot Morse-Penrose pseudo-inverse
   // and equivalent to merkel method
   ////////////////////////////////////////////


   printcr("calculate (J|U) matrix");
   STARTTIME(tbuff,ckstart);

   Matrix<complex<double> > JU(NFR,N3,"JU"); 
   JU = (Jmatrix|Up);
   Up.clear();
   Jmatrix.clear();
   //   Up.clear();
   STOPTIME(tbuff,ckstart);

   fname = "JU"+plasma_name+Nnn_str+Nmm_str+".R.out";
   data.resize() = real(JU);
   data.perline(1);
   data.textformat(text_nobraces);
   save(data,fname);
   fname = "JU"+plasma_name+Nnn_str+Nmm_str+".I.out";
   data = imag(JU);
   save(data,fname);
   

   printcr("calculate (Omega|Dinv)");
   STARTTIME(tbuff,ckstart);

   for (unsigned int kk=0; kk<N3; kk++) {
      double dtemp = Dinv[kk];
      for (unsigned int j =0; j<Nmodes; j++) {
	 FG_ortho(j,kk) = dtemp * FG_ortho(j,kk);
      }
   }

   STOPTIME(tbuff,ckstart);



   printcr("calculate Ginv matrix");
   STARTTIME(tbuff,ckstart);

   Matrix<complex<double> > Ginv(NFR,Nmodes,"Ginv"); 
   Ginv = (JU|adj(FG_ortho));

   FG_ortho.clear();
   JU.clear();
   STOPTIME(tbuff,ckstart);




   // calculate current
  
   cout << endl;
   cout<<"$ Calculate coil current ("<<NFR<<" x 1)"<<endl;

   STARTTIME(tbuff,ckstart);
   LAvector<complex<double> > IF(NFR,"IF");
   IF = (Ginv|FluxFP);

   STOPTIME(tbuff,ckstart);

   Ginv.clear();


   // ALL DONE, NOW SAVE TO FILES

   fname="NFGmodes"+plasma_name+methstr+Nnn_str+Nmm_str+N3_str+".out";
   save(N3,fname);

   massage(IF,1e-8);

   fname = "IF"+plasma_name+methstr+Nnn_str+Nmm_str+N3_str+".out";
   save_coefs(fname,CoefFileFormat_sincos,nnR,mmR,IF);

   return 0;
} // main()





