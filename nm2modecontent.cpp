/************************************************************************* 
 * 
 *   File Name    :  nm2modecontent.cpp
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  
 * 
 *
 *   SYNOPSIS     
 *     Stellarator COIL Design 
 *     This file finds the coil current from given plasma surface flux 
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
   enable_option(opt_Nnn);
   enable_option(opt_Nmm);
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

   LAvector <double> datavec("datavec");
   datavec.perline(1);
   datavec.textformat(text_nobraces);
   Matrix <double> data("data");
   data.perline(1);
   data.textformat(text_nobraces);
   Matrix <double> data2("data2");
   LAvector <double> datavec2("datavec2");

   string fname;

   // variables for measuring times
   struct tms tbuff;
   clock_t ckstart;

   // display COOLL mode
   cout << endl;
   display_execution_mode();
   cout << endl;



  
   ostringstream strmtemp;
   strmtemp <<Nnn;
   string Nnn_str(strmtemp.str());
   strmtemp.str("");
   strmtemp <<Nmm;
   string Nmm_str(strmtemp.str());
   strmtemp.str("");


   // Create Fourier Mode vectors
   LAvector<double> nn("nn");
   LAvector<double> mm("mm");
   unsigned int NF;
   bool mode00 = true;
   modevectors(NF,nn,mm,Nnn,Nmm,Nharm,Mharm,mode00);



   // these exclude the n=0,m=0 case
   LAvector<double> nnR("nnR");
   LAvector<double> mmR("mmR");
   unsigned int NFR;
   mode00 = false;
   modevectors(NFR,nnR,mmR,Nnn,Nmm,Nharm,Mharm,mode00);


 
   Matrix<complex<double> > Q(NFR,NFR,"Q");
   LAvector<complex<double> >  L(NFR,"L");

   cout<<"$ Load Omega matrix ("<<NFR<<"x"<<NFR<<")"<<endl;

   STARTTIME(tbuff,ckstart);

   // UR seems to always be effectively zero.  do we need it?

   // PRINT OUT FILENAMES!!!
   data.perline(1);
   data.textformat(text_nobraces);
   data.resize(NFR,NFR);
   fname = "omega." + omega_filename + ".PR.out";
   load(data,fname);
   
   data2.perline(1);
   data2.textformat(text_nobraces);
   data2.resize(NFR,NFR);
   fname = "omega." + omega_filename + ".PI.out";
   load(data2,fname);
   Q = mcomplex(data,data2);
   data.clear();
   data2.clear();

   datavec.perline(1);
   datavec.textformat(text_nobraces);
   datavec.resize(NFR);
   fname = "omega." + omega_filename + ".WR.out";
   load(datavec,fname);

   datavec2.perline(1);
   datavec2.textformat(text_nobraces);
   datavec2.resize(NFR);
   fname = "omega." + omega_filename + ".WI.out";
   load(datavec2,fname);
   L = vcomplex(datavec,datavec2);
   datavec.clear();
   datavec2.clear();

   STOPTIME(tbuff,ckstart);



   Matrix<complex<double> > P(NFR,NFR,"P");
   LAvector<complex<double> >  S(NFR,"S");

   cout<<"$ Load OmegaN matrix ("<<NFR<<"x"<<NFR<<")"<<endl;

   STARTTIME(tbuff,ckstart);

   data.perline(1);
   data.textformat(text_nobraces);
   data.resize(NFR,NFR);
   fname = "omegaN." + omega_filename + ".PR.out";
   load(data,fname);
   
   data2.perline(1);
   data2.textformat(text_nobraces);
   data2.resize(NFR,NFR);
   fname = "omegaN." + omega_filename + ".PI.out";
   load(data2,fname);
   P = mcomplex(data,data2);
   data.clear();
   data2.clear();

   datavec.perline(1);
   datavec.textformat(text_nobraces);
   datavec.resize(NFR);
   fname = "omegaN." + omega_filename + ".WR.out";
   load(datavec,fname);

   datavec2.perline(1);
   datavec2.textformat(text_nobraces);
   datavec2.resize(NFR);
   fname = "omegaN." + omega_filename + ".WI.out";
   load(datavec2,fname);
   S = vcomplex(datavec,datavec2);
   datavec.clear();
   datavec2.clear();

   STOPTIME(tbuff,ckstart);







  
   cout << endl;
   cout<<"$ Hermitian Transpose of Omega Eigenvector Matrix ("<<NFR<<"x"<<NFR<<")"<<endl;
     
   STARTTIME(tbuff,ckstart);

   Matrix<complex<double> > Qinv(NFR,NFR,"Qinv");
   Qinv = adj(Q);
     
   STOPTIME(tbuff,ckstart);
   Q.clear();

   cout << endl;
   cout<<"$ Calculate Transformation Matrix = (Q^(T*)|P) ("<<NFR<<"x"<<NFR<<")"<<endl;
     
   STARTTIME(tbuff,ckstart);

   Matrix<complex<double> > QtP(NFR,NFR,"QtP");
   QtP = (Qinv|P);
     
   STOPTIME(tbuff,ckstart);




   data.resize() = real(QtP);
   fname = "QtP.Nn="+Nnn_str+".Nm="+Nmm_str+".R.out";
   dispcr(fname);
   save(data,fname);
   data.resize() = imag(QtP);
   fname = "QtP.Nn="+Nnn_str+".Nm="+Nmm_str+".I.out";
   dispcr(fname);
   save(data,fname);


   return 0;
} // main()





