/************************************************************************* 
 * 
 *   File Name    :  
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  16 May 2003
 * 
 *
 *   SYNOPSIS     
 *
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


const double NEGLECT =  1e-12;




// Main Function for code

int main (int argc, char *argv[])
{
  disable_all_options();
  enable_option(opt_pf);
  enable_option(opt_Nphi);
  enable_option(opt_Ntheta);
  enable_option(opt_Nnn);
  enable_option(opt_Nmm);
 enable_option(opt_Nharm);
  enable_option(opt_Mharm);
 
  // parse command line input
 
  if (!parse_cmd(argc, argv))
     return 1;


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


  // Create Fourier Mode vectors
  Vector<double> nn("nn");
  Vector<double> mm("mm");
  unsigned int NF;
  bool mode00 = true;
  if ( (Nharm >1) ||(Mharm>1) )
     modevectors(NF,nn,mm,Nnn,Nmm,Nharm,Mharm,mode00);
  else
     modevectors(NF,nn,mm,Nnn,Nmm,mode00);

 // load the surface fourier coef's

  cout << "$ Loading SURFACE fourier coefficients from " << plasma_filename << endl;
  FourierSurface plasmafourier;
  if (load_fourier_surface(plasma_filename,plasmafourier))
    return 1;
  plasmafourier.RF().name("p.RF");
  plasmafourier.ZF().name("p.ZF");

  // print coef's
  //  printfouriercoefs(plasmafourier.nn(),plasmafourier.mm(),plasmafourier.RF(),plasmafourier.ZF(),10,18);
  


    


 
 
  // lay surface onto grid 
  
  Vector<p3vector<double> > X(Npts, "X");
  Vector<p3vector<double> > A(Npts, "A");

  cout << endl;
  cout <<"$ Mapping plasma surface fourier coefficients to "<<Ntheta<<" x "<<Nphi<<" (theta by phi) grid"<<endl;

  STARTTIME(tbuff,ckstart);

  expandsurface(X,A,plasmafourier,thetas,phis);

  STOPTIME(tbuff,ckstart);


  // create surface normals
  Vector<p3vector<double> > n(Npts, "n");
  for (unsigned int j =0; j<Npts; j++)
    n[j] = A[j] / norm(A[j]);


  Vector<double> r(Npts, "r");
  Vector<double> z(Npts, "z");
  for (unsigned int j =0; j<Npts; j++) {
    r[j] = sqrt(sqr(X[j].x()) + sqr(X[j].y()));
    z[j] = X[j].z();
  }
  dispcr(min(r));
  dispcr(max(r));

  ostringstream strm;

  // ALL DONE, NOW SAVE TO FILES
  p3vectorformat::textformat(text_nobraces);
  X.perline(1);
  X.textformat(text_nobraces);
  n.perline(1);
  n.textformat(text_nobraces);

  double x0 = 0.0;
  for (unsigned int j =0; j<Npts; j++) {
     if ( (X[j].y()==0.0) && (X[j].z()==0.0))
	x0 = X[j].x();
  }
  dispcr(x0);
  save(x0,"expandsurf.x0.out",6);  // 6 digit precision (5 after decimal)
  

  string plasma_rootname;
  string fn_ext;
  disect_filename(plasma_filename,plasma_rootname,fn_ext);

  strm.str("");
  strm << plasma_rootname << "_" << Ntheta<<"x"<<Nphi<<".out";
  print("Saving surface data to ");
  printcr(strm.str());
  save(X,strm.str());

  strm.str("");
  strm << plasma_rootname << "_normals" << Ntheta<<"x"<<Nphi<<".out";
  print("Saving surface normals to ");
  printcr(strm.str());
  save(n,strm.str());

  strm.str("");
  strm << plasma_rootname << "_Ntheta.out";
  print("Saving Ntheta to ");
  printcr(strm.str());
  save(Ntheta,strm.str());

  strm.str("");
  strm << plasma_rootname << "_Nphi.out";
  print("Saving Nphi to ");
  printcr(strm.str());
  save(Nphi,strm.str());
 

  r.perline(1);
  r.textformat(text_nobraces);
  strm.str("");
  strm << plasma_rootname << "_R_" << Ntheta<<"x"<<Nphi<<".out";
  print("Saving surface R data to ");
  printcr(strm.str());
  save(r,strm.str());


  z.perline(1);
  z.textformat(text_nobraces);
  strm.str("");
  strm << plasma_rootname << "_Z_" << Ntheta<<"x"<<Nphi<<".out";
  print("Saving surface Z data to ");
  printcr(strm.str());
  save(z,strm.str());


  thetas.perline(1);
  thetas.textformat(text_nobraces);
  strm.str("");
  strm << plasma_rootname << "_Thetas_" << Ntheta<<"x"<<Nphi<<".out";
  print("Saving Theta data to ");
  printcr(strm.str());
  save(thetas,strm.str());

  phis.perline(1);
  phis.textformat(text_nobraces);
  strm.str("");
  strm << plasma_rootname << "_Phis_" << Ntheta<<"x"<<Nphi<<".out";
  print("Saving phi data to ");
  printcr(strm.str());
  save(phis,strm.str());


  return 0;
} // main()





