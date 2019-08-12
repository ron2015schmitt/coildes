/************************************************************************* 
 * 
 *   File Name    :  
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  
 * 
 *
 *   SYNOPSIS     
 *     Calculate Bfield trajectory
 *
 **************************************************************************/


// Standard C libraries
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <ctime>

// Standard C++ libraries

#include <iostream>
#include <complex>

using namespace std;



// coil libraries
#include "coils.hpp"
#include "coilio.hpp"
#include "coils_cmdline.hpp"
#include "surface.hpp"

// field-line following code
#include "calcPoincare.hpp"

// defines bext, bplasma, bTotal
#include "bfield_ext.hpp"


const double NEGLECT =  1e-12;



// Main Function for code

int main (int argc, char *argv[])
{

  disable_all_options();
  enable_option(opt_NperCycle);
  enable_option(opt_Ncircuits);
  enable_option(opt_xyz);
  enable_option(opt_rzp);

  // parse command line input
 
  
   if (!parse_cmd(argc, argv))
     return 1;


  string fext = "out";
  //  string fext = plasma_extname;
 //  ios_base::fmtflags flags = ios_base::right | ios_base::scientific;
  string ftemp;
  p3vectorformat::textformat(text_nobraces);

  // variables for measuring times
  struct tms tbuff;
  clock_t ckstart;

  // display COOLL mode
  cout << endl;
  display_execution_mode();
  cout << endl;


  dispcr(NperCycle);
  dispcr(Ncircuits);


  // in cylindrical (R,Z,phi) coords
  // starting point should be on the plasma surface
  p3vector<double> startPosition(r0,z0,phi0);

  dispcr(startPosition);
  const unsigned int Npts = NperCycle*Ncircuits;

  LAvector<p3vector<double> > Bpath_rzp(Npts,"Bpath_rzp");
  Bpath_rzp = p3vector<double>(0,0,0);
  Bpath_rzp.textformat(text_nobraces);
  Bpath_rzp.perline(1);

  calcPoincare(Bpath_rzp,startPosition,Ncircuits,NperCycle);

//    LAvector<p3vector<double> > Bpath_xyz(Npts,"Bpath_xyz");
//    Bpath_xyz.textformat(text_nobraces);
//    Bpath_xyz.perline(1);

//    LAvector<p3vector<double> > B(Npts,"B");
//    B.textformat(text_nobraces);
//    B.perline(1);
//    LAvector<double> Bmag(Npts,"Bmag");

//    for (unsigned int i =0; i<Npts; i++) {
//           const double r =   Bpath_rzp[i][0];
//      const double z =   Bpath_rzp[i][1];
//      const double phi =   Bpath_rzp[i][2];

//      Bpath_xyz[i].x() = r*cos(phi);
//      Bpath_xyz[i].y() = r*sin(phi);
//      Bpath_xyz[i].z() = z;

//      bTotal(Bpath_xyz[i],B[i]);
//      Bmag[i] = norm(B[i]);
//    }

//    save(B,"B.xyz.out");
  
//    dispcr(min(Bmag));
//    dispcr(sum(Bmag)/Npts);
//    dispcr(max(Bmag));
  
  //  save(Bpath_xyz,"Bpath.xyz.out");

  save(NperCycle,"NperCycle.out");
  save(Ncircuits,"Ncircuits.out");
  printcr("saving data...");
  save(Bpath_rzp,"Bpath.rzp.out");
  
  return 0;
} // main()












