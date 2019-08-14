/************************************************************************* 
 * 
 *   File Name    : 
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  
 * 
 *
 *   SYNOPSIS     
 *
 * VERSION NOTES:
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
#include "coilfft.hpp"

// This is the file that defines the B field configuration



// Main Function for code

int main (int argc, char *argv[])
{

  disable_all_options();
  enable_option(opt_pf);

  // parse command line input
 
  if (!parse_cmd(argc, argv))
     return 1;


  // 
  //  ios_base::fmtflags flags = ios_base::right | ios_base::scientific;
  p3vectorformat::textformat(text_braces);

  Vector <double> datavec("datavec");
  datavec.textformat(text_braces);

  string fname;

  load(datavec,plasma_filename);

  string ext("out");
  disect_filename(plasma_filename,fname,ext);

  fname += ".nobraces.out";

  p3vectorformat::textformat(text_nobraces);
  datavec.textformat(text_nobraces);
  datavec.perline(1);

  save(datavec,fname);

  return 0;
} // main()





