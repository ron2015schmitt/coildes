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

#include "bfield_coils.hpp"

///////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////



const double NEGLECT =  1e-12;



// Main Function for code

int main (int argc, char *argv[])
{

  disable_all_options();
  enable_option(opt_cf);
  enable_option(opt_if);
  enable_option(opt_Nphi);
  enable_option(opt_Ntheta);
  enable_option(opt_Nnn);
  enable_option(opt_Nmm);
  enable_option(opt_xyz);
  enable_option(opt_rzp);
  enable_option(opt_NperCycle);
  enable_option(opt_Ncircuits);
  enable_option(opt_Itoroidal);
  enable_option(opt_Ipoloidal);
  enable_option(opt_Nharm);
  enable_option(opt_Mharm);
  enable_option(opt_append);
 
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

  ostringstream strmtemp;
   string fname;

  string fext = "out";
  //  string fext = plasma_extname;
 //  ios_base::fmtflags flags = ios_base::right | ios_base::scientific;
  string ftemp;
  p3vectorformat::textformat(text_nobraces);

  // variables for measuring times
//  struct tms tbuff;
//  clock_t ckstart;

  // display COOLL mode
  cout << endl;
  display_execution_mode();
  cout << endl;

  bool append_data;
  if (append_filename == "") 
     append_data = false;
  else
     append_data =true;


  setupcoils(coil_filename, current_filename, Itoroidal, Ipoloidal, 
	     Nnn, Nmm, Nphi, Ntheta, Nharm, Mharm);


  dispcr(Itoroidal);
  dispcr(Ipoloidal);
  dispcr(NperCycle);
  dispcr(Ncircuits);

  static int Ncircuits_planned = Ncircuits;  

  // in cylindrical (R,Z,phi) coords
  // starting point should be on the plasma surface
  p3vector<double> startPosition(r0,z0,phi0);

  unsigned int Np1 = 0;
  unsigned int Nc1 = 0;
  if (append_data) {
      //   const unsigned int LINESZ = 1024;
     string line;
     std::ifstream in;
     in.open(append_filename.data(),std::ios::in);
     if (!in) {
	cant_open_file(append_filename);
	return 1;
     }
     getline(in,line);
     //     dispcr(line);
     in.close();

     sscanf(line.data(),"%% NperCycle = %d; Ncircuits = %d",&Np1,&Nc1);

     cout<< "According to file header, '"<<append_filename<<"'" <<endl;
     cout << "has "<<Nc1<< " circuits of data at "<<Np1<<" datapoints per circuit."<<endl;
     if (Np1 != NperCycle) {
	cerr<< "ERROR: Number of datapoints per cycle at command line ("<<NperCycle<<") does not match that in file("<<Np1<<")."<<endl;
	print(myname);printcr(": aborting.");
	return (2);
     }

     if (Nc1 >= Ncircuits) {
	cout<< "Desired number of cycles ("<<Ncircuits << ") already calculated!"<<endl;
	return (0);
     }

     p3vector<double> PreviousStopPosition = startPosition;

     // get the last datapoint, which overrides command line position
     if (Nc1 > 0) {


	LAvector<p3vector<double> > Bpath_rzp(0,"Bpath_rzp");
	Bpath_rzp = p3vector<double>(0,0,0);
	Bpath_rzp.textformat(text_nobraces);
	Bpath_rzp.perline(1);
	if (load(Bpath_rzp,append_filename)) {
	   printcr("Above ERROR occurred in "+myname+".");
	   return 4;
	}

 	if (Bpath_rzp.size()!=(Np1*Nc1)) {
 	   cerr<< "ERROR: Number of datapoints in file ("<<Bpath_rzp.size()<<") does not match values given in file header, (NperCycle = "<<Np1<<") * (Ncircuits = "<<Nc1<<") = "<<(Np1*Nc1)<<"."<<endl;
 	   print(myname);printcr(": aborting.");
 	   return (5);
 	}
	PreviousStopPosition = Bpath_rzp[Bpath_rzp.size()-1];

	strmtemp <<"% NperCycle = "<<Np1 << "; Ncircuits = "<<Nc1<<";";
	string preamble_string(strmtemp.str());
	strmtemp.str("");
	//	dispcr(preamble_string);
	fname = "Bpath.fromcoils.temp.rzp.out";
	printcr("Copying existing data to:");
	dispcr(fname);
	Bpath_rzp.textformat(text_nobraces);
	Bpath_rzp.perline(1);
	save(Bpath_rzp,fname,std::ios::trunc,preamble_string);

// 	in.open(append_filename.data(),std::ios::in);
// 	if (!in) {
// 	   cant_open_file(append_filename);
// 	   return 1;
// 	}
// 	int cntr=0;
// 	while (in.good()) {
// 	   if (!(in >> PreviousStopPosition).fail())
// 	      cntr++;
// 	}
// 	in.close();
// 	if (cntr!=(Np1*Nc1)) {
// 	   cerr<< "ERROR: Number of datapoints in file ("<<cntr<<") does not match values given in file header, (NperCycle = "<<Np1<<") * (Ncircuits = "<<Nc1<<") = "<<(Np1*Nc1)<<"."<<endl;
// 	   print(myname);printcr(": aborting.");
// 	   return (3);
// 	}

     }
     dispcr(PreviousStopPosition);
     if ( ( PreviousStopPosition[0] == 0) && ( PreviousStopPosition[1] == 0) && ( PreviousStopPosition[2] == 0) ) {
	cerr<< "ERROR: Previous Stopping Position is all zeros."<<endl;
	print(myname);printcr(": aborting.");
	return (6);
	
     }
 
     cout<<"Data will be appended to '"<<append_filename<<"'."<<endl;
     p3vector<double> position(0,0,0);
     calc_1pt_Poincare(position,PreviousStopPosition);
     startPosition = position;

     Ncircuits_planned = Ncircuits - Nc1;
  }

  

  dispcr(startPosition);
  const unsigned int Npts = NperCycle*Ncircuits_planned;

  LAvector<p3vector<double> > Bpath_rzp(Npts,"Bpath_rzp");
  Bpath_rzp = p3vector<double>(0,0,0);
  Bpath_rzp.textformat(text_nobraces);
  Bpath_rzp.perline(1);

  cout << "Calculating "<<Ncircuits_planned<<" circuits of datapoints..."<<endl;
  int Ncircuits_actual = Ncircuits_planned;
  if (calcPoincare(Bpath_rzp,startPosition,Ncircuits_actual,NperCycle) > 1) {
  } 
  dispcr(Ncircuits_actual);

//   LAvector<p3vector<double> > Bpath_xyz(Npts,"Bpath_xyz");
//   Bpath_xyz.textformat(text_nobraces);
//   Bpath_xyz.perline(1);

//    LAvector<p3vector<double> > B(Npts,"B");
//    B.textformat(text_nobraces);
//    B.perline(1);
//    LAvector<p3vector<double> > Brzp(Npts,"Brzp");
//    Brzp.textformat(text_nobraces);
//    Brzp.perline(1);
//    LAvector<double> Bmag(Npts,"Bmag");
//    Bmag.textformat(text_nobraces);
//    Bmag.perline(1);

//    for (unsigned int i =0; i<Npts; i++) {
//      const double r =  Bpath_rzp[i][0];
//      const double z =  Bpath_rzp[i][1];
//      const double phi =  Bpath_rzp[i][2];
     
//      p3vector<double> X;
//      X.x() = r*cos(phi);
//      X.y() = r*sin(phi);
//      X.z() = z;


//      bTotal(X,B[i]);
//      Bmag[i] = norm(B[i]);
//      Brzp[i][0] = B[i].x()*cos(phi) + B[i].y()*sin(phi);
//      Brzp[i][1] = B[i].z();
//      Brzp[i][2] = -B[i].x()*sin(phi) + B[i].y()*cos(phi);
     
//      //     Bpath_xyz[i]= X;

//    }
  
//   dispcr(min(Bmag));
//   dispcr(sum(Bmag)/N);
//   dispcr(max(Bmag));
  
  //  save(Bpath_xyz,"Bpath.fromcoils.xyz.out");

  const unsigned int Npts_actual = NperCycle*Ncircuits_actual;
  if (Ncircuits_actual != Ncircuits_planned) {
     dispcr(Npts_actual);
     LAvector<p3vector<double> > Bpath_temp(Npts_actual,"Bpath_temp");
     for(unsigned int i = 0; i < Npts_actual; i++) {
	Bpath_temp[i] = Bpath_rzp[i];
     }
     Bpath_rzp.resize() = Bpath_temp;
  }



  strmtemp <<"% NperCycle = "<<NperCycle << "; Ncircuits = "<<Ncircuits_actual<<";";
  string preamble_string(strmtemp.str());
  strmtemp.str("");
  //  dispcr(preamble_string);
  Bpath_rzp.textformat(text_nobraces);
  Bpath_rzp.perline(1);
  fname = "Bpath.fromcoils.rzp.out";
  dispcr(fname);
  save(Bpath_rzp,fname,std::ios::trunc,preamble_string);
  
  if (append_data) {
     // load in previous data points
     Bpath_rzp.clear();
     fname = "Bpath.fromcoils.temp.rzp.out";
     if (load(Bpath_rzp,fname)) {
	printcr("Above ERROR occurred in "+myname+".");
	return 7;
     }

     // now save previous data points with new header
     strmtemp <<"% NperCycle = "<<NperCycle << "; Ncircuits = "<<(Ncircuits_actual+Nc1)<<";";
     string preamble_string(strmtemp.str());
     strmtemp.str("");
     Bpath_rzp.textformat(text_nobraces);
     Bpath_rzp.perline(1);
     save(Bpath_rzp,append_filename,std::ios::trunc,preamble_string);

     // load back new data points
     Bpath_rzp.clear();
     fname = "Bpath.fromcoils.rzp.out";
     if (load(Bpath_rzp,fname)) {
	printcr("Above ERROR occurred in "+myname+".");
	return 9;
     }

     // now append the new data to file
     Bpath_rzp.textformat(text_nobraces);
     Bpath_rzp.perline(1);
     save(Bpath_rzp,append_filename,std::ios::app);

  }

  //  save(NperCycle,"NperCycle.out");
  //  save(Ncircuits_actual+Nc1,"Ncircuits.out");
  
//   save(Bmag,"Bmag.intrajfromcoils.out");
//   save(B,"B.intrajfromcoils.out");
//   save(Brzp,"Brzp.intrajfromcoils.out");
  
  return 0;
} // main()





