/************************************************************************* 
 * 
 *   File Name    :  cmdline.cpp
 *   Platform     :  GNU C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  16 May 2003
 * 
 *
 *   SYNOPSIS     
 *     This code parses the input from the unix command line
 *
 **************************************************************************/


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>

//#define CMD_DEBUG 1

#include "coils_cmdline.hpp"


using namespace std;

std::string plasma_filename("surface_plasma.txt");
std::string plasma2_filename("surface2_plasma.txt");
std::string coil_filename("surface_coils.txt");
std::string current_filename("default_current.txt");
std::string current_filename2("");
std::string flux_filename("default_flux.txt");
std::string flux2_filename("default_flux2.txt");
std::string pert_flux_filename("pert_flux.txt");
std::string Bt_filename("BthetaF.txt");
std::string Bp_filename("BphiF.txt");
std::string resonantmode_filename("null.txt");
std::string omega_filename("");
std::string M_filename("");
std::string gamma_filename("");
std::string omegagamma_filename("omegagamma");
std::string lambda_filename("");
std::string JBt_filename("JBthetaF.txt");
std::string JBp_filename("JBphiF.txt");
std::string modes_filename("null.txt");
std::string myname("");
std::string FGo_filename("");
std::string append_filename("");

// DEFAULT INPUT FILENAME FOR PLASMA SURFACE
bool coilfilegiven = false;

// Number of Poloidal and Toroidal grid points
unsigned int Ntheta = 32;  // Poloidal
unsigned int Nphi = 32;    // Toroidal 

unsigned int NperCycle = 32;  
unsigned int Ncircuits = 32;  

unsigned int Nnn = 6;
unsigned int Nmm = 6;

unsigned int Nharm = 1;
unsigned int Mharm = 1;

unsigned int No = 0;
unsigned int N2 = 0;
unsigned int N3 = 0;
unsigned int Na = 0;
unsigned int Nb = 0;
unsigned int Npert = 0;
unsigned int Nprin = 0;
unsigned int Ncomb = 0;

double pos_x=0;
double pos_y=0;
double pos_z=0;

double r0 = 0; 
double z0 = 0;
double phi0 = 0;

double delta = 0.0;

double Itoroidal=0;
double Ipoloidal=0;

double relative_zero = 1e-8;

double alpha_theta = 0;
double alpha_phi = 0;


double perturbation_valueR = 0;
double perturbation_valueZ = 0;

bool scale_flux_basis = false;

bool iota_constraint = false;

double iota_given = 0;

double fluxshear_given  = 0;

double condnum = 0;

bool limitNmodes = false;

unsigned int NnnIF = 0;
unsigned int NmmIF = 0;


bool extradata_flag =false;

bool ignoreLargeModes = false;

const unsigned int MAX_FILES = 999;
unsigned int NUM_FILES = 0;

const unsigned int NUM_OPTIONS = 59;
const unsigned int MAX_PARMS_PER_OPTION = 3;


unsigned int NthetaI=0;  // Poloidal
unsigned int NphiI=0;    // Toroidal 

// Actual text for each command line option
char *opts[NUM_OPTIONS] = {
   "-pf", "-cf", "-Nphi", "-Ntheta", "-nmax", "-mmax", "-if", "-ff", "-No","-pff",
   "-xyz","-Nper","-Ntimes","-pf2","-del","-Itor","-Ipol", "-Btheta", "-Bphi", "-ff2",
   "-N2", "-Nfund", "-Mfund", "-zero", "-alpha_theta", "-alpha_phi", "-rmf", "-omega", "-Mf", "-N3", 
   "+sfx", "-gamma", "-pertvalR",  "-pertvalZ", "-omegagamma", "-Na", "-Nb", "-Npert","-Nprin","-Ncomb",
   "-lambda", "-JBtheta", "-JBphi", "-ic", "-iota", "-modesf", "-nmaxIF", "-mmaxIF", "-rzp", "-fluxshear",
   "-condnum", "-limitNmodes", "-FGof", "-append", "-NphiI", "-NthetaI","-extradata","-if2","-nohighmodes"
};




// you can also enable/disable options at runtime
bool opts_enable[NUM_OPTIONS] = {
   false,false,false,false,false,false,false,false,false,false,
   false,false,false,false,false,false,false,false,false,false, 
   false,false,false,false,false,false,false,false,false,false, 
   false,false,false,false,false,false,false,false,false,false,
   false,false,false,false,false,false,false,false,false,false,
   false,false,false,false,false,false,false,false,false
};


// number of parameters for each option
unsigned int opts_parms[NUM_OPTIONS] = {
   1,1,1,1,1,1,1,1,1,1,
   3,1,1,1,1,1,1,1,1,1,
   1,1,1,1,1,1,1,1,1,1,
   0,1,1,1,1,1,1,1,1,1,
   1,1,1,1,1,1,1,1,3,1,
   1,0,1,1,1,1,0,1,0
};




// text that describes each option
char *opts_help[NUM_OPTIONS] = {
   "<file> = File name of plasma surface fourier data (sin/cos series)\n\t",
   "<file> = File name for coil surface fourier data  (sin/cos series)\n\t",
   "<d> = Number of toroidal angle values to use for surface grid\n\t",
   "<d> = Number of polodial angle values to use for surface grid\n\t",
   "<d> = Number of the highest toroidal mode to use in fourier representation\n\t",
   "<d> = Number of the highest poloidal mode to use in fourier representation\n\t",
   "<file> = File name of coil current fourier data\n\t",
   "<file> = File name for plasma flux fourier data\n\t",
   "<d> = Number of eigevalues to keep in expansion of surface \n\tperturbation matrix (Gamma). If No=0, then all EVs will be kept\n\t",
   "<file> = File name for perturbed (approximated) plasma flux fourier data\n\t",
   "<x>,<y>,<z> = xyz coordinates\n\t",
   "<d> = Number of points to use per circuit\n\t",
   "<d> = Number of times to traverse the torus (in toroidal direction)\n\t",
   "<file> = File name of second plasma surface fourier data (sin/cos series)\n\t",
   "<d> = distance between plasma surface and coil surface\n\t",
   "<d> = toroidal current\n\t",
   "<d> = poloidal current\n\t",
   "<file> = File name for B_theta fourier data (sin/cos series)\n\t",
   "<file> = File name for B_phi fourier data (sin/cos series)\n\t",
   "<file> = File name for  flux fourier data\n\t",
   "<d> = Number of eigevalues to keep in expansion of inductance \n\tperturbation matrix (Omega). \n\tIf N2=0, then N2=No EVs will be kept\n\t",
   "<d> = fundamental toroidal mode (i.e. lowest non-zero mode) to use in fourier representation\n\tOnly harmonics of this number will be used in the expansion.\n\t",
   "<d> = fundamental poloidal mode (i.e. lowest non-zero mode) to use in fourier representation\n\tOnly harmonics of this number will be used in the expansion.\n\t",
   "<d> = 'relative zero'. Any mode which is this magnitude less than the strongest mode is set to zero.\n\t",
   "<d> = exponential damping coef for poloidal modes in the rho matrix\n\t",
   "<d> = exponential damping coef for toroidal modes in the rho matrix\n\t",
   "<file> = File name for resonant mode data\n\t",
   "<file> = Root File name to load omega matrices from.\n\tIf this string is empty, then omega is directly calculated instead.\n\t",
   "<file> = File name to load M matrix (inductance) from.\n\tIf this string is empty, then M is directly calculated instead.\n\t",
   "<d> = Number of eigevalues to keep in expansion of \n\tinductance matrix (M). If N3=0, then all EVs will be kept\n\t",
   "when this option is given, then each vector of teh flux basis is scaled by the corresponding singular value\n\t",
   "<file> = Root File name to load gamma matrices from.\n\tIf this string is empty, then gamma is directly calculated instead.\n\t",
   "<d> = size of perturbation in meters\n\t",
   "<d> = size of perturbation in meters\n\t",
   "<file> = Root File name to load omega-gamma matrices from.\n\t",
   "<d> = start number for OmegaGamma eigenvectors to keep\n\t",
   "<d> = finish number for OmegaGamma eigenvectors to keep\n\t",
   "<d> = number of the most sensitive perturbation modes to retain\n\t",
   "<d> = number of the pricipal modes to keep (i.e. largest magnitude)\n\t",
   "<d> = number of the combined modes (pert*prin) to keep (i.e. largest magnitude)\n\t",   "<file> = File name for lambda function spectrum\n\t",
   "<file> = File name for JB_theta fourier data (sin/cos series)\n\t",
   "<file> = File name for JB_phi fourier data (sin/cos series)\n\t",
   "<bool> = use iota constraint (true/false)\n\t",
   "<d> = iota\n\t",
   "<file> = File name that contains flux mode indicies to retain\n\t",
   "<d> = Number of the highest toroidal mode to use in fourier representation of I\n\t",
   "<d> = Number of the highest poloidal mode to use in fourier representation of I\n\t",
   "<r>,<z>,<p> = {r,z,phi} coordinates\n\t",
   "<d> = dpsi_t/dr = flux shear\n\t",
   "<d> = condition number used to define matrix singularity\n\t",
   "when this option is given, the number of Fourier-Green modes used is limited to number of magnetic modes specified.",
   "<file> = Root File name to load orthogonal FG functions and J and U matrices from.\n\tIf this string is empty, then it is directly calculated instead.\n\t",
   "<file> = File name to append data to\n\t",
   "<d> = Number of toroidal angle values to use for COIL surface grid\n\t",
   "<d> = Number of polodial angle values to use for COIL surface grid\n\t",
   "signifies that you want verbose data calculated and saved",
   "<file> = File name of 2nd set of coil current fourier data (optional)\n\t",
   "signifies that you want to ignore all flux modes whose mode number is outside of the Fourier grid used for the current.",
};


// text that for labeling the parameters of each option
char *opts_parm_text[NUM_OPTIONS][MAX_PARMS_PER_OPTION] = {
   {" <file>"},
   {" <file>"},
   {" <d>"},
   {" <d>"},
   {" <d>"},
   {" <d>"},
   {" <file>"},
   {" <file>"},
   {" <d>"},
   {" <file>"},
   {"<x>", " <y>", " <z>"},
   {" <d>"},
   {" <d>"},
   {" <file>"},
   {" <d>"},
   {" <d>"},
   {" <d>"},
   {" <file>"},
   {" <file>"},
   {" <file>"},
   {" <d>"},
   {" <d>"},
   {" <d>"},
   {" <d>"},
   {" <d>"},
   {" <d>"},
   {" <file>"},
   {" <file>"},
   {" <file>"},
   {" <d>"},
   {""},
   {" <file>"},
   {" <d>"},
   {" <d>"},
   {" <file>"},
   {" <d>"},
   {" <d>"},
   {" <d>"},
   {" <d>"},
   {" <d>"},
   {" <file>"},
   {" <file>"},
   {" <file>"},
   {" <bool>"},
   {" <d>"},
   {" <file>"},
   {" <d>"},
   {" <d>"},
   {"<r>", " <z>", " <p>"},
   {" <d>"},
   {" <d>"},
   {""},
   {" <file>"},
   {" <file>"},
   {" <d>"},
   {" <d>"},
   {""},
   {" <file>"},
   {""},
};


// Synopsis for the command
std::string synopsis("");




void disect_filename(const std::string fn,std::string& root,std::string& extension) {

      // extract head and tail of filenames
   uint j = fn.rfind(".");
   root = fn.substr(0,j);

   uint last = fn.length()-1;
   extension = fn.substr(j+1,last);

//    cout << "fn="<<fn;
//    cout << "j="<<j;
//    cout << "last="<<last;
//    cout << "root="<<root;
//    cout << "extension="<<fn;


}








//***********************************************************



/*************************************************************************
 *
 *   FUNCTION     : opt_func
 *   INPUTS       :  
 *     opt - option type
 *     parms - array of parameter strings that were input with the given option
 *   OUTPUTS      :  
 *     none
 *   I/O          :
 *    none
 *
 *   SYNOPSIS
 *
 *   This function is called by the command line parsing engine.
 *   It is called once for each option that is parsed. This function
 *   then performs the appropriate action for the option.
 *
 **************************************************************************/

void opts_func(t_opts opt, char *parms[])
{
 
#ifdef CMD_DEBUG
   {
      int debug;
      printf("found: %s",opts[opt]);
      for(debug=0; debug < opts_parms[opt]; debug++)
	 printf(" '%s'",parms[debug]);
      printf("\n");
   }
#endif

   string s;
   const int Np= opts_parms[opt];
   for(int j=0; j<Np; j++) {
      s = s + string(parms[j]);
      if (j <  Np-1)
	 s += " ";
   }
      
   std::istringstream strmline(s);

   switch(opt) {
   case opt_pf:
      plasma_filename = s;
      break;
   case opt_cf:
      coil_filename = s;
      coilfilegiven =true;
      break;
   case opt_Ntheta:
      strmline >> Ntheta;
      break;
   case opt_Nphi:
      strmline >> Nphi;
      break;
   case opt_Nnn:
      strmline >> Nnn;
      break;
   case opt_Nmm:
      strmline >> Nmm;
      break;
   case opt_if:
      current_filename = s;
      break;
   case opt_ff:
      flux_filename = s;
      break;
   case opt_No:
      strmline >> No;
      break;
   case opt_pff:
      pert_flux_filename = s;
      break;
   case opt_xyz:
      strmline >> pos_x >> pos_y >> pos_z;
      r0 = sqrt(pos_x*pos_x + pos_y*pos_y);
      z0 = pos_z;
      phi0 = atan2(pos_y,pos_x);
      break;
   case opt_NperCycle:
      strmline >> NperCycle;
      break;
   case opt_Ncircuits:
      strmline >> Ncircuits;
      break;
   case opt_pf2:
      plasma2_filename = s;
      break;
   case opt_del:
      strmline >> delta;
      break;
   case opt_Itoroidal:
      strmline >> Itoroidal;
      break;
   case opt_Ipoloidal:
      strmline >> Ipoloidal;
      break;
   case opt_Btf:
      Bt_filename = s;
      break;
   case opt_Bpf:
      Bp_filename = s;
      break;
   case opt_ff2:
      flux2_filename = s;
      break;
   case opt_N2:
      strmline >> N2;
      break;
   case opt_Nharm:
      strmline >> Nharm;
      break;
   case opt_Mharm:
      strmline >> Mharm;
      break;
   case opt_zero:
      strmline >> relative_zero;
      break;
   case opt_alpha_theta:
      strmline >> alpha_theta;
      break;
   case opt_alpha_phi:
      strmline >> alpha_phi;
      break;
   case opt_rmf:
      strmline >> resonantmode_filename;
      break;
   case opt_omegaf:
      strmline >> omega_filename;
      break;
   case opt_Mf:
      strmline >> M_filename;
      break;
   case opt_N3:
      strmline >> N3;
      break;
   case opt_sfx:
      scale_flux_basis = true;
      break;
   case opt_gammaf:
      strmline >> gamma_filename;
      break;
   case opt_pertvalR:
      strmline >>perturbation_valueR;
      break;
   case opt_pertvalZ:
      strmline >>perturbation_valueZ;
      break;
   case opt_ogf:
      strmline >> omegagamma_filename;
      break;
   case opt_Na:
      strmline >> Na;
      break;
   case opt_Nb:
      strmline >> Nb;
      break;
   case opt_Npert:
      strmline >> Npert;
      break;
   case opt_Nprin:
      strmline >> Nprin;
      break;
   case opt_Ncomb:
      strmline >> Ncomb;
      break;
   case opt_lambdaf:
      strmline >> lambda_filename;
      break;
   case opt_JBtf:
      JBt_filename = s;
      break;
   case opt_JBpf:
      JBp_filename = s;
      break;
   case opt_ic:
      if ( (s[0] == 't') || (s[0] == '1') ) 
	 iota_constraint = true;
      else if ( (s[0] == 'f') || (s[0] == '0') )
	 iota_constraint = false;   
      else
	 iota_constraint = false;   
      break;
   case opt_iota:
      strmline >> iota_given;
      break;
   case opt_modesf:
      modes_filename = s;
      break;
   case opt_NnnIF:
      strmline >> NnnIF;
      break;
   case opt_NmmIF:
      strmline >> NmmIF;
      break;
   case opt_rzp:
      strmline >> r0 >> z0 >> phi0;
      break;
   case opt_fluxshear:
      strmline >> fluxshear_given;
      break;
   case opt_condnum:
      strmline >> condnum;
      break;
   case opt_limitNmodes:
      limitNmodes = true;
      break;
   case opt_FGof:
      strmline >> FGo_filename;
      break;
   case opt_append:
      strmline >> append_filename;
      break;
   case opt_NthetaI:
      strmline >> NthetaI;
      break;
   case opt_NphiI:
      strmline >> NphiI;
      break;
   case opt_extradata:
      extradata_flag = true;
      break;
   case opt_if2:
      current_filename2 = s;
      break;
   case opt_ignoreLargeModes:
      ignoreLargeModes = true;
      break;
   default:
      std::cout<<"huh? at line# (BUG ENCOUNTERED)"<<__LINE__<<std::endl;
      break;

   } // switch(opt)

  
}




/*-------------------------COMMAND PARSING ENGINE-----------------------*/

/************************************************************************
 This code can parse out all your options and filenames from the command line.
 Options must start with either a '+' or a '-'.
************************************************************************/


void disable_option(t_opts opt) {
   opts_enable[opt] = false;
}

void enable_option(t_opts opt) {
   opts_enable[opt] = true;
}

void disable_all_options(void) {
   for(unsigned int opt = 0; opt<NUM_OPTIONS; opt++)
      opts_enable[opt] = false;
}


int opts_found[NUM_OPTIONS];


void help()
{

   cerr<<"------------------------------------------------------------------------------"<<endl; 
   cerr<<"Help for "<<myname<<" "<<endl;
   cerr<<endl<<"Usage:"<<endl;
   cerr<<" "<<myname<<" [-options]";
   if (NUM_FILES == MAX_FILES)
      cerr<<" file1 file2 ... fileN";
   else
      for(unsigned int i=0; i<NUM_FILES; i++)
	 cerr<<" file"<<(i+1);

   cerr<<endl;

   cerr<<endl<<"Synopsis:"<<endl<<synopsis<<endl;
    
   cerr<<"Options:"<<endl;
   cerr<<" -h,--help"<<endl<<"\tDisplays this help screen."<<endl;
   for (unsigned int i=0; i< NUM_OPTIONS; i++) {

      if (opts_enable[i]) {
	 cerr<<" "<<opts[i];
      
	 if (opts_parms[i] >1)
	    cerr<<" \"";

	 for(unsigned int j=0; j<opts_parms[i]; j++) {
	    if (j>0)
	       cerr<<",";
	    cerr<<(opts_parm_text[i][j]);
	 }
	 if (opts_parms[i] >1)
	    cerr<<"\"";


	 string preface;
	 if (opts_found[i])
	    preface = string("( you entered: ");
	 else
	    preface = string("( default is ");

	 cerr <<endl<<"\t"<< opts_help[i] << preface;

	 switch(i) {
	 case opt_pf:
	    cerr << plasma_filename;
	    break;
	 case opt_cf:
	    cerr << coil_filename;
	    break;
	 case opt_Ntheta:
	    cerr << Ntheta;
	    break;
	 case opt_Nphi:
	    cerr <<Nphi;
	    break;
	 case opt_Nnn:
	    cerr << Nnn;
	    break;
	 case opt_Nmm:
	    cerr <<Nmm;
	    break;
	 case opt_if:
	    cerr <<current_filename;
	    break;
	 case opt_ff:
	    cerr <<flux_filename;
	    break;
	 case opt_No:
	    cerr<<No;
	    break;
	 case opt_pff:
	    cerr << pert_flux_filename;
	    break;
	 case opt_xyz:
	    cerr << "x = "<< pos_x<<", y = "<<pos_y<<", z = "<<pos_z;
	    break;
	 case opt_NperCycle:
	    cerr<<NperCycle;
	    break;
	 case opt_Ncircuits:
	    cerr<<Ncircuits;
	    break;
	 case opt_pf2:
	    cerr  <<plasma2_filename;
	    break;
	 case opt_del:
	    cerr <<delta;
	    break;
	 case opt_Itoroidal:
	    cerr << Itoroidal;
	    break;
	 case opt_Ipoloidal:
	    cerr <<Ipoloidal;
	    break;
	 case opt_Btf:
	    cerr << Bt_filename;
	    break;
	 case opt_Bpf:
	    cerr <<  Bp_filename;
	    break;
	 case opt_ff2:
	    cerr << flux2_filename;
	    break;
	 case opt_N2:
	    cerr <<N2;
	    break;
	 case opt_Nharm:
	    cerr <<Nharm;
	    break;
	 case opt_Mharm:
	    cerr <<Mharm;
	    break;
	 case opt_zero:
	    cerr << relative_zero;
	    break;
	 case opt_alpha_theta:
	    cerr << alpha_theta;
	    break;
	 case opt_alpha_phi:
	    cerr << alpha_phi;
	    break;
	 case opt_rmf:
	    cerr << resonantmode_filename;
	    break;
	 case opt_omegaf:
	    cerr << omega_filename;
	    break;
	 case opt_Mf:
	    cerr<<  M_filename;
	    break;
	 case opt_N3:
	    cerr <<N3;
	    break;
	 case opt_sfx:
	    cerr << scale_flux_basis;
	    break;
	 case opt_gammaf:
	    cerr << gamma_filename;
	    break;
	 case opt_pertvalR:
	    cerr << perturbation_valueR;
	    break;
	 case opt_pertvalZ:
	    cerr << perturbation_valueZ;
	    break;
	 case opt_ogf:
	    cerr << omegagamma_filename;
	    break;
	 case opt_Na:
	    cerr <<Na;
	    break;
	 case opt_Nb:
	    cerr <<Nb;
	    break;
	 case opt_Npert:
	    cerr <<Npert;
	    break;
	 case opt_Nprin:
	    cerr <<Nprin;
	    break;
	 case opt_Ncomb:
	    cerr <<Ncomb;
	    break;
	 case opt_lambdaf:
	    cerr << lambda_filename;
	    break;
	 case opt_JBtf:
	    cerr << JBt_filename;
	    break;
	 case opt_JBpf:
	    cerr <<  JBp_filename;
	    break;
	 case opt_ic:
	    if (iota_constraint) 
	       cerr<< "true";
	    else
	       cerr<< "false";
	    break;
	 case opt_iota:
	    cerr << iota_given;
	    break;
	 case opt_modesf:
	    cerr << modes_filename;
	    break;
	 case opt_NnnIF:
	    cerr << NnnIF;
	    break;
	 case opt_NmmIF:
	    cerr <<NmmIF;
	    break;
	 case opt_rzp:
	    cerr << "r = "<< r0<<", z = "<<z0<<", p = "<<phi0;
	    break;
	 case opt_fluxshear:
	    cerr << fluxshear_given;
	    break;
	 case opt_condnum:
	    cerr << condnum;
	    break;
	 case opt_limitNmodes:
	    if (limitNmodes) 
	       cerr<< "true";
	    else
	       cerr<< "false";
	    break;
	 case opt_FGof:
	    cerr << FGo_filename;
            break;
	 case opt_append:
	    cerr << append_filename;
            break;
	 case opt_NthetaI:
	    cerr << NthetaI;
	    break;
	 case opt_NphiI:
	    cerr <<NphiI;
	    break;
	 case opt_extradata:
	    if (extradata_flag) 
	       cerr<< "true";
	    else
	       cerr<< "false";
	    break;
	 case opt_if2:
	    cerr <<current_filename2;
	    break;
	 case opt_ignoreLargeModes:
	    if (ignoreLargeModes) 
	       cerr<< "true";
	    else
	       cerr<< "false";
	    break;
	 default:
	    std::cout<<"bug? at line#"<<__LINE__<<" in "<< __FILE__<<std::endl;
	    break;
	 } // switch(opt)
	 cerr << " )"<<endl;
      }
   }
   cerr<<  endl<<"Notes:"<< endl<< "Options may be placed in any order." <<endl;

   cerr<<endl<<"------------------------------------------------------------------------------"<<endl; 

}




unsigned  int nfiles;
char *filenames[255];


bool parse_cmd (int argc, char *argv[])
{
   bool not_option;
   bool call_help =false;
   myname = std::string(argv[0]);
   nfiles = 0;


   for (unsigned int j=0; j< NUM_OPTIONS; j++)
      opts_found[j] = 0;

   for(unsigned int i=1; i< static_cast<unsigned int>(argc); i++){
      not_option = true;
      //    printf("argv[%d] = %s"<<endl,i,argv[i]);
      for(unsigned int j=0; j<NUM_OPTIONS; j++){  

	 if (opts_enable[j]) {

	    if (opts_parms[j] == 0){
	       /* check for  option without parameters */
	       if ( (strcmp(argv[i],opts[j])) == 0) {
		  //	    printf("no parms for option[%d]"<<endl<<"",i);
		  not_option = false;
		  opts_found[j] = 1;
		  opts_func((t_opts)j,NULL);
		  break;
	       }
	    } else
	       {
		  /* check for option that has parmeters*/
		  if ( ( strncmp(argv[i],opts[j],strlen(opts[j])) == 0) && 
		       ( strlen(opts[j]) == strlen(argv[i]) ) )  {
		     char stemp[80];
		     char *sptr = argv[i] + strlen(opts[j]);
		     char *sparms[5];

		     //	      printf("parms for option[%d] = %d "<<endl<<"",i,opts_parms[j]);

		     strncpy(stemp,sptr,78);
		     strcat(stemp,",");
		     sptr = stemp;

		     if ((sparms[0] = strtok(sptr,",")) == NULL) {
			sptr = argv[++i];
			if ((sparms[0] = strtok(sptr,",")) == NULL) {
			   cerr<<myname<<": Bad input (not enough parameters) for option "<<opts[j]<<endl;
			   return false;
			}
		     }
                     unsigned int k=0;
		     for(k=1;k< opts_parms[j] ; k++){
			if ((sparms[k] = strtok(NULL,",")) == NULL) {
			   cerr<<myname<<": Bad input (not enough parameters) for option "<<opts[j]<<endl;
			   return false;
			}
		     } 
		     if ((sparms[k] = strtok(NULL,",")) != NULL) {
			cerr<<myname<<": Bad input (too many parameters) for option "<<opts[j]<<endl;
			return false;
		     }

		     opts_found[j] = 1;
		     not_option = false;
		     opts_func((t_opts)j,sparms);
		     break;
		  } 
	       }  /* else (check for opts with parms) */
	 }
      } /*for(j=o) */

      //    printf("not_option[%d] = %d"<<endl<<"",i,not_option);
    
      if (not_option) {
	 if ( ( strcmp(argv[i],"-h") == 0 ) || ( strcmp(argv[i],"-help") == 0 ) || ( strcmp(argv[i],"--help") == 0 ) ) {
	    call_help =true;
	 } else if ( ( argv[i][0] == '-' ) || ( argv[i][0] == '+' )) {
	    cerr<<endl<< myname<<": Illegal option: "<< argv[i]<<endl<<endl;
	    return false;
	 } else {
	    filenames[nfiles++] = argv[i];
	 }
      }
    
   } /* for(i=0) */

  
  

#ifdef CMD_DEBUG
   {
      int debug;
      cerr<<"file list:";
      for(debug=0; debug < nfiles; debug++)
	 cerr<<" '"<<filenames[debug]<<"'";
      cerr<<endl;
   }
#endif

   if ((NUM_FILES != MAX_FILES) && (nfiles != NUM_FILES)){
      cerr<<endl<<myname<<": Expecting ("<< NUM_FILES<<") filenames on the command line"<<endl;
      return false;
   }



   // dynamic default values
   if (NthetaI==0)
     NthetaI = Ntheta;
   if (NphiI==0)
     NphiI= Nphi;


   if (call_help) {
      help();
      return false;
   }


   return true;
} // parse_cmd()


















   // create coil surface filename from plasma surface filename (if not given)
 //   if (!coilfilegiven) {
//       using namespace std;
//       string s = plasma_filename;
//       string substr = "plasma";
//       uint i = s.find(substr);
//       if (i==(string::npos)) {
// 	 std::cout <<"A" <<endl;
// 	 uint j = s.rfind(".");
// 	 string head = s.substr(0,j);
// 	 uint last = s.length()-1;
// 	 string tail = s.substr(j,last);
// 	 s = head + ".coils"+tail;
//       }else{
// 	 string head = s.substr(0,i);
// 	 uint last = s.length()-1;
// 	 uint Nsubstr = substr.length();
// 	 string tail = s.substr(i+Nsubstr,last);
// 	 s = head + "coils" + tail;
//       }
//       coil_filename = s;
    
