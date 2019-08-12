


/************************************************************************* 
 * 
 * Command Line Definitions
 *
 **************************************************************************/

#include <string>
#include <iostream>


// INPUT FILENAME FOR PLASMA & COIL SURFACES
extern std::string plasma_filename;
extern std::string plasma2_filename;
extern std::string coil_filename;
extern std::string current_filename;
extern std::string current_filename2;
extern std::string flux_filename;
extern std::string pert_flux_filename;
extern std::string Bt_filename;
extern std::string Bp_filename;
extern std::string flux2_filename;
extern std::string resonantmode_filename;
extern std::string omega_filename;
extern std::string gamma_filename;
extern std::string M_filename;
extern std::string omegagamma_filename;
extern std::string lambda_filename;
extern std::string JBt_filename;
extern std::string JBp_filename;
extern std::string modes_filename;
extern std::string myname;
extern std::string FGo_filename;
extern std::string append_filename;


extern unsigned int Ntheta;  // Poloidal
extern unsigned int Nphi;    // Toroidal 

extern unsigned int NperCycle;  
extern unsigned int Ncircuits;    

extern unsigned int No;    // eigenvalues to keep  in Gamma
extern unsigned int Nnn;
extern unsigned int Nmm;
extern double pos_x;
extern double pos_y;
extern double pos_z;
extern double delta;

extern double r0;
extern double z0;
extern double phi0;

extern double Itoroidal;
extern double Ipoloidal;

extern unsigned int N2;    // eigenvalues to keep in Omega
extern unsigned int N3;    // eigenvalues to keep in M

extern unsigned int Nharm;
extern unsigned int Mharm;

extern double relative_zero;
extern double alpha_theta;
extern double alpha_phi;

extern double perturbation_valueR;
extern double perturbation_valueZ;

extern bool scale_flux_basis;

extern unsigned int Na;
extern unsigned int Nb;

extern unsigned int Npert;
extern unsigned int Nprin;
extern unsigned int Ncomb;

extern bool iota_constraint;

extern double iota_given;
extern unsigned int NnnIF;
extern unsigned int NmmIF;

extern double fluxshear_given;
extern double condnum;
extern bool limitNmodes;

extern unsigned int NthetaI;  // Poloidal
extern unsigned int NphiI;    // Toroidal 


extern bool extradata_flag;

extern bool ignoreLargeModes;



// enumerated types for option
typedef enum {
   opt_pf=0, opt_cf, opt_Nphi, opt_Ntheta, opt_Nnn, opt_Nmm, opt_if, opt_ff, opt_No, opt_pff,
   opt_xyz, opt_NperCycle, opt_Ncircuits, opt_pf2, opt_del, opt_Itoroidal, opt_Ipoloidal, opt_Btf, opt_Bpf, opt_ff2, 
   opt_N2, opt_Nharm, opt_Mharm, opt_zero, opt_alpha_theta, opt_alpha_phi, opt_rmf, opt_omegaf, opt_Mf, opt_N3,
   opt_sfx, opt_gammaf, opt_pertvalR, opt_pertvalZ, opt_ogf, opt_Na, opt_Nb, opt_Npert, opt_Nprin, opt_Ncomb, 
   opt_lambdaf, opt_JBtf, opt_JBpf, opt_ic, opt_iota, opt_modesf, opt_NnnIF, opt_NmmIF, opt_rzp, opt_fluxshear,
   opt_condnum, opt_limitNmodes, opt_FGof, opt_append, opt_NphiI, opt_NthetaI, opt_extradata, opt_if2, opt_ignoreLargeModes
} t_opts;

void disable_option(t_opts opt);
void disable_all_options(void);
void enable_option(t_opts opt);

bool parse_cmd (int argc, char *argv[]);



extern std::string synopsis;


void disect_filename(const std::string fn,std::string& root,std::string& extension);

