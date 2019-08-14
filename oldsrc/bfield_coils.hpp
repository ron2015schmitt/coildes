
#include "matricks.hpp"



// bplasma is defined in bfield_plasma_p%.cpp
void bplasma(const  Matricks::p3vector<double>& X,  Matricks::p3vector<double>& B);


// defined in bfield_coils.cpp
void bCoils(const p3vector<double>& X, p3vector<double>& Bc);

void bTotal(const p3vector<double>& X, p3vector<double>& B) ;
 
void bTotalandbCoils(const p3vector<double>& X, p3vector<double>& B, p3vector<double>& Bc);


int setupcoils(std::string coil_filename, std::string current_filename, 
	       const double Itoroidal, const double Ipoloidal, 
	       unsigned int Nnn, unsigned int Nmm, 
	       unsigned int Nphi, unsigned int Ntheta, 
	       unsigned int Nfund, unsigned int Mfund);

