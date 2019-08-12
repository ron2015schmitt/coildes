



#ifndef _SURFACE_H
#define _SURFACE_H




void anglevectors(LAvector<double>& thetas, LAvector<double>& phis, 
		  const unsigned int Ntheta, const unsigned int  Nphi);


void modevectors(unsigned int& NF,  LAvector<double>& nn,  LAvector<double>& mm,
		 const unsigned int n_max, const unsigned int m_max,
		 const bool mode00 = true);

void modevectors(unsigned int& NF,  LAvector<double>& nn,  LAvector<double>& mm,
		 unsigned int& n_max, unsigned int& m_max,
		 const unsigned int N_harmonics, const unsigned int M_harmonics,
		 const bool mode00 = true);


// this function is depricated, use functions above
void modevectorsR(unsigned int& NF,  LAvector<double>& nn,  LAvector<double>& mm,
		 const unsigned int n_max, const unsigned int m_max);


void make_mode_map(unsigned int& NF_A, const unsigned int& Nnn_A, const unsigned int& Nmm_A,
	      const unsigned int Ndelta_A, const unsigned int Mdelta_A,
	      const bool mode00_A,
	      unsigned int& NF_B, const unsigned int& Nnn_B, const unsigned int& Nmm_B,
	      const unsigned int Ndelta_B, const unsigned int Mdelta_B,
	      const bool mode00_B,
	      LAvector<unsigned int>& Fmap);

void fseries(const LAvector<double>& nn,  const LAvector<double>& mm,
	     const LAvector<double>& thetas, const LAvector<double>& phis,
	     Matrix<complex<double> >& f);


void remove_mode00(const LAvector<complex<double> >& vF, LAvector<complex<double> >& vFR);
void remove_mode00(const Matrix<complex<double> >& fs, Matrix<complex<double> >& fsR);

double calcAspectRatio(const FourierSurface& fsurface); 


void expandsurfaceandbases
(
 LAvector<p3vector<double> >&X, LAvector<p3vector<double> >&dA_dthetadphi,
 LAvector<p3vector<double> >& dx_dr, LAvector<p3vector<double> >& dx_dtheta, LAvector<p3vector<double> >& dx_dphi, 
 LAvector<p3vector<double> >& grad_r, LAvector<p3vector<double> >& grad_theta, LAvector<p3vector<double> >& grad_phi, 
 const FourierSurface& fsurface,
 const LAvector<double>& thetas, const LAvector<double>& phis 
 );



void expandsurface(
 LAvector<p3vector<double> >&X, LAvector<p3vector<double> >&dA_dthetadphi,
 const FourierSurface& fsurface,
 const LAvector<double>& thetas, const LAvector<double>& phis 
 );


void transformsurface(
 const LAvector<p3vector<double> >&X,
 FourierSurface& fsurface,
 const Matrix<complex<double> >& fs,
 const LAvector<double>& nn, const LAvector<double>& mm
 );


void expandfunction(
 LAvector<complex<double> >&func,
 const LAvector<complex<double> > &funcF,
 const Matrix<complex<double> >& f_ortho
 );

void expandfunction(
 LAvector<double>&funcREAL,
 const LAvector<complex<double> > &funcF,
 const Matrix<complex<double> >& f_ortho
 );

void transformfunction(
 const LAvector<double>&func,
 LAvector<complex<double> > &funcF,
 const Matrix<complex<double> >& f_ortho
 );


void completetheseries(const LAvector<double>& nn, const LAvector<double>& mm,
		       LAvector<complex<double> >& xF);



enum CondenseMode {CondenseMode_zeros_only, CondenseMode_zeros_conjugates};


void condensetheseries(const LAvector<double>& nn1, const LAvector<double>& mm1, const LAvector<complex<double> >& xF1, 
		       Matrix<double>& Aout,
		       const CondenseMode mode, const double cutoff = numeric_limits<double>::epsilon() );



int convert_sincos2exp(const LAvector<double>& ncosin, const LAvector<double>& mcosin,
		       const LAvector<double>& xCOS, const LAvector<double>& xSIN,
		       LAvector<double>& nn, LAvector<double>& mm,
		       LAvector<complex<double> >& xF, const bool giveAllWarnings = true);

int convert_exp2sincos(LAvector<double>& ncosin, LAvector<double>& mcosin,
		       LAvector<double>& xCOS,  LAvector<double>& xSIN,
		       const LAvector<double>& nn, const LAvector<double>& mm,
		       const LAvector<complex<double> >& xF);


#endif
