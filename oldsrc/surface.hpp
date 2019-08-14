



#ifndef _SURFACE_H
#define _SURFACE_H




void anglevectors(Vector<double>& thetas, Vector<double>& phis, 
		  const unsigned int Ntheta, const unsigned int  Nphi);


void modevectors(unsigned int& NF,  Vector<double>& nn,  Vector<double>& mm,
		 const unsigned int n_max, const unsigned int m_max,
		 const bool mode00 = true);

void modevectors(unsigned int& NF,  Vector<double>& nn,  Vector<double>& mm,
		 unsigned int& n_max, unsigned int& m_max,
		 const unsigned int N_harmonics, const unsigned int M_harmonics,
		 const bool mode00 = true);


// this function is depricated, use functions above
void modevectorsR(unsigned int& NF,  Vector<double>& nn,  Vector<double>& mm,
		 const unsigned int n_max, const unsigned int m_max);


void make_mode_map(unsigned int& NF_A, const unsigned int& Nnn_A, const unsigned int& Nmm_A,
	      const unsigned int Ndelta_A, const unsigned int Mdelta_A,
	      const bool mode00_A,
	      unsigned int& NF_B, const unsigned int& Nnn_B, const unsigned int& Nmm_B,
	      const unsigned int Ndelta_B, const unsigned int Mdelta_B,
	      const bool mode00_B,
	      Vector<unsigned int>& Fmap);

void fseries(const Vector<double>& nn,  const Vector<double>& mm,
	     const Vector<double>& thetas, const Vector<double>& phis,
	     Matrix<complex<double> >& f);


void remove_mode00(const Vector<complex<double> >& vF, Vector<complex<double> >& vFR);
void remove_mode00(const Matrix<complex<double> >& fs, Matrix<complex<double> >& fsR);

double calcAspectRatio(const FourierSurface& fsurface); 


void expandsurfaceandbases
(
 Vector<p3vector<double> >&X, Vector<p3vector<double> >&dA_dthetadphi,
 Vector<p3vector<double> >& dx_dr, Vector<p3vector<double> >& dx_dtheta, Vector<p3vector<double> >& dx_dphi, 
 Vector<p3vector<double> >& grad_r, Vector<p3vector<double> >& grad_theta, Vector<p3vector<double> >& grad_phi, 
 const FourierSurface& fsurface,
 const Vector<double>& thetas, const Vector<double>& phis 
 );



void expandsurface(
 Vector<p3vector<double> >&X, Vector<p3vector<double> >&dA_dthetadphi,
 const FourierSurface& fsurface,
 const Vector<double>& thetas, const Vector<double>& phis 
 );


void transformsurface(
 const Vector<p3vector<double> >&X,
 FourierSurface& fsurface,
 const Matrix<complex<double> >& fs,
 const Vector<double>& nn, const Vector<double>& mm
 );


void expandfunction(
 Vector<complex<double> >&func,
 const Vector<complex<double> > &funcF,
 const Matrix<complex<double> >& f_ortho
 );

void expandfunction(
 Vector<double>&funcREAL,
 const Vector<complex<double> > &funcF,
 const Matrix<complex<double> >& f_ortho
 );

void transformfunction(
 const Vector<double>&func,
 Vector<complex<double> > &funcF,
 const Matrix<complex<double> >& f_ortho
 );


void completetheseries(const Vector<double>& nn, const Vector<double>& mm,
		       Vector<complex<double> >& xF);



enum CondenseMode {CondenseMode_zeros_only, CondenseMode_zeros_conjugates};


void condensetheseries(const Vector<double>& nn1, const Vector<double>& mm1, const Vector<complex<double> >& xF1, 
		       Matrix<double>& Aout,
		       const CondenseMode mode, const double cutoff = numeric_limits<double>::epsilon() );



int convert_sincos2exp(const Vector<double>& ncosin, const Vector<double>& mcosin,
		       const Vector<double>& xCOS, const Vector<double>& xSIN,
		       Vector<double>& nn, Vector<double>& mm,
		       Vector<complex<double> >& xF, const bool giveAllWarnings = true);

int convert_exp2sincos(Vector<double>& ncosin, Vector<double>& mcosin,
		       Vector<double>& xCOS,  Vector<double>& xSIN,
		       const Vector<double>& nn, const Vector<double>& mm,
		       const Vector<complex<double> >& xF);


#endif
