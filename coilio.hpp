




#include <string>

int load_fourier_surface(const std::string& fname, FourierSurface& fsurface);

int save_fourier_surface(const std::string& fname, const FourierSurface& fsurface);

enum CoefFileFormat {CoefFileFormat_sincos,CoefFileFormat_complexexp,CoefFileFormat_sincos_RHC2LHC};

int load_coefs(const std::string fname, CoefFileFormat format, 
	       const LAvector<double>& nn, const LAvector<double>& mm,
	       LAvector<complex<double> >& xF, const bool giveAllWarnings =true);

int save_coefs(const std::string fname, CoefFileFormat format, 
	       const LAvector<double>& nn, const LAvector<double>& mm,
	       const LAvector<complex<double> >& xF, const bool rankorder = false);


int load_garabedian_surface(const std::string& fname, FourierSurface& fsurface);

int load_garabedian_current(const std::string& fname, FourierSurface& fsurface);
