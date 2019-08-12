/************************************************************************* 
 * 
 *   File Name    :  coilio.cpp
 *   Platform     :  Red Hat LINUX 
 *   Author       :  Ron Schmitt
 *   Date         :  16 May 2003
 * 
 *
 *   SYNOPSIS     
 *    Loads Fourier coefficients from a file
 *
 **************************************************************************/


// Standard C++ Libraries

#include <fstream>
#include <iostream>

using namespace std;

#include "coils.hpp"
#include "coilio.hpp"



// For loading surface coef's
int load_fourier_surface(const std::string& fname, FourierSurface& fsurface)
{
  

  COOLL::Matrix<double> Mtemp("load_fourier_surface::Mtemp");
  Mtemp.textformat(text_nobraces);
  if (COOLL::load(Mtemp,fname)) {
    return 1;
  }

  const unsigned int Nmodes = Mtemp.Nrows();

  fsurface.size(Nmodes);

  fsurface.nn() = Mtemp.col(0);
  fsurface.mm() = Mtemp.col(1);
  fsurface.RF() = Mtemp.col(2);
  fsurface.ZF() = Mtemp.col(3);
  
  return 0;
}

// For saving surface coef's
int save_fourier_surface(const std::string& fname, const FourierSurface& fsurface)
{
  
  const unsigned int w = 19;
  const unsigned int precision = 10;

  std::ostringstream stream;
  stream << "%"<<setw(w-1)<< "n"<<setw(w)<<"m";
  stream <<setw(w)<<"COS"<<setw(w)<<"SIN";

  const unsigned int Nmodes = fsurface.nn().size();
  COOLL::Matrix<double> Mtemp(Nmodes,4,"save_fourier_surface::Mtemp");


  Mtemp.col(0) = fsurface.nn();
  Mtemp.col(1) = fsurface.mm();
  Mtemp.col(2) = fsurface.RF();
  Mtemp.col(3) = fsurface.ZF();

  Mtemp.width(w);
  Mtemp.textformat(text_nobraces);
  if (  save(Mtemp,fname,precision,stream.str()) ) {
    return 1;
  }
  
  return 0;
}





int save_coefs(const std::string fname, CoefFileFormat format, 
	       const LAvector<double>& nn, const LAvector<double>& mm,
	       const LAvector<complex<double> >& xF, const bool rankorder ) {

  const unsigned int w = 19;
  const unsigned int precision = 10;


  COOLL::Matrix<double> Atemp("save_coefs::Atemp");

  std::ostringstream stream;


  switch (format) {
  case CoefFileFormat_sincos:
    {
      stream << "%"<<setw(w-1)<< "n"<<setw(w)<<"m";
      stream <<setw(w)<<"COS"<<setw(w)<<"SIN";
      LAvector<double> ncosin;
      LAvector<double> mcosin;
      LAvector<double> xCOS;
      LAvector<double> xSIN;
      convert_exp2sincos(ncosin,mcosin,xCOS,xSIN,nn,mm,xF);
      const unsigned int Nmodes = xCOS.size();
      Atemp.resize(Nmodes,4);
      Atemp.col(0) = ncosin;
      Atemp.col(1) = mcosin;
      Atemp.col(2) = xCOS;
      Atemp.col(3) = xSIN;
    }      
    break;
  case CoefFileFormat_complexexp:
    {
      stream << "%"<<setw(w-1)<< "n"<<setw(w)<<"m";
      stream <<setw(w)<<"REAL"<<setw(w)<<"IMAG";
      condensetheseries(nn,mm,xF,Atemp,CondenseMode_zeros_conjugates);
    }
    break;
  default:
    break;
  } //switch


  const unsigned int Na = Atemp.Nrows();
  COOLL::Matrix<double> Btemp(Na,4,"save_coefs::Btemp");

  if (rankorder) {
    LAvector<double> magsqr = sqr(vcast<double>(Atemp.col(2))) + sqr(vcast<double>(Atemp.col(3)));
    LAvector<unsigned int> indvec(Na,"save_coefs::indvec");
    indvec = sortwind(magsqr);
    for (unsigned int r=0; r<Na; r++)
      for (unsigned int c=0; c<4; c++)
	Btemp(r,c) = Atemp(indvec[Na-1-r],c);
  } else {
    Btemp = Atemp;
  }

		
  Atemp.resize(0,0);

  Btemp.textformat(text_nobraces);
  Btemp.width(w);
  if (  save(Btemp,fname,precision,stream.str())) {
    std::cerr<<"ERROR: could not save to "<<fname<<std::endl;    
    return 1;
  }


  return 0;
}



// nn and mm must be preset using "modevectors()"

int load_coefs(const std::string fname, CoefFileFormat format, 
	       const LAvector<double>& nn, const LAvector<double>& mm,
	       LAvector<complex<double> >& xF, const bool giveAllWarnings ) {

  COOLL::Matrix<double> Mtemp("load_coefs::Mtemp");
  Mtemp.textformat(text_nobraces);
  if (COOLL::load(Mtemp,fname)) {
    std::cerr<<"ERROR: could not open "<<fname<<std::endl;    
    return 1;
  }
  const unsigned int Nmodes = Mtemp.Nrows();

  switch (format) {
  case CoefFileFormat_sincos:
    {
      COOLL::LAvector<double> ncopy = nn;
      COOLL::LAvector<double> mcopy = mm;
	
      COOLL::LAvector<double> nFILE(Nmodes);
      nFILE = Mtemp.col(0);
      COOLL::LAvector<double> mFILE(Nmodes);
      mFILE = Mtemp.col(1);
      COOLL::LAvector<double> xcos(Nmodes);
      xcos = Mtemp.col(2);
      COOLL::LAvector<double> xsin(Nmodes);
      xsin = Mtemp.col(3);

      convert_sincos2exp(nFILE,mFILE,xcos,xsin,ncopy,mcopy,xF,giveAllWarnings);      
    }      
    break;
  case CoefFileFormat_sincos_RHC2LHC:
    {
       // this works for converting RightHandCoordinate system to LeftHandCoordinate system
       // AND vice versa!
      LAvector<double> ncopy = nn;
      LAvector<double> mcopy = mm;
	
      LAvector<double> nFILE(Nmodes);
      nFILE = Mtemp.col(0);
      LAvector<double> mFILE(Nmodes);
      mFILE = Mtemp.col(1);
      LAvector<double> xcos(Nmodes);
      xcos = Mtemp.col(2);
      LAvector<double> xsin(Nmodes);
      xsin = Mtemp.col(3);


      xcos=cos(mFILE*PI)*xcos;
      xsin=cos(mFILE*PI)*xsin;
      mFILE=-mFILE;

      convert_sincos2exp(nFILE,mFILE,xcos,xsin,ncopy,mcopy,xF,giveAllWarnings);      
    }      
    break;
  case CoefFileFormat_complexexp:
    {
      bool warned_already = false;
      const unsigned int NF = nn.size();
      for(unsigned int j = 0; j<Nmodes; j++) {
	int n=int(Mtemp(j,0));
	int m=int(Mtemp(j,1));
	std::complex<double> A = std::complex<double>(Mtemp(j,2),Mtemp(j,3));
	bool found = false;
	
	for(unsigned int k = 0; k<NF; k++) {
	  if ( (int(nn[k])==n) && (int(mm[k])==m) ){
	    xF[k] = A;
	    found=true;
	  }
	}
	if (!found) {
	   if (giveAllWarnings || !warned_already) {
	      int nmax = int(max(abs(nn)));
	      int mmax = int(max(abs(mm)));
	      std::cerr<<"WARNING: mode (n="<<n<<",m="<<m<<") in file '"<<fname<<"' was not included due to"<<std::endl;
	      std::cerr<<"       maximum number of modes (nmax="<<nmax<<",mmax="<<mmax<<") that "<<std::endl;
	      std::cerr<<"       were specified at command line OR"<<std::endl;
	      std::cerr<<"       harmonic numbers  (Nharm="<<",Mharm="<<") that "<<std::endl;
	      std::cerr<<"       were specified at command line"<<std::endl;
              if (!giveAllWarnings)
                std::cerr<<"       FURTHER WARNINGS WILL BE SUPPRESSED"<<std::endl;
              warned_already = true;
	   }
	}
      }
      // insert coef's for all conjugate modes (as needed)
      completetheseries(nn,mm,xF);
    }
    break;
  default:
    break;
  } //switch

  return 0;
}



// For loading surface coef's, formatted according to garabedian's paper
// data is m n delta

int load_garabedian_surface(const std::string& fname, FourierSurface& fsurface)
{


  COOLL::Matrix<double> Mtemp("load_garabedian_surface::Mtemp");
  Mtemp.textformat(text_nobraces);
  if (COOLL::load(Mtemp,fname)) {
    return 1;
  }

  //  dispcr(Mtemp);
  const unsigned int Nmodes = Mtemp.Nrows();

  COOLL::LAvector<double> nn(Nmodes,"load_garabedian_surface::nn");
  COOLL::LAvector<double> mm(Nmodes,"load_garabedian_surface::mm");
  COOLL::LAvector<double> del(Nmodes,"load_garabedian_surface::del");

  mm = Mtemp.col(0);
  nn = Mtemp.col(1);
  del = Mtemp.col(2);


  // convert from garabedian notation to my notation
  // second reversal is due to fact that u is defined by garabedian
  // as counter clockwise direction, whereas I define as clockwise
  // direction so as to form right-handed coord. system
  mm = -(1-mm);
  nn = 2*nn;

  //  dispcr(mm);
  //dispcr(nn);
  //dispcr(del);

  const unsigned int Nmax = static_cast<unsigned int>(max(abs(nn)));
  const unsigned int Mmax = static_cast<unsigned int>(max(abs(mm)));
  const unsigned int NM2 = (Nmax+1)*(2*Mmax+1);
  COOLL::LAvector<double> nn2(NM2,"load_garabedian_surface::nn2");
  COOLL::LAvector<double> mm2(NM2,"load_garabedian_surface::mm2");
  COOLL::LAvector<double> RF2(NM2,"load_garabedian_surface::RF2");
  COOLL::LAvector<double> ZF2(NM2,"load_garabedian_surface::ZF2");
  nn2 = 0;
  mm2 = 0;
  RF2 = 0;
  ZF2 = 0;

  //  dispcr(Nmax);
  //  dispcr(Mmax);
  //  dispcr(NM2);

  unsigned int k =0;
  for( unsigned int n=0; n <= Nmax; n++){
    for( int m = -int(Mmax); m <= int(Mmax); m++, k++){
      //      disp(n);disp(m);dispcr(k);
      nn2[k] = double(n);
      mm2[k] = double(m);
    }
  }

  //  dispcr(nn2);
  //  dispcr(mm2);


  for( unsigned int i=0; i < Nmodes; i++){
    const  int n = static_cast<int>(nn[i]);
    const  int m = static_cast<int>(mm[i]);
    if ( (n==0) && (m==0) ) {
      unsigned int i2 = Mmax + m + (2*Mmax+1)*abs(n);
      RF2[i2] += del[i];
      ZF2[i2] = 0;
    } else if ( (n==0) && (m>0) ) {
      unsigned int i2 = Mmax + m + (2*Mmax+1)*abs(n);
      RF2[i2] += del[i];
      ZF2[i2] += del[i];
    } else if ( (n==0) && (m<0) ) {
      unsigned int i2 = Mmax + abs(m) + (2*Mmax+1)*abs(n);
      RF2[i2] += del[i];
      ZF2[i2] += -del[i];
    } else if (n >= 1) {
      unsigned int i2 = Mmax + m + (2*Mmax+1)*n;
      RF2[i2] += del[i];
      ZF2[i2] += del[i];
    } else if (n<=-1)  {
      unsigned int i2 = Mmax + -m + (2*Mmax+1)*abs(n);
      RF2[i2] += +del[i];
      ZF2[i2] += -del[i];
    } 
  }


//   dispcr(nn2);
//   dispcr(mm2);
//   dispcr(RF2);
//   dispcr(ZF2);

  unsigned int Nmodes2 = 0;
  for( unsigned int i=0; i < NM2; i++){
    if ( (RF2[i]!=0) || (ZF2[i]!=0) )
      Nmodes2++;
  }
  

//   dispcr(Nmodes2);

  fsurface.nn().resize(Nmodes2);
  fsurface.mm().resize(Nmodes2);
  fsurface.RF().resize(Nmodes2);
  fsurface.ZF().resize(Nmodes2);

  unsigned int i2=0;
  for( unsigned int i=0; i < NM2; i++){
    if ( (RF2[i]!=0) || (ZF2[i]!=0) ){
      fsurface.nn()[i2] = nn2[i];
      fsurface.mm()[i2] = mm2[i];
      fsurface.RF()[i2] = RF2[i];
      fsurface.ZF()[i2] = ZF2[i];
      i2++;
    }
  }

  return 0;
}


// For loading surface coef's, formatted according to garabedian's paper
// data is m n delta

int load_garabedian_current(const std::string& fname, FourierSurface& fsurface)
{


  COOLL::Matrix<double> Mtemp("load_garabedian_current::Mtemp");
  Mtemp.textformat(text_nobraces);
  if (COOLL::load(Mtemp,fname)) {
    return 1;
  }

  //  dispcr(Mtemp);
  const unsigned int Nmodes = Mtemp.Nrows();

  COOLL::LAvector<double> nn(Nmodes,"load_garabedian_current::nn");
  COOLL::LAvector<double> mm(Nmodes,"load_garabedian_current::mm");
  COOLL::LAvector<double> kappa(Nmodes,"load_garabedian_current::del");

  mm = Mtemp.col(0);
  nn = Mtemp.col(1);
  kappa = Mtemp.col(2);


  // convert from garabedian notation to my notation
  mm = -mm;
  nn = -2*nn;
  kappa = kappa;

  unsigned int Nmodes2 = 0;
  for( unsigned int i=0; i < Nmodes; i++){
    if ( (kappa[i]!=0) )
      Nmodes2++;
  }
  

//   dispcr(Nmodes2);

  fsurface.nn().resize(Nmodes2);
  fsurface.mm().resize(Nmodes2);
  fsurface.RF().resize(Nmodes2);
  fsurface.ZF().resize(Nmodes2);

  unsigned int i2=0;
  for( unsigned int i=0; i < Nmodes; i++){
    if ( (kappa[i]!=0) ) {
      fsurface.nn()[i2] = nn[i];
      fsurface.mm()[i2] = mm[i];
      fsurface.RF()[i2] = 0;
      fsurface.ZF()[i2] = kappa[i];
      i2++;
    }
  }

  return 0;
}
