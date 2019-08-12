
#include "coils.hpp"

// massage in place--real numbers
void massage(LAvector<double>& x, const double relative_small)
{
  
  double maxval = max(abs(x));
  double SMALL = maxval*relative_small;
  const unsigned int N = x.size();
  for (int i = 0; i<int(N); i++){
    if (abs(x[i]) < SMALL)
      x[i] = 0;
  }
}
void massage(Matrix<double>& x, const double relative_small)
{
  
  double maxval = max(abs(x));
  double SMALL = maxval*relative_small;
  const unsigned int N = x.size();
  for (int i = 0; i<int(N); i++){
    if (abs(x(i)) < SMALL)
      x(i) = 0;
  }
}


void massage(FourierSurface& f, const double relative_small)
{
  
  double maxval = 0;
  LAvector<double>& RF = f.RF();
  LAvector<double>& ZF = f.ZF();
  const unsigned int N = RF.size();

  //find max value
  for (int i = 0; i<int(N); i++){
    if (abs(RF[i]) > maxval)
      maxval = abs(RF[i]);
    if (abs(ZF[i]) > maxval)
      maxval = abs(ZF[i]);
  }

  double SMALL = maxval*relative_small;

  dispcr(maxval);
  dispcr(SMALL);


  // null out small values
  for (int i = 0; i<int(N); i++){
    if (abs(RF[i]) < SMALL)
      RF[i] = 0;
    if (abs(ZF[i]) < SMALL)
      ZF[i] = 0;
  }
  dispcr(RF);
  dispcr(ZF);
}




// massage in place--complex numbers
void massage(LAvector<complex<double> >& x, const double relative_small)
{
  
  double maxreal = max(abs(real(x)));
  double maximag = max(abs(imag(x)));
  double maxboth = max(maxreal,maximag);
  double SMALL = maxboth*relative_small;
  const unsigned int N = x.size();
  double rl,im;
  for (int i = 0; i<int(N); i++){
    rl = real(x[i]);
    if (abs(rl) < SMALL)
      rl = 0;
    im = imag(x[i]);
    if (abs(im) < SMALL)
      im = 0;
    x[i] = complex<double>(rl,im);
  }
}

void massage(Matrix<complex<double> >& x, const double relative_small)
{
  
  double maxreal = max(abs(real(x)));
  double maximag = max(abs(imag(x)));
  double maxboth = max(maxreal,maximag);
  double SMALL = maxboth*relative_small;
  const unsigned int N = x.Nrows();
  const unsigned int M = x.Ncols();
  double rl,im;
  for (int i = 0; i<int(N); i++){
    for (int j = 0; j<int(M); j++){
      rl = real(x(i,j));
      if (abs(rl) < SMALL)
	 rl = 0;
      im = imag(x(i,j));
      if (abs(im) < SMALL)
	 im = 0;
      x(i,j) = complex<double>(rl,im);
    }
  }
}




// massage in place--complex numbers
void massage_absolute(LAvector<complex<double> >& x, const double absolute_small)
{
  
  double SMALL = absolute_small;
  const unsigned int N = x.size();
  double rl,im;
  for (int i = 0; i<int(N); i++){
    rl = real(x[i]);
    if (abs(rl) < SMALL)
      rl = 0;
    im = imag(x[i]);
    if (abs(im) < SMALL)
      im = 0;
    x[i] = complex<double>(rl,im);
  }
}

void massage_absolute(Matrix<complex<double> >& x, const double absolute_small)
{
  
  double SMALL = absolute_small;
  const unsigned int N = x.Nrows();
  const unsigned int M = x.Ncols();
  double rl,im;
  for (int i = 0; i<int(N); i++){
    for (int j = 0; j<int(M); j++){
      rl = real(x(i,j));
      if (abs(rl) < SMALL)
	 rl = 0;
      im = imag(x(i,j));
      if (abs(im) < SMALL)
	 im = 0;
      x(i,j) = complex<double>(rl,im);
    }
  }
}






double gridspacing(const LAvector<p3vector<double> >& XX) {

  double delx = XX[1].x() - XX[0].x();
  double dely = XX[1].y() - XX[0].y();
  double delz = XX[1].z() - XX[0].z();
  return sqrt(sqr(delx)+sqr(dely)+sqr(delz));

}



double coil2plasmaspacing(const LAvector<p3vector<double> >& XX, const LAvector<p3vector<double> >& XXcoil) {

  double delx = XXcoil[0].x() - XX[0].x();
  double dely = XXcoil[0].y() - XX[0].y();
  double delz = XXcoil[0].z() - XX[0].z();
  return sqrt(sqr(delx)+sqr(dely)+sqr(delz));

}



void printfouriercoefs(const LAvector<double>& nn,const LAvector<double>& mm, const LAvector<double>& v1) {
  const unsigned int NF = nn.size();
  const char* s = v1.debugtxt().c_str();

  cout.setf(ios::right);
  cout.setf(ios::showpos);
  cout.flush();
  printf("%6s%6s%16s\n","n","m",s);
  for (unsigned int i = 0; i< NF; i++){
    if ( (v1[i]!=0.0) ) {
      cout.precision(0);
      cout.unsetf(ios::showpoint);
      cout << setw(6) <<nn[i];
      cout << setw(6) <<mm[i];
      cout.precision(5);
      cout.setf(ios::showpoint);
      cout << setw(16) <<v1[i];
      cout <<endl;
      cout.flush();
    }
  }
  cout << "\n";
  cout.flush();
}


void printfouriercoefs(const LAvector<double>& nn, const LAvector<double>& mm, 
		       const LAvector<double>&  v1,
		       const LAvector<double>&  v2,
		       const unsigned int precision,
		       const unsigned int fieldwidth)
{
  const unsigned int NF = nn.size();

  cout.setf(ios::right);
  cout.setf(ios::showpos);
  cout.flush();

  cout << setw(6) <<"n";
  cout << setw(6) <<"m";
  cout << setw(fieldwidth)<< v1.debugtxt();
  cout << setw(fieldwidth)<< v2.debugtxt();

  cout <<endl;

  for (unsigned int i = 0; i< NF; i++){
    double sum = 0;
    sum = abs(v1[i])+abs(v2[i]);
    if ( (sum!=0.0) ) {
      cout.precision(0);
      cout.unsetf(ios::showpoint);
      cout <<setw(6) <<nn[i];
      cout << setw(6) <<mm[i];
      cout.precision(precision);
      cout.setf(ios::showpoint);
      cout << setw(fieldwidth) << v1[i];
      cout << setw(fieldwidth) << v2[i];
      cout <<endl;
      cout.flush();
    }
  }
  cout << "\n";
  cout.flush();
}





void printfouriercoefs(const LAvector<double>& nn, const LAvector<double>& mm, 
		       const LAvector<double>&  v1,
		       const LAvector<double>&  v2,
		       const LAvector<double>&  v3,
		       const LAvector<double>&  v4,
		       const unsigned int fieldwidth) {
  const unsigned int NF = nn.size();

  cout.setf(ios::right);
  cout.setf(ios::showpos);
  cout.flush();

  cout << setw(6) <<"n";
  cout << setw(6) <<"m";
  cout << setw(fieldwidth)<< v1.debugtxt();
  cout << setw(fieldwidth)<< v2.debugtxt();
  cout << setw(fieldwidth)<< v3.debugtxt();
  cout << setw(fieldwidth)<< v4.debugtxt();

  cout <<endl;

  for (unsigned int i = 0; i< NF; i++){
    double sum = 0;
    sum = abs(v1[i])+abs(v2[i])+abs(v3[i])+abs(v4[i]);
    if ( (sum!=0.0) ) {
      cout.precision(0);
      cout.unsetf(ios::showpoint);
      cout <<setw(6) <<nn[i];
      cout << setw(6) <<mm[i];
      cout.precision(5);
      cout.setf(ios::showpoint);
      cout << setw(fieldwidth) << v1[i];
      cout << setw(fieldwidth) << v2[i];
      cout << setw(fieldwidth) << v3[i];
      cout << setw(fieldwidth) << v4[i];
      cout <<endl;
      cout.flush();
    }
  }
  cout << "\n";
  cout.flush();
}

