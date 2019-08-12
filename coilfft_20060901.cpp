/************************************************************************* 
 * 
 *   File Name    :  coilfwd2.cpp
 *   Platform     :  gnu C++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  16 May 2003
 * 
 *
 *   SYNOPSIS     
 *     This file finds the plasma surface flux for a given coil current.
 *
 **************************************************************************/



#include "coils.hpp"
#include "surface.hpp"
#include "coilfft.hpp"



using namespace std;
using namespace COOLL;



void null_small_values(fftw_complex *x, const unsigned int N, const double relatively_small) {

  double maxval = 0;

  for (unsigned int i = 0; i<N; i++){
    double rl = abs(x[i][0]);
    maxval = max(maxval,rl);
    double im = abs(x[i][1]);
    maxval = max(maxval,im);
  }

  const double small = maxval*relatively_small;

  for (int i = 0; i<int(N); i++){
    double rl = x[i][0];
    if (abs(rl) < small)
      x[i][0] = 0.0;
    double im = x[i][1];
    if (abs(im) < small)
      x[i][1] = 0.0;
  }
}



//use N1 = Nphi, N2=Ntheta
void testfunc(const LAvector<complex<double> >& v, LAvector<complex<double> >& vF, 
	   const unsigned int N1, const unsigned int N2,
	   const unsigned int maxNF1, const unsigned int maxNF2,
	   const unsigned int N1delta, const unsigned int N2delta,
	   const double neglect, const double scale)
{  

  if (maxNF1 > (N1/2)) {
    cerr<<"bad value given for maxNF1 ("<<maxNF1<<")"<<endl;
    return;
  }
  if (maxNF2 > (N2/2)) {
    cerr<<"bad value given for maxNF2 ("<<maxNF2<<")"<<endl;
    return;
  }

  const unsigned int Npts=N1*N2;
  const unsigned int NF=Npts;

  fftw_complex *x1 = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * (Npts)));

   for (unsigned int i=0; i<Npts; i++) {
    x1[i][0] = v[i].real();
    x1[i][1] = v[i].imag();
  }
 
  fftw_complex *x2 = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * (NF)));

  {
    fftw_plan p;
    const int rank =1;
    const int n = N1;
    const int howmany = N2;
    const int istride =1;
    const int idist =N1;
    const int ostride =1;
    const int odist =N1;
    p = fftw_plan_many_dft(rank,&n,howmany,x1,NULL,istride,idist,x2,NULL,ostride,odist,-1,FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  }
  fftw_free(x1);


  null_small_values(x2,NF,neglect);
  

  fftw_complex *x3 = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * Npts));

  {
    fftw_plan p;
    const int rank =1;
    const int n = N2;
    const int howmany = N1;
    const int istride =N1;
    const int idist =1;
    const int ostride =1;
    const int odist =N2;
    p = fftw_plan_many_dft(rank,&n,howmany,x2,NULL,istride,idist,x3,NULL,ostride,odist,-1,FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  }
  fftw_free(x2);


  null_small_values(x3,Npts,neglect);


  for (unsigned int n=0; n<Npts; n++) {
     const double C = scale/double(Npts);
     double tempR = C*x3[n][0];
     double tempI = C*x3[n][1];
     vF[n] = std::complex<double>(tempR,tempI);
  }


  fftw_free(x3);

}




//use N1 = Nphi, N2=Ntheta
void fft2d(const LAvector<double>& v, LAvector<complex<double> >& vF, 
	   const unsigned int N1, const unsigned int N2,
	   const unsigned int maxNF1, const unsigned int maxNF2,
	   const unsigned int N1delta, const unsigned int N2delta,
	   const double neglect, const double scale)
{  

  if (maxNF1 > (N1/2)) {
    cerr<<"bad value given for maxNF1 ("<<maxNF1<<")"<<endl;
    return;
  }
  if (maxNF2 > (N2/2)) {
    cerr<<"bad value given for maxNF2 ("<<maxNF2<<")"<<endl;
    return;
  }

  const unsigned int Npts=N1*N2;
  double *x1 = static_cast<double *>(fftw_malloc(sizeof(double) * Npts));
  for (unsigned int i=0; i<Npts; i++) {
    x1[i] = v[i];
  }
 
  const unsigned int NF=(N1/2+1)*(N2);
  fftw_complex *x2 = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * (NF)));

  {
    fftw_plan p;
    const int rank =1;
    const int n = N1;
    const int howmany = N2;
    const int istride =1;
    const int idist =N1;
    const int ostride =1;
    const int odist =N1/2+1;
    p = fftw_plan_many_dft_r2c(rank,&n,howmany,x1,NULL,istride,idist,x2,NULL,ostride,odist,FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  }
  fftw_free(x1);


  null_small_values(x2,NF,neglect);
  

  const unsigned int N_x3 = (N2)*(maxNF1+1);
  fftw_complex *x3 = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * N_x3));

  {
    fftw_plan p;
    const int rank =1;
    const int n = N2;
    const int howmany = maxNF1+1;
    const int istride =N1/2+1;
    const int idist =1;
    const int ostride =1;
    const int odist =N2;
    p = fftw_plan_many_dft(rank,&n,howmany,x2,NULL,istride,idist,x3,NULL,ostride,odist,-1,FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  }
  fftw_free(x2);


  null_small_values(x3,N_x3,neglect);


  unsigned int k =0;
  for( int n=-maxNF1; n <= int(maxNF1); n+=N1delta){
    for( int m = -maxNF2; m <= int(maxNF2); m+=N2delta){
      const double C= scale/double(Npts);
      double tempR=0;
      double tempI=0;
      int j=0;
      
      if (n>=0) {
	if (m>=0)
	  j = n*N2+m;
	else
	  j = n*N2+m+N2;
	tempR = C*x3[j][0];
	tempI = C*x3[j][1];
      } else{
	if (m>0)
	  j = abs(n)*N2-m+N2;
	else
	  j = abs(n)*N2-m;
	tempR = C*x3[j][0];
	tempI = -C*x3[j][1];
      }
      vF[k++] = std::complex<double>(tempR,tempI);
    }
  }


  fftw_free(x3);

}



//use N1 = Nphi, N2=Ntheta
void fft2d(const LAvector<complex<double> >& v, LAvector<complex<double> >& vF, 
	   const unsigned int N1, const unsigned int N2,
	   const unsigned int maxNF1, const unsigned int maxNF2,
	   const unsigned int N1delta, const unsigned int N2delta,
	   const double neglect, const double scale)
{  

  if (maxNF1 > (N1/2)) {
    cerr<<"bad value given for maxNF1 ("<<maxNF1<<")"<<endl;
    return;
  }
  if (maxNF2 > (N2/2)) {
    cerr<<"bad value given for maxNF2 ("<<maxNF2<<")"<<endl;
    return;
  }

  const unsigned int Npts=N1*N2;
  const unsigned int NF=Npts;

  fftw_complex *x1 = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * (Npts)));

   for (unsigned int i=0; i<Npts; i++) {
    x1[i][0] = v[i].real();
    x1[i][1] = v[i].imag();
  }
 
  fftw_complex *x2 = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * (NF)));

  {
    fftw_plan p;
    const int rank =1;
    const int n = N1;
    const int howmany = N2;
    const int istride =1;
    const int idist =N1;
    const int ostride =1;
    const int odist =N1;
    p = fftw_plan_many_dft(rank,&n,howmany,x1,NULL,istride,idist,x2,NULL,ostride,odist,-1,FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  }
  fftw_free(x1);


  null_small_values(x2,NF,neglect);
  

  const unsigned int N_x3 = (N2)*(maxNF1+1);
  fftw_complex *x3 = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * N_x3));

  {
    fftw_plan p;
    const int rank =1;
    const int n = N2;
    const int howmany = maxNF1+1;
    const int istride =N1;
    const int idist =1;
    const int ostride =1;
    const int odist =N2;
    p = fftw_plan_many_dft(rank,&n,howmany,x2,NULL,istride,idist,x3,NULL,ostride,odist,-1,FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  }
  fftw_free(x2);


  null_small_values(x3,N_x3,neglect);


  unsigned int k =0;
  for( int n=-maxNF1; n <= int(maxNF1); n+=N1delta){
    for( int m = -maxNF2; m <= int(maxNF2); m+=N2delta){
      const double C= scale/double(Npts);
      double tempR=0;
      double tempI=0;
      int j=0;
      
      if (n>=0) {
	if (m>=0)
	  j = n*N2+m;
	else
	  j = n*N2+m+N2;
	tempR = C*x3[j][0];
	tempI = C*x3[j][1];
      } else{
	if (m>0)
	  j = abs(n)*N2-m+N2;
	else
	  j = abs(n)*N2-m;
	tempR = C*x3[j][0];
	tempI = C*x3[j][1];
      }
      vF[k++] = std::complex<double>(tempR,tempI);
    }
  }


  fftw_free(x3);

}





//use N1 = Nphi, N2=Ntheta
void ifft2d(LAvector<complex<double> >& v, const LAvector<complex<double> >& vF, 
	   const unsigned int N1, const unsigned int N2,
	   const unsigned int maxNF1, const unsigned int maxNF2,
	   const unsigned int N1delta, const unsigned int N2delta,
	   const double neglect, const double scale)
{  

  if (maxNF1 > (N1/2)) {
    cerr<<"bad value given for maxNF1 ("<<maxNF1<<")"<<endl;
    return;
  }
  if (maxNF2 > (N2/2)) {
    cerr<<"bad value given for maxNF2 ("<<maxNF2<<")"<<endl;
    return;
  }

  const unsigned int Npts=N1*N2;
  //  const unsigned int NF=Npts;

  fftw_complex *x1 = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * (Npts)));

  for( int n=0; n <= int(Npts); n++) {
     x1[n][0] = 0;
     x1[n][1] = 0;
  }

  unsigned int k =0;
  for( int n=-maxNF1; n <= int(maxNF1); n+=N1delta){
    for( int m = -maxNF2; m <= int(maxNF2); m+=N2delta){
       //      const double C= scale/double(Npts);
      const double C= 1.0;
      const double tempR=vF[k].real();
      const double tempI=vF[k].imag();
      k++;
      int j1=0;
      int j2=0;
      
      if (n>=0) 
	 j1 = N2*(n);
      else
	 j1 = N2*(n+N1);
      if (m>=0)
	 j2 = m;
      else
	 j2 = m+N2;

      const int j = j1 + j2;

      x1[j][0] = C * tempR;
      x1[j][1] = C * tempI;
    }
  }

  for( int n=0; n <= int(Npts); n++) {
     v[n] = std::complex<double>(x1[n][0],x1[n][1]);
  }
  fftw_free(x1);return;


  fftw_complex *x2 = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * (Npts)));
  {
    fftw_plan p;
    const int rank =1;
    const int n = N1;
    const int howmany = N2;
    const int istride =1;
    const int idist =N1;
    const int ostride =1;
    const int odist =N1;
    p = fftw_plan_many_dft(rank,&n,howmany,x1,NULL,istride,idist,x2,NULL,ostride,odist,1,FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  }
  fftw_free(x1);


  null_small_values(x2,Npts,neglect);

  fftw_complex *x3 = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * Npts));

  {
    fftw_plan p;
    const int rank =1;
    const int n = N2;
    const int howmany = N1;
    const int istride =N1;
    const int idist =1;
    const int ostride =1;
    const int odist =N2;
    p = fftw_plan_many_dft(rank,&n,howmany,x2,NULL,istride,idist,x3,NULL,ostride,odist,1,FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  }
  fftw_free(x2);


  null_small_values(x3,Npts,neglect);


  for( int n=0; n <= int(Npts); n++) {
     v[n] = std::complex<double>(x3[n][0],x3[n][1]);
  }

  fftw_free(x3);


}














void fft2d(const double * const v, std::complex<double> * const vF,
	   const unsigned int N1, const unsigned int N2,
	   const unsigned int maxNF1, const unsigned int maxNF2,
	   const unsigned int v_stride, const unsigned int vF_stride,
	   const int sign,const double neglect)
{  

  if (maxNF1 > (N1/2)) {
    cerr<<"bad value given for maxNF1 ("<<maxNF1<<")"<<endl;
    return;
  }
  if (maxNF2 > (N2/2)) {
    cerr<<"bad value given for maxNF2 ("<<maxNF2<<")"<<endl;
    return;
  }

  const unsigned int Npts=N1*N2;
  double *x1 = static_cast<double *>(fftw_malloc(sizeof(double) * Npts));
  for (unsigned int i=0; i<Npts; i++) {
    x1[i] = v[i*v_stride];
  }
 
  const unsigned int NF=(N1/2+1)*(N2);
  fftw_complex *x2 = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * (NF)));

  {
    fftw_plan p;
    const int rank =1;
    const int n = N1;
    const int howmany = N2;
    const int istride =1;
    const int idist =N1;
    const int ostride =1;
    const int odist =N1/2+1;
    p = fftw_plan_many_dft_r2c(rank,&n,howmany,x1,NULL,istride,idist,x2,NULL,ostride,odist,FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  }
  fftw_free(x1);


  null_small_values(x2,NF,neglect);

  //conjugate x2 if sign is positive
  if (sign>0) {
    for (unsigned int k=0; k<NF; k++) {
      x2[k][1] = - x2[k][1];
    }
  }

  const unsigned int N_x3 = (N2)*(maxNF1+1);
  fftw_complex *x3 = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * N_x3));

  {
    fftw_plan p;
    const int rank =1;
    const int n = N2;
    const int howmany = maxNF1+1;
    const int istride =N1/2+1;
    const int idist =1;
    const int ostride =1;
    const int odist =N2;
    p = fftw_plan_many_dft(rank,&n,howmany,x2,NULL,istride,idist,x3,NULL,ostride,odist,sign,FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  }
  fftw_free(x2);

  null_small_values(x3,N_x3,neglect);


  unsigned int k =0;
  for( int n=-maxNF1; n <= int(maxNF1); n++){
    for( int m = -maxNF2; m <= int(maxNF2); m++){
      const double C= 1/double(Npts);
      double tempR=0;
      double tempI=0;
      int j=0;
      
      if (n>=0) {
	if (m>=0)
	  j = n*N2+m;
	else
	  j = n*N2+m+N2;
	tempR = C*x3[j][0];
	tempI = C*x3[j][1];
      } else{
	if (m>0)
	  j = abs(n)*N2-m+N2;
	else
	  j = abs(n)*N2-m;
	tempR = C*x3[j][0];
	tempI = -C*x3[j][1];
      }
      vF[k++*vF_stride] = std::complex<double>(tempR,tempI);
    }
  }

  fftw_free(x3);

}




void fft2d(const std::complex<double> * const v, std::complex<double> * const vF,
	   const unsigned int N1, const unsigned int N2,
	   const unsigned int maxNF1, const unsigned int maxNF2,
	   const unsigned int v_stride, const unsigned int vF_stride,
	   const int sign,const double neglect)
{  

  if (maxNF1 > (N1/2)) {
    cerr<<"bad value given for maxNF1 ("<<maxNF1<<")"<<endl;
    return;
  }
  if (maxNF2 > (N2/2)) {
    cerr<<"bad value given for maxNF2 ("<<maxNF2<<")"<<endl;
    return;
  }

  const unsigned int Npts=N1*N2;
  fftw_complex *x1 = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * Npts));
  for (unsigned int i=0; i<Npts; i++) {
    x1[i][0] = v[i*v_stride].real();
    x1[i][1] = v[i*v_stride].imag();
  }
 
  const unsigned int NF=N1*(N2);
  fftw_complex *x2 = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * (NF)));

  {
    fftw_plan p;
    const int rank =1;
    const int n = N1;
    const int howmany = N2;
    const int istride = 1;
    const int idist = N1;
    const int ostride = 1;
    const int odist = N1;
    p = fftw_plan_many_dft(rank,&n,howmany,x1,NULL,istride,idist,x2,NULL,ostride,odist,sign,FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  }
  fftw_free(x1);


  null_small_values(x2,NF,neglect);

  const unsigned int NF1 = min(2*maxNF1+1,N1);

  const unsigned int N_x3 =  (N2)*NF1;
  fftw_complex *x3 = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) *N_x3));

  {
    fftw_plan p;
    const int rank =1;
    const int n = N2;
    const int howmany = maxNF1+1;
    const int istride =N1;
    const int idist =1;
    const int ostride =1;
    const int odist =N2;
    p = fftw_plan_many_dft(rank,&n,howmany,x2,NULL,istride,idist,x3,NULL,ostride,odist,sign,FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  }

  if (maxNF1>0) {
    fftw_plan p;
    const int rank =1;
    const int n = N2;
    const int howmany = NF1-(maxNF1+1);
    const int istride =N1;
    const int idist =1;
    const int ostride =1;
    const int odist =N2;
    p = fftw_plan_many_dft(rank,&n,howmany,&(x2[N1-howmany]),NULL,istride,idist,&(x3[N2*(maxNF1+1)]),NULL,ostride,odist,sign,FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  }
  fftw_free(x2);

  null_small_values(x3,N_x3,neglect);

  unsigned int k =0;
  for( int n=-maxNF1; n <= int(maxNF1); n++){
    for( int m = -maxNF2; m <= int(maxNF2); m++){
      const double C= 1/double(Npts);
      double tempR=0;
      double tempI=0;
      int j=0;
      
      if (n>=0) {
	if (m>=0)
	  j = n*N2 + m;
	else
	  j = n*N2 + m+N2;
      } else{
	if (m>=0)
	  j = (n+NF1)*N2 + m;
	else
	  j = (n+NF1)*N2 + m+N2;
      }
      tempR = C*x3[j][0];
      tempI = C*x3[j][1];
      vF[k++*vF_stride] = std::complex<double>(tempR,tempI);
    }
  }

  fftw_free(x3);

}







// note that this function clobbers M
void fft_of_M(Matrix<double>& M, Matrix<complex<double> >& MF, 
	      const unsigned int Nphi, const unsigned int Ntheta,
	      const unsigned int Nnn, const unsigned int Nmm,
	      const unsigned int Ndelta, const unsigned int Mdelta,
	      const double neglect)
{

   unsigned int NF = (2*Nnn+1)*(2*Nmm+1);
  const unsigned int Npts = Nphi*Ntheta;
  const double PI = 3.141592653589793;

  Matrix<complex<double> > MxF(Npts,NF,"fft_of_M::MxF");
  
  for (unsigned int i=0; i<Npts; i++) 
    fft2d( &(M(i,0)), &(MxF(i,0)), Nphi, Ntheta, Nnn, Nmm, 1, 1, +1, neglect);

  // free up memory used by M
  M.resize(0,0);


  Matrix<complex<double> > MFtemp(NF,NF,"MFtemp"); 

  for (unsigned int k=0; k<NF; k++) 
    fft2d( &(MxF(0,k)), &(MFtemp(0,k)), Nphi, Ntheta, Nnn, Nmm, NF, NF, -1, neglect);


  // free up memory used by M
  MxF.resize(0,0);


  // remove n=0,m=0 mode and unwanted modes  and scale output

   unsigned int NFrow;
  LAvector<unsigned int> Fmap_row;
  make_mode_map(NFrow,Nnn,Nmm,Ndelta,Mdelta,false,NF,Nnn,Nmm,1,1,true,Fmap_row);
   unsigned int NFcol;
  LAvector<unsigned int> Fmap_col;
  make_mode_map(NFcol,Nnn,Nmm,Ndelta,Mdelta,false,NF,Nnn,Nmm,1,1,true,Fmap_col);

  MF.resize(NFrow,NFcol); 
  for (unsigned int r =0; r<NFrow; r++) {
    for (unsigned int c =0; c<NFcol; c++) {
       const double C2= 4.0*PI*PI;
       MF(r,c) = C2*MFtemp(Fmap_row[r], Fmap_col[c]);
    }
  }

  return;
}












// note that this function clobbers M
void half_fft_of_M(Matrix<double>& M, Matrix<complex<double> >& MF, 
	      const unsigned int Nphi, const unsigned int Ntheta,
	      const unsigned int Nnn, const unsigned int Nmm,
		   const unsigned int Ndelta, const unsigned int Mdelta,
	      const double neglect)
{

  const unsigned int Npts = Nphi*Ntheta;
  unsigned int NF = (2*Nnn+1)*(2*Nmm+1);
  const double PI = 3.141592653589793;
  Matrix<complex<double> > MFtemp(Npts,NF,"fft_of_M::MxF");

  for (unsigned int i=0; i<Npts; i++) 
    fft2d( &(M(i,0)), &(MFtemp(i,0)), Nphi, Ntheta, Nnn, Nmm, 1, 1, +1, neglect);

  // free up memory used by M
  M.resize(0,0);


  LAvector<unsigned int> Fmap_col;
  unsigned int NFcol = 0;

  make_mode_map(NFcol,Nnn,Nmm,Ndelta,Mdelta,false,NF,Nnn,Nmm,1,1,true,Fmap_col);


  MF.resize(Npts,NFcol); 
  for (unsigned int r =0; r<Npts; r++) {
    for (unsigned int c =0; c<NFcol; c++) {
       const double C2= 2.0*PI;
       MF(r,c) = C2*MFtemp(r, Fmap_col[c]);
    }
  }

  return;
}


// note that this function clobbers Bdotgradf
void omega_half_fft(Matrix<complex<double> >& Bdotgradf, Matrix<complex<double> >& Omega, 
		    const unsigned int Nphi, const unsigned int Ntheta,
		    const unsigned int Nnn, const unsigned int Nmm,
		    const unsigned int Ndelta, const unsigned int Mdelta,
		    const double neglect)
{

  const unsigned int Npts = Nphi*Ntheta;
 unsigned int NF = (2*Nnn+1)*(2*Nmm+1);
  const double PI = 3.141592653589793;
  //  Matrix<complex<double> > Ftemp1(NF-1,Npts,"Omega_half_fft::Ftemp1");
  Matrix<complex<double> > Ftemp2(NF-1,NF,"Omega_half_fft::Ftemp2");
  const unsigned int NF2 = Bdotgradf.Nrows();

  dispcr(NF);
  dispcr(NF2);
//   for (unsigned int i=0; i<Npts; i++) {
//     for (unsigned int k =0; k<NF2; k++) {
//       Ftemp1(k,i) = -conj(Bdotgradf(i,k));
//       //      Bdotgradf(i,c) = -conj(Bdotgradf(i,c));
//     }
//   }


  for (unsigned int c=0; c<NF2; c++) 
    fft2d( &(Bdotgradf(c,0)), &(Ftemp2(c,0)), Nphi, Ntheta, Nnn, Nmm, 1, 1, +1, neglect);
  

  // free up memory
  Bdotgradf.resize(0,0);
  //  Ftemp1.resize(0,0);



  // remove n=0,m=0 mode and unwanted modes  and scale output

  unsigned int NFcol;
  LAvector<unsigned int> Fmap_col;
  make_mode_map(NFcol,Nnn,Nmm,Ndelta,Mdelta,false,NF,Nnn,Nmm,1,1,true,Fmap_col);

  dispcr(NFcol);

  //  Omega.resize(NFrow,NFcol); 
  for (unsigned int r =0; r<NF2; r++) {
    for (unsigned int c =0; c<NFcol; c++) {
       const double C2= 2*PI;
      Omega(c,r) = C2*Ftemp2(r, Fmap_col[c]);
    }
  }


  return;
}












