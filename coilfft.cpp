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
      double rl = abs(x[i][0]);
      double im = abs(x[i][1]);
      if (sqrt(rl*rl+im*im) < small){
	 x[i][0] = 0.0;
	 x[i][1] = 0.0;
      }

      if (rl < (im*relatively_small))
	 x[i][0] = 0.0;
      if (im < (rl*relatively_small))
	 x[i][1] = 0.0;
   }
}




////////////////////////////////////////////////////////////////////////////
//  2D FFT of a REAL VECTOR
////////////////////////////////////////////////////////////////////////////


//use N1 = Nphi, N2=Ntheta
void fft2d(const LAvector<double>& v, LAvector<complex<double> >& vF, 
	   const unsigned int N1, const unsigned int N2,
	   const unsigned int maxNF1, const unsigned int maxNF2,
	   const unsigned int N1delta, const unsigned int N2delta,
	   const double neglect, const double scale, const bool mode00)
{  
   const int sign = -1;
   //   printcr("fft of real LAvector");


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

   //conjugate x2 if sign is positive
   if (sign>0) {
      for (unsigned int k=0; k<NF; k++) {
	 x2[k][1] = - x2[k][1];
      }
   }

//  const unsigned int NF1 = min(2*maxNF1+1,N1);
//  const unsigned int NF2 = N2;
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
	    if (m<=0)
	       j = abs(n)*N2-m;
	    else
	       j = abs(n)*N2-m+N2;
	    tempR = C*x3[j][0];
	    tempI = -C*x3[j][1];
	 }
	 if ((!mode00) && (n == 0) && (m == 0)) 
	    ;
	 else 
	    vF[k++] = std::complex<double>(tempR,tempI);
      }
   }

   fftw_free(x3);


}


////////////////////////////////////////////////////////////////////////////
//  2D FFT of a COMPLEX VECTOR
////////////////////////////////////////////////////////////////////////////
//use N1 = Nphi, N2=Ntheta
void fft2d(const LAvector<complex<double> >& v, LAvector<complex<double> >& vF, 
	   const unsigned int N1, const unsigned int N2,
	   const unsigned int maxNF1, const unsigned int maxNF2,
	   const unsigned int N1delta, const unsigned int N2delta,
	   const double neglect, const double scale, const bool mode00)
{  

   const int sign = -1;
   //   printcr("fft of complex LAvector");


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
    x1[i][0] = v[i].real();
    x1[i][1] = v[i].imag();
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
  const unsigned int NF2 = N2;

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
  for( int n=-maxNF1; n <= int(maxNF1); n+=N1delta){
    for( int m = -maxNF2; m <= int(maxNF2); m+=N2delta){
       const double C= scale/double(Npts);
       int j1=0;
       int j2=0;
      
       if (n>=0) 
	  j1 = NF2*(n);
       else
	  j1 = NF2*(n+NF1);
       if (m>=0)
	  j2 = m;
       else
	  j2 = m+NF2;
       
       int j = j1 + j2;
       
       double tempR = C*x3[j][0];
       double tempI = C*x3[j][1];
       
       if ((!mode00) && (n == 0) && (m == 0)) 
	  ;
       else 
	  vF[k++] = std::complex<double>(tempR,tempI);
    }
  }

  fftw_free(x3);




}


//////////////////////////////////////////////////////////////////


// void testfunc(const LAvector<complex<double> >& v, LAvector<complex<double> >& vF, 
// 	   const unsigned int N1, const unsigned int N2,
// 	   const unsigned int maxNF1, const unsigned int maxNF2,
// 	   const unsigned int N1delta, const unsigned int N2delta,
// 	   const double neglect, const double scale, const bool mode00)
// {  

//    const int sign = -1;
//    //   printcr("fft of complex LAvector");


//   if (maxNF1 > (N1/2)) {
//     cerr<<"bad value given for maxNF1 ("<<maxNF1<<")"<<endl;
//     return;
//   }
//   if (maxNF2 > (N2/2)) {
//     cerr<<"bad value given for maxNF2 ("<<maxNF2<<")"<<endl;
//     return;
//   }

//   const unsigned int Npts=N1*N2;
//   fftw_complex *x1 = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * Npts));
//   for (unsigned int i=0; i<Npts; i++) {
//     x1[i][0] = v[i].real();
//     x1[i][1] = v[i].imag();
//   }
 
//   const unsigned int NF=N1*(N2);
//   fftw_complex *x2 = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * (NF)));

//   {
//     fftw_plan p;
//     const int rank =1;
//     const int n = N1;
//     const int howmany = N2;
//     const int istride = 1;
//     const int idist = N1;
//     const int ostride = 1;
//     const int odist = N1;
//     p = fftw_plan_many_dft(rank,&n,howmany,x1,NULL,istride,idist,x2,NULL,ostride,odist,sign,FFTW_ESTIMATE);
//     fftw_execute(p);
//     fftw_destroy_plan(p);
//   }
//   fftw_free(x1);


//   null_small_values(x2,NF,neglect);

//   /////////////
//   LAvector <double> datavec("datavec");
//   datavec.perline(1);
//   datavec.textformat(text_nobraces);

//   LAvector<complex<double> > v2(NF,"v2");
//   for (unsigned int i=0; i<Npts; i++) {
//      v2[i] = std::complex<double>(x2[i][0],x2[i][1]);
//   }  
//   datavec.resize() = real(v2);
//   datavec.perline(1);
//   save(datavec,"v2.R.out");
//   datavec = imag(v2);
//   save(datavec,"v2.I.out");
//   ///////////////

//   const unsigned int NF1 = N1;
//   const unsigned int NF2 = N2;

//   fftw_complex *x3 = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * NF));

//   {
//     fftw_plan p;
//     const int rank =1;
//     const int n = N2;
//     const int howmany = N1;
//     const int istride =N1;
//     const int idist =1;
//     const int ostride =1;
//     const int odist =N2;
//     p = fftw_plan_many_dft(rank,&n,howmany,x2,NULL,istride,idist,x3,NULL,ostride,odist,sign,FFTW_ESTIMATE);
//     fftw_execute(p);
//     fftw_destroy_plan(p);
//   }


//   for (unsigned int i=0; i<Npts; i++) {
//     x2[i][0] = 0;
//     x2[i][1] = 0;
//   }

//   {
//     fftw_plan p;
//     const int rank =1;
//     const int n = N2;
//     const int howmany = N1;
//     const int ostride =N1;
//     const int odist =1;
//     const int istride =1;
//     const int idist =N2;
//     p = fftw_plan_many_dft(rank,&n,howmany,x3,NULL,istride,idist,x2,NULL,ostride,odist,+1,FFTW_ESTIMATE);
//     fftw_execute(p);
//     fftw_destroy_plan(p);
//   }


//   ////////////////////////
//   for (unsigned int i=0; i<Npts; i++) {
//      v2[i] = std::complex<double>(x2[i][0],x2[i][1])/double(N2);
//   }  
//   datavec.resize() = real(v2);
//   datavec.perline(1);
//   save(datavec,"v2B.R.out");
//   datavec = imag(v2);
//   save(datavec,"v2B.I.out");
//   ////////////////////

//   fftw_free(x2);
//   null_small_values(x3,NF,neglect);

//   unsigned int k =0;
//   for( int n=-maxNF1; n <= int(maxNF1); n+=N1delta){
//     for( int m = -maxNF2; m <= int(maxNF2); m+=N2delta){
//        const double C= scale/double(Npts);
//        int j1=0;
//        int j2=0;
      
//        if (n>=0) 
// 	  j1 = NF2*(n);
//        else
// 	  j1 = NF2*(n+NF1);
//        if (m>=0)
// 	  j2 = m;
//        else
// 	  j2 = m+NF2;
       
//        int j = j1 + j2;
       
//        double tempR = C*x3[j][0];
//        double tempI = C*x3[j][1];
       
//        if ((!mode00) && (n == 0) && (m == 0)) 
// 	  ;
//        else 
// 	  vF[k++] = std::complex<double>(tempR,tempI);
//     }
//   }

//   fftw_free(x3);




// }










////////////////////////////////////////////////////////////////////////////
//  2D inverse FFT to a COMPLEX VECTOR
////////////////////////////////////////////////////////////////////////////
//use N1 = Nphi, N2=Ntheta
void ifft2d(LAvector<complex<double> >& v, const LAvector<complex<double> >& vF, 
	    const unsigned int N1, const unsigned int N2,
	    const unsigned int maxNF1, const unsigned int maxNF2,
	    const unsigned int N1delta, const unsigned int N2delta,
	    const double neglect, const double scale, const bool mode00)
{  



  const unsigned int Npts=N1*N2;
  const unsigned int NF=N1*(N2);

   const int sign = +1;
   //   printcr("ifft of complex LAvector");


  if (maxNF1 > (N1/2)) {
    cerr<<"bad value given for maxNF1 ("<<maxNF1<<")"<<endl;
    return;
  }
  if (maxNF2 > (N2/2)) {
    cerr<<"bad value given for maxNF2 ("<<maxNF2<<")"<<endl;
    return;
  }


  fftw_complex *x1 = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * Npts));
  for (unsigned int i=0; i<Npts; i++) {
    x1[i][0] = 0;
    x1[i][1] = 0;
  }

  const unsigned int NF1 = N1;
  const unsigned int NF2 = N2;

  unsigned int k =0;
  for( int n=-maxNF1; n <= int(maxNF1); n+=N1delta){
    for( int m = -maxNF2; m <= int(maxNF2); m+=N2delta){
       int j1=0;
       int j2=0;
      
       if (n>=0) 
	  j1 = NF2*(n);
       else
	  j1 = NF2*(n+NF1);
       if (m>=0)
	  j2 = m;
       else
	  j2 = m+NF2;
       
       int j = j1 + j2;
       
       if ((!mode00) && (n == 0) && (m == 0)) {
	  x1[j][0] = 0;
	  x1[j][1] = 0;
       } else {
	  x1[j][0] = vF[k].real();
	  x1[j][1] = vF[k].imag();
	  k++;
       }

    }
  }



  fftw_complex *x2 = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * (NF)));
  {
    fftw_plan p;
    const int rank =1;
    const int n = N2;
    const int howmany = N1;
    const int ostride =N1;
    const int odist =1;
    const int istride =1;
    const int idist =N2;
    p = fftw_plan_many_dft(rank,&n,howmany,x1,NULL,istride,idist,x2,NULL,ostride,odist,sign,FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  }
  fftw_free(x1);

  null_small_values(x2,NF,neglect);

  fftw_complex *x3 = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * (NF)));

  {
    fftw_plan p;
    const int rank =1;
    const int n = N1;
    const int howmany = N2;
    const int ostride = 1;
    const int odist = N1;
    const int istride = 1;
    const int idist = N1;
    p = fftw_plan_many_dft(rank,&n,howmany,x2,NULL,istride,idist,x3,NULL,ostride,odist,sign,FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  }
  fftw_free(x2);

  null_small_values(x3,NF,neglect);


  for (unsigned int i=0; i<Npts; i++) {
     v[i] = std::complex<double>(x3[i][0],x3[i][1])*scale;
  }  

  fftw_free(x3);

}







////////////////////////////////////////////////////////////////////////////
//  2D inverse FFT to a REAL VECTOR
// This function just throws away the imaginary part of the vector
// after it has been transformed to ordinary space, assuming that
// the user has passed a fourier series vF that obeys the
// symmetry exhibited by teh fourier series of real vectors
////////////////////////////////////////////////////////////////////////////
//use N1 = Nphi, N2=Ntheta
void ifft2d(LAvector<double>& v, const LAvector<complex<double> >& vF, 
	    const unsigned int N1, const unsigned int N2,
	    const unsigned int maxNF1, const unsigned int maxNF2,
	    const unsigned int N1delta, const unsigned int N2delta,
	    const double neglect, const double scale, const bool mode00)
{  



  const unsigned int Npts=N1*N2;
  const unsigned int NF=N1*(N2);

   const int sign = +1;
   //   printcr("ifft of complex LAvector");


  if (maxNF1 > (N1/2)) {
    cerr<<"bad value given for maxNF1 ("<<maxNF1<<")"<<endl;
    return;
  }
  if (maxNF2 > (N2/2)) {
    cerr<<"bad value given for maxNF2 ("<<maxNF2<<")"<<endl;
    return;
  }


  fftw_complex *x1 = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * Npts));
  for (unsigned int i=0; i<Npts; i++) {
    x1[i][0] = 0;
    x1[i][1] = 0;
  }

  const unsigned int NF1 = N1;
  const unsigned int NF2 = N2;

  unsigned int k =0;
  for( int n=-maxNF1; n <= int(maxNF1); n+=N1delta){
    for( int m = -maxNF2; m <= int(maxNF2); m+=N2delta){
       int j1=0;
       int j2=0;
      
       if (n>=0) 
	  j1 = NF2*(n);
       else
	  j1 = NF2*(n+NF1);
       if (m>=0)
	  j2 = m;
       else
	  j2 = m+NF2;
       
       int j = j1 + j2;
       
       if ((!mode00) && (n == 0) && (m == 0)) {
	  x1[j][0] = 0;
	  x1[j][1] = 0;
       } else {
	  x1[j][0] = vF[k].real();
	  x1[j][1] = vF[k].imag();
	  k++;
       }

    }
  }



  fftw_complex *x2 = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * (NF)));
  {
    fftw_plan p;
    const int rank =1;
    const int n = N2;
    const int howmany = N1;
    const int ostride =N1;
    const int odist =1;
    const int istride =1;
    const int idist =N2;
    p = fftw_plan_many_dft(rank,&n,howmany,x1,NULL,istride,idist,x2,NULL,ostride,odist,sign,FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  }
  fftw_free(x1);

  null_small_values(x2,NF,neglect);

  fftw_complex *x3 = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * (NF)));

  {
    fftw_plan p;
    const int rank =1;
    const int n = N1;
    const int howmany = N2;
    const int ostride = 1;
    const int odist = N1;
    const int istride = 1;
    const int idist = N1;
    p = fftw_plan_many_dft(rank,&n,howmany,x2,NULL,istride,idist,x3,NULL,ostride,odist,sign,FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  }
  fftw_free(x2);

  null_small_values(x3,NF,neglect);


  for (unsigned int i=0; i<Npts; i++) {
     v[i] = scale*x3[i][0];
  }  

  fftw_free(x3);

}








///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////








void fft2d(const double * const v, std::complex<double> * const vF,
	   const unsigned int N1, const unsigned int N2,
	   const unsigned int maxNF1, const unsigned int maxNF2,
	   const unsigned int v_stride, const unsigned int vF_stride,
	   const int sign,const double neglect)
{  


   //   printcr("fft2d real* complex*");

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

//  const unsigned int NF1 = min(2*maxNF1+1,N1);
//  const unsigned int NF2 = N2;
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
	    if (m<=0)
	       j = abs(n)*N2-m;
	    else
	       j = abs(n)*N2-m+N2;
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

   //   printcr("fft2d complex* complex*");

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
  const unsigned int NF2 = N2;

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
       int j1=0;
       int j2=0;
      
       if (n>=0) 
	  j1 = NF2*(n);
       else
	  j1 = NF2*(n+NF1);
       if (m>=0)
	  j2 = m;
       else
	  j2 = m+NF2;
       
       int j = j1 + j2;
       
       double tempR = C*x3[j][0];
       double tempI = C*x3[j][1];
       
       vF[k++*vF_stride] = std::complex<double>(tempR,tempI);
    }
  }

  fftw_free(x3);

}







////////////////////////////////////////////////////



// note that this function clobbers M
// this version allows different parameters for plasma and coil
void fft_of_M(Matrix<double>& M, Matrix<complex<double> >& MF, 
	      const unsigned int Nphi_P, const unsigned int Ntheta_P,
	      const unsigned int Nnn_P, const unsigned int Nmm_P,
	      const unsigned int Nphi_C, const unsigned int Ntheta_C,
	      const unsigned int Nnn_C, const unsigned int Nmm_C,
	      const unsigned int Ndelta, const unsigned int Mdelta,
	      const double neglect)
{

   MF.resize(0,0);

   unsigned int NF_P = (2*Nnn_P+1)*(2*Nmm_P+1);
   unsigned int NF_C = (2*Nnn_C+1)*(2*Nmm_C+1);
   const unsigned int Npts_P = Nphi_P*Ntheta_P;
   //  const unsigned int Npts_C = Nphi_C*Ntheta_C;
   const double PI = 3.141592653589793;

   Matrix<complex<double> > MxF(Npts_P,NF_C,"fft_of_M::MxF");

   // Fourier transform coil side at each point on plasma surface, creating the Fourier-Green functions
	  
   for (unsigned int i=0; i<Npts_P; i++) 
      fft2d( &(M(i,0)), &(MxF(i,0)), Nphi_C, Ntheta_C, Nnn_C, Nmm_C, 1, 1, +1, neglect);

   // free up memory used by M
   M.resize(0,0);

   Matrix<complex<double> > MFtemp(NF_P,NF_C,"MFtemp"); 

   // Now transform each Fourier-Green function

   for (unsigned int k=0; k<NF_C; k++) 
      fft2d( &(MxF(0,k)), &(MFtemp(0,k)), Nphi_P, Ntheta_P, Nnn_P, Nmm_P, NF_C, NF_C, -1, neglect);


   // free up memory used by MxF
   MxF.resize(0,0);


   // remove n=0,m=0 mode and unwanted modes and scale output

   unsigned int NFrow;
   LAvector<unsigned int> Fmap_row;
   make_mode_map(NFrow,Nnn_P,Nmm_P,Ndelta,Mdelta,false,NF_P,Nnn_P,Nmm_P,1,1,true,Fmap_row);
   unsigned int NFcol;
   LAvector<unsigned int> Fmap_col;
   make_mode_map(NFcol,Nnn_C,Nmm_C,Ndelta,Mdelta,false,NF_C,Nnn_C,Nmm_C,1,1,true,Fmap_col);

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
void fft_of_M(Matrix<double>& M, Matrix<complex<double> >& MF, 
	      const unsigned int Nphi, const unsigned int Ntheta,
	      const unsigned int Nnn, const unsigned int Nmm,
	      const unsigned int Ndelta, const unsigned int Mdelta,
	      const double neglect)
{

   MF.resize(0,0);

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

   MF.resize(0,0);

   const unsigned int Npts_C = Nphi*Ntheta;
   const unsigned int Npts_P = M.Nrows();
   unsigned int NF = (2*Nnn+1)*(2*Nmm+1);
   const double PI = 3.141592653589793;
   Matrix<complex<double> > MFtemp(Npts_P,NF,"half_fft_of_M::MFtemp");

   for (unsigned int i=0; i<Npts_P; i++) 
      fft2d( &(M(i,0)), &(MFtemp(i,0)), Nphi, Ntheta, Nnn, Nmm, 1, 1, +1, neglect);

   // free up memory used by M
   //   M.resize(0,0);


   LAvector<unsigned int> Fmap_col;
   unsigned int NFcol = 0;


   // remove (0,0) mode, non-harmonics of fundamental, etc.
   make_mode_map(NFcol,Nnn,Nmm,Ndelta,Mdelta,false,NF,Nnn,Nmm,1,1,true,Fmap_col);


   MF.resize(Npts_P,NFcol); 
   for (unsigned int r =0; r<Npts_P; r++) {
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

   dispcr(Omega.Nrows());
   dispcr(Omega.Ncols());
   Omega.resize(0,0);
//   const unsigned int Npts = Nphi*Ntheta;
   unsigned int NF = (2*Nnn+1)*(2*Nmm+1);
   const double PI = 3.141592653589793;
   //  Matrix<complex<double> > Ftemp1(NF-1,Npts,"Omega_half_fft::Ftemp1");
   const unsigned int NF2 = Bdotgradf.Nrows();
   Matrix<complex<double> > Ftemp2(NF2,NF,"Omega_half_fft::Ftemp2");

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

   Omega.resize(NFcol,NF2); 
   dispcr(Omega.Nrows());
   dispcr(Omega.Ncols());
   for (unsigned int r =0; r<NF2; r++) {
      for (unsigned int c =0; c<NFcol; c++) {
	 const double C2= 2*PI;
	 Omega(c,r) = C2*Ftemp2(r, Fmap_col[c]);
      }
   }


   return;
}




void omega_half_fft2(Matrix<complex<double> >& FuncMatF, Matrix<complex<double> >& Omega, 
		    const unsigned int Nphi, const unsigned int Ntheta,
		    const unsigned int Nnn, const unsigned int Nmm,
		    const unsigned int Ndelta, const unsigned int Mdelta,
		    const double neglect)
{



dispcr(Nphi);
dispcr(Ntheta);
dispcr(Nnn);
dispcr(Nmm);
dispcr(Ndelta);
dispcr(Mdelta);

   Omega.resize(0,0);
//   const unsigned int Npts = Nphi*Ntheta;
   unsigned int NF = (2*Nnn+1)*(2*Nmm+1);
   const double PI = 3.141592653589793;
   const unsigned int NFc = FuncMatF.Ncols();

   Matrix<complex<double> > Ftemp(NF,NFc,"Omega_half_fft::Ftemp");


   for (unsigned int c=0; c<NFc; c++) 
      fft2d( &(FuncMatF(0,c)), &(Ftemp(0,c)), Nphi, Ntheta, Nnn, Nmm, FuncMatF.Ncols(), Ftemp.Ncols(), -1, neglect);
  

   // free up memory
   FuncMatF.resize(0,0);



   // remove n=0,m=0 mode and unwanted modes  and scale output

   unsigned int NFrow;
   LAvector<unsigned int> Fmap_row;
   make_mode_map(NFrow,Nnn,Nmm,Ndelta,Mdelta,false,NF,Nnn,Nmm,1,1,true,Fmap_row);

   dispcr(NFrow);

   Omega.resize(NFrow,NFc); 
   dispcr(Omega.Nrows());
   dispcr(Omega.Ncols());
   for (unsigned int r =0; r<NFrow; r++) {
      for (unsigned int c =0; c<NFc; c++) {
	 const double C2= 2*PI;
	 Omega(r,c) = C2*Ftemp(Fmap_row[r],c);
      }
   }


   return;
}












