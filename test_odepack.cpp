

// Standard C libraries
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <tgmath.h>

// Standard C++ libraries

#include <iostream>
#include <complex>

using namespace std;


// LSODE (FEX, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
//        ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JEX, MF)
// from ODEPAC

typedef int FINT;



//      SUBROUTINE DLSODE (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
//                       ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
extern "C" {
  void dlsode_(
	       void(*)(const int &, const double&, const double [], double []),
	       const int & NEQ,
	       double *Y,
	       double &T,
	       double &TOUT,
	       int &ITOL,
	       double &RTOL,
	       double &ATOL,
	       int &ITASK,
	       int &ISTATE,
	       int &IOPT,
	       double *RWORK,
	       const int &LRW,
	       int *IWORK,
	       const int &LIW,
	       void(*)(const int&, const double&, const double [], const int&, const int&, double*, const int&),
	       int &MF
	       );
}


void f1 (const int& neq, const double& t, const double y[], double ydot[]) {
  ydot[0] = y[1];
  ydot[1] = 3.0L * (1.0L - y[0]*y[0]) * y[1] - y[0];
  return;
}

//SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
//   INTEGER  NEQ, ML, MU, NROWPD
//   DOUBLE PRECISION  T, Y(*), PD(NROWPD,*)
//
//  Fortran is Column major 
// A(1,1) A(2,1) A(3,1) A(1,2) A(2,2) A(3,2) A(1,3) A(2,3) A(3,3)
//  C arrays in row-major order:
// A[0][0] A[0][1] A[0][2] A[1][0] A[1][1] A[1][2] A[2][0] A[2][1] A[2][2]

void jac1 (const int& neq, const double& t, const double y[], const int& ml, const int& mu, double* pd, const int& nrowpd) {
  //  printf("neq=%d, t=%f, x=y[0]=%f,  xdot=y[1]=%f, ml=%d, mu=%d, nrowpd=%d\n",neq,t,y[0],y[1],ml,mu,nrowpd);

  // A(1,1) 
  pd[0] = 0.0L;
  // A(2,1)
  pd[1] = -6.0L * y[0]*y[1] - 1.0L;
  // A(1,2)
  pd[2] = 1.0L;
  // A(2,2)
  pd[3] = 3.0L * (1.0L - y[0]*y[0]);  
  
  return;
}




int main (int argc, char *argv[])
{
  printf("\n");
  printf("C++ version of Demonstration program for DLSODE package\n");

  const int neq = 2;
  const int liw = 64;
  const int lrw = 1024;

  int  i, iopar, iopt, iout, istate, itask, itol;
  int leniw, lenrw, mband, meth, mf, miter;
  int ml, mu, nerr, nfe, nfea, nje, nout, nqu, nst;
  int lout = 6;

  double atol, er, erm, ero, hu, rtol, t, tnext;
  const double tstart = 1.39283880203L;
  const double dtout = 2.214773875L;

  double y[neq];
  int iwork[liw];
  double rwork[lrw];
	  
  nerr = 0;

  // First problem
  printf("Problem 1:  Van der Pol oscillator:\n");
  nout = 4;
  printf("  xdotdot - 3*(1 - x**2)*ẋ + x = 0,\n");
  printf("                x(0) = 2, ẋ(0) = 0,\n");
  printf(" neq = %3d\n",neq);
  printf(" itol = %3d    rtol = % 10.1f    atol = % 10.1f\n", itol,  rtol, atol);


  //             NST    IWORK(11)  Number of steps taken for the problem so
  //                               far.
  //             NFE    IWORK(12)  Number of f evaluations for the problem
  //                               so far.
  //             NJE    IWORK(13)  Number of Jacobian evaluations (and of
  //                               matrix LU decompositions) for the problem
  //                               so far.
  //             NQU    IWORK(14)  Method order last used (successfully).
  //             LENRW  IWORK(17)  Length of RWORK actually required.  This
  //                               is defined on normal returns and on an
  //                               illegal input return for insufficient
  //                               storage.
  //             LENIW  IWORK(18)  Length of IWORK actually required.  This
  //                               is defined on normal returns and on an
  //                               illegal input return for insufficient
  //                               storage.



  
  for(meth = 1; meth<=2; meth++) {
    for(miter = 0; miter<=3; miter++) {
      mf = 10*meth + miter;
      printf("Solution with mf = %d \n", mf);
      printf("%15s %15s %15s %7s %4s %4s %4s %4s %6s %6s \n", "t","x","xdot","istate","nst","nfe","nje","nqu","lenrw","leniw");
      y[0] = 2.0L;
      y[1] = 0.0L;
      itask = 1;
      istate = 1;
      t = 0.0L;
      tnext = tstart;
      ero = 0.0L;
      itol = 1;
      rtol = 0.0e0L;
      atol = 1.0e-6L;
      iopt = 0;

      for(iout = 1; iout<=nout; iout++) {
	dlsode_(f1,neq,y,t,tnext,itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,jac1,mf);
	nst = iwork[11-1];
	nfe = iwork[12-1];
	nje = iwork[13-1];
	nqu = iwork[14-1];
	lenrw = iwork[17-1];
	leniw = iwork[18-1];
	printf("% 15.5e % 15.5e % 15.5e %7d %4d %4d %4d %4d %6d %6d\n", t, y[1-1], y[2-1], istate, nst, nfe, nje, nqu, lenrw, leniw);
	if (istate < 0) {
	  continue;
	}
	if (iout % 2 == 0) {
	  er = abs(y[1-1])/atol;
	  ero = max(ero,er);
	  if (er > 1000.0L) {
	    printf(" Warning: error exceeds 1000 * tolerance.  (er=%f)\n", er);
	    nerr++;
	  }
	}
	tnext = tnext + dtout;
      } // iout


      if (istate <  0) {
	nerr = nerr + 1;
      }
      nfea = nfe;
      if (miter == 2) {
	nfea = nfe - neq*nje;
      }
      if (miter == 3) {
	nfea = nfe - nje;
      }

      printf("Final statistics for this run:\n");
      printf("  rwork size = %3d   iwork size = %3d \n", lenrw, leniw);
      printf("  number of steps = %3d\n", nst);
      printf("  number of f-s   = %3d\n", nfe);
      printf("  (excluding J-s) = %3d\n", nfea);
      printf("  number of J-s   = %3d\n", nje);
      printf("  error overrun = % 10.2f \n", ero);

      
    } // miter
  } // meth
}
