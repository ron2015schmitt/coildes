#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>

#include "matricks.hpp"
#include "flsode.hpp"
#include "calcPoincare.hpp"



using namespace std;
using namespace Matricks;


typedef struct {
   double phisection;
   double phishape;
} poincareParamsF;






// LSODE (FEX, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
//        ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JEX, MF)
// from ODEPACK, from LLNL (Lawerence Livermore Nation Labs)

extern "C" {void lsode_(void(*)(const FINT &,const double&,const double*,double*),
			FINT & NEQ,
			double *Y,
			double &T,
			double &TOUT,
			FINT &ITOL,
			double &RTOL,
			double &ATOL,
			FINT &ITASK,
			FINT &ISTATE,
			FINT &IOPT,
			double *RWORK,
			FINT &LRW,
			FINT *IWORK,
			FINT &LIW,
			FINT &JEX,
			FINT &MF, FINT &DEBUG);}



//positions and return valuer are given as (r,z,phi) (cylindrical coords)

bool calcPoincare(Vector<p3vector<double> >& position,
		  const p3vector<double>& startPositionRZP,
		  int& numTransits,
		  const int pointsPerTransit)
{

   bool debug = false;
   // debug = true;
  
   if (debug) {
      cout << "calcPoincare started..." << endl;
      dispcr(startPositionRZP);
      dispcr(numTransits);
      dispcr(pointsPerTransit);
      cr();
   }

   // some constants
//   const double pi =          3.141592653589793;
//   const long double piL =    3.1415926535897932385L;
   const long double twopiL = 6.2831853071795864769L;
   const long double circuitL = twopiL;

   FINT neq = 2;
   FINT mf = 10;

   // test to see if this is better
   //mf = 22;


   FINT lrw;
   FINT liw;
   switch (mf) {
   case 10:
      lrw = 20+16*neq;
      liw = 20;
      break;
   case 21:
   case 22:
      lrw = 22 + 9*neq + neq*neq;
      liw = 20 +neq;
      break;
   case 24:
   case 25:
      break;
   }
   FINT jac=0;
   double aTol = 1e-10;
   double rTol = 0;
   FINT iTol =1;
   FINT iTask=1;
   FINT iState=1;
   double y[neq];
   double rWork[lrw];
   FINT iWork[liw];
   FINT lsode_debug = 1;
   double t;
   double tout;

   // aliases
   double & phi = t;
   double & phiNext = tout;
   double & r = y[0];
   double & z = y[1];
  
   // debugging outputs from LSODE
   const double & hu = rWork[11-1];
   const double & hcur = rWork[12-1];
   const double & tcur = rWork[13-1];
   const double & tolsf = rWork[14-1];
   const FINT & nst = iWork[11-1];
   const FINT & nfe = iWork[12-1];
   const FINT & nje = iWork[13-1];
   const FINT & nqu = iWork[14-1];
   const FINT & nqcur = iWork[15-1];
   const FINT & imxer = iWork[16-1];
   const FINT & lenrw = iWork[17-1];
   const FINT & leniw = iWork[18-1];

   // optional inputs to LSODE (set iOpt=1 to use)
   FINT iOpt = 0;
   double & h0 = rWork[5-1];
   double & hmax = rWork[6-1];
   double & hmin = rWork[7-1];
   FINT & maxord = iWork[5-1];
   FINT & mxstep = iWork[6-1];
   FINT & mxhnil = iWork[7-1];

   /////////////////////////////////////////////
   // debugging.............
   iOpt =0;
   h0 = 5e-4;
   hmax = 5e-4;
   hmin =0;
   mxstep = 300;
   mxhnil = 10;
   /////////////////////////////////////////////

  

   int count = 0;
   bool fail = false;
   const unsigned int Npts = numTransits*pointsPerTransit;


   // display counters
   unsigned int dispcount = 0;
   unsigned int dispcnt = 0;
   unsigned int dispNpts = numTransits*pointsPerTransit;
   unsigned int dispCheckPoint = static_cast<unsigned int>(dispNpts*0.01);
   if (dispCheckPoint==0)
      dispCheckPoint = 1;

   bool abort_loop = false;

   printcr("Starting Trajectory Calculation.");
   printcr("Type ?<Enter> for interactive options.");

   for(int i = 0; i < numTransits; i++) {
      for(int j = 0; j < pointsPerTransit; j++) {

	 // check for user input from keyboard
	 streamsize size = cin.rdbuf()->in_avail();
	 if (size>0) {
	    string instring;
	    cin >> instring;
	    if ( instring == "abort") {
	       abort_loop =true;
	       printcr("*** 'abort' command received!");
	       printcr("This Trajectory will be aborted when present circuit is completed,");
	       cout << " for a total of "<<i+1<<" of "<<numTransits<<" circuits."<<endl;
	       cout << " Need to compute "<<(pointsPerTransit-j)<<" more points (out of "<<(pointsPerTransit)<<") to finish present circuit"<<endl;
	       dispNpts = (pointsPerTransit-j);
	       dispcount = 0;
	       dispcnt = 0;
	       dispCheckPoint = static_cast<unsigned int>(dispNpts*0.01);
	       if (dispCheckPoint==0)
		  dispCheckPoint = 1;
	       printcr("% Progress Counter has been reset.");
	    }else {
	       printcr("Interactive options. Type code word then <Enter>.");
	       printcr("Valid Code words:");
	       printcr(" 'abort' - abort code run early and save results.");
	       //	       double tcompleted = double(i)+double(j)/double(pointsPerTransit);
	       double pctcompleted = 100*(double(j)+double(i)*double(pointsPerTransit))/double(Npts);
	       cout << i<<" "<<j<<"/"<<pointsPerTransit <<" Trajectories completed ("<<pctcompleted<<"%)." <<endl;
	    }
	 }

	 // % display counter
	 dispcnt++;
	 if (++dispcount == dispCheckPoint) {
	    print(Matricks::round(double(dispcnt)/dispNpts*100));cout <<" %"<<endl;
	    dispcount =0;
	 }

	 if (count == 0) {
	    // initialize the first time through
	    position[count++] = startPositionRZP;
	    r = startPositionRZP[0];
	    z = startPositionRZP[1];
	    phi = startPositionRZP[2];
	    continue;
	 }

	 const double phi_wrapped = double(j)/double(pointsPerTransit)*double(circuitL);
	 phiNext =  phi_wrapped + double(circuitL*double(i));


	 if (debug) {
	    cr();
	    printcr("------------------------------------------------------------------------");
	    disp(i);disp(j);cr();
	    cout << "INPUT to lsode"<<endl;
	    disp(phi);disp(phiNext);cr();
	    disp(r);disp(z);cr();
	    disp(iState);cr();
	    disp(neq);disp(iTol);disp(rTol);disp(iTask);disp(iOpt);cr();
	    disp(lrw);disp(liw);disp(jac);disp(mf);cr();
	    if (iOpt) {
	       disp(h0);disp(hmax);disp(hmin);cr();
	       disp(maxord);disp(mxstep);disp(mxhnil);cr();
	    }
	    cout.flush();
	 }   

   
	 lsode_(flsode,neq,y,t,tout,iTol,rTol,aTol,iTask,iState,iOpt,rWork,lrw,iWork,liw,jac,mf,lsode_debug);
	  
	  

	 if (debug |(iState!=2)) {
	    cout << "OUTPUT from lsode"<<endl;
	    disp(phi);disp(phiNext);cr();
	    disp(r);disp(z);cr();
	    disp(iState);disp(hu);disp(hcur);cr();
	    disp(tcur);disp(tolsf);cr();
	    disp(nst);disp(nfe);disp(nje);disp(nqu);disp(nqcur);cr();
	    disp(imxer);disp(lenrw);disp(leniw);cr();
	    cout.flush();
	    cout << "OUTPUT from flsode (at r,z,phi given above)"<<endl;    
	    double dy_dt[neq];
	    double & dr_dt = dy_dt[0];
	    double & dz_dt = dy_dt[1];
	    flsode(neq, t, y, dy_dt);
	    disp(dr_dt);disp(dz_dt);cr();
	    double xx,yy,zz,Bx,By,Bz,Br,Bphi;
	    disp(r);disp(z);cr();
	    flsode_details(neq, t, y, dy_dt, xx,yy,zz,Bx,By,Bz,Br,Bphi);
	    disp(dr_dt);disp(dz_dt);cr();
	    disp(xx);disp(yy);disp(zz);cr();
	    disp(Bx);disp(By);disp(Bz);cr();
	    disp(Br);disp(Bz);disp(Bphi);cr();
	    double Bmag = sqrt(Bx*Bx+By*By+Bz*Bz);
	    dispcr(Bmag);

	 }
	 if(iState != 2) {
	    fail = true;
	    cout << "calcPoincare failed!"<<" LSODE error code = "<<iState<<"\n\n";
	    disp(i);disp(j);cr();
	    numTransits = i + 1;
	    return fail;
	 }
	 // note that we don't use "phi" because we need
	 // phi to be exactly the same on each circuit for
	 // sorting purposes in post-processing
	 p3vector<double> temp(r,z,phi_wrapped);
	 position[count++] = temp;
      }
      if (abort_loop) {
	 int temp = i + 1;
	 cout<<"Field line trajectory aborted after "<<temp<<" circuits instead of "<<numTransits<< "circuits." <<endl;
	 numTransits = temp;
	 dispcr(numTransits);
	 dispcr(Npts);
	 dispcr(count);
	 return false;
      }
   }
   dispcr(Npts);
   dispcr(count);
   return fail;
}












bool calc_1pt_Poincare(p3vector<double> & finalposition,
		       const p3vector<double>& startPositionRZP)
{

   bool debug = false;
   // debug = true;
  
   if (debug) {
      cout << "calcPoincare started..." << endl;
      dispcr(startPositionRZP);
      cr();
   }

   // some constants
//   const double pi =          3.141592653589793;
//   const long double piL =    3.1415926535897932385L;
   const long double twopiL = 6.2831853071795864769L;
   const long double circuitL = twopiL;
   p3vector<double> position;

   FINT neq = 2;
   FINT mf = 10;

   // test to see if this is better
   //mf = 22;


   FINT lrw;
   FINT liw;
   switch (mf) {
   case 10:
      lrw = 20+16*neq;
      liw = 20;
      break;
   case 21:
   case 22:
      lrw = 22 + 9*neq + neq*neq;
      liw = 20 +neq;
      break;
   case 24:
   case 25:
      break;
   }
   FINT jac=0;
   double aTol = 1e-10;
   double rTol = 0;
   FINT iTol =1;
   FINT iTask=1;
   FINT iState=1;
   double y[neq];
   double rWork[lrw];
   FINT iWork[liw];
   FINT lsode_debug = 1;
   double t;
   double tout;

   // aliases
   double & phi = t;
   double & phiNext = tout;
   double & r = y[0];
   double & z = y[1];

   double phi_wrapped = phiNext;  
   // debugging outputs from LSODE
   const double & hu = rWork[11-1];
   const double & hcur = rWork[12-1];
   const double & tcur = rWork[13-1];
   const double & tolsf = rWork[14-1];
   const FINT & nst = iWork[11-1];
   const FINT & nfe = iWork[12-1];
   const FINT & nje = iWork[13-1];
   const FINT & nqu = iWork[14-1];
   const FINT & nqcur = iWork[15-1];
   const FINT & imxer = iWork[16-1];
   const FINT & lenrw = iWork[17-1];
   const FINT & leniw = iWork[18-1];

   // optional inputs to LSODE (set iOpt=1 to use)
   FINT iOpt = 0;
   double & h0 = rWork[5-1];
   double & hmax = rWork[6-1];
   double & hmin = rWork[7-1];
   FINT & maxord = iWork[5-1];
   FINT & mxstep = iWork[6-1];
   FINT & mxhnil = iWork[7-1];

   /////////////////////////////////////////////
   // debugging.............
   iOpt =0;
   h0 = 5e-4;
   hmax = 5e-4;
   hmin =0;
   mxstep = 300;
   mxhnil = 10;
   /////////////////////////////////////////////

  

//   int count = 0;
   bool fail = false;
//   const unsigned int Npts = 1;



   // initialize the first time through
   r = startPositionRZP[0];
   z = startPositionRZP[1];
   phi = startPositionRZP[2];

   phiNext = finalposition[2];
   phi_wrapped = phiNext;
   while (phi_wrapped > circuitL) 
      phi_wrapped -= circuitL;  

   while (phiNext < phi) 
      phiNext += circuitL;  
  
   disp(phi);disp(phiNext);disp(phi_wrapped);
   cr();
  
   if (debug) {
      cr();
      printcr("------------------------------------------------------------------------");
      cout << "INPUT to lsode"<<endl;
      disp(phi);disp(phiNext);cr();
      disp(r);disp(z);cr();
      disp(iState);cr();
      disp(neq);disp(iTol);disp(rTol);disp(iTask);disp(iOpt);cr();
      disp(lrw);disp(liw);disp(jac);disp(mf);cr();
      if (iOpt) {
	 disp(h0);disp(hmax);disp(hmin);cr();
	 disp(maxord);disp(mxstep);disp(mxhnil);cr();
      }
      cout.flush();
   }   
  
   
   lsode_(flsode,neq,y,t,tout,iTol,rTol,aTol,iTask,iState,iOpt,rWork,lrw,iWork,liw,jac,mf,lsode_debug);
	  
	  

   if (debug |(iState!=2)) {
      cout << "OUTPUT from lsode"<<endl;
      disp(phi);disp(phiNext);cr();
      disp(r);disp(z);cr();
      disp(iState);disp(hu);disp(hcur);cr();
      disp(tcur);disp(tolsf);cr();
      disp(nst);disp(nfe);disp(nje);disp(nqu);disp(nqcur);cr();
      disp(imxer);disp(lenrw);disp(leniw);cr();
      cout.flush();
      cout << "OUTPUT from flsode (at r,z,phi given above)"<<endl;    
      double dy_dt[neq];
      double & dr_dt = dy_dt[0];
      double & dz_dt = dy_dt[1];
      flsode(neq, t, y, dy_dt);
      disp(dr_dt);disp(dz_dt);cr();
      double xx,yy,zz,Bx,By,Bz,Br,Bphi;
      disp(r);disp(z);cr();
      flsode_details(neq, t, y, dy_dt, xx,yy,zz,Bx,By,Bz,Br,Bphi);
      disp(dr_dt);disp(dz_dt);cr();
      disp(xx);disp(yy);disp(zz);cr();
      disp(Bx);disp(By);disp(Bz);cr();
      disp(Br);disp(Bz);disp(Bphi);cr();
      double Bmag = sqrt(Bx*Bx+By*By+Bz*Bz);
      dispcr(Bmag);
     
   }
   if(iState != 2) {
      fail = true;
      cout << "calcPoincare failed!"<<" LSODE error code = "<<iState<<"\n\n";
      return fail;
   }

   p3vector<double> temp(r,z,phi_wrapped);
   finalposition = temp;
  
   return false;
}
