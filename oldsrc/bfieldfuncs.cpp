#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>

#include "matricks.hpp"
#include "coils.hpp"
#include "bfieldfuncs.hpp"

using namespace std;
using namespace Matricks;



const double pi = 4*atan2(1.0,1.0);
const double muOverPi = 4e-7;


//calculate b field of a straight wire along z-axis

void bfieldwire(p3vector<double>& B, const p3vector<double>& X, const double I) {
  const double Bo  = mu0/(2.0*PI)*I;
  
  const double r = sqrt(sqr(X.x()) + sqr(X.y()));
  const double Bphi = Bo / r;
  
  const double phi = atan2(X.y(),X.x());
  B.x() = -Bphi*sin(phi);
  B.y() = +Bphi*cos(phi);
  B.z() = 0.0;
}
 

//calculate b field of a straight wire in z direction at position x0,y0

void bfieldwire_general(Matricks::p3vector<double>& B, const Matricks::p3vector<double>& X,const double x0, const double y0, const double I) {

  const double Bo  = mu0/(2.0*PI)*I;
  const Matricks::p3vector<double> deltaX(X.x()-x0, X.y()-y0, X.z());
  const double r = sqrt(sqr(deltaX.x()) + sqr(deltaX.y()));
  const double Bphi = Bo / r;
  
  const double phi = atan2(deltaX.y(),deltaX.x());
  B.x() = -Bphi*sin(phi);
  B.y() = +Bphi*cos(phi);
  B.z() = 0.0;

}


// calculate b field of a loop of current


// implementation of auxiliary functions

inline double delK(const double k2)
{
  const double t = 1. - k2;
  if(t < 1e-9)
    {
      //cout.precision(15);
      cout << "Warning - loss of accuracy in delK(k2=" << k2 << ")\n";
      cout << "t = 1 - k2 =" << t << endl;
    }
  const double delK = (((0.032024666*t+0.054544409)*t+0.097932891)*t+
	  1.3862944)-(((0.010944912*t+0.060118519)*t+
		       0.12475074)*t+0.5)*log(t);
  //cout << "delK = " << delK << endl;
  return delK;
}


inline double delE(const double k2)
{
  const double t = 1.-k2;
  if(t < 1e-9)
    {
      //cout.precision(15);
      cout << "Warning - loss of accuracy in delE(k2=" << k2 << ")\n";
      cout << "t = 1 - k2 =" << t << endl;
    }
  const double delE = ((0.040905094*t+0.085099193)*t+0.44479204)*t+
    1.0-(((0.01382999*t+0.08150224)*t+
	  0.24969795)*t)*log(t);
  //cout << "delE = " << delE << endl;
  return delE;
}

 
void bLocal(const p3vector<double>& position,
	    const double & coilRadius,
	    p3vector<double>& bLocal)
{
 
  const double x = position[0];
  const double y = position[1];
  const double z = position[2];
  
  const double rho = sqrt(x*x+y*y);
  
  const double alpha = rho/coilRadius;
  const double beta = z/coilRadius;
  const double q = (1.+alpha)*(1.+alpha)+beta*beta;
  const double kVal = 4*alpha/q;
  //  cout.precision(20);
  //cout << "alpha,beta,q,kVal = " << alpha << ", " << beta << ", " << q << ", " << kVal << endl;
  const double alpha2 = alpha*alpha;
  const double beta2 = beta*beta;
  const double coilRadius2 = coilRadius*coilRadius;
  if(kVal>1e-7)
    {
      const double gamma = z/rho;
      const double factor = muOverPi/(2.0*coilRadius*sqrt(q));
      const double ellipE = delE(kVal);
      const double ellipK = delK(kVal);
      bLocal[2] = factor*(ellipE*(1.-alpha2-beta2)/(q-4*alpha)+ellipK);
      const double bR = factor*gamma*(ellipE*(1+alpha2+beta2)/(q-4*alpha)-ellipK);
      //cout << "ellipE, ellipK,bZ,bR = " << ellipE << ", " << ellipK << ", " 
      //	   << bLocal[2] << ", " << bR << endl;
      bLocal[0] = bR*x/rho;
      bLocal[1] = bR*y/rho;
    }
  else
    {
      const double z2 = z*z;
      if(rho>1e-100)
	{
	  const double length = coilRadius2+z2;
	  const double bR = muOverPi*pi*3*z*rho/(4*length*length*sqrt(1+beta2)) + 
	    5*coilRadius*muOverPi*pi*z*(11*coilRadius2-3*z2)*rho*rho*rho/
	    (8*pow(coilRadius2+z2,4)*sqrt(1+beta2));
	  bLocal[0] = bR*x/rho;
	  bLocal[1] = bR*y/rho;
	}
      else
	{
	  bLocal[0] = 0;
	  bLocal[1] = 0;
	}
      bLocal[2] = coilRadius*muOverPi*pi/(2*(coilRadius2+z2)*sqrt(1+beta2)) + 
	3*coilRadius*muOverPi*pi*(coilRadius2-4*z2)*rho*rho/
	(8*pow(coilRadius2+z2,3)*sqrt(1+beta2));
      
    }

}

 




void bfieldcoil(const p3vector<double>& coilPosition,
		const p3vector<double>& coilOrient,
		const p3vector<double>& coilPlaneOrient,
		const double& coilRadius,
		const p3vector<double>& fieldPosition,
		p3vector<double>& bField,
		const double I)
{
//  const double coilHeight = 0;
//  const double coilWidth = 0;
//  const int numWindingsInR = 1;
//  const int numWindingsInZ = 1;

  p3vector<double> pZero("pZero");
  
  // compute vector difference between field and coil positions.
  pZero = fieldPosition - coilPosition;
  
  // Use the normalized coil orientation vector (normal to the coil plane) as the
  // z-axis;
  p3vector<double> unitZ( "unitZ");
  unitZ = coilOrient/norm(coilOrient);

  // Use the normalized coil plane orientation vector as the x-axis.
  p3vector<double> unitX("unitX");
  unitX = coilPlaneOrient/norm(coilPlaneOrient);

  //Check to see that unitX and unitZ are perp. as expected.
  if(abs(unitX|unitZ)>10e-10) {
    cout << "unitX and unitZ are not perp!!\n";
    dispcr(unitX);
    dispcr(unitZ);
  }

  p3vector<double> unitY("unitY");
  unitY[0] = unitZ[1]*unitX[2] - unitZ[2]*unitX[1];
  unitY[1] = unitZ[2]*unitX[0] - unitZ[0]*unitX[2];
  unitY[2] = unitZ[0]*unitX[1] - unitZ[1]*unitX[0];

  // now normalize unitY
  unitY = unitY/norm(unitY);
  

  // now calculate the projection of pZero onto the new unit vectors;
  // these are now the coordinates in the local coordinate system.
  p3vector<double> local("local");
  local[0] = unitX|pZero;
  local[1] = unitY|pZero;
  local[2] = unitZ|pZero;

  p3vector<double> localBField("localBField");

  localBField = 0;

  bLocal(local,coilRadius,localBField);	
  
  // Now because the unitX and unitZ vectors are orthonormal, transformation fo the B-field 
  // vector back to the lab frame is simple.  Note: the local y coordinate of B is zero by construction.
  for(int i =0; i < 3; i++) {
    bField[i] = localBField[0]*unitX[i] + localBField[1]*unitY[i] + localBField[2]*unitZ[i];
    bField[i] = I * bField[i];
  }

}






