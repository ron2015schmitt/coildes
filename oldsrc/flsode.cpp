#include "flsode.hpp"
#include "matricks.hpp"
#include <math.h>
#include <stdlib.h>
#include <iostream>

// assumes bTotal(X,B) function exists somewhere. bTotal calculates the
// desired B field at the given coord (X).
extern void bTotal(const Matricks::p3vector<double>& X, Matricks::p3vector<double>& B);


inline int dummy(const int a)
{
    return a;
}

void flsode(const int &neq, const double &t, const double y[], double dy_dt[])
{

  using namespace Matricks;
  using namespace std;
  dummy(neq);
  const double & phi = t;
  const double & r = y[0];
  const double & z = y[1];
  double & dr_dt = dy_dt[0];
  double & dz_dt = dy_dt[1];

  const double xx = r*cos(phi);
  const double yy = r*sin(phi);

  const p3vector<double> fieldPos(xx,yy,z);
  p3vector<double> bField;
  bTotal(fieldPos,bField);
  

  // these are vector calculus components of B
  // i.e. with respect to the orthonormal basis
  // (ie. neither covariant, nor contravariant form)
  const double Br = bField.x()*cos(phi)+bField.y()*sin(phi);
  const double Bz = bField.z();
  const double Bphi = -bField.x()*sin(phi)+bField.y()*cos(phi);

  dr_dt = r*Br/Bphi;
  dz_dt = r*Bz/Bphi;
  

  return;
}




void flsode_details(const int &neq, const double &t, const double y[], double dy_dt[], 
		    double& x0, double& y0, double& z0,
		    double& Bx0, double& By0, double& Bz0,
		    double& Br0, double& Bphi0) {
  using namespace Matricks;
  using namespace std;

  const double & phi = t;
  const double & r = y[0];
  const double & z = y[1];
  double & dr_dt = dy_dt[0];
  double & dz_dt = dy_dt[1];
  dummy(neq);
  disp(t);disp(y[0]);dispcr(y[1]);
  disp(phi);disp(r);dispcr(z);

  const double xx = r*cos(phi);
  const double yy = r*sin(phi);

  const p3vector<double> fieldPos(xx,yy,z);
  p3vector<double> bField;
  bTotal(fieldPos,bField);
  

  // these are vector calculus components of B
  // i.e. with respect to the orthonormal basis
  // (ie. neither covariant, nor contravariant form)
  const double Br = bField.x()*cos(phi)+bField.y()*sin(phi);
  const double Bz = bField.z();
  const double Bphi = -bField.x()*sin(phi)+bField.y()*cos(phi);

  dr_dt = r*Br/Bphi;
  dz_dt = r*Bz/Bphi;
  

  // *** flsode_detail code ***
  x0=xx;
  y0=yy;
  z0=z;
  Bx0=bField.x();
  By0=bField.y();
  Br0 = Br;
  Bz0 = Bz;
  Bphi0 =Bphi;
  // *** flsode_detail code ***

  return;
}
