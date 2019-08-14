

#ifndef CREATESURFACE_H
#define CREATESURFACE_H

void createsurface(const Vector<p3vector<double> >& XX, const Vector<p3vector<double> >& AA,
		   Vector<p3vector<double> >& XXnew,  Vector<p3vector<double> >& AAnew,
		   const unsigned int Ntheta, const unsigned int Nphi, 
		   const double dtheta, const double dphi,
		   const double scale);



#endif
