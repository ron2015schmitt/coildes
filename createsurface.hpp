

#ifndef CREATESURFACE_H
#define CREATESURFACE_H

void createsurface(const LAvector<p3vector<double> >& XX, const LAvector<p3vector<double> >& AA,
		   LAvector<p3vector<double> >& XXnew,  LAvector<p3vector<double> >& AAnew,
		   const unsigned int Ntheta, const unsigned int Nphi, 
		   const double dtheta, const double dphi,
		   const double scale);



#endif
