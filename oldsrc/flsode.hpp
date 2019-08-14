typedef int FINT;

 
void flsode(const FINT &neq, const double &t, const double y[], double dy_dt[]);

void flsode_details(const int &neq, const double &t, const double y[], double dy_dt[], 
		    double& x0, double& y0, double& z0,
		    double& Bx0, double& By0, double& Bz0,
		    double& Br0, double& Bphi0);





