
#include "matricks.hpp"

void bfieldwire(Matricks::p3vector<double>& B, const Matricks::p3vector<double>& X, const double I);

void bfieldwire_general(Matricks::p3vector<double>& B, const Matricks::p3vector<double>& X,const double x0, const double y0, const double I) ;

void bfieldcoil(const Matricks::p3vector<double>& coilPosition,
		const Matricks::p3vector<double>& coilOrient,
		const Matricks::p3vector<double>& coilPlaneOrient,
		const double& coilRadius,
		const Matricks::p3vector<double>& fieldPosition,
		Matricks::p3vector<double>& bField, const double I);
