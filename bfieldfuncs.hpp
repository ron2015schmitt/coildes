
#include "cooll.hpp"

void bfieldwire(COOLL::p3vector<double>& B, const COOLL::p3vector<double>& X, const double I);

void bfieldwire_general(COOLL::p3vector<double>& B, const COOLL::p3vector<double>& X,const double x0, const double y0, const double I) ;

void bfieldcoil(const COOLL::p3vector<double>& coilPosition,
		const COOLL::p3vector<double>& coilOrient,
		const COOLL::p3vector<double>& coilPlaneOrient,
		const double& coilRadius,
		const COOLL::p3vector<double>& fieldPosition,
		COOLL::p3vector<double>& bField, const double I);
