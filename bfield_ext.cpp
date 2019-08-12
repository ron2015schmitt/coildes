#include "bfield_ext.hpp"


void bTotal(const  COOLL::p3vector<double>& X,  COOLL::p3vector<double>& B) {
   using namespace COOLL;

  p3vector<double> Bp;
  bplasma(X,Bp);

  p3vector<double> Bext;
  bext(X,Bext);
 
  B = Bext + Bp;

}
