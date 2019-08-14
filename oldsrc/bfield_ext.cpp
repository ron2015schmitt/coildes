#include "bfield_ext.hpp"


void bTotal(const  Matricks::p3vector<double>& X,  Matricks::p3vector<double>& B) {
   using namespace Matricks;

  p3vector<double> Bp;
  bplasma(X,Bp);

  p3vector<double> Bext;
  bext(X,Bext);
 
  B = Bext + Bp;

}
