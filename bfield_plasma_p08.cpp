#include "matricks.hpp"
#include "bfieldfuncs.hpp"

// plasma p08 HSX

void bplasma(const  Matricks::p3vector<double>& X,  Matricks::p3vector<double>& B) {

   using namespace Matricks;
   B= X;  // just to get rid of compiler warning
    B = p3vector<double>(0,0,0);
   
}


