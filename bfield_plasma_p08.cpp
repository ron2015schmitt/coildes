#include "cooll.hpp"
#include "bfieldfuncs.hpp"

// plasma p08 HSX

void bplasma(const  COOLL::p3vector<double>& X,  COOLL::p3vector<double>& B) {

   using namespace COOLL;
   B= X;  // just to get rid of compiler warning
    B = p3vector<double>(0,0,0);
   
}


