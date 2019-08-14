#include "matricks.hpp"


// bTotal is defined in bfield_ext.cpp
void bTotal(const  Matricks::p3vector<double>& X,  Matricks::p3vector<double>& B);

// bext is defined in bfield_ext_%.cpp
void bext(const Matricks::p3vector<double>& X, Matricks::p3vector<double>& Bext);

// bplasma is defined in bfield_plasma_p%.cpp
void bplasma(const  Matricks::p3vector<double>& X,  Matricks::p3vector<double>& B);

