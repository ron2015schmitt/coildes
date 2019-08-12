/*************************************************************************
*
*   File Name    :  matlabio.hpp
*   Platform     :  Red Hat LINUX
*   Author       :  Ron Schmitt
*   Date         :  February 3, 2003
*
*
*   SYNOPSIS
*     This file contains functions for loading (and saving) matrices 
*     and vectors from (and to) Matlab .mat files.
*
*   Change History:
*
**************************************************************************/


#ifndef MATLABIO_H
#define MATLABIO_H

// You must create typedefs for the following types:
// Dmatrix
// Dvector
// You can use a specific math library or just use STL classes

#include "matricks.hpp"
using namespace Matricks;

typedef Vector<double> Dvector;

typedef Matrix<double> Dmatrix;
/////////////////////////////////////


int LoadMatrixFromMATFile(char *file, char *varname, Dmatrix** Ap);

int LoadVectorFromMATFile(char *file, char *varname, Dvector** vp);

int SaveMatrixToMATFile(char *file, char *varname, Dmatrix& Ap);
int SaveVectorToMATFile(char *file, char *varname, Dvector& vp);

int PrintMATFileDirectory(char *file);



#endif
