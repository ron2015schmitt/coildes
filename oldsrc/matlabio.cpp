/*************************************************************************
 *
 *   File Name    :  matlabio.cpp
 *   Platform     :  Red Hat LINUX, g++ compiler
 *   Author       :  Ron Schmitt
 *   Date         :  February 3, 2003
 *
 *
 *   SYNOPSIS
 *     This file contains functions for loading (and saving) MTL matrices 
 *     and vectors from (and to) Matlab .mat files.
 *
 *   Change History:
 *
 **************************************************************************/


#include <iostream>
#include <complex>

using namespace std;

//MATLAB .MAT file utility header
// Provided by mathworks
// you must link as 
#include "mat.h"

#include "matlabio.hpp"


/////////////////////////////////////////////////////////////////////////
//                          Type Definitions
/////////////////////////////////////////////////////////////////////////







/*************************************************************************
 *
 *   FUNCTION     : LoadMatrixFromMATFile 
 *   INPUTS       :  
 *     file - string containing filename
 *     varname - string containing variable name
 *   OUTPUTS      :  
 *   I/O          :
      A - Matrix where data is to be copied to.
      *
      *   SYNOPSIS
      *
      *   This function loads the variable 'varname' from the matlab .mat file
      *   name 'file' into the Matrix 'A'.
      *
      *   Matrix must be double precision real and dense.
      *
      *  NOTES:  could expand this function by letting type of matrix be any
      *          matrix and then load in appropriate type at runtime.
      *
      **************************************************************************/

int LoadMatrixFromMATFile(char *file, char *varname,  Dmatrix** Ap)
{

  MATFile *pmat;
  mxArray *pa;


  // open the matlab data file
  
  pmat = matOpen(file, "r");
  if (pmat == NULL) {
    printf("Error opening file %s\n", file);
    return(1);
  }


  // get pointer for specified variable

  pa = matGetVariable(pmat, varname);
  if (pa == NULL) {
    printf("Error reading variable '%s' in file '%s'\n",varname, file);
    return(1);
  }


  // extract some of the basic info about this variable
  // and check that the data to be sure it is of appropriate type

  int classID = mxGetClassID(pa);
  int dim = mxGetNumberOfDimensions(pa);
  bool iscomplex=mxIsComplex(pa);


 
  if (classID != mxDOUBLE_CLASS) {
    printf("Error in file %s: Data of variable '%s' is not of type double precision (%d).\n", file,varname,classID);
    return(1);
  }

  if (dim != 2) {
    printf("Error in file %s: Data of variable '%s' is not a matrix (%d).\n",file,varname,dim);
    return(1);
  }

  if (iscomplex) {
    printf("Error in file %s: Data of variable '%s' is complex when it should be real.\n", file,varname);
    return(1);
  }
 
  // get pointer to actual data
  int Nrows= mxGetM(pa);
  int Ncols= mxGetN(pa);
  double *data=(double*)mxGetData(pa);

  // constructor for matrix
  *Ap = new Dmatrix(Nrows,Ncols);
  

  // copy data. note matlab stores data in column major order
  int ii=0;
  for(int k = 0; k<Ncols ; k++)
    for(int j = 0; j<Nrows ; j++) 
      (**Ap)(j,k)=data[ii++];
 
  // deallocate the memory for the array

  mxDestroyArray(pa);
  
  // close the file

  if (matClose(pmat) != 0) {
    printf("Error closing file %s\n",file);
    return(1);
  }

  // print some stuff

  printf("successfully loaded variable '%s': ", varname);
  printf("%d x %d matrix\n",Nrows,Ncols);

  return(0);
}







/*************************************************************************
 *
 *   FUNCTION     : LoadVectorFromFile
 *   INPUTS       :  
 *     file - string containing filename
 *     varname - string containing variable name
 *   OUTPUTS      :  
 *   I/O          :
 *      A - vector where data is to be copied to.
 *
 *
 *   SYNOPSIS
 *
 *   This function loads the variable 'varname' from the matlab .mat file
 *   name 'file' into the vector 'A'.
 *
 *   Vector must be of doubles and real.
 *
 **************************************************************************/


int LoadVectorFromMATFile(char *file, char *varname, Dvector** vp)
{

  MATFile *pmat;
  mxArray *pa;
  int N;

  // open the matlab data file
  
  pmat = matOpen(file, "r");
  if (pmat == NULL) {
    printf("Error opening file %s\n", file);
    return(1);
  }


  // get pointer for specified variable

  pa = matGetVariable(pmat, varname);
  if (pa == NULL) {
    printf("Error reading variable '%s' in file '%s'\n",varname, file);
    return(1);
  }


  // extract some of the basic info about this variable
  // and check that the data to be sure it is of appropriate type

  int classID = mxGetClassID(pa);
  int dim = mxGetNumberOfDimensions(pa);
  bool iscomplex=mxIsComplex(pa);


 
  if (classID != mxDOUBLE_CLASS) {
    printf("Error in file %s: Data of variable '%s' is not of type double precision (%d).\n", file,varname,classID);
    return(1);
  }

  if (dim != 2) {
    printf("Error in file %s: Variable '%s' is not a matrix (%d).\n",file,varname,dim);
    return(1);
  }

  if (iscomplex) {
    printf("Error in file %s: Data of variable '%s' is complex when it should be real.\n", file,varname);
    return(1);
  }
 
  // get pointer to actual data
  int Nrows= mxGetM(pa);
  int Ncols= mxGetN(pa);
  double *data=(double*)mxGetData(pa);

  // Note that a matlab vector may be of row or column type

  if ((Nrows > 1) && (Ncols > 1)) {
    printf("Error in file %s: Variable '%s' is not a vector (%d x %d).\n",file,varname,Nrows,Ncols);
    return(1);
  } else {
    N=Nrows*Ncols;
  }

  // constructor for vector
  *vp=new Dvector(N);
 
  // copy data 
  for( int kk=0;kk<N;kk++) 
    (**vp)[kk]=data[kk];
 

  // deallocate the memory for the array

  mxDestroyArray(pa);
  
  // close the file

  if (matClose(pmat) != 0) {
    printf("Error closing file %s\n",file);
    return(1);
  }

  // print some stuff

  printf("successfully loaded variable '%s': ", varname);
  printf("%d element Vector\n",N);


  return(0);
}







int SaveMatrixToMATFile(char *file, char *varname, Dmatrix& A)
{

  MATFile *pmat;
  mxArray *pa;
  int N = A.Nrows();
  int M = A.Ncols();
  int NM = N*M;
  

  pmat = matOpen(file, "u");  // use "u" to update an existing file
  if (pmat == NULL) {
    pmat = matOpen(file, "w");  // use "w" to create a new file
    if (pmat == NULL) {
      printf("Error creating file %s\n", file);
      printf("(Do you have write permission in this directory?)\n");
      return(EXIT_FAILURE);
    }
  }

  pa = mxCreateDoubleMatrix(N,M,mxREAL);
  if (pa == NULL) {
      printf("%s : Out of memory on line %d\n", __FILE__, __LINE__);
      printf("Unable to create mxArray.\n");
      return(EXIT_FAILURE);
  }

  double* data;
  data = new double[NM];
  int ii=0;
  for(int k = 0; k<M ; k++)
    for(int j = 0; j<N ; j++) 
      data[ii++]=A(j,k);  // matlab stores matrices in column order

  //  for( ii = 0; ii<NM ; ii++)
  //cout << data[ii] <<endl;

  memcpy((void *)(mxGetPr(pa)), (void *)data, sizeof(*data)*NM);
  

  int status = matPutVariableAsGlobal(pmat,varname , pa);
  if (status != 0) {
      printf("Error using matPutVariableAsGlobal\n");
      return(EXIT_FAILURE);
  } 
  
  
  /* clean up */
  mxDestroyArray(pa);


  if (matClose(pmat) != 0) {
    printf("Error closing file %s\n",file);
    return(EXIT_FAILURE);
  }



  printf("successfully saved variable '%s': ", varname);
  printf("%d x %d matrix\n",N,M);



  return(0);
}












int SaveVectorToMATFile(char *file, char *varname, Dvector& vp)
{

  MATFile *pmat;
  mxArray *pa;
  int N = vp.size();

  

  pmat = matOpen(file, "u");  // use "u" to update an existing file
  if (pmat == NULL) {
    pmat = matOpen(file, "w");  // use "w" to create a new file
    if (pmat == NULL) {
      printf("Error creating file %s\n", file);
      printf("(Do you have write permission in this directory?)\n");
      return(EXIT_FAILURE);
    }
  }

  int dims[2];
  dims[0]=N;
  dims[1]=1;

  pa = mxCreateNumericArray(2,dims, mxDOUBLE_CLASS, mxREAL);
  if (pa == NULL) {
      printf("%s : Out of memory on line %d\n", __FILE__, __LINE__);
      printf("Unable to create mxArray.\n");
      return(EXIT_FAILURE);
  }

  double* data;
  data = new double[N];
  for(int i = 0; i<N ; i++)
    data[i]=vp[i];
  
  memcpy((void *)(mxGetPr(pa)), (void *)data, sizeof(*data)*N);
  

  int status = matPutVariableAsGlobal(pmat,varname , pa);
  if (status != 0) {
      printf("Error using matPutVariableAsGlobal\n");
      return(EXIT_FAILURE);
  } 
  
  
  /* clean up */
  mxDestroyArray(pa);


  if (matClose(pmat) != 0) {
    printf("Error closing file %s\n",file);
    return(EXIT_FAILURE);
  }


  printf("successfully saved variable '%s': ", varname);
  printf("%d element Vector\n",N);




  return(0);
}










char classtypes[][25]={ "mxUNKNOWN_CLASS",
        "mxCELL_CLASS",
        "mxSTRUCT_CLASS",
        "mxLOGICAL_CLASS",
        "mxCHAR_CLASS",
        "mxSPARSE_CLASS",
        "mxDOUBLE_CLASS",
        "mxSINGLE_CLASS",
        "mxINT8_CLASS",
        "mxUINT8_CLASS",
        "mxINT16_CLASS",
        "mxUINT16_CLASS",
        "mxINT32_CLASS",
        "mxUINT32_CLASS",
        "mxINT64_CLASS",
        "mxUINT64_CLASS",
        "mxFUNCTION_CLASS",
        "mxOPAQUE_CLASS",
        "mxOBJECT_CLASS"};





int PrintMATFileDirectory(char *file)
{

  MATFile *pmat;

  // open the matlab data file
  
  pmat = matOpen(file, "r");
  if (pmat == NULL) {
    printf("Error opening file %s\n", file);
    return(1);
  }

  char    **dir;
  int     ndir;



  /*
   * get directory of MAT-file
   */
  dir = matGetDir(pmat, &ndir);
  if (dir == NULL) {
    printf("Error reading directory of file %s\n", file);
    return(1);
  } else {
    printf("Directory of %s:\n", file);
    mxArray *pa;

    for (int i=0; i < ndir; i++) {
      // get pointer for specified variable

      pa = matGetVariable(pmat, dir[i]);
      if (pa == NULL) {
	printf("Error reading variable '%s' in file '%s'\n",dir[i], file);
	return(1);
      }
      // extract some of the basic info about this variable
      // and check that the data to be sure it is of appropriate type

      int classID = mxGetClassID(pa);
      int dim = mxGetNumberOfDimensions(pa);
      bool iscomplex=mxIsComplex(pa);
      char * scplx;
      if (iscomplex)
	scplx = "(Complex)";
      else
	scplx = "";
      printf("%s: %s%s, %d dimensions ",dir[i],classtypes[classID],scplx,dim);
      const int*  Np = mxGetDimensions(pa);
      for (int k = 0; k<dim; k++) {
	if (k==0)
	  cout << "[";
	cout << Np[k];
	if (k != (dim-1))
	  cout << " x ";
	else
	  cout << "]";

      }
      cout << endl;
    }
  }
  mxFree(dir);








  // close the file

  if (matClose(pmat) != 0) {
    printf("Error closing file %s\n",file);
    return(1);
  }

  return(0);

}
