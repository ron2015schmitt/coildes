
#include <string>
#include <iostream>
using namespace std;


#include "matlabio.hpp"
#include "matdumpcmdline.hpp"


// Main Function for code

int main (int argc, char *argv[])
{
  vn=0; mn=0;
  printdir =false;

  // parse command line input
 
  parse_cmd(argc, argv);


  char* fn=filenames[0];
  if (printdir)
    PrintMATFileDirectory(fn);


  Dvector* vp=NULL;
  Dvector** vpp = &vp;
  for (int i = vn; i>0; i--){
    LoadVectorFromMATFile(fn,vnames[i-1],vpp);
    Dvector vtemp = *vp;
    cout<< vnames[i-1] <<"=" << vtemp <<endl;
  }

  Dmatrix* Ap=NULL;
  Dmatrix** App = &Ap;
  for (int i = mn; i>0; i--){
    LoadMatrixFromMATFile(fn,mnames[i-1],App);
    Dmatrix mtemp = *Ap;
    cout<< mnames[i-1] <<"=" << mtemp <<endl;
  }
}
