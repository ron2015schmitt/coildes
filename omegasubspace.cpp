
#include "omegasubspace.hpp"


void omegasubspace(LAvector<p3vector<double> >& grad_r, Matrix<complex<double> >& fs, 
	      LAvector<double>& nn, LAvector<double>& mm,
	      Matrix<complex<double> >& fsR2,
	      LAvector<double>& nnR2, LAvector<double>& mmR2){

  const unsigned int NF=nn.size();
  const unsigned int Npts=grad_r.size();

  LAvector<double> grad_r_mag(Npts,"grad_r_mag");

  for (unsigned int i=0; i<Npts; i++)
    grad_r_mag[i] = norm(grad_r[i]);


  LAvector<complex<double> > grmfs(NF,"grmfs");
  grmfs = sqr((2*PI))/Npts*(grad_r_mag|fs);


  const double grmgrm  = sqr((2*PI))/Npts*(grad_r_mag|grad_r_mag);

  Matrix<complex<double> > fsTEMP(Npts,NF,"fsTEMP");
  for (unsigned int i=0; i<Npts; i++)
    for (unsigned int k=0; k<NF; k++)
      fsTEMP(i,k)=fs(i,k) - (grmfs[k]/grmgrm) * grad_r_mag[i];  

  LAvector<complex<double> > grmfsTEMP(NF,"grmfsTEMP");
  grmfsTEMP = sqr((2*PI))/Npts*(grad_r_mag|fsTEMP);


  // re-orthogonalize foruier series
  // remove any eigenvectors that are close to zero in length
  dispcr(NF);
  unsigned int NFnew=0;
  LAvector<unsigned int> modeindex(NF,"modeindex");
  for (unsigned int k=0; k<NF; k++) {
    LAvector<complex<double> > vk(Npts,"vk");
    vk=fsTEMP.col(k);  
    double magnitude = (2*PI)/sqrt(double(Npts))*(norm(vk));
    //    disp(k);dispcr(magnitude);
    if (magnitude > MODEMIN) {
      modeindex[NFnew] = k;
      NFnew++;
    }else {
      cout << "Mode (m="<<mm[k]<<",n="<<nn[k]<<") has been removed"<<endl;
    }
  }
  dispcr(NFnew);

  nnR2.resize(NFnew);
  mmR2.resize(NFnew);

  fsR2.resize(Npts,NFnew);
  for (unsigned int k=0; k<NFnew; k++) {
    mmR2[k] = mm[modeindex[k]];
    nnR2[k] = nn[modeindex[k]];
    LAvector<complex<double> > vold(Npts,"vold");
    vold=fsTEMP.col(modeindex[k]);  
    LAvector<complex<double> > v(Npts,"v");
    v=vold;  
    for (unsigned int j=0; j<k; j++) {
      LAvector<complex<double> > vj(Npts,"vj");
      vj=fsR2.col(j);  
      complex<double> ratio = (conj(vj)|vold)/norm(vj);
      v = v - ratio * vj;
    }
    fsR2.col(k) = v;
  }
}
