

#include "coils.hpp"
#include "coilio.hpp"
#include "surface.hpp"
#include "bfieldfuncs.hpp"
#include "bfield_coils.hpp"
#include "coilfft.hpp"




// at first I thought the intel compiler did not like this, so had to change to pointers
// I FOUND THAT THE PROBLEM WAS ELSEWHERE SO I SHOULD BE ABLE TO
// CHANGE THIS BACK TO WHAT IT WAS

// LAvector<p3vector<double> > coilPos("coilPos");
// LAvector<p3vector<double> > coilOrient("coilOrient");
// LAvector<p3vector<double> > coilPlane("coilPlane");
// LAvector<double> coilcurrent("coilcurrent");
// LAvector<double> coilRadius("coilRadius");
// LAvector<p3vector<double> > Acoil("setupcoils::Acoil");
//LAvector<p3vector<double> > jcoil("setupcoils::jcoil");

LAvector<p3vector<double> > *coilPos;
LAvector<p3vector<double> > *coilOrient;
LAvector<p3vector<double> > *coilPlane;
LAvector<double> *coilcurrent;
LAvector<double> *coilRadius;
LAvector<p3vector<double> > *Acoil;
LAvector<p3vector<double> > *jcoil;








// this function assumes that both the constant currents (Itoroidal and Iooloidal) and 
// the kappa current potential are included in jcoil
void bCoils(const p3vector<double>& X, p3vector<double>& B) {
   const unsigned int Npts = (*coilPos).size();

   B=0.0;
   for (unsigned int i = 0; i<Npts; i++) {
     
      const p3vector<double> R = X - (*coilPos)[i];
      //const p3vector<double> Ac = (*Acoil)[i];
      const double Rnorm= norm(R);
      const double r3 = 1/(Rnorm*Rnorm*Rnorm);

      
      const p3vector<double> Btemp =  cross((*jcoil)[i],R) * r3;
      
      B = B+Btemp;
      //    disp(Btemp);  dispcr(B);

   }
   //  dispcr(B);
}                                                                                



void bTotal(const p3vector<double>& X, p3vector<double>& B) {
                                                                                
   p3vector<double> Bc;
   bCoils(X,Bc);
 
   p3vector<double> Bp;
   bplasma(X,Bp);
 
   //  dispcr(Bc);dispcr(Bp);
   B = Bc + Bp;

    
}
 
void bTotalandbCoils(const p3vector<double>& X, p3vector<double>& B, p3vector<double>& Bc) {
                                                                                
   bCoils(X,Bc);

   p3vector<double> Bp;
   bplasma(X,Bp);
 

   B = Bc + Bp;
    
}





int setupcoils(std::string coil_filename, std::string current_filename, 
	       const double Itoroidal, const double Ipoloidal, 
	       unsigned int Nnn, unsigned int Nmm, 
	       unsigned int Nphi, unsigned int Ntheta, 
	       unsigned int Nfund, unsigned int Mfund)
{
      
   // variables for measuring times
//  struct tms tbuff;
//   clock_t ckstart;

   LAvector<double> nn("setupcoils::nn");
   LAvector<double> mm("setupcoils::mm");
   unsigned int NF;
  bool mode00 = true;
  modevectors(NF,nn,mm,Nnn,Nmm,Nfund,Mfund,mode00);

   const unsigned int Npts = Ntheta*Nphi;
   LAvector<double> thetas(Npts,"setupcoils::thetas");
   LAvector<double> phis(Npts,"setupcoils::phis");
   anglevectors(thetas, phis, Ntheta, Nphi);


//   cout << endl;
//   cout<<"$ Generate fourier series matrix ("<<Npts<<" x "<<NF<<")"<<endl;
    
//   STARTTIME(tbuff,ckstart);
    
//   Matrix<complex<double> > fs(Npts,NF,"setupcoils::fs");
//   fseries(nn,mm,thetas,phis,fs);
//   STOPTIME(tbuff,ckstart);

   // load the coil surface fourier coef's
                                                                                                  
   cout << "$ Loading COIL SURFACE fourier coefficients from " << coil_filename << endl;
   FourierSurface coilfourier;
   if (load_fourier_surface(coil_filename,coilfourier))
      return 2;
                                                                                                  
   coilfourier.RF().name("c.RF");
   coilfourier.ZF().name("c.ZF");
                                                                                                  
   // print coef's
   //   printfouriercoefs(coilfourier.nn(),coilfourier.mm(),coilfourier.RF(),coilfourier.ZF(),10,18);


   // load the coil CURRENT fourier coef's
   // at some point, add code so that user can select type of coef's from command line
 
   cout <<endl<< "$ Loading COIL CURRENT fourier coefficients from " << current_filename << endl;
   LAvector<complex<double> > IF(NF,"IF");
   if (load_coefs(current_filename,CoefFileFormat_sincos,nn,mm,IF))
      return 3;
//   LAvector<double> IFreal(NF,"IFreal");
//   LAvector<double> IFimag(NF,"IFimag");
 //  IFreal = real(IF);
//   IFimag = imag(IF);
 
   // print coef's
   //   printfouriercoefs(nn,mm,IFreal,IFimag,10,18);
 




//    coilPos.resize(Npts);
//    coilOrient.resize(Npts);
//    coilRadius.resize(Npts);
//    coilPlane.resize(Npts);
//    coilcurrent.resize(Npts);
//    Acoil.resize(Npts);
//    jcoil.resize(Npts);

   coilPos = new LAvector<p3vector<double> >(Npts,"coilPos");
   coilOrient = new LAvector<p3vector<double> >(Npts,"coilOrient") ;
   coilPlane = new LAvector<p3vector<double> >(Npts,"coilPlane");
   coilcurrent = new LAvector<double>(Npts,"coilcurrent");
   coilRadius = new LAvector<double>(Npts,"coilRadius"); 
   Acoil= new LAvector<p3vector<double> >(Npts,"Acoil"); 
   jcoil = new LAvector<p3vector<double> >(Npts,"jcoil"); 


 //  LAvector<double>  cc(Npts,"coilcurrent");


printcr("ifft of current potential");
//   expandfunction((*coilcurrent),IF,fs);
    ifft2d(*coilcurrent,IF,Nphi,Ntheta,Nnn,Nmm,Nfund,Mfund,1e-10,1/(2*PI),mode00);


 
   LAvector<p3vector<double> > dx_dr(Npts, "setupcoils::dx_dr");
   LAvector<p3vector<double> > dx_dtheta(Npts,"setupcoils::dx_dtheta");
   LAvector<p3vector<double> > dx_dphi(Npts,"setupcoils::dx_dphi");
   LAvector<p3vector<double> > grad_r(Npts,"setupcoils::grad_r");
   LAvector<p3vector<double> > grad_theta(Npts,"setupcoils::grad_theta");
   LAvector<p3vector<double> > grad_phi(Npts,"setupcoils::grad_phi");

   cout << endl;
   cout <<"$ Mapping coil surface fourier coefficients to "<<Ntheta<<" x "<<Nphi<<" (theta by phi) grid"<<endl;


   expandsurfaceandbases(*coilPos,*Acoil,dx_dr,dx_dtheta,dx_dphi,grad_r,grad_theta,grad_phi,coilfourier,thetas,phis);


   //     dx_dr.textformat(text_nobraces);
   //     dx_dr.perline(1);
   //     dx_dtheta.textformat(text_nobraces);
   //     dx_dtheta.perline(1);
   //     dx_dphi.textformat(text_nobraces);
   //     dx_dphi.perline(1);

   //     save(dx_dr,"dx_dr.out");
   //     save(dx_dtheta,"dx_dtheta.out");
   //     save(dx_dphi,"dx_dphi.out");



  printcr("calc geometrics for current potential");


   for (unsigned int i = 0; i<Npts; i++) {
      const double dtdp = (2*PI)/double(Ntheta) * (2*PI)/double(Nphi);
      (*Acoil)[i] = (*Acoil)[i] * dtdp;
      double Anorm = norm((*Acoil)[i]);
      (*coilOrient)[i] = (*Acoil)[i] / Anorm;
      (*coilRadius)[i] = sqrt(Anorm / PI);
      (*coilPlane)[i] = dx_dtheta[i] / norm(dx_dtheta[i]);
   }


   // calculate current density fourier coef's
   // ** these are actually j_theta*Jacobian and j_phi*Jacobian **
    
  printcr("calc current denisty componets in Fourier space");

   LAvector<complex<double> > j_thetaF(NF,"j_thetaF");
   LAvector<complex<double> > j_phiF(NF,"j_phiF");
   j_thetaF = complex<double>(0,1)*nn*IF;
   j_phiF = complex<double>(0,-1)*mm*IF;

   LAvector<double> j_theta(Npts,"j_theta");
   LAvector<double> j_phi(Npts,"j_phi");

   
//   expandfunction(j_theta,j_thetaF,fs);
//   expandfunction(j_phi,j_phiF,fs);

  printcr("ifft of  current denisty components");

    ifft2d(j_theta,j_thetaF,Nphi,Ntheta,Nnn,Nmm,Nfund,Mfund,1e-10,1/(2*PI),mode00);
    ifft2d(j_phi,j_phiF,Nphi,Ntheta,Nnn,Nmm,Nfund,Mfund,1e-10,1/(2*PI),mode00);


  printcr("put it all together");

   const double dphi_by_dtheta = 2.0*PI* 2.0*PI/double(Npts);

   // jcoil is actually j*Jacobian*mu0div4pi * dphi_by_dtheta
   for (unsigned int i = 0; i<Npts; i++) {
      //      double Jinv = dot(grad_r[i],cross(grad_theta[i],grad_phi[i]));
      j_theta[i] =  (Ipoloidal/(2*PI) + j_theta[i]);
      j_phi[i] =  (Itoroidal/(2*PI) + j_phi[i]);
      (*jcoil)[i] = j_theta[i] * dx_dtheta[i]   +    j_phi[i] * dx_dphi[i];
      (*jcoil)[i] = mu0div4pi * dphi_by_dtheta * (*jcoil)[i];
   }

   return 0;
}






