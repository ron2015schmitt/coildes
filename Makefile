

SHELL=/bin/bash


DIR_MATRICKS = ./matricks
INC_MATRICKS = -I $(DIR_MATRICKS) 
LIB_MATRICKS = -L$(DIR_MATRICKS) -lmatricks
LIBFILE_MATRICKS = $(DIR_MATRICKS)/libmatricks.a


INCLUDES = $(INC_MATRICKS)
LIBS =  $(LIB_MATRICKS) 




# C++ compiler
CPPOPT = 
CPPC = g++
CPPFLAGS = $(CPPOPT) $(INCLUDES)

# linker
LNKOPT =
LDFLAGS = $(LNKOPT)  
LNK = g++

# FORTRAN COMPILER
FC = gfortran



#notes:  -v option on g++/gcc is useful for 
# use this variable to set gcc/g++ compile options from command line
# just write 'CLOPTS=yada yada yada' at end of command line
# example for generating gprof output:
#   make coils 'COPTS=-pg' 'LOPTS=-pg'

#LNKOPT=-no-ipo  #-cxxlib-gcc #-i-static -static #-nodefaultlibs

#INC_FFTW = -I $(FFTWDIR)/include 

#LIB_LAPACK= -llapack 
#LIB_FFTW = -L$(FFTWDIR)/lib -lfftw3 


#LIB_LSODE= -L$(LSODEDIR) -llsode
#LSODEDIR = /home/rs2015/lsode
#FFTWDIR = /home/rs2015/fftw-3.0.1


# do we need -lircmt?
#LIB_FORTRAN = -lfio -lblas -lf77math -lf90math -L/usr/local/intel/fc/9.0/lib -lifcore -lifcoremt -lcxa -lcprts -lunwind -lirc  -lsvml -lompstub -limf  -lcxaguard -lifport -lguide -lguide_stats 
#LIB_FORTRAN = -L/usr/local/intel/fc/9.0/lib -lifcore


#LIBS_C =  -L/usr/local/intel/cc/9.0/lib $(LIB_COOLL) 
#LIBS_FIELD = $(LIB_LSODE)  $(LIB_NAG) $(LIB_FORTRAN)



#EXEC_MAIN = Bnormal2current flux2current current2Bnormal current2flux scoild
#EXEC_MAIN_FFT = Bnormal2current_fft flux2current_fft current2Bnormal_fft current2flux_fft scoild_fft scoild_fft_II scoild_fft_IIb scoild_fft_III scoild_fft_IV scoild_fft_V scoild_fft_VI scoild_fft_VII scoild_fft_VIII scoild_fft_IX scoild_fft_X scoild_fft_XI scoild_fft_Xa scoild_fft_XIa  scoild_fft_XII scoild_fft_XIIa scoild_fft_XIII scoild_fft_XIIIa  scoild_fft_XIV  scoild_fft_XV
#EXEC_BTRAJ =  calcBtraj calcBtrajfromCoils
#EXEC_BONSURF = calcBonsurf calcBonsurffromCoils 
#EXEC_AUXILIARY =  mkcoilsurf  expandsurf surfF findpertsurface deltasurfs deltaflux percentdelta reduceflux rankorder findpertsurface_fft expandfunc calcIpol onlyharms expandI cullF cutcoils iotafromIF projectBw divfuncs calcGammaOmega calcOmegaN calcOmegaN_fft cart2contravariantF cart2contravariantJF b2lambda lambda2b calcOmegas calcMs calcJ testfft
#EXEC_OTHER = garabedian2sincos garabedian2sincosII  garabedian_current2sincos
#EXEC_TEMP = temp

#EXEC = $(EXEC_MAIN) $(EXEC_MAIN_FFT)  $(EXEC_AUXILIARY) $(EXEC_OTHER) 

EXEC = 

TESTS = 

CPROGS = 


#prevent any default rules from being used
.SUFFIXES:

# don't delete .o files
.PRECIOUS: %.o calcBtraj_% calcBtrajfromCoils_% calcBonsurf_% calcBonsurffromCoils_% 

all: $(EXEC)


tests: $(TESTS)

clean:
	@command rm -f *.o
	@command rm -f $(EXEC)
	@command rm -f $(TESTS)
	@command rm -f core.*
	@command rm calcBtraj_* calcBtrajfromCoils_* calcBonsurf_* calcBonsurffromCoils_* 
def:
	echo "nothing done"


$(LIBFILE_MATRICKS):
	cd matricks && ./configure


%.o: %.cpp coils.hpp $(LIBFILE_MATRICKS)
ifdef CAREFUL
	$(CPPC) $(CPPFLAGS) -D"MATRICKS_CAREFUL=1" -c $*.cpp -o $@
else
	$(CPPC) $(CPPFLAGS) -c $*.cpp -o $@
endif
 

%: %.o 
	$(CPPC) $(LDFLAGS) $*.o -o $@ $(LIBS) 

calc_%: calcBtraj_% calcBtrajfromCoils_% calcBonsurf_% calcBonsurffromCoils_% 
	@echo

current2Bnormal: current2Bnormal.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)

current2flux: current2flux.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)

mkcoilsurf: mkcoilsurf.o coils.o coilio.o coils_cmdline.o surface.o createsurface.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)

expandsurf: expandsurf.o coils.o coilio.o coils_cmdline.o surface.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)

expandfunc: expandfunc.o coils.o coilio.o coils_cmdline.o surface.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)

surfF: surfF.o coils.o coilio.o coils_cmdline.o surface.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)

deltasurfs: deltasurfs.o coils.o coilio.o coils_cmdline.o surface.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)

percentdelta: percentdelta.o coils.o coilio.o coils_cmdline.o surface.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)

calcBonsurf_%: bfield_plasma_%.o bfield_ext_%.o bfield_ext.o calcBonsurf.o coils.o coilio.o coils_cmdline.o surface.o bfieldfuncs.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)

calcBonsurffromCoils_%: bfield_plasma_%.o bfield_coils.o calcBonsurffromCoils.o coils.o coilio.o coils_cmdline.o surface.o bfieldfuncs.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)

deltaflux: deltaflux.o coils.o coilio.o coils_cmdline.o surface.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK) 


garabedian2sincos: garabedian2sincos.o coils.o coilio.o coils_cmdline.o surface.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)

garabedian2sincosII: garabedian2sincosII.o coils.o coilio.o coils_cmdline.o surface.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)

garabedian_current2sincos: garabedian_current2sincos.o coils.o coilio.o coils_cmdline.o surface.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)

rankorder: rankorder.o coils.o coilio.o coils_cmdline.o surface.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)

calcIpol: calcIpol.o coils.o coilio.o coils_cmdline.o surface.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)

onlyharms: onlyharms.o coils.o coilio.o coils_cmdline.o surface.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)

expandI: expandI.o coils.o coilio.o coils_cmdline.o surface.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)

cullF: cullF.o coils.o coilio.o coils_cmdline.o surface.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)

cutcoils: cutcoils.o coils.o coilio.o coils_cmdline.o surface.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)

calcGammaOmega: calcGammaOmega.o coils.o coilio.o coils_cmdline.o surface.o gammamatrix.o omegamatrix.o gradfvector.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)

calcExactGammaOmega_%: bfield_plasma_%.o bfield_coils.o calcExactGammaOmega.o coils.o coilio.o coils_cmdline.o surface.o bfieldfuncs.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)

cart2contravariantF: cart2contravariantF.o coils.o coilio.o coils_cmdline.o surface.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)

cart2contravariantJF: cart2contravariantJF.o coils.o coilio.o coils_cmdline.o surface.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)

divfuncs: divfuncs.o coils.o coilio.o coils_cmdline.o surface.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)


lambda2b: lambda2b.o coils.o coilio.o coils_cmdline.o surface.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)



# These codes require LAPACK
Bnormal2current: Bnormal2current.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C) $(LIB_LAPACK) 

flux2current: flux2current.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C) $(LIB_LAPACK) 

scoild: scoild.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o gammamatrix.o omegamatrix.o gradfvector.o omegasubspace.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK) 

findpertsurface: findpertsurface.o coils.o coilio.o coils_cmdline.o surface.o omegamatrix.o gradfvector.o omegasubspace.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK) 

reduceflux: reduceflux.o coils.o coilio.o coils_cmdline.o surface.o gammamatrix.o omegamatrix.o gradfvector.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK) 

flux2current0: flux2current0.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C) $(LIB_LAPACK) 


calcOmegaN: calcOmegaN.o coils.o coilio.o coils_cmdline.o surface.o omegamatrix.o gradfvector.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C) $(LIB_LAPACK) 

projectBw: projectBw.o coils.o coilio.o coils_cmdline.o surface.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C) $(LIB_LAPACK) 


# These codes require FFTW
current2Bnormal_fft: current2Bnormal_fft.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o coilfft.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C) $(LIB_FFTW)

current2flux_fft: current2flux_fft.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o  coilfft.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C) $(LIB_FFTW) 

iotafromIF: iotafromIF.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o  coilfft.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C) $(LIB_FFTW) 

b2lambda: b2lambda.o coils.o coilio.o coils_cmdline.o surface.o  coilfft.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C) $(LIB_FFTW) 

calcOmegas: calcOmegas.o coils.o coilio.o coils_cmdline.o surface.o  coilfft.o omegamatrix.o omegamatrix_fft.o gradfvector.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C) $(LIB_FFTW) 

calcMs: calcMs.o coils.o coilio.o coils_cmdline.o surface.o  coilfft.o inductancematrix.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C) $(LIB_FFTW) 

calcJ: calcJ.o coils.o coilio.o coils_cmdline.o surface.o  coilfft.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C) $(LIB_FFTW) 

testfft: testfft.o coils.o coilio.o coils_cmdline.o surface.o  coilfft.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C) $(LIB_FFTW) 


# These codes require FFTW & LAPACK

calcOmegaN_fft: calcOmegaN_fft.o coils.o coilio.o coils_cmdline.o surface.o omegamatrix_fft.o gradfvector.o coilfft.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C) $(LIB_LAPACK)   $(LIB_FFTW)

Bnormal2current_fft: Bnormal2current_fft.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o  coilfft.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C) $(LIB_LAPACK)  $(LIB_FFTW)

flux2current_fft: flux2current_fft.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o  coilfft.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C) $(LIB_LAPACK)  $(LIB_FFTW)
 
scoild_fft: scoild_fft.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o gammamatrix.o omegamatrix.o gradfvector.o omegasubspace.o  coilfft.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK)  $(LIB_FFTW)

scoild_fft_II: scoild_fft_II.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o gammamatrix.o omegamatrix.o gradfvector.o omegasubspace.o  coilfft.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK)  $(LIB_FFTW)

findpertsurface_fft: findpertsurface_fft.o coils.o coilio.o coils_cmdline.o surface.o omegamatrix_fft.o gradfvector.o omegasubspace.o coilfft.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK) $(LIB_FFTW)

scoild_fft_III: scoild_fft_III.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o gammamatrix.o omegamatrix.o gradfvector.o omegasubspace.o  rhomatrix.o coilfft.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK)  $(LIB_FFTW)

scoild_fft_IV: scoild_fft_IV.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o gammamatrix.o omegamatrix.o gradfvector.o omegasubspace.o  rhomatrix.o coilfft.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK)  $(LIB_FFTW)

scoild_fft_IIb: scoild_fft_IIb.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o gammamatrix.o omegamatrix.o gradfvector.o omegasubspace.o  coilfft.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK)  $(LIB_FFTW)

scoild_fft_V: scoild_fft_V.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o gammamatrix.o omegamatrix.o gradfvector.o omegasubspace.o  rhomatrix.o coilfft.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK)  $(LIB_FFTW)

scoild_fft_VI: scoild_fft_VI.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o gammamatrix.o omegamatrix.o gradfvector.o omegasubspace.o  rhomatrix.o coilfft.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK)  $(LIB_FFTW)

scoild_fft_VII: scoild_fft_VII.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o gammamatrix.o omegamatrix.o gradfvector.o omegasubspace.o  rhomatrix.o coilfft.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK)  $(LIB_FFTW)

scoild_fft_VIII: scoild_fft_VIII.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o gammamatrix.o omegamatrix.o gradfvector.o omegasubspace.o  rhomatrix.o coilfft.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK)  $(LIB_FFTW)

scoild_fft_IX: scoild_fft_IX.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o gammamatrix.o omegamatrix.o gradfvector.o omegasubspace.o  rhomatrix.o coilfft.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK)  $(LIB_FFTW)

scoild_fft_X: scoild_fft_X.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o omegamatrix_fft.o  coilfft.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK)  $(LIB_FFTW)

scoild_fft_Xa: scoild_fft_Xa.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o omegamatrix_fft.o  coilfft.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK)  $(LIB_FFTW)

scoild_fft_XI: scoild_fft_XI.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o omegamatrix_fft.o  coilfft.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK)  $(LIB_FFTW)

scoild_fft_XIa: scoild_fft_XIa.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o omegamatrix_fft.o  coilfft.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK)  $(LIB_FFTW)

scoild_fft_XII: scoild_fft_XII.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o omegamatrix_fft.o  coilfft.o gradfvector.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK)  $(LIB_FFTW)

scoild_fft_XIIa: scoild_fft_XIIa.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o omegamatrix_fft.o  coilfft.o gradfvector.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK)  $(LIB_FFTW)

scoild_fft_XIII: scoild_fft_XIII.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o omegamatrix_fft.o  coilfft.o gradfvector.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK)  $(LIB_FFTW)

scoild_fft_XIIIa: scoild_fft_XIIIa.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o omegamatrix_fft.o  coilfft.o gradfvector.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK)  $(LIB_FFTW)

scoild_fft_XIV: scoild_fft_XIV.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o omegamatrix.o  coilfft.o gradfvector.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK)  $(LIB_FFTW)

scoild_fft_XV: scoild_fft_XV.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o omegamatrix.o  coilfft.o gradfvector.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK)  $(LIB_FFTW)






# These codes require the use of FORTRAN FIELD LINE CODE

calcBtraj_%:  bfield_plasma_%.o bfield_ext_%.o bfield_ext.o calcBtraj.o coils.o coilio.o coils_cmdline.o surface.o bfieldfuncs.o calcPoincare.o flsode.o
	$(LNK) $(LDFLAGS) $^ -o $@  $(LIBS_C) $(LIBS_FIELD)

calcBtrajfromCoils_%: bfield_plasma_%.o bfield_coils.o calcBtrajfromCoils.o coils.o coilio.o coils_cmdline.o surface.o bfieldfuncs.o calcPoincare.o flsode.o 
	$(LNK) $(LDFLAGS) $^ -o $@  $(LIBS_C) $(LIBS_FIELD)


#temp executable for trying stuff out
temp: temp.o coils.o coilio.o coils_cmdline.o surface.o 
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)



scoild_fft_II_test: scoild_fft_II_test.o coils.o coilio.o coils_cmdline.o surface.o inductancematrix.o gammamatrix.o omegamatrix.o gradfvector.o omegasubspace.o  coilfft.o
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C)  $(LIB_LAPACK)  $(LIB_FFTW)




