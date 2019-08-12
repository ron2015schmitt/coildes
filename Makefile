

SHELL=/bin/bash
COOLLDIR = /home/rs2015/intel_versions/cool_1.91beta
LSODEDIR = /home/rs2015/lsode
FFTWDIR = /home/rs2015/fftw-3.0.1


#notes:  -v option on g++/gcc is useful for 

# use this variable to set gcc/g++ compile options from command line
# just write 'CLOPTS=yada yada yada' at end of command line
# example for generating gprof output:
#   make coils 'COPTS=-pg' 'LOPTS=-pg'

COPTS = 
LOPTS =

#-O1    optimize for maximum speed, but disable some optimizations which
#       increase code size for a small speed benefit.
#-O2    enable optimizations (DEFAULT)
#-O3    enable -O2 plus more aggressive optimizations that may not improve
#       performance for all programs
#-O0    disable optimizations
#-O     same as -O2
#-Os    enable speed optimizations, but disable some optimizations which
#       increase code size for small speed benefit
#-ax<codes> generate code specialized for processors specified by <codes>
#           while also generating generic IA-32 code.  <codes> includes
#           one or more of the following characters:
#    K  Intel Pentium III and compatible Intel processors
#    W  Intel Pentium 4 and compatible Intel processors
#    N  Intel Pentium 4 and compatible Intel processors.  Enables new
#       optimizations in addition to Intel processor-specific optimizations
#    P  Intel Pentium 4 processors with SSE3 extensions
#    B  Intel Pentium M and compatible Intel processors
#-x<codes>  generate specialized code to run exclusively on processors
#           indicated by <codes> as described above.
#
#Enable and specify the scope of Interprocedural (IP) Optimizations:
#-ip     enable single-file IP optimizations (within files)
#-ipo[n] enable multi-file IP optimizations (between files)
#-ipo-c  generate a multi-file object file (ipo_out.o)
#-ip-no-inlining    disable full and partial inlining (requires -ip or -ipo)
#-ip-no-pinlining   disable partial inlining (requires -ip or -ipo)
#-ipo-separate      create one object file for every source file
#                   (overrides -ipo[n])
#-w                 disable all warnings
#-w<n>              control diagnostics:
#   n=0               display errors (same as -w)
#   n=1               display warnings and errors (DEFAULT)
#   n=2               display remarks, warnings, and errors
#-static        prevents linking with shared libraries
#-rcd          rounding mode to enable fast float-to-int conversions

CXXOPT = -w -xW #-O3 

#use either (CXX -no-ipo and LNK -no-ipo) or (CXX -ipo-c and LNK -ipo)
CXX = /usr/local/intel/cc/9.0/bin/icpc -no-ipo 


#-i-dynamic     link Intel provided libraries dynamically
#-i-static      link Intel provided libraries statically
#-dynamic-linker<file>
#               select dynamic linker other than the default
#-no-cpprt      do not link in C++ runtime libraries
#-nodefaultlibs do not use standard libraries when linking
#-nostartfiles  do not use standard startup files when linking
#-nostdlib      do not use standard libraries and startup files when linking
#-static        prevents linking with shared libraries
#-shared        produce a shared object
#-static-libcxa link Intel libcxa C++ library statically
#-shared-libcxa link Intel libcxa C++ library dynamically, overrides the default
#               behavior when -static is used
#-cxxlib-<mode> tell the compiler which C++ run-time libraries to use
#               gcc[=dir] - link using C++ run-time libraries provided with gcc
#                           (default on systems running gcc 3.2 or above)
#                           dir is an optional top-level location for the gcc
#                           binaries and libraries
#               icc       - link using C++ run-time libraries provided by Intel
#                           (default on systems running a gcc version lower
#                            than 3.2)
LNKOPT=-no-ipo  #-cxxlib-gcc #-i-static -static #-nodefaultlibs
LNK = /usr/local/intel/cc/9.0/bin/icpc 


INC_FFTW = -I $(FFTWDIR)/include 
INC_COOLL = -I $(COOLLDIR) 
INCLUDES = $(INC_COOLL) $(INC_FFTW)

# intel lapack
LIB_LAPACK =  -L/usr/local/intel/mkl/8.0/lib/32 -lmkl_lapack -lmkl_ia32 -lguide -lpthread
# standard lapack
#LIB_LAPACK= -llapack 
LIB_FFTW = -L$(FFTWDIR)/lib -lfftw3 

LIB_COOLL = -L$(COOLLDIR) -lcooll

LIB_LSODE= -L$(LSODEDIR) -llsode

LIB_NAG= -L/usr/local/lib -lnag 

# do we need -lircmt?
#LIB_FORTRAN = -lfio -lblas -lf77math -lf90math -L/usr/local/intel/fc/9.0/lib -lifcore -lifcoremt -lcxa -lcprts -lunwind -lirc  -lsvml -lompstub -limf  -lcxaguard -lifport -lguide -lguide_stats 
LIB_FORTRAN = -L/usr/local/intel/fc/9.0/lib -lifcore


LIBS =  $(LIB_COOLL) 
LIBS_C =  -L/usr/local/intel/cc/9.0/lib $(LIB_COOLL) 
LIBS_FIELD = $(LIB_LSODE)  $(LIB_NAG) $(LIB_FORTRAN)


CXXFLAGS =  $(CXXOPT)  $(COPTS) $(INCLUDES)

LDFLAGS =  $(LOPTS) $(LNKOPT)  

INC = 



EXEC_MAIN = Bnormal2current flux2current current2Bnormal current2flux scoild
EXEC_MAIN_FFT = Bnormal2current_fft flux2current_fft current2Bnormal_fft current2flux_fft scoild_fft scoild_fft_II scoild_fft_IIb scoild_fft_III scoild_fft_IV scoild_fft_V scoild_fft_VI scoild_fft_VII scoild_fft_VIII scoild_fft_IX scoild_fft_X scoild_fft_XI scoild_fft_Xa scoild_fft_XIa  scoild_fft_XII scoild_fft_XIIa scoild_fft_XIII scoild_fft_XIIIa  scoild_fft_XIV  scoild_fft_XV
#EXEC_BTRAJ =  calcBtraj calcBtrajfromCoils
#EXEC_BONSURF = calcBonsurf calcBonsurffromCoils 
EXEC_AUXILIARY =  mkcoilsurf  expandsurf surfF findpertsurface deltasurfs deltaflux percentdelta reduceflux rankorder findpertsurface_fft expandfunc calcIpol onlyharms expandI cullF cutcoils iotafromIF projectBw divfuncs calcGammaOmega calcOmegaN calcOmegaN_fft cart2contravariantF cart2contravariantJF b2lambda lambda2b calcOmegas calcMs calcJ testfft
EXEC_OTHER = garabedian2sincos garabedian2sincosII  garabedian_current2sincos
EXEC_TEMP = temp

EXEC = $(EXEC_MAIN) $(EXEC_MAIN_FFT)  $(EXEC_AUXILIARY) $(EXEC_OTHER) 


TESTS = 

CPROGS = 


#prevent any default rules from being used
.SUFFIXES:

# don't delete .o files
.PRECIOUS: %.o calcBtraj_% calcBtrajfromCoils_% calcBonsurf_% calcBonsurffromCoils_% 

all: $(EXEC)




tests: $(TEST)

clean:
	@command rm -f *.o
	@command rm -f $(EXEC)
	@command rm -f $(TESTS)
	@command rm -f core.*
	@command rm calcBtraj_* calcBtrajfromCoils_* calcBonsurf_* calcBonsurffromCoils_* 
def:
	echo "nothing done"




%.o: %.cpp coils.hpp 
ifdef CAREFUL
	$(CXX) $(CXXFLAGS) -D"COOLL_CAREFUL=1" -c $*.cpp -o $@
else
	$(CXX) $(CXXFLAGS) -c $*.cpp -o $@
endif
 

%: %.o 
	$(CXX) $(LDFLAGS) $*.o -o $@ $(LIBS_C) 

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




