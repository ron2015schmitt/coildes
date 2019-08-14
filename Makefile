


SHELL:=/bin/bash
RED:='\033[0;31m'
NC:='\033[0m' # No Color
BOLDON:='\e[1m'
BOLDOFF:='\e[0m'



# this resides in the bin subdirectory
COILDES_DIR:=$(shell coildesdir)
ifndef COILDES_DIR
$(error PATH is not set properly. Refer to coildes configuration instructions.)
endif


NAME_MATRICKS = matricks
DIR_MATRICKS = $(COILDES_DIR)/matricks
INC_MATRICKS = -I $(DIR_MATRICKS) 
LIB_MATRICKS = -L$(DIR_MATRICKS) -l$(NAME_MATRICKS)
LIBFILE_MATRICKS = $(DIR_MATRICKS)/lib$(NAME_MATRICKS).a


NAME_ODEPACK = odepack
DIR_ODEPACK = $(COILDES_DIR)/odepack/src
INC_ODEPACK = -I $(DIR_ODEPACK) 
LIB_ODEPACK = -L$(DIR_ODEPACK) -l$(NAME_ODEPACK)
LIBFILE_ODEPACK = $(DIR_ODEPACK)/lib$(NAME_ODEPACK).a
FILE_DEMO_LSODE = $(DIR_ODEPACK)/demos/LSODE.o



INCLUDES = $(INC_MATRICKS) $(INC_ODEPACK)
LIBS =  $(LIB_MATRICKS) $(LIB_ODEPACK)


# C++ compiler
CPPOPT = -finline-functions -finline-limit=750 -O0
CPPC = g++
CPPFLAGS = $(CPPOPT) $(INCLUDES)

# linker
LNKOPT =
LDFLAGS = $(LNKOPT)
LNKOPT_FORTRAN = -lgfortran
LNK = g++

# FORTRAN COMPILER
FC = gfortran



#INC_FFTW = -I $(FFTWDIR)/include 
#LIB_LAPACK= -llapack 
#LIB_FFTW = -L$(FFTWDIR)/lib -lfftw3 
#FFTWDIR = /home/rs2015/fftw-3.0.1
#LIBS_FIELD = $(LIB_LSODE)  $(LIB_NAG) $(LIB_FORTRAN)


#prevent any default rules from being used
.SUFFIXES:

# don't delete .o files
.PRECIOUS: %.o




%.o: %.f 
	$(FC) -c $*.f -o $*.o


%.o: %.cpp 
ifdef CAREFUL
	$(CPPC) $(CPPFLAGS) -D"MATRICKS_CAREFUL=1" -c $*.cpp -o $@
else
	$(CPPC) $(CPPFLAGS) -c $*.cpp -o $@
endif


%: %.o 
	$(CPPC) $(LDFLAGS) $*.o -o $@ $(LIBS) 

###########################################################
#      matricks submodule build
###########################################################
$(LIBFILE_MATRICKS): 
	cd matricks && ./configure

###########################################################
#      ODEPACK submodule build
###########################################################
$(LIBFILE_ODEPACK):  $(DIR_ODEPACK)/opkda1.o  $(DIR_ODEPACK)/opkda2.o  $(DIR_ODEPACK)/opkdmain.o  
	@echo "building ODEPACK library..."
	cd $(DIR_ODEPACK) && ar rUuv libodepack.a $(^F)


###########################################################
# EXEC_ODEPACK
###########################################################


test_odepack: test_odepack.o  $(LIBFILE_ODEPACK)
	$(LNK) $(LDFLAGS) $^ -o $@ $(LIBS_C) $(LNKOPT_FORTRAN)

EXEC_ODEPACK += test_odepack
EXEC += $(EXEC_ODEPACK)
TESTS += test_odepack





###########################################################
# EXEC
###########################################################

all: $(EXEC)

tests: $(TESTS)

clean:
	@command rm -f *.o
	@command rm -f *~
	@command rm -f $(EXEC)
	@command rm -f $(TESTS)
	@command rm -f core.*

fullpurge:
	@cd $(COILDES_DIR) && make clean
	@cd $(COILDES_DIR)/src && make clean
	@cd $(COILDES_DIR)/bin && make clean
	@cd $(COILDES_DIR)/bin && make clean
	@cd matricks && ./deconfigure
	@cd $(DIR_ODEPACK) && make clean

def:
	echo "nothing done"
