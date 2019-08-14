


SHELL:=/bin/bash
RED:='\033[0;31m'
CYAN='\033[0;36m'
NC:='\033[0m' # No Color
BOLDON:='\e[1m'
BOLDOFF:='\e[0m'


define echovar
 echo -e $(BOLDON)"${1}"$(BOLDOFF)" = ${${1}}"
endef

define title
 echo -e $(CYAN)"${1}"$(NC)
endef

define hr
 echo -e $(BOLDON)"-------------------------------------------------------------------------------"$(BOLDOFF)
endef


###########################################################
#      Directory Structure
###########################################################

# this resides in the bin subdirectory
DIR_COILDES:=$(shell coildesdir)
ifndef DIR_COILDES
$(error PATH is not set properly. Refer to coildes configuration instructions.)
endif

DIR_SRC = $(DIR_COILDES)/src
DIR_BIN = $(DIR_COILDES)/bin



###########################################################
#      MATRICKS submodule 
###########################################################

NAME_MATRICKS = matricks
DIR_MATRICKS = $(DIR_COILDES)/matricks
INC_MATRICKS = -I $(DIR_MATRICKS) 
LIB_MATRICKS = -L$(DIR_MATRICKS) -l$(NAME_MATRICKS)
LIBFILE_MATRICKS = $(DIR_MATRICKS)/lib$(NAME_MATRICKS).a

INCLUDES += $(INC_MATRICKS)
LIBS += $(LIB_MATRICKS)

###########################################################
#      ODEPACK submodule 
# the odepack repo has no Makefiles, so we copy source
# and build in our src directory
###########################################################
NAME_ODEPACK = odepack
DIR_ODEPACK = $(DIR_COILDES)/odepack/src
INC_ODEPACK =  
BLD_ODEPACK = $(DIR_COILDES)/src
LIB_ODEPACK = -L$(BLD_ODEPACK) -l$(NAME_ODEPACK)
LIBFILE_ODEPACK = $(BLD_ODEPACK)/lib$(NAME_ODEPACK).a

INCLUDES += $(INC_ODEPACK)
LIBS += $(LIB_ODEPACK)


#################################################################
#  COMPILERS
#
# Note from
# https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html
# "If you use multiple -O options, with or without level
#  numbers, the last such option is the one that is effective."
#
# COPTS provided at the command line will supercede the following
#       That is the values below will be written over.
# CADD provided at the command line will append to values below
##################################################################
# C++ compiler
CPPC = g++
COPTS = -finline-functions -finline-limit=750 -O0
CADD ?= 
CFLAGS = $(COPTS) $(CADD) $(INCLUDES)


# FORTRAN COMPILER
FC = gfortran
LNKOPT_FORTRAN = -lgfortran

# linker
override LOPTS +=
LFLAGS = $(LOPTS)
LNK = g++



###########################################################
#  TARGETS
###########################################################

#prevent any default rules from being used
.SUFFIXES:

# don't delete .o files
.PRECIOUS: %.o




%.o: %.f 
	$(FC) -c $*.f -o $*.o


%.o: %.cpp 
ifdef CAREFUL
	$(CPPC) $(CFLAGS) -D"MATRICKS_CAREFUL=1" -c $*.cpp -o $@
else
	$(CPPC) $(CFLAGS) -c $*.cpp -o $@
endif


%: %.o 
	$(CPPC) $(LFLAGS) $*.o -o $@ $(LIBS) 

###########################################################
#      matricks submodule build
###########################################################
$(LIBFILE_MATRICKS): 
	cd matricks && ./configure

info_matricks:
	@echo
	@$(call hr)
	@$(call echovar,DIR_MATRICKS)
	@$(call echovar,LIB_MATRICKS)
	@$(call echovar,INC_MATRICKS)
	@$(call echovar,LIBFILE_MATRICKS)
	@$(call hr)
	@echo


###########################################################
#      ODEPACK submodule build
###########################################################
REQ_ODEPACK = opkda1 opkda2 opkdmain
REQ_ODEPACK_O = $(addsuffix .o, $(REQ_ODEPACK))
REQ_ODEPACK_F = $(addsuffix .f, $(REQ_ODEPACK))

# copy the .f files to *our* source directory since ODEPACK provide no makefile
$(addprefix $(BLD_ODEPACK)/, $(REQ_ODEPACK_F) ):
	cp $(DIR_ODEPACK)/$(@F) $(BLD_ODEPACK)/


$(LIBFILE_ODEPACK): $(addprefix $(BLD_ODEPACK)/, $(REQ_ODEPACK_O) )
	echo "building ODEPACK library in $(BLD_ODEPACK) ..."
	cd $(BLD_ODEPACK) && ar rUuv lib$(NAME_ODEPACK).a $(^F)


info_odepack:
	@echo
	@$(call hr)
	@$(call title,"ODEPACK")
	@$(call echovar,DIR_ODEPACK)
	@$(call echovar,LIB_ODEPACK)
	@$(call echovar,INC_ODEPACK)
	@$(call echovar,BLD_ODEPACK)
	@$(call echovar,LIBFILE_ODEPACK)
	@$(call echovar,REQ_ODEPACK)
	@$(call echovar,REQ_ODEPACK_F)
	@$(call echovar,REQ_ODEPACK_O)
	@echo $(BLD_ODEPACK)/opkda1.f
	@$(call hr)
	@echo

###########################################################
#      General Stuff
###########################################################

info:
	@echo
	@$(call hr)
	@$(call title,"Directories")
	@$(call echovar,DIR_COILDES)
	@$(call echovar,DIR_SRC)
	@$(call echovar,DIR_BIN)
	@echo
	@$(call title,"Includes")
	@$(call echovar,INCLUDES)
	@echo
	@$(call title,"Libraries")
	@$(call echovar,LIBS)
	@echo
	@$(call title,"C++ Compiler")
	@$(call echovar,CPPC)
	@$(call echovar,COPTS)
	@$(call echovar,CFLAGS)
	@echo
	@$(call title,"Fortran")
	@$(call echovar,FC)
	@$(call echovar,LOPTS_FORTRAN)
	@echo
	@$(call title,"Linker")
	@$(call echovar,LNK)
	@$(call echovar,LOPTS)
	@$(call echovar,LFLAGS)
	@$(call hr)
	@echo

clean:
	@command rm -f *.o
	@command rm -f *.a
	@command rm -f *~
	@command rm -f $(EXEC)
	@command rm -f $(TESTS)
	@command rm -f core.*

cleanall:
	@cd $(DIR_COILDES) && make clean
	@cd $(DIR_SRC) && make clean
	@cd $(DIR_SRC) && rm -fr $(REQ_ODEPACK_F)
	@cd $(DIR_BIN) && make clean


fullpurge: cleanall
#	@make --makefile $(DIR_COILDES)/Makefile --directory $(DIR_ODEPACK) clean
	@echo "run $(DIR_MATRICKS)/deconfigure ..."
	@cd $(DIR_MATRICKS) && ./deconfigure


all:
	@cd $(DIR_SRC) && make exec tests

def:
	echo "nothing done"
