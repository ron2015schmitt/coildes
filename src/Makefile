
include ../Makefile

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

exec: $(EXEC)

tests: $(TESTS)
