# include the environment, compilation, and linking flags

include ../Make_flags
include ../Link_flags

#------------------------------------------------


SOBJS= driver.o COMMON/shared_utils.o
COMMON_INCLUDE = -I./COMMON -I../PRIMMESRC/COMMONSRC
override INCLUDE += $(COMMON_INCLUDE)
LIBDIRS += -L../

USE_NATIVE    ?= yes
USE_PETSC     ?= $(if $(findstring undefined,$(origin PETSC_DIR)),no,yes)
USE_PARASAILS ?= $(if $(findstring undefined,$(origin PARASAILS_LIB_DIR)),no,yes)
USE_MPI       ?= $(if $(findstring mpi,$(CC)),yes,no)

ifeq ($(USE_MPI), yes)
  DEFINES += -DUSE_MPI
endif

ifeq ($(USE_NATIVE), yes)
  DEFINES += -DUSE_NATIVE
  SOBJS += COMMON/csr.o COMMON/mat.o COMMON/ssrcsr.o COMMON/mmio.o
  SOBJSdouble += COMMON/ilut.o COMMON/amux.o
  SOBJSdoublecomplex += COMMON/zilut.o COMMON/zamux.o
endif

ifeq ($(USE_PARASAILS), yes)
  ifneq ($(USE_MPI), yes)
    $(error "PARASAILS needs MPI")
  endif
  DEFINES += -DUSE_PARASAILS
  SOBJS += COMMON/parasailsw.o COMMON/csr.o COMMON/ssrcsr.o COMMON/mmio.o
  override INCLUDE += -I$(PARASAILS_INCLUDE_DIR)
  LIBDIRS += -L$(PARASAILS_LIB_DIR)
  LIBS += -lParaSails
endif

ifeq ($(USE_PETSC), yes)
  include ${PETSC_DIR}/lib/petsc/conf/variables
  DEFINES += -DUSE_PETSC
  SOBJS += COMMON/petscw.o COMMON/mmio.o
  override INCLUDE += $(PETSC_CCPPFLAGS)
  FINCLUDE += $(PETSC_FCPPFLAGS)
  LIBDIRS += $(PETSC_C_SH_LIB_PATH)
  LIBS += $(PETSC_LIB)
endif


OBJS = $(sort $(SOBJS))
OBJSdouble = $(sort $(SOBJSdouble)) $(patsubst %.o,%double.o,$(OBJS))
OBJSdoublecomplex = $(sort $(SOBJSdoublecomplex)) $(patsubst %.o,%doublecomplex.o,$(OBJS))

.PHONY: clean veryclean

primme_double: $(OBJSdouble) ../libprimme.a 
	$(CLDR) -o primme_double $(OBJSdouble) $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

primme_doublecomplex: $(OBJSdoublecomplex) ../libprimme.a 
	$(CLDR) -o primme_doublecomplex $(OBJSdoublecomplex) $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

ex_dseq: ex_dseq.o ../libprimme.a 
	$(CLDR) -o ex_dseq ex_dseq.o $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

ex_zseq: ex_zseq.o ../libprimme.a 
	$(CLDR) -o ex_zseq ex_zseq.o $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

ex_petsc: ex_petsc.o ../libprimme.a 
	$(CLDR) -o ex_petsc ex_petsc.o $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

ex_dseqf77: ex_dseqf77.o ../libprimme.a 
	$(FLDR) -o ex_dseqf77 ex_dseqf77.o $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

ex_zseqf77: ex_zseqf77.o ../libprimme.a 
	$(FLDR) -o ex_zseqf77 ex_zseqf77.o $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

ex_petscf77: ex_petscf77.o ../libprimme.a 
	$(FLDR) -o ex_petscf77 ex_petscf77.o $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

ex_petscf77ptr: ex_petscf77ptr.o ../libprimme.a 
	$(FLDR) -o ex_petscf77ptr ex_petscf77ptr.o $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

%double.o: %.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c $< -o $@

%doublecomplex.o: %.c
	$(CC) $(CFLAGS) $(DEFINES) -DUSE_DOUBLECOMPLEX $(INCLUDE) -c $< -o $@

.c.o:
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c $< -o $@

.f.o .F.o:
	$(F77) $(FFLAGS) $(FINCLUDE) $(COMMON_INCLUDE) -c $< -o $@

all: primme_double primme_doublecomplex

ifeq ($(USE_PETSC), yes)
simple_examples: ex_dseq ex_zseq ex_petsc ex_dseqf77 ex_zseqf77 ex_petscf77 ex_petscf77ptr
else
simple_examples: ex_dseq ex_zseq ex_dseqf77 ex_zseqf77
endif

ifeq ($(USE_MPI), yes)
  MPIRUN ?= mpirun -np 4
endif

all_tests_double: primme_double
	@echo "Please wait, this could take a while...";\
	ok="0";for i in tests/test_??? ; do \
		echo "********** Test $$i double ***********"; \
		$(MPIRUN) ./primme_double $$i || ok="1"; \
	done > tests.log 2>&1;\
	if test $$ok -eq 0 ; then \
		echo "All tests passed!"; \
	else\
		cat tests.log;\
		echo "Some tests fail. Please consider to send us the file";\
		echo "TEST/tests.log if the software doesn't work as expected.";\
		exit 1;\
	fi

all_tests_doublecomplex: primme_doublecomplex
	@echo "Please wait, this could take a while...";\
	ok="0";	for i in tests/test_??? ; do \
		echo "********** Test $$i complex ***********"; \
		sed 's/sol_.../&.complex/' $$i > ._test00;\
		$(MPIRUN) ./primme_doublecomplex ._test00 || ok="1"; \
		rm ._test00;\
	done > tests.log 2>&1;\
	if test $$ok -eq 0 ; then \
		echo "All tests passed!"; \
	else\
		cat tests.log;\
		echo "Some tests fail. Please consider to send us the file";\
		echo "TEST/tests.log if the software doesn't work as expected.";\
		exit 1;\
	fi


all_tests_double_save: primme_double
	@for i in tests/test_??? ; do \
		sed 's/checkXFile/saveXFile/' $$i > ._test00;\
		echo "********** Test $$i double ***********"; \
		./primme_double ._test00 || ok="1"; \
	done > tests.log 2>&1;\
	rm -f ._test00

all_tests_doublecomplex_save: primme_doublecomplex
	@for i in tests/test_??? ; do \
		sed 's/checkXFile/saveXFile/;s/sol_.../&.complex/' $$i > ._test00;\
		echo "********** Test $$i complex ***********"; \
		./primme_doublecomplex ._test00 || ok="1"; \
	done > tests.log 2>&1;\
	rm -f ._test00

test_double: ex_dseq
	@./ex_dseq > tests.log 2>&1 && grep -q "10 eigenpairs converged" tests.log \
	    && echo "Test passed!" || \
	       (cat tests.log; echo "Something went wrong. Please consider to send us the file";\
                echo "TEST/tests.log if the software doesn't work as expected.";\
                exit 1;)

test_doublecomplex: ex_zseq
	@./ex_zseq > tests.log 2>&1 && grep -q "10 eigenpairs converged" tests.log \
	    && echo "Test passed!" || \
	       (cat tests.log; echo "Something went wrong. Please consider to send us the file";\
                echo "TEST/tests.log if the software doesn't work as expected.";\
                exit 1;)

test: test_double test_doublecomplex

clean:
	@rm -f $(OBJSdouble) $(OBJSdoublecomplex) *.o

veryclean: clean
	@rm -f primme_double primme_doublecomplex seqf77_dprimme seqf77_zprimme ex_dseq ex_zseq ex_petsc


COMMON/csr.c: COMMON/csr.h COMMON/mmio.h
COMMON/csr.h: COMMON/num.h
COMMON/mat.c: COMMON/native.h
COMMON/mmio.c: COMMON/mmio.h
COMMON/native.h: COMMON/csr.h
COMMON/parasailsw.c: COMMON/parasailsw.h COMMON/csr.h
COMMON/parasailsw.h: COMMON/csr.h
COMMON/petscw.c: COMMON/petscw.h COMMON/mmio.h
COMMON/petscw.h: COMMON/num.h
COMMON/shared_utils.c: COMMON/shared_utils.h
COMMON/driver.c: COMMON/shared_utils.h COMMON/native.h COMMON/parasailsw.h COMMON/petscw.h
