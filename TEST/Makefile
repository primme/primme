# include the environment, compilation, and linking flags

include ../Make_flags
include ../Link_flags

#------------------------------------------------


SOBJS= COMMON/shared_utils.o COMMON/ioandtest.o
COMMON_INCLUDE = -I./COMMON -I../PRIMMESRC/COMMONSRC -I../PRIMMESRC/SVDS/COMMONSRC
override INCLUDE += $(COMMON_INCLUDE)
LIBDIRS += -L../

USE_NATIVE    ?= yes
USE_PETSC     ?= $(if $(findstring undefined,$(origin PETSC_DIR)),no,yes)
USE_PARASAILS ?= $(if $(findstring undefined,$(origin PARASAILS_LIB_DIR)),no,yes)
USE_MPI       ?= $(if $(findstring mpi,$(CC)),yes,no)
USE_RSB       ?= $(if $(findstring undefined,$(origin LIBRSB_LIB_DIR)),no,yes)

ifeq ($(USE_MPI), yes)
  DEFINES += -DUSE_MPI
endif

ifeq ($(USE_NATIVE), yes)
  DEFINES += -DUSE_NATIVE
  SOBJS += COMMON/csr.o COMMON/mat.o COMMON/ssrcsr.o COMMON/mmio.o
  SOBJSdouble += COMMON/ilut.o COMMON/matvec.o
  SOBJSdoublecomplex += COMMON/zilut.o COMMON/zmatvec.o
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

ifeq ($(USE_RSB), yes)
  DEFINES += -DUSE_RSB
  SOBJS += COMMON/rsbw.o
  override INCLUDE += -I$(LIBRSB_INCLUDE_DIR)
  LIBDIRS += -L$(LIBRSB_LIB_DIR)
  override LIBS += -lrsb
  override LDFLAGS += -fopenmp
  override CFLAGS += -fopenmp
endif


OBJS = $(sort $(SOBJS))
OBJSdouble = $(sort $(SOBJSdouble)) $(patsubst %.o,%double.o,$(OBJS))
OBJSdoublecomplex = $(sort $(SOBJSdoublecomplex)) $(patsubst %.o,%doublecomplex.o,$(OBJS))

.PHONY: clean veryclean

primme_double: $(OBJSdouble) driverdouble.o ../libprimme.a Makefile
	$(CLDR) -o primme_double $(OBJSdouble) driverdouble.o $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

primme_doublecomplex: $(OBJSdoublecomplex) driverdoublecomplex.o ../libprimme.a Makefile
	$(CLDR) -o primme_doublecomplex $(OBJSdoublecomplex) driverdoublecomplex.o $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

ex_dseq ex_zseq ex_petsc exsvd_dseq: %: %.o ../libprimme.a Makefile
	$(CLDR) -o $@ $@.o $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

ex_dseqf77 ex_zseqf77 ex_petscf77 ex_petscf77ptr: % : %.o ../libprimme.a Makefile
	$(FLDR) -o $@ $@.o $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

primmesvds_double: $(OBJSdouble) driversvdsdouble.o ../libprimme.a Makefile
	$(CLDR) -o primmesvds_double $(OBJSdouble) driversvdsdouble.o $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

primmesvds_doublecomplex: $(OBJSdoublecomplex) driversvdsdoublecomplex.o ../libprimme.a Makefile
	$(CLDR) -o primmesvds_doublecomplex $(OBJSdoublecomplex) driversvdsdoublecomplex.o $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 


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
	ok="0";for i in tests/test_0?? ; do \
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
	@for i in tests/test_0?? ; do \
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
COMMON/ioandtest.c: COMMON/num.h COMMON/ioandtest.h
COMMON/driver.c: COMMON/shared_utils.h COMMON/native.h COMMON/parasailsw.h COMMON/petscw.h

ex_dseq.c: ex_c.m4
	m4 ex_c.m4 > $@

ex_zseq.c: ex_c.m4
	m4 -D USE_COMPLEX -D ADVANCED ex_c.m4 > $@

ex_petsc.c: ex_c.m4
	m4 -D USE_PETSC ex_c.m4 > $@

ex_dseqf77.f: ex_f77.m4
	m4 ex_f77.m4 > $@

ex_zseqf77.f: ex_f77.m4
	m4 -D USE_COMPLEX ex_f77.m4 > $@

ex_petscf77.F: ex_f77.m4
	m4 -D USE_PETSC ex_f77.m4 > $@

ex_petscf77ptr.F: ex_f77.m4
	m4 -D USE_PETSC -D USE_POINTER ex_f77.m4 > $@

exsvd_dseq.c: exsvd_c.m4
	m4 exsvd_c.m4 > $@

exsvd_dseq0.c: exsvd_c.m4
	m4 -D ADVANCED_HYBRID exsvd_c.m4 > $@
