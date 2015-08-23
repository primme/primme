# include the environment, compilation, and linking flags

include ../Make_flags
include ../Link_flags

#------------------------------------------------


SOBJS= driver.o COMMON/shared_utils.o
INCLUDE += -I./COMMON -I../PRIMMESRC/COMMONSRC
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
  INCLUDE += -I$(PARASAILS_INCLUDE_DIR)
  LIBDIRS += -L$(PARASAILS_LIB_DIR)
  LIBS += -lParaSails
endif

ifeq ($(USE_PETSC), yes)
  include ${PETSC_DIR}/lib/petsc/conf/variables
  DEFINES += -DUSE_PETSC
  SOBJS += COMMON/petscw.o COMMON/mmio.o
  INCLUDE += $(PETSC_CCPPFLAGS)
  LIBDIRS += $(PETSC_C_SH_LIB_PATH)
  LIBS += $(PETSC_LIB)
endif


OBJS = $(sort $(SOBJS))
OBJSdouble = $(sort $(SOBJSdouble)) $(patsubst %.o,%double.o,$(OBJS))
OBJSdoublecomplex = $(sort $(SOBJSdoublecomplex)) $(patsubst %.o,%doublecomplex.o,$(OBJS))

.PHONY: clean veryclean

primme_double: $(OBJSdouble) ../libprimme.a Makefile
	$(CLDR) -o primme_double $(OBJSdouble) $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

primme_doublecomplex: $(OBJSdoublecomplex) ../libprimme.a Makefile
	$(CLDR) -o primme_doublecomplex $(OBJSdoublecomplex) $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

ex_dseq: ex_dseq.o ../libprimme.a Makefile
	$(CLDR) -o ex_dseq ex_dseq.o $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

%double.o: %.c
	$(CC) $(HIGH_OPT) $(CFLAGS) $(DEFINES) $(INCLUDE) -c $< -o $@

%doublecomplex.o: %.c
	$(CC) $(HIGH_OPT) $(CFLAGS) $(DEFINES) -DUSE_DOUBLECOMPLEX $(INCLUDE) -c $< -o $@

.c.o:
	$(CC) $(HIGH_OPT) $(CFLAGS) $(DEFINES) $(INCLUDE) -c $< -o $@

.f.o:
	$(F77) $(HIGH_OPT) $(FFLAGS) $(INCLUDE) -c $< -o $@

all: primme_double primme_doublecomplex

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
		echo "DTEST/tests.log if the software doesn't work as expected.";\
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
		echo "DTEST/tests.log if the software doesn't work as expected.";\
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

test_double: primme_double
	./primme_double tests/test_003

test_doublecomplex: primme_doublecomplex
	sed 's/sol_.../\0.complex/' tests/test_003 > ._test00;\
	./primme_doublecomplex ._test00;\
	rm -f ._test00

test: test_double test_doublecomplex


seqf77_dprimme_FOBJS = driver_f77seq.o COMMON/ilut.o COMMON/amux.o
seqf77_dprimme: $(seqf77_dprimme_FOBJS) ../libprimme.a
	$(FLDR) -o seqf77_dprimme $(seqf77_dprimme_FOBJS) $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

seqf77_zprimme_FOBJS = driver_f77seqz.o COMMON/zamux.o
seqf77_zprimme: $(seqf77_zprimme_FOBJS) ../libprimme.a
	$(FLDR) -o seqf77_zprimme $(seqf77_zprimme_FOBJS) $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

clean:
	@rm -f $(OBJSdouble) $(OBJSdoublecomplex) $(seqf77_dprimme_FOBJS) $(seqf77_zprimme_FOBJS)

veryclean: clean
	@rm -f primme_double primme_doublecomplex seqf77_dprimme seqf77_zprimme


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
