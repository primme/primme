# include the environment, compilation, and linking flags

include ../Make_flags
include ../Link_flags

#------------------------------------------------


SOBJS= COMMON/shared_utils.o COMMON/ioandtest.o
COMMON_INCLUDE = -I./COMMON -I../include
override INCLUDE += $(COMMON_INCLUDE)
LIBDIRS += -L../lib

USE_NATIVE    ?= yes
USE_PETSC     ?= $(if $(findstring undefined,$(origin PETSC_DIR)),no,yes)
USE_PARASAILS ?= $(if $(findstring undefined,$(origin PARASAILS_LIB_DIR)),no,yes)
USE_MPI       ?= $(if $(findstring mpi,$(CC)),yes,$(USE_PETSC))
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

primme_double: $(OBJSdouble) driverdouble.o
	$(CLDR) -o primme_double $(OBJSdouble) driverdouble.o $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

primme_doublecomplex: $(OBJSdoublecomplex) driverdoublecomplex.o
	$(CLDR) -o primme_doublecomplex $(OBJSdoublecomplex) driverdoublecomplex.o $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

primmesvds_double: $(OBJSdouble) driversvdsdouble.o
	$(CLDR) -o primmesvds_double $(OBJSdouble) driversvdsdouble.o $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

primmesvds_doublecomplex: $(OBJSdoublecomplex) driversvdsdoublecomplex.o
	$(CLDR) -o primmesvds_doublecomplex $(OBJSdoublecomplex) driversvdsdoublecomplex.o $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

%double.o: %.c
	$(CC) $(CFLAGS) $(DEFINES) -DUSE_DOUBLE $(INCLUDE) -c $< -o $@

%doublecomplex.o: %.c
	$(CC) $(CFLAGS) $(DEFINES) -DUSE_DOUBLECOMPLEX $(INCLUDE) -c $< -o $@

.f.o .F.o:
	$(F77) $(FFLAGS) $(FINCLUDE) $(COMMON_INCLUDE) -c $< -o $@

all: drivers

drivers: primme_double primme_doublecomplex primmesvds_double primmesvds_doublecomplex

primme_double primme_doublecomplex primmesvds_double primmesvds_doublecomplex: ../lib/libprimme.a

ifeq ($(USE_MPI), yes)
  MPIRUN ?= mpirun -np 4
endif

TESTS_double = tests/test_0?? tests/testi-*.F
TESTS_doublecomplex = tests/test_0?? tests/test_1?? tests/testi-*.F
TESTSsvds_double = tests/test_2??
TESTSsvds_doublecomplex = tests/test_2??
EXTRA =

all_tests: all_tests_double all_tests_doublecomplex all_testssvds_double all_testssvds_doublecomplex
all_tests_save: all_tests_double_save all_tests_doublecomplex_save all_testssvds_double_save all_testssvds_doublecomplex_save

all_tests_double all_tests_doublecomplex \
all_testssvds_double all_testssvds_doublecomplex: all_tests%: primme% tests_primme_interface
	@echo "Please wait, this could take a while...";\
	ok="0";for i in $(TESTS$*) ; do \
		echo "********** Test $$i $* ***********"; \
		sed 's/sol_[^ ]*/&$*/' $$i > ._test00;\
		echo "$(EXTRA)" >> ._test00;\
		$(MPIRUN) ./primme$* ._test00 || ok="1"; \
	done > tests.log 2>&1;\
	if test $$ok -eq 0 ; then \
		awk 'BEGIN{c=i=mv=t=0} /^Iterations/{i+=$$3;c++} /^Matvecs/{mv+=$$3} /^Wallclock/{t+=$$4} END{printf("Tests: %d Iterations: %d Matvecs: %d Time: %d s\n", c,i,mv,t)}' tests.log; \
		echo "All tests passed!"; \
	else\
		cat tests.log;\
		echo "Some tests fail. Please consider to send us the file";\
		echo "tests/tests.log if the software doesn't work as expected.";\
		exit 1;\
	fi

all_tests_double_save all_tests_doublecomplex_save \
all_testssvds_double_save all_testssvds_doublecomplex_save: all_tests%_save: primme% tests_primme_interface
	@echo "Please wait, this could take a while...";\
	ok="0";for i in $(TESTS$*) ; do \
		echo "********** Test $$i $* ***********"; \
		sed 's/sol_[^ ]*/&$*/; s/checkXFile/saveXFile/' $$i > ._test00;\
		$(MPIRUN) ./primme$* ._test00 || ok="1"; \
	done > tests.log 2>&1;\
	if test $$ok -eq 0 ; then \
		echo "All tests passed!"; \
	else\
		cat tests.log;\
		echo "Some tests fail. Please consider to send us the file";\
		echo "tests/tests.log if the software doesn't work as expected.";\
		exit 1;\
	fi

T_methods = DEFAULT_METHOD DYNAMIC DEFAULT_MIN_TIME DEFAULT_MIN_MATVECS Arnoldi GD_plusK GD_Olsen_plusK JD_Olsen_plusK JDQR JDQMR JDQMR_ETol STEEPEST_DESCENT LOBPCG_OrthoBasis LOBPCG_OrthoBasis_Window 
T_sizes = 0 1 2 3 4 5 6 7 10 100

tests_primme_interface: $(patsubst %,laplace%.mtx,$(T_sizes))
	@rm -f tests/testi-*; \
	for method in $(T_methods); do \
	for n in $(T_sizes); do \
	for nevals in 0 1 2 3 4 5 6 15 100; do \
	[ $$nevals -gt $$n ] && continue || true; \
	for target in primme_smallest primme_largest primme_closest_abs primme_closest_geq; do \
	case "$$target-$$n-$$nevals-$$method" in \
		primme_closest_leq-[123]-*) continue ;; \
		primme_closest_leq-[4567]-1-*) ;; \
		primme_closest_leq-[4567]-*) continue ;; \
		primme_closest_leq-10-[12]-*) ;; \
		primme_closest_leq-10-*) continue ;; \
		primme_closest_geq-4-4-*) continue ;; \
		primme_closest_geq-5-5-*) continue ;; \
		primme_closest_geq-6-6-*) continue ;; \
		primme_closest_geq-7-7-*) continue ;; \
		primme_closest_[gl]eq-100-*-LOBPCG*) continue ;; \
	esac; \
	for proj in primme_proj_RR primme_proj_refined ; do \
		case "$$target-$$proj-$$method" in \
			primme_closest*-primme_proj_RR-LOBPCG*) continue ;; \
			primme_closest*-STEEPEST_DESCENT*) continue ;; \
			primme_closest*-Arnoldi*) continue ;; \
			primme_closest*-GD*) continue ;; \
			primme_closest*) ;; \
			*-primme_proj_RR*) ;; \
			*) continue ;; \
		esac; \
		f="tests/testi-$${n}-$${method}-$${nevals}-$${target}-$${proj}"; \
		echo "driver.matrixFile = laplace$$n.mtx" > $$f.F; \
		echo "driver.checkXFile = tests/sol_testi-$${n}-$${nevals}-$${target}" >> $$f.F; \
		echo "driver.PrecChoice = noprecond" >> $$f.F; \
		echo "primme.numEvals = $$nevals" >> $$f.F; \
		echo "primme.eps = 1e-6" >> $$f.F; \
		echo "primme.numTargetShifts = 1" >> $$f.F; \
		echo "primme.targetShifts  = 0.5" >> $$f.F; \
		echo "primme.target = $$target" >> $$f.F; \
		echo "primme.projection.projection = $$proj" >> $$f.F; \
		echo "primme.maxMatvecs = 50000" >> $$f.F; \
		echo "method = PRIMME_$$method" >> $$f.F; \
	done; done; done; done; done

laplace%.mtx:
	@echo "%%MatrixMarket matrix coordinate real symmetric" > $@; \
	echo "$* $* $$(($*>0?2*$*-1:0))" >> $@; \
	for i in `seq 1 $*`; do \
		echo "$$i $$i 2.0" >> $@; \
		[  $$i -eq $* ] || echo "$$i $$((i+1)) -1.0" >> $@; \
	done

clean:
	@rm -f $(OBJSdouble) $(OBJSdoublecomplex) *.o tests.log tests/*.F $(patsubst %,laplace%.mtx,$(T_sizes)) ._test00

veryclean: clean
	@rm -f primme_double primme_doublecomplex primmesvds_double primmesvds_doublecomplex


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

.PHONY: clean veryclean all drivers examples
