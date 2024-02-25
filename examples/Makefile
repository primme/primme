# include the environment, compilation, and linking flags

include ../Make_flags
include ../Link_flags

#------------------------------------------------


INCLUDE := -I../include
FINCLUDE := -I../include
LIBS := ../lib/libprimme.a $(LIBS)

EXAMPLES_C = ex_eigs_dseq ex_eigs_zseq ex_eigs_zseq_normal ex_svds_dseq ex_svds_zseq
EXAMPLES_CXX = ex_eigs_zseqxx ex_svds_zseqxx
EXAMPLES_F = ex_eigs_dseqf77 ex_eigs_zseqf77 ex_svds_dseqf77 ex_svds_zseqf77 ex_eigs_dseqf90 ex_svds_dseqf90

USE_PETSC     ?= $(if $(findstring undefined,$(origin PETSC_DIR)),no,yes)
USE_MPI       ?= $(if $(findstring mpi,$(CC)),yes,$(USE_PETSC))

ifeq ($(USE_PETSC), yes)
  include ${PETSC_DIR}/lib/petsc/conf/variables
  INCLUDE += $(PETSC_CCPPFLAGS)
  FINCLUDE += $(PETSC_FCPPFLAGS)
  LIBDIRS += $(PETSC_C_SH_LIB_PATH)
  LIBS += $(PETSC_LIB)
  EXAMPLES_C += ex_eigs_petsc ex_svds_petsc ex_eigs_mpi
  EXAMPLES_F += ex_eigs_petscf77 ex_eigs_petscf77ptr ex_svds_petscf77 ex_svds_petscf77ptr
  MPIRUN ?= mpirun -np 4
else ifeq ($(USE_MPI),yes)
  EXAMPLES_C +=  ex_eigs_mpi
  MPIRUN ?= mpirun -np 4
  CC ?= mpicc
  CLDR ?= mpicc
endif

ifeq ($(PRIMME_WITH_MAGMA),yes)
  EXAMPLES_C += ex_eigs_dmagma ex_svds_dmagma
  INCLUDE += -I$(MAGMADIR)/include -I$(MAGMADIR)/sparse/include -I$(CUDADIR)/include
endif
ifeq ($(PRIMME_WITH_CUDA),yes)
  EXAMPLES_C += ex_eigs_dcublas ex_eigs_zcublas ex_eigs_zcublas0
endif
ifeq ($(PRIMME_WITH_HIPBLAS),yes)
  EXAMPLES_C += ex_eigs_dhipblas ex_eigs_zhipblas ex_eigs_zhipblas0
endif


$(EXAMPLES_C): %: %.o
	$(CLDR) -o $@ $@.o $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS)

$(EXAMPLES_CXX): %: %.o
	$(CXX) -o $@ $@.o $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS)

$(EXAMPLES_F): % : %.o
	$(FLDR) -o $@ $@.o $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c $< -o $@

%.o: %.cxx
	$(CXX) $(CXXFLAGS) $(DEFINES) $(INCLUDE) -c $< -o $@

.f.o .F.o:
	$(F77) $(FFLAGS) $(FINCLUDE) -c $< -o $@

%.o: %.f90
	$(FC) $(FFLAGS) $(FINCLUDE) -c $< -o $@

all: examples

examples: examples_C examples_CXX examples_F
examples_C: $(EXAMPLES_C)
examples_CXX: $(EXAMPLES_CXX)
examples_F: $(EXAMPLES_F)

test_examples: test_examples_C test_examples_F test_examples_CXX
test_examples_C test_examples_F test_examples_CXX: test_examples% : examples%
	@ok="0";for ex in $(EXAMPLES$*); do \
		echo "=========== Executing ./$$ex"; \
	        case $$ex in \
		*petsc*) ${MPIRUN} ./$$ex < /dev/null || ok="1";; \
		*)                 ./$$ex || ok="1"; \
	        esac; \
	done < /dev/null > tests.log 2>&1; \
	if test $$ok -eq 0 ; then \
		echo "All tests passed!"; \
	else\
		cat tests.log;\
		echo "Some tests fail. Please consider to send us the file";\
		echo "examples/tests.log if the software doesn't work as expected.";\
		exit 1;\
	fi

clean:
	@rm -f *.o tests.log

veryclean: clean
	@rm -f $(EXAMPLES_C) $(EXAMPLES_F) $(EXAMPLES_CXX)

ifeq ($(UNAME), FreeBSD)
  M4 ?= m4 -g
else
  M4 ?= m4
endif

ex_mine.c: ex1_c.m4
	$(M4) ex1_c.m4 > $@

ex_eigs_dseq.c: ex_c.m4
	$(M4) ex_c.m4 > $@

ex_eigs_zseq.c: ex_c.m4
	$(M4) -D USE_COMPLEX -D ADVANCED ex_c.m4 > $@

ex_eigs_zseqxx.cxx: ex_c.m4
	$(M4) -D USE_COMPLEX -D USE_COMPLEX_CXX -D ADVANCED ex_c.m4 > $@

ex_eigs_petsc.c: ex_c.m4
	$(M4) -D USE_PETSC ex_c.m4 > $@

ex_eigs_dseqf77.f: ex_f77.m4
	$(M4) ex_f77.m4 > $@

ex_eigs_zseqf77.f: ex_f77.m4
	$(M4) -D USE_COMPLEX ex_f77.m4 > $@

ex_eigs_petscf77.F: ex_f77.m4
	$(M4) -D USE_PETSC ex_f77.m4 > $@

ex_eigs_petscf77ptr.F: ex_f77.m4
	$(M4) -D USE_PETSC -D USE_POINTER ex_f77.m4 > $@

ex_svds_dseq.c: exsvd_c.m4
	$(M4) exsvd_c.m4 > $@

ex_svds_zseq.c: exsvd_c.m4
	$(M4) -D USE_COMPLEX -D ADVANCED_HYBRID exsvd_c.m4 > $@

ex_svds_zseqxx.cxx: exsvd_c.m4
	$(M4) -D USE_COMPLEX -D USE_COMPLEX_CXX -D ADVANCED_HYBRID exsvd_c.m4 > $@

ex_svds_petsc.c: exsvd_c.m4
	$(M4) -D USE_PETSC exsvd_c.m4 > $@

ex_svds_dseqf77.f: exsvd_f77.m4
	$(M4) exsvd_f77.m4 > $@

ex_svds_zseqf77.f: exsvd_f77.m4
	$(M4) -D USE_COMPLEX exsvd_f77.m4 > $@

ex_svds_petscf77.F: exsvd_f77.m4
	$(M4) -D USE_PETSC exsvd_f77.m4 > $@

ex_svds_petscf77ptr.F: exsvd_f77.m4
	$(M4) -D USE_PETSC -D USE_POINTER exsvd_f77.m4 > $@

.PHONY: clean veryclean all examples
