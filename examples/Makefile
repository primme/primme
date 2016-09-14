# include the environment, compilation, and linking flags

include ../Make_flags
include ../Link_flags

#------------------------------------------------


override INCLUDE += -I../include
override FINCLUDE += -I../include
LIBDIRS += -L../lib

EXAMPLES_C = ex_eigs_dseq ex_eigs_zseq ex_svds_dseq ex_svds_zseq
EXAMPLES_CXX = ex_eigs_zseqxx ex_svds_zseqxx
EXAMPLES_F77 = ex_eigs_dseqf77 ex_eigs_zseqf77 ex_svds_dseqf77 ex_svds_zseqf77

USE_PETSC     ?= $(if $(findstring undefined,$(origin PETSC_DIR)),no,yes)

ifeq ($(USE_PETSC), yes)
  include ${PETSC_DIR}/lib/petsc/conf/variables
  override INCLUDE += $(PETSC_CCPPFLAGS)
  override FINCLUDE += $(PETSC_FCPPFLAGS)
  LIBDIRS += $(PETSC_C_SH_LIB_PATH)
  LIBS += $(PETSC_LIB)
  EXAMPLES_C += ex_eigs_petsc ex_svds_petsc
  EXAMPLES_F77 += ex_eigs_petscf77 ex_eigs_petscf77ptr ex_svds_petscf77 ex_svds_petscf77ptr
endif

$(EXAMPLES_C): %: %.o
	$(CLDR) -o $@ $@.o $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

$(EXAMPLES_CXX): %: %.o
	$(CXX) -o $@ $@.o $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

$(EXAMPLES_F77): % : %.o
	$(FLDR) -o $@ $@.o $(LIBDIRS) $(INCLUDE) $(LIBS) $(LDFLAGS) 

.c.o:
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c $< -o $@

%.o: %.cxx
	$(CXX) $(CXXFLAGS) $(DEFINES) $(INCLUDE) -c $< -o $@

.f.o .F.o:
	$(F77) $(FFLAGS) $(FINCLUDE) $(COMMON_INCLUDE) -c $< -o $@

all: examples

examples: examples_C examples_CXX examples_F77
examples_C: $(EXAMPLES_C)
examples_CXX: $(EXAMPLES_CXX)
examples_F77: $(EXAMPLES_F77)

test_examples: test_examples_C test_examples_F77 test_examples_CXX
test_examples_C test_examples_F77 test_examples_CXX: test_examples% : examples%
	@ok="0";for ex in $(EXAMPLES$*); do \
		echo "=========== Executing ./$$ex"; \
		./$$ex || ok="1"; \
	done > tests.log 2>&1; \
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
	@rm -f $(EXAMPLES_C) $(EXAMPLES_F77) $(EXAMPLES_CXX)

ex_eigs_dseq.c: ex_c.m4
	m4 ex_c.m4 > $@

ex_eigs_zseq.c: ex_c.m4
	m4 -D USE_COMPLEX -D ADVANCED ex_c.m4 > $@

ex_eigs_zseqxx.cxx: ex_c.m4
	m4 -D USE_COMPLEX -D USE_COMPLEX_CXX -D ADVANCED ex_c.m4 > $@

ex_eigs_petsc.c: ex_c.m4
	m4 -D USE_PETSC ex_c.m4 > $@

ex_eigs_dseqf77.f: ex_f77.m4
	m4 ex_f77.m4 > $@

ex_eigs_zseqf77.f: ex_f77.m4
	m4 -D USE_COMPLEX ex_f77.m4 > $@

ex_eigs_petscf77.F: ex_f77.m4
	m4 -D USE_PETSC ex_f77.m4 > $@

ex_eigs_petscf77ptr.F: ex_f77.m4
	m4 -D USE_PETSC -D USE_POINTER ex_f77.m4 > $@

ex_svds_dseq.c: exsvd_c.m4
	m4 exsvd_c.m4 > $@

ex_svds_zseq.c: exsvd_c.m4
	m4 -D USE_COMPLEX -D ADVANCED_HYBRID exsvd_c.m4 > $@

ex_svds_zseqxx.cxx: exsvd_c.m4
	m4 -D USE_COMPLEX -D USE_COMPLEX_CXX -D ADVANCED_HYBRID exsvd_c.m4 > $@

ex_svds_petsc.c: exsvd_c.m4
	m4 -D USE_PETSC exsvd_c.m4 > $@

ex_svds_dseqf77.f: exsvd_f77.m4
	m4 exsvd_f77.m4 > $@

ex_svds_zseqf77.f: exsvd_f77.m4
	m4 -D USE_COMPLEX exsvd_f77.m4 > $@

ex_svds_petscf77.F: exsvd_f77.m4
	m4 -D USE_PETSC exsvd_f77.m4 > $@

ex_svds_petscf77ptr.F: exsvd_f77.m4
	m4 -D USE_PETSC -D USE_POINTER exsvd_f77.m4 > $@

.PHONY: clean veryclean all examples
