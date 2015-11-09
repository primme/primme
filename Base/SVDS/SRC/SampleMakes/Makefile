#include Makefile.opts

COMPILE_C   = $(CC) $(CFLAGS) $(INCLUDEMPI) $(INCLUDESCALAPACK) -c
COMPILE_F77 = g77 -c -O3


PSGESVD = test_psgesvd
PDGESVD = test_pdgesvd
PCGESVD = test_pcgesvd
PZGESVD = test_pzgesvd

all: test_psgesvd test_pdgesvd test_pcgesvd test_pzgesvd

test_psgesvd.c: test_pgesvd.base.c
	sed -f arithmetics/cs.sed $< > $@

test_pcgesvd.c: test_pgesvd.base.c
	sed -f arithmetics/cc.sed $< > $@

test_pdgesvd.c: test_pgesvd.base.c
	sed -f arithmetics/cd.sed $< > $@

test_pzgesvd.c: test_pgesvd.base.c
	sed -f arithmetics/cz.sed $< > $@

%.o: %.c
	$(COMPILE_C)   $*.c -o $@
SOURCES_C = test_psgesvd.c test_pdgesvd.c test_pcgesvd.c test_pzgesvd.c 
OBJECTS_C = $(SOURCES_C:%.c=%.o)

psgesvd.o: src-scalapack/psgesvd.f
	$(COMPILE_F77) src-scalapack/psgesvd.f -o psgesvd.o

pcgesvd.o: src-scalapack/pcgesvd.f
	$(COMPILE_F77) src-scalapack/pcgesvd.f -o pcgesvd.o

pdgesvd.o: src-scalapack/pdgesvd.f
	$(COMPILE_F77) src-scalapack/pdgesvd.f -o pdgesvd.o

pzgesvd.o: src-scalapack/pzgesvd.f
	$(COMPILE_F77) src-scalapack/pzgesvd.f -o pzgesvd.o

OBJECTS_F77 = psgesvd.o pdgesvd.o pcgesvd.o pzgesvd.o

LINK = $(F77) $(LDFLAGS)
LIBS = $(LIBSCALAPACK) $(LIBBLACS) $(LIBBLAS) $(LIBMPI)

test_psgesvd: test_psgesvd.o psgesvd.o
	$(LINK) test_psgesvd.o psgesvd.o $(LIBS) -o $@

test_pdgesvd: test_pdgesvd.o pdgesvd.o
	$(LINK) test_pdgesvd.o pdgesvd.o $(LIBS) -o $@

test_pcgesvd: test_pcgesvd.o pcgesvd.o
	$(LINK) test_pcgesvd.o pcgesvd.o $(LIBS) -o $@

test_pzgesvd: test_pzgesvd.o pzgesvd.o
	$(LINK) test_pzgesvd.o pzgesvd.o $(LIBS) -o $@

testing: test_psgesvd test_pdgesvd test_pcgesvd test_pzgesvd
	cp svd-small.dat svd.dat
	$(MPIRUN) -np 8 ./test_psgesvd -f >  svd.out
	$(MPIRUN) -np 8 ./test_pdgesvd -f >> svd.out
	$(MPIRUN) -np 8 ./test_pcgesvd -f >> svd.out
	$(MPIRUN) -np 8 ./test_pzgesvd -f >> svd.out

clean:
	rm -f  $(PSGESVD) $(PDGESVD) $(PCGESVD) $(PZGESVD) $(OBJECTS_F77) $(OBJECTS_C) $(SOURCES_C)

