#-----------------------------------------------------------------
# make 
#   lib	      makes the libprimme.a library (both double and complex)
#   libd      makes the libdprimme.a library (only double)
#   libz      makes the libzprimme.a library (only complex)
#   seqs      makes sequential test programs in DTEST and ZTEST
#   pars      makes parallel test programs (only in DTEST currently)
#   depends   builds the dependencies needed for making the test programs
#   		ddepends_seq  (double precision sequential)
#   		ddepends_par  (double precision parallel)
#   		zdepends_seq  (Complex sequential)
#   clean     removes all .o a.out and core files
#   all       makes main library and all test programs
#-----------------------------------------------------------------
include Make_flags

# Build the libprimme library, as well as all TESTS
all: lib primme_double primme_complex

# The following are sequential executables of the test programs
seqs: primme_double primme_complex seqf77_dprimme seqf77_zprimme

# The following builds the parallel executables
pars: primme_double primme_complex

.PHONY: lib libd libz clean backup

#------------------------ Libraries ------------------------------
lib:
	@(\
	cd $(SRC) ;\
	echo "-----------------------------------"; \
	echo "     Making the PRIMME library     "; \
	echo "  Includes both double and complex "; \
	echo "-----------------------------------"; \
	make lib;\
	)

libd:
	@(\
	cd $(SRC) ;\
	echo "---------------------------------------"; \
	echo " Making the double-only PRIMME library "; \
	echo "---------------------------------------"; \
	make libd;\
	)

libz:
	@(\
	cd $(SRC) ;\
	echo "----------------------------------------"; \
	echo " Making the complex-only PRIMME library "; \
	echo "----------------------------------------"; \
	make libz;\
	)
#------------------- DOUBLE TEST PROGRAMS ------------------------
primme_double:
	@(\
	cd $(TESTDIR) ;\
	echo "----------------------------------------------"; \
	echo " Making the executable for the double C       "; \
	echo "----------------------------------------------"; \
	make primme_double; \
	)

seqf77_dprimme:
	@(\
	cd $(TESTDIR) ;\
	echo "------------------------------------------------"; \
	echo " Making the executable for double sequential f77 "; \
	echo "------------------------------------------------"; \
	make seqf77_dprimme; \
	)

#------------------- COMPLEX TEST PROGRAMS ----------------------
primme_complex:
	@(\
	cd $(TESTDIR) ;\
	echo "----------------------------------------------"; \
	echo " Making the executable for the complex C      "; \
	echo "----------------------------------------------"; \
	make primme_doublecomplex; \
	)

seqf77_zprimme:
	@(\
	cd $(TESTDIR) ;\
	echo "------------------------------------------------"; \
	echo " Making executable for complex sequential f77 "; \
	echo "------------------------------------------------"; \
	make seqf77_zprimme; \
	)

clean: 
	@(\
	echo "--------------------------------------------------"; \
	echo " Cleaning .o, a.out, cores, from all directories "; \
	echo "--------------------------------------------------"; \
	cd $(SRC) ;\
	make -f Makefile clean;\
	rm -f libprimme.a;\
	echo " From Test directories";\
	cd $(TESTDIR);\
	make clean;\
	echo "--------------------------------------------------"; \
	)

test:
	@(\
	cd $(TESTDIR) ;\
	echo "------------------------------------------------"; \
	echo " Test double sequential C                       "; \
	echo "------------------------------------------------"; \
	make test_double; \
	)


check_style:
	( grep $$'\t' -R . --include='*.[chf]' && echo "Please don't use tabs!" ) || true

backup: 
	@(\
	cd $(TOP)/../ ;\
	tar cvf back_$(shell date "+%m_%d_%y_%H_%M").tar PRIMME;\
	gzip back_*.tar;\
	)

distribution:
	@(\
	cd $(TOP)/../ ;\
	tar cvf primme_v1.11.tar \
	PRIMME/PRIMMESRC PRIMME/DTEST  PRIMME/ZTEST  \
	PRIMME/readme.txt PRIMME/COPYING.txt \
	PRIMME/readme.html PRIMME/doc.pdf PRIMME/primmestyle.css \
	PRIMME/Abstract_stathopoulos.pdf PRIMME/Paper_stathopoulos.pdf \
	PRIMME/Make_flags PRIMME/Link_flags PRIMME/makefile;\
	gzip primme_v1.11.tar;\
	)

