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
all: lib depends seqs pars 

# The following are sequential executables of the test programs
seqs: seq_dprimme seq_zprimme seqf77_dprimme seqf77_zprimme

# The following builds the parallel executables
pars: par_dprimme 

# Build the dependencies for .c and .h files in the test directories
depends: ddepends_seq zdepends_seq ddepends_par 

.PHONY: lib libd libz clean backup ddepends_seq ddepends_par zdepends_seq 

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
par_dprimme:
	@(\
	cd $(DTESTDIR) ;\
	echo "----------------------------------------------"; \
	echo " Making the executable for the double parallel C"; \
	echo "----------------------------------------------"; \
	make -f Makefile_par par_dprimme; \
	)

seq_dprimme:
	@(\
	cd $(DTESTDIR) ;\
	echo "------------------------------------------------"; \
	echo " Making the executable for double sequential C"; \
	echo "------------------------------------------------"; \
	make -f Makefile_seq seq_dprimme; \
	)

seqf77_dprimme:
	@(\
	cd $(DTESTDIR) ;\
	echo "------------------------------------------------"; \
	echo " Making the executable for double sequential f77 "; \
	echo "------------------------------------------------"; \
	make -f Makefile_f77seq seqf77_dprimme; \
	)

#------------------- COMMPLEX TEST PROGRAMS ----------------------
seq_zprimme:
	@(\
	cd $(ZTESTDIR) ;\
	echo "------------------------------------------------"; \
	echo " Making the executable for complex sequential C"; \
	echo "------------------------------------------------"; \
	make -f Makefile_seq seq_zprimme; \
	)

seqf77_zprimme:
	@(\
	cd $(ZTESTDIR) ;\
	echo "------------------------------------------------"; \
	echo " Making executable for complex sequential f77 "; \
	echo "------------------------------------------------"; \
	make -f Makefile_f77seq seqf77_zprimme; \
	)

#------------------------- Test dependencies -----------------------
ddepends_seq: 
	@(\
	cd $(DTESTDIR);\
	rm -f ddependencies_seq;\
	echo "----------------------------------"; \
	echo " Creating double seq dependencies "; \
	echo "----------------------------------"; \
	make -f Makefile_seq ddependencies_seq;\
	)

ddepends_par: 
	@(\
	cd $(DTESTDIR);\
	rm -f ddependencies_par;\
	echo "----------------------------------"; \
	echo " Creating double par dependencies "; \
	echo "----------------------------------"; \
	make -f Makefile_par ddependencies_par;\
	)

zdepends_seq: 
	@(\
	cd $(ZTESTDIR);\
	rm -f zdependencies_seq;\
	echo "-----------------------------------"; \
	echo " Creating complex seq dependencies "; \
	echo "-----------------------------------"; \
	make -f Makefile_seq zdependencies_seq;\
	)

clean: 
	@(\
	echo "--------------------------------------------------"; \
	echo " Cleaning .o, a.out, cores, from all directories "; \
	echo "--------------------------------------------------"; \
	cd $(SRC) ;\
	make -f Makefile clean;\
	cd $(DTESTDIR) ; \
	make -f Makefile_par clean;\
	make -f Makefile_seq clean;\
	make -f Makefile_f77seq clean;\
	cd $(ZTESTDIR) ; \
	make -f Makefile_seq clean;\
	make -f Makefile_f77seq clean;\
	)

backup: 
	@(\
	cd $(TOP)/../ ;\
	tar cvf back_$(shell date "+%m_%d_%y_%H_%M").tar PRIMME;\
	gzip back_*.tar;\
	)

distribution:
	@(\
	cd $(TOP)/../ ;\
	tar cvf primme_v1.1.tar \
	PRIMME/PRIMMESRC PRIMME/DTEST  PRIMME/ZTEST  \
	PRIMME/readme.txt PRIMME/COPYING.txt \
	PRIMME/Make_flags PRIMME/Link_flags PRIMME/makefile;\
	gzip primme_v1.1.tar;\
	)

