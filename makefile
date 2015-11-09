#-----------------------------------------------------------------
# make 
#   lib	      makes the libprimme.a library (both double and complex)
#   libd      makes the libdprimme.a library (only double)
#   libz      makes the libzprimme.a library (only complex)
#   clean     removes all .o a.out and core files
#   test      build and execute a simple example
#-----------------------------------------------------------------
include Make_flags

.PHONY: lib libd libz clean backup test

#------------------------ Libraries ------------------------------
lib:
	@(\
	cd PRIMMESRC;\
	echo "-----------------------------------"; \
	echo "     Making the PRIMME library     "; \
	echo "  Includes both double and complex "; \
	echo "-----------------------------------"; \
	make lib;\
	)

libd:
	@(\
	cd PRIMMESRC;\
	echo "---------------------------------------"; \
	echo " Making the double-only PRIMME library "; \
	echo "---------------------------------------"; \
	make libd;\
	)

libz:
	@(\
	cd PRIMMESRC;\
	echo "----------------------------------------"; \
	echo " Making the complex-only PRIMME library "; \
	echo "----------------------------------------"; \
	make libz;\
	)
clean: 
	@(\
	echo "--------------------------------------------------"; \
	echo " Cleaning .o, a.out, cores, from all directories "; \
	echo "--------------------------------------------------"; \
	make -C PRIMMESRC clean;\
	echo " From Test directories";\
	make -C TEST clean;\
	echo "--------------------------------------------------"; \
	)

test:
	@(\
	cd TEST ;\
	echo "------------------------------------------------"; \
	echo " Test double sequential C                       "; \
	echo "------------------------------------------------"; \
	make test_double; \
	)


check_style:
	( grep '	' -R . --include='*.[chfmF]' && echo "Please don't use tabs!" ) || true

backup: 
	@(\
	cd .. ;\
	tar cvf back_$(shell date "+%m_%d_%y_%H_%M").tar PRIMME;\
	gzip back_*.tar;\
	)

distribution:
	@(\
	cd .. ;\
	tar cvf primme_v1.11.tar \
	PRIMME/PRIMMESRC PRIMME/DTEST  PRIMME/ZTEST  \
	PRIMME/readme.txt PRIMME/COPYING.txt \
	PRIMME/readme.html PRIMME/doc.pdf PRIMME/primmestyle.css \
	PRIMME/Abstract_stathopoulos.pdf PRIMME/Paper_stathopoulos.pdf \
	PRIMME/Make_flags PRIMME/Link_flags PRIMME/makefile;\
	gzip primme_v1.11.tar;\
	)

