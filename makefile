#-----------------------------------------------------------------
# make 
#   lib	      makes the libprimme.a library (both double and complex)
#   libd      makes the libdprimme.a library (only double)
#   libz      makes the libzprimme.a library (only complex)
#   clean     removes all .o a.out and core files
#   test      build and execute a simple example
#-----------------------------------------------------------------
include Make_flags

.PHONY: lib libd libz clean test

#------------------------ Libraries ------------------------------
# Making the PRIMME library
# Includes float, double and complex counterparts
lib:
	@make -C src ../lib/$(LIBRARY)

clean: 
	@make -C src clean;\
	make -C TEST clean

test:
	@\
	echo "------------------------------------------------"; \
	echo " Test double sequential C                       "; \
	echo "------------------------------------------------"; \
	make -C TEST test_double USE_PETSC=no

all_test:
	@make -C TEST all_tests USE_PETSC=no

check_style:
	( grep '	' -R . --include='*.[chfmF]' && echo "Please don't use tabs!" ) || true

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

