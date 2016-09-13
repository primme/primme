#-----------------------------------------------------------------
# make 
#   lib	      makes the libprimme.a library (both double and complex)
#   libd      makes the libdprimme.a library (only double)
#   libz      makes the libzprimme.a library (only complex)
#   clean     removes all .o a.out and core files
#   test      build and execute a simple example
#-----------------------------------------------------------------
include Make_flags

.PHONY: lib clean test all_tests check_style

#------------------------ Libraries ------------------------------
# Making the PRIMME library
# Includes float, double and complex counterparts
lib:
	@make -C src ../lib/$(LIBRARY)

clean: 
	@make -C src clean

test:
	@\
	echo "------------------------------------------------"; \
	echo " Test C examples                                "; \
	echo "------------------------------------------------"; \
	make -C examples test_examples_C USE_PETSC=no

all_test:
	@make -C examples test_examples;\
	make -C tests all_tests

check_style:
	( grep '	' -R . --include='*.[chfmF]' && echo "Please don't use tabs!" ) || true

tags:
	@make -C src ../tags
