#-----------------------------------------------------------------
# make 
#   lib	      makes the libprimme.a library
#   solib     makes the libprimme.so library
#   clean     removes all *.o files
#   test      build and execute simple examples
#-----------------------------------------------------------------
include Make_flags

.PHONY: lib clean test all_tests check_style

#------------------------ Libraries ------------------------------
# Making the PRIMME library
# Includes float, double and complex counterparts
lib:
	@make -C src ../lib/$(LIBRARY)

solib:
	@make -C src ../lib/$(SOLIBRARY)

clean: 
	@make -C src clean

test:
	@\
	echo "------------------------------------------------"; \
	echo " Test C examples                                "; \
	echo "------------------------------------------------"; \
	make -C examples veryclean test_examples_C USE_PETSC=no

all_tests:
	@make -C examples veryclean test_examples;\
	make -C tests veryclean all_tests

check_style:
	( grep '	' -R . --include='*.[chfmF]' && echo "Please don't use tabs!" ) || true

tags:
	@make -C src ../tags

.NOTPARALLEL:
