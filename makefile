#-----------------------------------------------------------------
# make 
#   lib	      makes the libprimme.a library
#   solib     makes the libprimme.so library
#   matlab    make libprimme.a compatible with MATLAB and the
#             module for MATLAB
#   octave    make libprimme.a and the Octave module
#   python    make libprimme.a and the python module
#   clean     removes all *.o files
#   test      build and execute simple examples
#-----------------------------------------------------------------
include Make_flags

.PHONY: lib clean test all_tests check_style matlab octave python tags

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

matlab:
	@make clean lib CFLAGS="${CFLAGS} -DPRIMME_BLASINT_SIZE=64"
	@make -C Matlab matlab

octave:
	@make clean lib
	@make -C Matlab octave

python:
	@make clean lib
	@make -C Python

check_style:
	( grep '	' -R . --include='*.[chfmF]' && echo "Please don't use tabs!" ) || true

tags:
	@make -C src ../tags

.NOTPARALLEL:
