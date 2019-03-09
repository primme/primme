#--------------------------------------------------------------
# Makefile for building PRIMME MEX
#-------------------------------------------------------------

include ../Make_flags

# Name of MEX source file
MEXNAME = primme_mex

# PRIMME includes and linking
INCLUDE = -I../include
LIBS = -L../lib -lprimme -lm

ifeq ($(PRIMME_WITH_MAGMA),yes)
INCLUDE += -I$(MAGMADIR)/include
LIBS += -L$(MAGMADIR)/lib -lmagma  \
        -L$(CUDADIR)/lib64 -lcublas -lcudart
endif

#------------------------------------------------------------
# MATLAB configuration

# MATLAB MEX compiler
MEX ?= mex
ifneq ($(findstring MEX,$(shell $(MEX) -version 2>&1 || true)),MEX)
   MEX_ERROR = $(error Invalid MATLAB`s mex, $(MEX). Set MEX=/path/to/Matlab/bin/mex)
endif

# Linux
# MEX = /usr/local/matlab/bin/mex
#
# Mac OS X
# MEX = /Applications/MATLAB_R2014a.app/bin/mex
# MEX = /Applications/MATLAB_R2015a.app/bin/mex

# Extra compiler flags 
MEXFLAGS += -O -largeArrayDims
# For debug:
# MEXFLAGS += -g -largeArrayDims
MW_NVCC_PATH ?= $(CUDADIR)/bin

# Extra libraries to link
MEXLIBS += -lmwlapack -lmwblas

# MATLAB executable app (only for testing)
MATLAB ?= matlab

#------------------------------------------------------------
# Octave configuration

# Octave compiler
MKOCTFILE ?= mkoctfile

# Extra compiler flags 
# OCTFLAGS += -g

# Octave binary (only for testing)
OCTAVE ?= octave --no-gui -H -q

all: matlab

matlab:
	$(MEX_ERROR)
	$(MEX) $(MEXFLAGS) $(LIBS) $(MEXLIBS) $(INCLUDE) -output $(MEXNAME) $(MEXNAME).cpp

matlab-cuda:
	$(MEX_ERROR)
	$(MATLAB) -nodesktop -nodisplay -nojvm -r "setenv('MW_NVCC_PATH','$(MW_NVCC_PATH)'); try;mexcuda $(MEXFLAGS) -DUSE_GPUARRAY $(LIBS) $(MEXLIBS) $(INCLUDE) -output $(MEXNAME) $(MEXNAME).cpp; catch me,disp(getReport(me));end, exit" < /dev/null

octave:
	$(MKOCTFILE) --mex $(OCTFLAGS) $(LIBS) $(INCLUDE) $(MEXNAME).cpp

test: test_matlab

test_matlab:
	@trap 'echo MATLAB was killed' 9;\
	rm -f test.log;\
	nohup $(MATLAB) -nodesktop -nodisplay -r "try,test_primme,catch me,disp(getReport(me));end;exit" -logfile test.log > /dev/null; \
	if [ -s test.log ] && grep -q 'Success' test.log; then \
		echo "Test passed!" ;\
	else\
		cat test.log || true ;\
		echo ;\
		echo "Something went wrong or the MATLAB output was not captured.";\
		echo "Please make sure PRIMME is compiled with -DPRIMME_BLASINT_SIZE=64.";\
		echo "Consider to send us the file 'test.log'" ;\
		echo "if the software doesn't work as expected." ;\
		exit 1;\
	fi ;\

test_octave:
	@trap 'echo Octave was killed' 9;\
	rm -f test.log;\
	nohup $(OCTAVE) test_primme.m > test.log; \
	if [ -s test.log ] && grep -q 'Success' test.log; then \
		echo "Test passed!" ;\
	else\
		cat test.log || true ;\
		echo ;\
		echo "Something went wrong or the Octave output was not captured.";\
		echo "Probably PRIMME shouldn't be compiled with -DPRIMME_BLASINT_SIZE=64.";\
		echo "Consider to send us the file 'test.log'" ;\
		echo "if the software doesn't work as expected." ;\
		exit 1;\
	fi ;\

clean:
	rm -f $(MEXNAME).[^c]*

.PHONY: all matlab octave test test_matlab test_octave clean 
