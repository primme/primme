#--------------------------------------------------------------
# Makefile for building PRIMME MEX
#-------------------------------------------------------------

#------------------------------------------------------------
# Name of MEX source file
#------------------------------------------------------------
MEXNAME = PRIMME_mex

#------------------------------------------------------------
# MATLAB MEX compiler
#------------------------------------------------------------
MEX ?= mex
# Linux
# MEX = /usr/local/matlab/bin/mex
#------------------------------------------------------------
# Mac OS X
# MEX = /Applications/MATLAB_R2014a.app/bin/mex
# MEX = /Applications/MATLAB_R2015a.app/bin/mex
#------------------------------------------------------------

#------------------------------------------------------------
# MATLAB binary (only for testing)
#------------------------------------------------------------
MATLAB ?= matlab

# Modify -O as -g to enter debug model
MEXFLAGS += -O -largeArrayDims

#------------------------------------------------------------
# Link to PRIMME and MATLAB BLAS/LAPACK
#------------------------------------------------------------
# Set TOP to the path where PRIMME is installed
TOP = ../

LIBS = -L../lib -lprimme -lm  -lmwlapack -lmwblas
INCLUDE = -I../include

all:
	$(MEX) $(MEXFLAGS) $(LIBS) $(INCLUDE) -output $(MEXNAME) $(MEXNAME).c

test:
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

clean:
	rm -f $(MEXNAME).[^c]*

.PHONY: all test clean
