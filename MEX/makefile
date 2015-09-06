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

LIBS = -L$(TOP) -lprimme -lm  -lmwlapack -lmwblas
INCLUDE = -I$(TOP)/PRIMMESRC/COMMONSRC

all:
	$(MEX) $(MEXFLAGS) $(LIBS) $(INCLUDE) -output $(MEXNAME) $(MEXNAME).c

test:
	$(MATLAB) -nodesktop -nodisplay -r "E=primme_eigs(diag(1:50),6);assert(norm(E-(50:-1:45)')<1e-9);exit" &> test.log < /dev/null
	@if `grep -q 'Assertion failed' test.log` ; then \
	   cat test.log;\
	       echo;\
	       echo "Something went wrong. Please consider to send us the file 'test.log'";\
	       echo "if the software doesn't work as expected.";\
	       exit 1;\
	else\
	  echo "Test passed!";\
	fi

clean:
	rm -f $(MEXNAME).[^c]*

.PHONY: all test clean
