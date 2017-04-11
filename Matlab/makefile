#--------------------------------------------------------------
# Makefile for building PRIMME MEX
#-------------------------------------------------------------

# Name of MEX source file
MEXNAME = primme_mex

# PRIMME includes and linking
INCLUDE = -I../include
LIBS = -L../lib -lprimme -lm

#------------------------------------------------------------
# MATLAB configuration

# MATLAB MEX compiler
MEX ?= mex
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
	$(MEX) $(MEXFLAGS) $(LIBS) $(MEXLIBS) $(INCLUDE) -output $(MEXNAME) $(MEXNAME).cpp

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
