#--------------------------------------------------------------
#	Makefile for building PRIMME MEX
#-------------------------------------------------------------

#------------------------------------------------------------
#	Name of MEX source file
#------------------------------------------------------------
  MEXNAME = PRIMME_mex

#------------------------------------------------------------
#   Modify suffix name of MEX binary file for different OS
#------------------------------------------------------------
# LINUX SUSE 64 BIT
  MEXSUFFIX = mexa64
#------------------------------------------------------------
# MAC OXS 64 BIT
# MEXSUFFIX = mexmaci64
#------------------------------------------------------------

#------------------------------------------------------------
#    Modify home directory of MATLAB and show path of mex
#------------------------------------------------------------
# LINUX SUSE 64 BIT
  MATLABHOME = /usr/local/matlab
  MEX = $(MATLABHOME)/bin/mex
#------------------------------------------------------------
# MAC OXS 64 BIT
# MATLABHOME = /Applications/MATLAB_R2014a.app
# MEX = $(MATLABHOME)/bin/mex
#------------------------------------------------------------

# Modify compiler users want to use, or use default ones.
# For MATLAB 2012, the recommended compilers are gcc 4.4.x
# but any older compilers should also been supported.
# The details are refered to MATLAB supported compilers.
CXX = /usr/bin/gcc
#CXX = /home/scratch/lfwu/usr/local/gcc4.4.6/bin/gcc 
#FC = gfortran
#CC = gcc
CFLAGS= -fPIC -m64

# Modify -O as -g to entre debug model
MEXFLAGS = -cxx CC='$(CXX)' CXX='$(CXX)' LD='$(CXX)'
MEXFLAGS += -v -O -largeArrayDims

#------------------------------------------------------------
#	Link to PRIMME liabray, math library, MATLAB BLAS
#       and LAPACK libraries  
#------------------------------------------------------------
# Modify this to show the path of the installed PRIMME directory
#TOP = /Users/teddy/Dropbox/Lingfei/MyWork/LWERAS/PRIMME_MEX/PRIMME
TOP = ../
 
LIBS = -lprimme -lm  -lmwlapack -lmwblas
LIB = -L$(TOP)

INCLUDE = -I$(MATLABHOME)/extern/include/ -I$(TOP)/PRIMMESRC/COMMONSRC

$(MEXNAME).$(MEXSUFFIX): $(MEXNAME).o
	$(MEX) $(MEXFLAGS) $(LIB) $(LIBS) -output $(MEXNAME) $^

$(MEXNAME).o: $(MEXNAME).c
	$(CC) $(CFLAGS) $(INCLUDE) -c $^

clean:
	rm -f *.o; rm -f *.mex*
