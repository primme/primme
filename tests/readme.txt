
This directory has the test programs for PRIMME.
It should have at least:

- Makefile             makefile to build the programs.
- driver.c             eigenvalue driver;
                       reads MTX (http://math.nist.gov/MatrixMarket/mmio-c.html)
                       and PETSc (http://www.mcs.anl.gov/petsc/) matrices,
                       and can be configured to run
                       in parallel or sequential, double or double complex,
                       w/ or w/o preconditioning, and with a variety of
                       preconditioners. For simpler examples see below.
- driversvds.c         singular value driver;
                       similar features than driver.c
- COMMON/              with source used by driver.c and driversvds.c.
    csr.h, csr.c       routines for matrices CSR
    mmio.h, mmio.c     MatrixMarket IO routines.
    native.h, mat.c    wrapper for CSR matrix and sequential ILUT.
    num.h              constants
    parasailsw.h, .c   wrapper for ParaSails matrix and preconditioner.
    petscw.h, .c       wrapper for PETSc matrices and preconditioners.
    shared_utils.h, .c IO routines for primme_params and driver options.
    ssrcsr.c           routine to convert from Sym Sparse Row to CSR (from Sparskit).
    amux.f             routine for CSR matrix-vector product (from Sparskit).
    ilut.f             routine for sequential ILUT (from Sparskit).
    zamux.f            routine for complex CSR matrix-vector product (from Sparskit).
    zilut.f            routine for complex sequential ILUT (from Sparskit).
- DriverConf           example of driver configuration file used by the driver.
- MinConf, LeanConf,
  FullConf             examples of PRIMME configuration file used by the driver.
- LUNDA.mtx            matrix used for testing and in DriverConf as an example.
- tests/               configuration files for testing purpose.

The Makefile can perform the next actions:

make primme_double          build eigenvalue driver in double.
make primme_doublecomplex     "     "      "        in complex double.
make primmesvds_double      build singular value driver in double.
make primmesvds_doublecomplex     "     "      "            in complex double.
make all_tests              test all configurations in "tests"
make clean                  remove object files.
make veryclean              remove object and program files.


* Compile driver with PETSc

First set PETSC_DIR and PETSC_ARCH to valid values for your PETSc installation.
If PETSc is using MPI (probably it is), run

  make primme_double USE_PETSC=yes

If not (this is rare), run

   make primme_double USE_PETSC=yes USE_MPI=no


* Compile driver with ParaSails (https://computation.llnl.gov/casc/parasails/)

Set PARASAILS_INCLUDE_DIR and PARASAILS_LIB_DIR to the corresponding include path
and library path. Then run  

  make primme_double USE_PARASAILS=yes USE_MPI=yes

        --------------------------------------------------------------
	The comments in the sample drivers show how to run executables
        --------------------------------------------------------------
--------------------------------------------------------------------------------
 Note: these sample drivers are meant to be expository and not high performance
--------------------------------------------------------------------------------
