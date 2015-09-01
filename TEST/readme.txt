
This directory has the test programs as well as usage examples of PRIMME.
It should have:

- Makefile             makefile to build the programs.
- driver.c             driver that can read MTX and PETSc binary matrices.
- COMMON/              with source used by driver.c.
    csr.h, csr.c       routines for matrices CSR
    mmio.h, mmio.c     MatrixMarket IO routines.
    native.h, mat.c    wrapper for CSR matrix and sequential ILUT.
    num.h              constants
    parasailsw.h, .c   wrapper for Parasails matrix and preconditioner.
    petscw.h, .c       wrapper for PETSc matrices and preconditioners.
    shared_utils.h, .c IO routines for primme_params and driver options.
    ssrcsr.c           routine to convert from Sym Sparse Row to CSR.
    amux.f             routine for CSR matrix-vector product.
    ilut.f             routine for sequential ILUT.
    zamux.f            routine for complex CSR matrix-vector product.
    zilut.f            routine for complex sequential ILUT.
- DriverConf           example of driver configuration file used by the driver.
- MinConf, LeanConf,
  FullConf             examples of PRIMME configuration file used by the driver.
- LUNDA.mtx            matrix used for testing and in DriverConf as an example.
- tests/               configuration files for testing purpose.
- ex_dseq.c            example of sequential program calling PRIMME.
- ex zseq.c            example of sequential complex program.
- ex_petsc.c           example of PETSc program.

The Makefile can perform the next actions:

make primme_double          build driver in double floats.
make primme_doublecomplex     "     "    in complex double floats.
make ex_dseq                build example in C
make ex_zseq                  "     "
make ex_petsc                 "     "
make ex_dseqf77             build example in Fortran
make ex_zseqf77               "     "
make ex_petscf77              "     "
make ex_petscf77ptr           "     "
make test                   build and execute a simple example in double and complex.
make all_tests_double       test all configurations in "tests" for doubles.
make all_tests_doublecomplex  "   "          "      "     "    for complex.
make clean                  remove object files.
make veryclean              remove object and program files.


* Compile driver with PETSc

First set PETSC_DIR and PETSC_ARCH to valid values for your PETSc installation.
If PETSc is using MPI (probably it is), run

  make primme_double USE_PETSC=yes

If not (really? quite unusual), run

   make primme_double USE_PETSC=yes USE_MPI=no


* Compile driver with Parasails

Set PARASAILS_INCLUDE_DIR and PARASAILS_LUB_DIR to the corresponding include path
and library path. Then run  

  make primme_double USE_PARASAILS=yes USE_MPI=yes

