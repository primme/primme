
This directory has the test programs as well as usage examples of PRIMME.
It should have at least:

- Makefile               makefile to build the programs.
- readme.txt             this file
- ex_eigs_dseq.c         eigenvalue sequential example in C using double
- ex_eigs_zseq.c                            "    "          using double complex
- ex_eigs_dseqf77.f      eigenvalue sequential example in F77 using double
- ex_eigs_zseqf77.f                         "    "            using double complex
- ex_eigs_zseqxx.cxx     eigenvalue sequential example in C++ using double complex
- ex_eigs_petsc.c        eigenvalue PETSc example in C
- ex_eigs_petscf77.F                        "    "     in F77
- ex_eigs_petscf77ptr.F                     "    "            using pointers
- ex_svds_dseq.c         singular value sequential example in C using double
- ex_svds_zseq.c                            "    "              using double complex
- ex_svds_dseqf77.f      singular value sequential example in F77 using double
- ex_svds_zseqf77.f                         "    "                using double complex
- ex_svds_zseqxx.cxx     singular value sequential example in C++ using double complex
- ex_svds_petsc.c        singular value PETSc example in C
- ex_svds_petscf77.F                        "    "    in F77
- ex_svds_petscf77ptr.F                     "    "           using pointers

The Makefile can perform the next actions:

make examples               build all the examples
make ex_eigs_dseq           build that example; replace by the name of example to build
make test_examples          build and execute all the examples
make clean                  remove object files.
make veryclean              remove object and program files.


* Compile examples with PETSc

First set PETSC_DIR and PETSC_ARCH to valid values for your PETSc installation.
Execute the action, e.g.:

  make ex_eigs_petsc USE_PETSC=yes

--------------------------------------------------------------------------------
 Note: these examples are meant to be expository and not high performance
--------------------------------------------------------------------------------

* Note for developers:

  The examples are generated from the *.m4 files. Please modify the corresponding
  m4 file instead of the examples.
