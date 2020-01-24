
Notes for Developers
====================

::

    The information in this document is relevant for modifying and maintaining PRIMME.
    Please find the documentation about using PRIMME in the user manual.

The intention is to have a library:

- interoperable and easy to integrate into other projects,
- with error and resources management, and
- maintainable and extendable.

The high-performance community is fond of C/C++, and software based on that has an advantage for the first bullet. We are not going to consider other languages like D, nim, and Rust, although they may make a strong case for the rest of the bullets.
Also, the evaluation of the last bullet is very subjective, and the high variability introduced by the human factor makes the task of taking rational decisions about languages and programming strategies cumbersome.

Numerical codes, like PRIMME, are mostly made of functions that admit book-case abstraction regarding:

- the number system parametrized by the precision (half, single, double, quadruple) and the complexity (real, complex); and
- architecture/platform (host/accelerator, shared and distributed memory).

Even rudimentary forms of metaprogramming such as parametric and ad hoc polymorphism can help experienced developers to write and maintain these functions for all the variants in less time. Besides, I am not down to reimplement every function 32 times.

The previous and current state of the framework
-------------------------------------------

In version 1 of PRIMME, a script translated the generic code into the specialized version for every precision. A few memory allocations are made at the beginning of the execution of the main routine, and internal functions passed pointers to workspaces, making error and memory management simple.

In versions 2, the source code uses macros and conditional compilation for implementing polymorphism. The tool ctemplate_ generates macros associated with the internal functions that implicitly change the function's name for every type, following a similar scheme as Fortran 77 functions are named.  See section func_names_. The tool also generates forward definitions for all the function variants. Error and memory management are handled similarly as PETSc. Dynamic memory is used instead of workspaces, but the memory was not freed in many cases if an error happens.

In version 3, which is the current version, the source files include themselves several times instead of calling the compiler several times for the different types. It is introduced support for half precision and GPU. More macros are added to call for a specific type and also for the CPU version of a function. CHKERR_ tracks memory allocations and notifies if an allocation was not freed.

The terrifying aspect of the framework is that the source code communicates with ctemplate_ using macros. See an example of this in the definition of TEMPLATE_PLEASE_. It is powerful, but not straightforward to debug and modify. If I did my job good enough, you don't need to know how that works.

Code style
----------

The C source code is intended to be space indented with size 3 spaces and a maximum row size of 80 characters (yes, the first developer was born before the 80s and initiated in F77, and I'm also old enough not to consider that a completely obsolete practice). A useful setup for vim can be::

    :set ts=3 sw=3 expandtab colorcolumn=80

Opening braces for functions and control flow statements are kept in the same line. if-else statements are styled either in Stroustrup or `} else {`. The style for breaking long expressions is free, although most of them are formatted by vim or by the style described in the file `.clang-format`. See ClangFormat_ for integration in several editors. I'm not sure if preprocessor statements follow a style.

Language standard
-----------------

The C source should honor C99 (don't ask me why) and should compile with popular C++ compilers without adding flags. A C++ compiler can handle most of the C99 standard. A problematic exception is complex types. The C99 `complex` (a macro for the C99 keyword `_Complex`) is not part of any C++ standard. The code explicitly deals with complex types, as explained in section complex_.

To check if the code follows the standard, compilers usually have special flags for that. The next actions should pass without compiler warnings or errors for gcc and clang::

    make clean lib CFLAGS='-std=c99 -Wpedantic -Wall'
    make clean lib CC=g++

Structure of the source file
----------------------------

Most of the files start with a long comment that includes the license, the name of the file, and a brief description of the purpose. 

The `.c` source files requiring the automatic generation of headers and auto-specialization should follow the following guidelines in including headers.
Note that all generated headers are listed in the Makefile macros AUTOMATED_HEADERS_LINALG_ and AUTOMATED_HEADERS_OTHERS_. The generated headers are named after the source file name.

#. The macro `THIS_FILE` should be defined before any `include`. 
#. The headres `common.h` and `template.h` should be included next, or other headers that include them, for instance `numerical.h`.
#. All generated headers should be included under the conditional compilation `#ifndef CHECK_TEMPLATE`

Types X/H/SCALAR and REAL
-------------------------

Be aware that every file is compiled 16 times, one for every type, which are half, float, double, quad, and their complex versions, and with GPU. Every time only one of the next macros is defined `USE_HALF`, `USE_FLOAT`, `USE_DOUBLE`, `USE_QUAD`, `USE_HALFCOMPLEX`, `USE_FLOATCOMPLEX`, ..., `USE_HALF_MAGMA`, ..., `USE_HALFCOMPLEX_MAGMA`, ... The rest of the macros used in PRIMME are defined depending on which of the previous macros is defined.

The following macros define types:

- XSCALAR_ and XREAL_ define the type of pointers passed as eigenvalues, eigenvectors and norms on the `primme` call, and also used in primme_params callbacks, such as matrixMatvec. For instance, when USE_FLOAT is defined, `sprimme` is compiled, and XREAL and XSCALAR_ define `float`. When USE_FLOATCOMPLEX is defined, `cprimme` is compiled, and XREAL is `float` and `XSCALAR` is `PRIMME_COMPLEX_FLOAT`.

- REAL_ and SCALAR_ define the types for PRIMME internal vectors and matrices with dimension as large as the problem matrix. In GPU configurations the pointer is on device memory. Because of that pointer arithmetic is allowed, but not accessing unless the type variant is not GPU. For instance, the dense matrices storing the basis of eigenvectors V, and A times V, the matrix W, are of type SCALAR.

- HREAL_ and HSCALAR_ define the types for vectors and matrices with dimensions of the rank of the basis. These pointers are on CPU and can be accessed. For instance, the projected matrix :math:`H = V^T A V` and its eigenvectors are of type HSCALAR, and the eigenvalues and the residual norms are of type HREAL. The precision of XSCALAR, SCALAR and HSCALAR, and the corresponding REAL versions are the same, except for half precision, where XSCALAR and SCALAR are half precision and HSCALAR is single precision (float).

.. _func_names :

Internal functions name and calling convention
----------------------------------------------

Public internal functions whose prototype depends on SCALAR/REAL or their variants should have TEMPLATE_PLEASE_ on their definitions and the function name should end in _Sprimme. For instance::

    TEMPLATE_PLEASE int dummy_Sprimme(SCALAR *v, PRIMME_INT n, primme_context ctx) {
       ...
    }

When calling the function, please match the argument types with the suffix of the function's name. For instance::

    SCALAR *x;
    dummy_Sprimme(x, n, ctx);
    REAL *x;
    dummy_Rprimme(x, n, ctx);
    HSCALAR *x;
    dummy_SHprimme(x, n, ctx);
    HREAL *x;
    dummy_RHprimme(x, n, ctx);
    XSCALAR *x;
    dummy_SXprimme(x, n, ctx);
    XREAL *x;
    dummy_RXprimme(x, n, ctx);

For calling the function for a particular precision, add the letter h, s, d or q for the half, single, double and quad precision as follows::

    SCALAR *x; void *x_half;
    Num_copy_astype_Sprimme(x, n, x_half, primme_op_half, ctx); // copy x into x_half as half precision
    dummy_Shprimme(x_half, n, ctx);
    HSCALAR *x; void *x_single;
    Num_copy_astype_Sprimme(x, n, x_single, primme_op_float, ctx); // copy x into x_single as single precision
    dummy_SHsprimme(x_single, n, ctx);
  

Conditional compilation and type inspection
-------------------------------------------

The following macros are defined to be used in conditional compilation, the `#if` statements:

- USE_COMPLEX_, only defined for complex variants, USE_HALFCOMPLEX, USE_FLOATCOMPLEX, ...
- USE_HOST_, only defined for CPU variants, USE_HALF, USE_FLOAT, ...
- USE_MAGMA_, only defined for GPU variants, USE_HALF_MAGMA, USE_FLOAT_MAGMA, ...
- SUPPORTED_TYPE_, only defined for supported variants, for instance, it is only defined for USE_HALF if the user defines the macro PRIMME_WITH_HALF, and it is only defined for USE_FLOAT_MAGMA when the user defines PRIMME_WITH_MAGMA.
- SUPPORTED_HALF_TYPE_, only defined for variants whose half version is supported.

The next macros return a member of the enum primme_op_datatype, primme_op_half, primme_op_float, primme_op_double or primme_op_quad. They cannot be used in conditional compilations:

- PRIMME_OP_SCALAR and PRIMME_OP_REAL: the precision for SCALAR and REAL
- PRIMME_OP_XSCALAR and PRIMME_OP_XREAL: the precision for XSCALAR and XREAL
- PRIMME_OP_HSCALAR and PRIMME_OP_HREAL: the precision for HSCALAR and HREAL

The macro MACHINE_EPSILON_ has the machine epsilon of the type for SCALAR and REAL.

.. _complex :

Complex, half, quad
-------------------

The header `include/primme.h` defines the complex types for half, single, double and quad, named PRIMME_COMPLEX_HALF, PRIMME_COMPLEX_FLOAT, PRIMME_COMPLEX_DOUBLE, PRIMME_COMPLEX_QUAD. Use the following macros for expressions with XSCALAR/SCALAR/HSCALAR type:

- REAL_PART_ (A): the real part of A
- IMAGINARY_PART_ (A): the imaginary part of A
- ABS_ (A): the absolute value of A
- CONJ_ (A): the complex conjugate of A

No C or C++ standard requires to support half quadruple precision, and neither their complex versions. Quadruple and complex quadruple are fully supported in gcc and clang. Half precision is supported by gcc for architectures with native arithmetic support (see gccHalf_). Clang supports a storage type __fp16, and the arithmetic is done by promoting the value to single precision. For some reason, std::complex<__fp16> does not work. So PRIMME defines a set of macros that implement complex arithmetic in that case by promoting the half complex values to float complex. For the following definitions `A` is SCALAR and `B` is HSCALAR, which should have support for complex arithmetic.

- SET_ZERO_ (A)       : set A = 0
- SET_COMPLEX_ (A, B) : set A = B
- TO_COMPLEX_ (A)     : cast A to HSCALAR
- PLUS_EQUAL_ (A, B)  : set A += B
- MULT_EQUAL_ (A, B)  : set A `*=` B

Memory and error management
---------------------------

Recent versions of PRIMME are using dynamic memory to manage the memory. In general, the use of dynamic memory simplifies the code by not having to take care of providing enough working space for all subsequent calls. The small drawback of dynamic memory is to mingle with error management. The goal is to avoid writing specific code to free allocated memory in case of an error happening in the body of a function.

By default, calls to PRIMME internal functions should be made under an error checker macro, CHKERR_, CHKERRM_ or CHKERRA_, if the function returns an error code. These macros expect the variable `ctx`, which is a `struct` with information about the allocations besides other things. Consider the next function::

    TEMPLATE_PLEASE int dummy_Sprimme(SCALAR *v, PRIMME_INT n, primme_context ctx) {
        SCALAR *x;
        CHKERR(Num_malloc_Sprimme(n, &x, ctx));
        CHKERR(Num_copy_Sprimme(n, v, 1, x, 1, ctx));
        CHKERR(Num_free_Sprimme(x, ctx));
        return 0;
    }

If Num_malloc_Sprimme or Num_copy_Sprimme or Num_free_Sprimme return a nonzero value, the function dummy_Sprimme immediately returns that value to the caller. If the error happens in Num_copy_Sprimme, the CHKERR in which the dummy_Sprimme call is enclosed frees the array allocated by Num_malloc_Sprimme. A function must notify if some allocations are not going to be freed on purpose after the function finishing with no error, by calling `Mem_keep_frame(ctx)`.

WTF is this!? Why not using C++
-------------------------------

You're right! We don't have much of an excuse for not using C++, a language as well-established and multiplatform as C, and with support for polymorphism and RAII and exceptions. The advantages of that support would be to have a cleaner code without conditional compilations and fewer macros, and clearer error messages than the ones that C gives involving macros. However, there are a few drawbacks that can be worked out. The minor issues include that error messages involving std::complex can be hard to read, and debugging C++ functions is slightly more tedious. The most pressing issue is that to remove most the C macros, we need the more advanced, and recent, parts of the C++ standard, such as partial template specialization and if-constexpr, C++ 14 and C++ 17 respectively. Additionally, one can implement a reverse communication interface for PRIMME using coroutines, currently implemented in Boost, and likely part of C++ 20.


.. _ClangFormat : https://clang.llvm.org/docs/ClangFormat.html
.. _gccHalf : https://gcc.gnu.org/onlinedocs/gcc/Half-Precision.html
.. _PETSc : http://www.mcs.anl.gov/petsc/
.. _SCALAR : https://github.com/primme/primme/blob/master/src/include/common.h#L133
.. _XSCALAR : https://github.com/primme/primme/blob/master/src/include/common.h#L134
.. _HSCALAR : https://github.com/primme/primme/blob/master/src/include/common.h#L135
.. _REAL : https://github.com/primme/primme/blob/master/src/include/common.h#L141
.. _XREAL : https://github.com/primme/primme/blob/master/src/include/common.h#L142
.. _HREAL : https://github.com/primme/primme/blob/master/src/include/common.h#L143
.. _HREAL : https://github.com/primme/primme/blob/master/src/include/common.h#L143
.. _ctemplate : https://github.com/primme/primme/blob/master/src/tools/ctemplate
.. _TEMPLATE_PLEASE : https://github.com/primme/primme/blob/master/src/include/template.h#L266
.. _USE_HOST : https://github.com/primme/primme/blob/master/src/include/template.h#L85
.. _USE_MAGMA : https://github.com/primme/primme/blob/master/src/include/template.h#L91
.. _USE_REAL : https://github.com/primme/primme/blob/master/src/include/template.h#L106
.. _USE_COMPLEX : https://github.com/primme/primme/blob/master/src/include/template.h#L111
.. _SUPPORTED_HALF_TYPE : https://github.com/primme/primme/blob/master/src/include/template.h#L157
.. _SUPPORTED_TYPE : https://github.com/primme/primme/blob/master/src/include/template.h#L167
.. _REAL_PART : https://github.com/primme/primme/blob/master/src/include/template.h#L178
.. _IMAGINARY_PART : https://github.com/primme/primme/blob/master/src/include/template.h#L179
.. _ABS : https://github.com/primme/primme/blob/master/src/include/template.h#L180
.. _CONJ : https://github.com/primme/primme/blob/master/src/include/template.h#L181
.. _SET_ZERO : https://github.com/primme/primme/blob/master/src/include/template.h#L209
.. _SET_COMPLEX : https://github.com/primme/primme/blob/master/src/include/template.h#L210
.. _TO_COMPLEX : https://github.com/primme/primme/blob/master/src/include/template.h#L211
.. _PLUS_EQUAL : https://github.com/primme/primme/blob/master/src/include/template.h#L216
.. _MULT_EQUAL : https://github.com/primme/primme/blob/master/src/include/template.h#L217
.. _MACHINE_PRECISON : https://github.com/primme/primme/blob/master/src/include/common.h#L141
.. _CHKERR : https://github.com/primme/primme/blob/master/src/include/common.h#L447
.. _CHKERRM : https://github.com/primme/primme/blob/master/src/include/common.h#L479
.. _CHKERRA : https://github.com/primme/primme/blob/master/src/include/common.h#L507
.. _Num_dot_Sprimme : https://github.com/primme/primme/blob/master/src/include/blaslapack.h#L1700
.. _Num_dot_Rprimme : https://github.com/primme/primme/blob/master/src/include/blaslapack.h#L1703
.. _AUTOMATED_HEADERS_LINALG : https://github.com/primme/primme/blob/master/src/Makefile#L51
.. _AUTOMATED_HEADERS_OTHERS : https://github.com/primme/primme/blob/master/src/Makefile#L55

