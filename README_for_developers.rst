
Notes for Developers
====================

.. warning::

    The information in this document is relevant for modifying and maintaining PRIMME.
    Please find the documentation about using PRIMME in the user manual.

The intention is to have a software:

- with great interoperability and easy to integrate in other projects,
- basic error and resources management,
- with parametric and ad hoc polymorphic code, and
- avoid forward declarations.

Languages like D and nim are good candidates, but the high-performance community is very conservative, and a software based on C/C++ has an advantage for the first bullet.

The importance of the last two bullets is a matter of taste. I think that machines with gigabytes of memory and processing giga-Hertz operations can in few seconds specialize the code for eight different number precisions and check in every function call that the function exists and the type of the arguments matches for some function declared at least in the same file, for a project of less than 40,000 lines of code, most of them limited to 80 characters.

Previous and current state of the framework
-------------------------------------------

LAPACK influenced the first versions of PRIMME: a script translated the generic code into the specialized version for every precision, and PRIMME internal functions passed pointers to workspaces.

In versions 2, the source code uses macros and conditional compilations for implementing polymorphism, in a very similar way to PETSc. The macros names SCALAR and REAL are copied from PETSc. Although PETSc code is polymorphic, the framework is built for a single configuration (for instance, double precision complex or just real single precision) and does not allow to compile and call functions in the framework for other types. PRIMME supports that by changing the name of the symbol for every type, in a similar way as Fortran functions are named. The tool ctemplates parses the code and generates a macro with the same name as the function's name but with a value that changes with the type. For instance the macro Num_dot_Sprimme has the value Num_dot_sprimme for real single and Num_dot_zprimme for complex double. The tool ctemplates generates also other macros for calling the real version of the function, for instance, Num_dot_Rprimme. The tool ctemplates also generates the forward definitions for all the function variants. Error and memory management are also copied from PETSc. Dynamic memory is used instead of workspaces, but memory was not free in many cases if an error happens.

In version 3 the source files self-include themselves several times instead of being compiled several times for the different types. Also more macros are added to call for a specific type and also for the CPU version of a function, in order to support GPUs, and half and quad precision. Future versions may use C11 _Generic instead of macros. Also this version tracks memory allocations and notifies if an allocation was not freed, and free allocations in case of error automatically.

The most terrifying aspect of the framework is that the source code communicates with ctemplates using macros. See src/include/templates.h. It is very powerful, but not straightforward to debug and modify. If I did my job good enough, you don't need to know how that works.

Code style
----------

The C source code is intended to be space indented with size 3 spaces and a maximum row size of 80 characters (yes, the first developer was born before the 80s and was initiated in F77, and I'm also old enough to not consider that a totally obsolete practice). A useful setup for vim can be::

    :set ts=3 sw=3 expandtab colorcolumn=80

Opening braces for functions and control flow statements are kept in the same line. if-else statements are styled either in Stroustrup or "} else {". The style for breaking long expressions is free, although most of them are formatted by vim or by the style described in the file `.clang-format`. See ClangFormat_ for integration in several editors. I'm not sure if preprocessor statements follow a style...

Language standard
-----------------

The C source should honor C99 strictly (don't ask me why) and should compile with popular C++ compilers without adding flags. Most C code can be handled by a C++ compiler. A problematic exception is complex types. The C99 `complex` (a macro for the C99 keyword `_Complex`) is not part of any C++ standard. The code explicitly deals with complex types, as explained in section complex_.

To check if the code follows the standard, compilers have usually special flags for that. The next actions should pass without compiler warnings or errors for gcc and clang::

    make clean lib CFLAGS='-std=c99 -Wpedantic -Wall'
    make clean lib CC=g++


Types X/H/SCALAR and REAL
-------------------------

Be aware that every file is compiled 16 times, one for every type, which are half, float, double, quad, and their complex versions, and with GPU. Every time only one of the next macros is defined USE_HALF, USE_FLOAT, USE_DOUBLE, USE_QUAD, USE_HALFCOMPLEX, USE_FLOATCOMPLEX, ..., USE_HALF_MAGMA, ..., USE_HALFCOMPLEX_MAGMA, ... The headers src/include/common.h and src/include/template.h define the values of the rest of the macros used in PRIMME.

The next macros define types:

- XSCALAR and XREAL define the type of pointers passed as eigenvalues, eigenvectors and norms on the `primme` call, and also used in primme_params callbacks, such as matrixMatvec. For instance, when USE_FLOAT is defined, `sprimme` is compiled, and XREAL and XSCALAR define `float`. When USE_FLOATCOMPLEX is defined, `cprimme` is compiled, and XREAL is `float` and `XSCALAR` is `PRIMME_COMPLEX_FLOAT`.

- REAL and SCALAR define the types for PRIMME internal vectors and matrices with dimension as large as the problem matrix. In GPU configurations the pointer is on device memory. Because of that pointer arithmetic is allowed, but not accessing unless the type variant is not GPU. For instance the dense matrices storing the basis of eigenvectors V, and A times V, the matrix W, are of type SCALAR.

- HREAL and HSCALAR define the types for vectors and matrices with dimensions of the rank of the basis. These pointers are on CPU and can be accessed. For instance, the projected matrix :math:`H = V^T A V` and its eigenvectors are of type HSCALAR, and the eigenvalues and the residual norms are of type HREAL. The precision of XSCALAR, SCALAR and HSCALAR, and the corresponding REAL versions are the same, except for half precision, where XSCALAR and SCALAR are half precision and HSCALAR is single precision (float).

Internal functions name and calling convention
----------------------------------------------

Public internal functions whose prototype depends on SCALAR/REAL or their variants, should have TEMPLATE_PLEASE on their definitions and the function name should end in _Sprimme. For instance::

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

The next macros are defined to be used in conditional compilation, the #if statements:

- USE_COMPLEX, only defined for complex variants, USE_HALFCOMPLEX, USE_FLOATCOMPLEX, ...
- USE_HOST, only defined for CPU variants, USE_HALF, USE_FLOAT, ...
- USE_MAGMA, only defined for GPU variants, USE_HALF_MAGMA, USE_FLOAT_MAGMA, ...
- SUPPORTED_TYPE, only defined for supported variants, for instance it is only defined for USE_HALF if the user defines the macro PRIMME_WITH_HALF, and it is only defined for USE_FLOAT_MAGMA when the user defines PRIMME_WITH_MAGMA.
- SUPPORTED_HALF_TYPE, only defined for variants whose half version is supported.

The next macros return a member of the enum primme_op_datatype, primme_op_half, primme_op_float, primme_op_double or primme_op_quad. They cannot be used in conditional compilations:

- PRIMME_OP_SCALAR and PRIMME_OP_REAL: the precision for SCALAR and REAL
- PRIMME_OP_XSCALAR and PRIMME_OP_XREAL: the precision for XSCALAR and XREAL
- PRIMME_OP_HSCALAR and PRIMME_OP_HREAL: the precision for HSCALAR and HREAL

The macro MACHINE_EPSILON has the machine epsilon of the type for SCALAR and REAL.

.. _complex :

Complex, half, quad
-------------------

The header :file:`include/primme.h` defines the complex types for half, single, double and quad, named PRIMME_COMPLEX_HALF, PRIMME_COMPLEX_FLOAT, PRIMME_COMPLEX_DOUBLE, PRIMME_COMPLEX_QUAD. Use the next macros for expressions with XSCALAR/SCALAR/HSCALAR type:

- REAL_PART(A): the real part of A
- IMAGINARY_PART(A): the imaginary part of A
- ABS(A): the absolute value of A
- CONJ(A): the complex conjugate of A

No C or C++ standard requires to support half quadruple precision, and neither their complex versions. Quadruple and complex quadruple are fully supported in gcc and clang. Half precision is supported by gcc for architectures with native arithmetic support (see gccHalf_). Clang supports a storage type __fp16, and the arithmetic is done by promoting the value to single precision. For some reason, std::complex<__fp16> does not work. So PRIMME defines a set of macros that implement complex arithmetic in that case by promoting the half complex values to float complex. For the next definitions `A` is SCALAR and `B` is HSCALAR, which should have support for complex arithmetic.

- SET_ZERO(A)       : set A = 0
- SET_COMPLEX(A, B) : set A = B
- TO_COMPLEX(A)     : cast A to HSCALAR
- PLUS_EQUAL(A, B)  : set A += B
- MULT_EQUAL(A, B)  : set A `*=` B

Memory and error management
---------------------------

Recent versions of PRIMME are using dynamic memory to manage the memory. In general the use of dynamic memory simplifies the code by not having to take care of providing enough working space for all subsequent calls. The small drawback of dynamic memory is to mingle with error management. The goal is to avoid writing specific code to free allocated memory in case of an error happening in the body of a function.

By default, calls to PRIMME internal functions should be made under an error checker macro, CHKERR, CHKERRM or CHKERRA, if the function returns an error code. Also these macros expects the variable ctx, which is a struct with information about the allocations besides other things. Consider the next function::

    TEMPLATE_PLEASE int dummy_Sprimme(SCALAR *v, PRIMME_INT n, primme_context ctx) {
        SCALAR *x;
        CHKERR(Num_malloc_Sprimme(n, &x, ctx));
        CHKERR(Num_copy_Sprimme(n, v, 1, x, 1, ctx));
        CHKERR(Num_free_Sprimme(x, ctx));
        return 0;
    }

If Num_malloc_Sprimme or Num_copy_Sprimme or Num_free_Sprimme return a nonzero value, the function dummy_Sprimme immediately returns that value to the caller. If the error happens in Num_copy_Sprimme, the array allocated by Num_malloc_Sprimme is freed by CHKERR in which the call to dummy_Sprimme is in. A function must notify if some allocations are not going to be freed on purpose after the function finished with no error, by calling `Mem_keep_frame(ctx)`.

WTF is this!? Why not using C++
-------------------------------

You're right! We don't have nowadays a technical justification for not using C++, a language as well-established and multiplatform as C, and with support for polymorphism and RAII and exceptions. The advantages of that support would be to have a cleaner code without conditional compilations and most of the macros, and clearer error messages than the ones that C gives involving macros. There are few drawbacks that can be worked out. Error messages involving std::complex can be hard to read. The most pressing issue is that some parts of the code will be much nicer with partial template specialization and if-constexpr, but these features are only available in recent C++ standard, C++ 14 and C++ 17 respectively. Also one can implement a reverse communication interface for PRIMME using coroutines, currently implemented in Boost, and probably part of C++ 20.


.. _ClangFormat : https://clang.llvm.org/docs/ClangFormat.html
.. _gccHalf: https://gcc.gnu.org/onlinedocs/gcc/Half-Precision.html
