# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

# ----------------------------------------------------------------------------
# If you submit this package back to Spack as a pull request,
# please first remove this boilerplate and all FIXME comments.
#
# This is a template package file for Spack.  We've put "FIXME"
# next to all the things you'll want to change. Once you've handled
# them, you can save this file and test your package like this:
#
#     spack install primme
#
# You can edit this file again by typing:
#
#     spack edit primme
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------

from spack import *


class Primme(MakefilePackage):
    """"PRIMME, pronounced as prime, is a high-performance library for computing a few eigenvalues/eigenvectors, and
    singular values/vectors. PRIMME is especially optimized for large, difficult problems. Real symmetric and complex
    Hermitian problems, standard 𝐴𝑥 = 𝜆𝑥 and generalized 𝐴𝑥 = 𝜆𝐵𝑥, are supported. Besides standard eigenvalue
    problems with a normal matrix are supported. It can find largest, smallest, or interior singular/eigenvalues, and can
    use preconditioning to accelerate convergence. PRIMME is written in C99, but complete interfaces are provided for
    Fortran, MATLAB, Python, and R."""

    homepage = "http://www.cs.wm.edu/~andreas/software/"
    url      = "https://github.com/primme/primme/archive/refs/tags/v3.2.tar.gz"


    version('3.2', sha256='8ff242a356cea465c9728a26cb6e0487712d9ae51050a362de487e3b13a2fe9b')

    def edit(self, spec, prefix):
        env['PREFIX'] = prefix
