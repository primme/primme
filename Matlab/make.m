function make
%MAKE Make PRIMME's Matlab module

   % Get the path of this directory
   [p,~,~] = fileparts(mfilename('fullpath'));

   % Save the current path and changed the working path to here
   w = pwd;
   cd(p);

   % Compile PRIMME
   try
      if exist('OCTAVE_VERSION')
         mex -O -I../include -I../src/include -DF77UNDERSCORE ../src/eigs/auxiliary_eigs.c ../src/eigs/auxiliary_eigs_normal.c ../src/eigs/convergence.c ../src/eigs/correction.c ../src/eigs/factorize.c ../src/eigs/init.c ../src/eigs/inner_solve.c ../src/eigs/main_iter.c ../src/eigs/ortho.c ../src/eigs/primme_c.c ../src/eigs/primme_f77.c ../src/eigs/primme_interface.c ../src/eigs/restart.c ../src/eigs/solve_projection.c ../src/eigs/update_projection.c ../src/eigs/update_W.c ../src/linalg/auxiliary.c ../src/linalg/blaslapack.c ../src/linalg/magma_wrapper.c ../src/linalg/memman.c ../src/linalg/wtime.c ../src/svds/primme_svds_c.c ../src/svds/primme_svds_f77.c ../src/svds/primme_svds_interface.c primme_mex.cpp  -largeArrayDims -lm -output primme_mex
      else
         mex -O -I../include -I../src/include -DF77UNDERSCORE -DPRIMME_BLASINT_SIZE=64 ../src/eigs/auxiliary_eigs.c ../src/eigs/auxiliary_eigs_normal.c ../src/eigs/convergence.c ../src/eigs/correction.c ../src/eigs/factorize.c ../src/eigs/init.c ../src/eigs/inner_solve.c ../src/eigs/main_iter.c ../src/eigs/ortho.c ../src/eigs/primme_c.c ../src/eigs/primme_f77.c ../src/eigs/primme_interface.c ../src/eigs/restart.c ../src/eigs/solve_projection.c ../src/eigs/update_projection.c ../src/eigs/update_W.c ../src/linalg/auxiliary.c ../src/linalg/blaslapack.c ../src/linalg/magma_wrapper.c ../src/linalg/memman.c ../src/linalg/wtime.c ../src/svds/primme_svds_c.c ../src/svds/primme_svds_f77.c ../src/svds/primme_svds_interface.c primme_mex.cpp  -largeArrayDims -lmwlapack -lmwblas -lm -output primme_mex
      end
   catch me
      cd(w);
      rethrow(me);
   end

   % Change back the current path
   cd(w)
