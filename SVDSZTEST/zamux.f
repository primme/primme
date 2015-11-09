c zamux is part of the SPARSKIT library by Y. Saad at the University of Minnesota
c For the full SPARSKIT library visit:
c                          www.cs.umn.edu/~saad/
c
c----------------------------------------------------------------------c
      subroutine zamux (n, x, y, a,ja,ia) 
      complex*16  x(*), y(*), a(*) 
      integer n, ja(*), ia(*)
c-----------------------------------------------------------------------
c         A times a vector
c----------------------------------------------------------------------- 
c multiplies a matrix by a vector using the dot product form
c Matrix A is stored in compressed sparse row storage.
c
c on entry:
c----------
c n     = row dimension of A
c x     = complex array of length equal to the column dimension of
c         the A matrix.
c a, ja,
c    ia = input matrix in compressed sparse row format.
c
c on return:
c-----------
c y     = complex array of length n, containing the product y=Ax
c
c-----------------------------------------------------------------------
c local variables
c
      complex*16 t
      integer i, k
c-----------------------------------------------------------------------
      do 100 i = 1,n
c
c     compute the inner product of row i with vector x
c 
         t = DCMPLX(0.0d0,0.0d0)
         do 99 k=ia(i), ia(i+1)-1 
            t = t + a(k)*x(ja(k))
 99      continue
c
c     store result in y(i) 
c
         y(i) = t
 100  continue
c
      return
c---------end-of-zamux---------------------------------------------------
c-----------------------------------------------------------------------
      end
