c amux is part of the SPARSKIT library by Y. Saad at the University of Minnesota
c For the full SPARSKIT library visit:
c                          www.cs.umn.edu/~saad/
c
c-----------------------------------------------------------------------
      subroutine atmux (n, x, y, a, ja, ia)
      real*8 x(*), y(*), a(*) 
      integer n, ia(*), ja(*)
c-----------------------------------------------------------------------
c         transp( A ) times a vector
c----------------------------------------------------------------------- 
c multiplies the transpose of a matrix by a vector when the original
c matrix is stored in compressed sparse row storage. Can also be
c viewed as the product of a matrix by a vector when the original
c matrix is stored in the compressed sparse column format.
c-----------------------------------------------------------------------
c
c on entry:
c----------
c n     = row dimension of A
c x     = real array of length equal to the column dimension of
c         the A matrix.
c a, ja,
c    ia = input matrix in compressed sparse row format.
c
c on return:
c-----------
c y     = real array of length n, containing the product y=transp(A)*x
c
c-----------------------------------------------------------------------
c     local variables 
c
      integer i, k 
c-----------------------------------------------------------------------
c
c     zero out output vector
c 
      do 1 i=1,n
         y(i) = 0.0
 1    continue
c
c loop over the rows
c
      do 100 i = 1,n
         do 99 k=ia(i), ia(i+1)-1 
            y(ja(k)) = y(ja(k)) + x(i)*a(k)
 99      continue
 100  continue
c
      return
c-------------end-of-atmux---------------------------------------------- 
c-----------------------------------------------------------------------
      end
