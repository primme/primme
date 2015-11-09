      subroutine atmuxr (m, n, x, y, a, ja, ia)
      real*8 x(*), y(*), a(*) 
      integer m, n, ia(*), ja(*)
c-----------------------------------------------------------------------
c         transp( A ) times a vector, A can be rectangular
c----------------------------------------------------------------------- 
c See also atmux.  The essential difference is how the solution vector
c is initially zeroed.  If using this to multiply rectangular CSC 
c matrices by a vector, m number of rows, n is number of columns.
c-----------------------------------------------------------------------
c
c on entry:
c----------
c m     = column dimension of A
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
      do 1 i=1,m
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
c-------------end-of-atmuxr--------------------------------------------- 
c-----------------------------------------------------------------------
      end
