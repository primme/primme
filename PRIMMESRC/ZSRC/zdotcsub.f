c-----------------------------------------------------------------------
      subroutine zdotcsub( zdotcr, n, x, incx, y, incy ) 
      double complex zdotcr, x(*), y(*)
      integer n, incx, incy
c-----------------------------------------------------------------------
c     interface for zdotc
c----------------------------------------------------------------------- 
      external zdotc
      double complex zdotc
c
      zdotcr = zdotc( n, x, incx, y, incy )
c
      return
c-----end-of-zdotcsub---------------------------------------------------
c-----------------------------------------------------------------------
      end
