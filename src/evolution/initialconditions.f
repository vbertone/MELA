************************************************************************
*
*     initialconditions.f:
*
*     These routines return PDFs/FFs at the intitial scale in N space.
*
************************************************************************
      subroutine electronPDFsn(N, npdf)
*
      implicit none
**
*     Input Variables
*
      double complex N
**
*     Internal Variables
*
      integer ipdf
**
*     Output Variables
*
      double complex npdf(-3:3)
*
*     Initialize PDFs to zero
*
      do ipdf = -3,3
         npdf(ipdf) = (0d0,0d0)
      enddo
*
      return
      end
*
************************************************************************
      subroutine electronFFsn(N, nff)
*
      implicit none
**
*     Input Variables
*
      double complex N
**
*     Internal Variables
*
      integer ipdf
**
*     Output Variables
*
      double complex nff(-3:3)
*
*     Initialize FFs to zero
*
      do ipdf = -3,3
         nff(ipdf) = (0d0,0d0)
      enddo
*
      return
      end
