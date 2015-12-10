************************************************************************
*
*     initCouplings.f:
*
*     Initialize the values of alphas at the initial scale Q20 and at 
*     the heavy quark thresholds.
*
************************************************************************
      subroutine initCouplings
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/massthrs.h"
**
*     Internal Variables
*
      integer i
      double precision eps
      double complex aQCD
      double complex q2thm(4:6)
      parameter(eps=1d-7)
*
      do i=4,6
         q2thm(i) = q2th(i) - eps
      enddo
*
      ASC   = AQCD(q2th(4))
      ASCM  = AQCD(q2thm(4))
*
      ASB   = AQCD(q2th(5))
      ASBM  = AQCD(q2thm(5))
*
      AST   = AQCD(q2th(6))
      ASTM  = AQCD(q2thm(6))
*
      return
      end
