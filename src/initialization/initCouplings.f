************************************************************************
*
*     initCouplings.f:
*
*     Initialize the values of alpha at the initial scale Q20 and at 
*     the lepton thresholds.
*
************************************************************************
      subroutine initCouplings
*
      implicit none
*
      include "../commons/alpha.h"
      include "../commons/massthrs.h"
**
*     Internal Variables
*
      integer i
      double precision eps
      double complex aQED
      double complex q2thm(3)
      parameter(eps=1d-7)
*
      do i=1,3
         q2thm(i) = q2th(i) * ( 1d0 - eps )
      enddo
*
      AE  = AQED(q2th(1))
      AEM = AQED(q2thm(1))
*
      AM  = AQED(q2th(2))
      AMM = AQED(q2thm(2))
*
      AT  = AQED(q2th(3))
      ATM = AQED(q2thm(3))
*
      return
      end
