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
      double precision aQED
*
      do i=1,9
         ath(i) = aQED(q2th(i))
      enddo
*
      return
      end
************************************************************************
      subroutine Geta0(a0)
      implicit none
      include "../commons/alpha.h"
      double precision a0
      a0 = ath(1)
      return
      end
      
