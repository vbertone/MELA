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
      AE  = AQED(q2th(1))
      AU  = AQED(q2th(2))
      AD  = AQED(q2th(3))
      AS  = AQED(q2th(4))
      AM  = AQED(q2th(5))
      AC  = AQED(q2th(6))
      AT  = AQED(q2th(7))
      AB  = AQED(q2th(8))
      ATP = AQED(q2th(9))
*
      return
      end
