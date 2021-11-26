************************************************************************
*
*     SetThresholds.f:
*
*     Sets the value of the thresholds
*
************************************************************************
      subroutine SetThresholds(me,mu,md,ms,mm,mc,mt,mb,mtp)
*
      implicit none
*
      include "../commons/massthrs.h"
**
*     Input Variables
*
      double precision me,mu,md,ms,mm,mc,mt,mb,mtp
*
      Q2TH(1) = ME**2
      Q2TH(2) = MU**2
      Q2TH(3) = MD**2
      Q2TH(4) = MS**2
      Q2TH(5) = MM**2
      Q2TH(6) = MC**2
      Q2TH(7) = MT**2
      Q2TH(8) = MB**2
      Q2TH(9) = MTP**2
*
      return
      end
************************************************************************
      subroutine GetThresholds(q2thrs)
      implicit none
      include "../commons/massthrs.h"      
      double precision q2thrs(9)
      q2thrs=q2th
      return
      end
