************************************************************************
*
*     SetAlpha.f:
*
*     Sets the value of alpha at a give reference scale
*
************************************************************************
      subroutine SetAlpha(ain,Qin)
*
      implicit none
*
      include "../commons/alpha.h"
**
*     Input Variables
*
      double precision ain, Qin
*
      aref = ain
      Q2ref = Qin**2
*
      return
      end
