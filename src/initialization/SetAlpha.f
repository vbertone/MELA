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
************************************************************************      
      subroutine GetAlphaRef(arefout)
      implicit none
      include "../commons/alpha.h"      
      double precision arefout
      arefout = aref
      return
      end
************************************************************************      
      subroutine GetAlphaQref(Qrefout)
      implicit none
      include "../commons/alpha.h"      
      double precision Qrefout
      Qrefout = dsqrt(Q2ref)
      return
      end
************************************************************************            
      subroutine GetAlphaFix(aemfixout)
      implicit none
      include "../commons/alpha.h"
      logical aemfixout
      aemfixout = aemfix
      return
      end
************************************************************************                  
      subroutine SetAlphaFix(aemfixin)
      implicit none
      include "../commons/alpha.h"
      logical aemfixin
      aemfix = aemfixin
      return
      end
      
