************************************************************************
*
*     NStructureFunctions.f:
*
*     Routine that returns the evolved N-space distributions.
*
************************************************************************
      subroutine NStructureFunctions(N,nQ,Q,SFN)
*
      implicit none
*
      include "../commons/renfacscales.h"
      include "../commons/ns.h"
      include "../commons/ipt.h"
**
*     Input Variables
*
      integer nQ
      double complex Q(100)
      double complex N
**
*     Internal Variables
*
      integer iptbkp
      integer isf,ifl,id,iQ
      double complex muF(100)
      double complex xfevN(13)
      double complex CF(3,6,13)
**
*     Output Variables
*
      double complex SFN(3,6)
*
*     Call evolved N-space distributions in the factorization scale.
*     (All scales rescaled by KFAQ but the fist one)
*
      muF(1) = Q(1)
      do iQ=2,nQ
         muF(iQ) = kfacQ * Q(iQ)
      enddo
*
      if(ns.eq."FFNS")then
         iptbkp = ipt
         ipt    = max(0,ipt-1)
         call NDistributions(N,nQ,muF,xfevN)
         ipt    = iptbkp
      else
         call NDistributions(N,nQ,muF,xfevN)
      endif
*
*     Call ceofficient functions at the final scale
*
      call NCoefficientFunctions(N,Q(nQ),CF)
*
*     Compute N-space structure functions
*
      do isf=1,3
         do ifl=1,6
            SFN(isf,ifl) = (0d0,0d0)
            do id=1,13
               SFN(isf,ifl) = SFN(isf,ifl) + CF(isf,ifl,id) * xfevN(id)
            enddo
          enddo
      enddo
*
      return
      end
