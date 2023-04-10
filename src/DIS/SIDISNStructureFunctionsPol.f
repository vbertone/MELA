************************************************************************
*
*     SIDISNStructureFunctionsPol.f:
*
*     Routine that returns the semi-inclusive polarised structure
*     functions in N-space.
*
************************************************************************
      subroutine SIDISNStructureFunctionsPol(N, M, Q0, Q, SFN)
*
      implicit none
*
      include "../commons/evol.h"
      include "../commons/renfacscales.h"
      include "../commons/alphas.h"
      include "../commons/ipt.h"
      include "../commons/pol.h"
**
*     Input Variables
*
      double complex Q0, Q
      double complex N, M
**
*     Internal Variables
*
      integer i
      double complex muF
      double complex Qv(2)
      double complex xfevN(13), xfphN(-6:6)
      double complex xdevN(13), xdphN(-6:6)
      double complex G11QQ_SIDIS, G11GQ_SIDIS, G11QG_SIDIS
      double complex G1qq, G1gq, G1qg
      double complex Lqq, Lgq, Lqg
      double complex ch
      character*5 evolbkp
**
*     Output Variables
*
      double complex SFN
*
*     Call evolved N-space distributions in the factorization scale.
*     (All scales rescaled by KFAQ but the fist one)
*
      muF = kfacQ * Q
      Qv(1) = Q0
      Qv(2) = muF
*
*     Backup evolution mode
*
      evolbkp = evol
*
*     Call PDFs at the final scale
*
      pol  = "ON"
      evol = "SPACE"
      call NDistributions(N, 2, Qv, xfevN)
      call evln2lhac(xfevN, xfphN)
*
*     Call FFs at the final scale and rotate them into the physical
*     basis
*
      pol  = "OFF"
      evol = "TIME"
      call NDistributions(M, 2, Qv, xdevN)
      call evln2lhac(xdevN, xdphN)
*
*     Restore evolution mode
*
      evol = evolbkp
*
*     Call coefficient functions
*
      G1qq = (1d0, 0d0)
      G1gq = (0d0, 0d0)
      G1qg = (0d0, 0d0)
      if (ipt.ge.1) then
         G1qq = G1qq + asDIS * G11QQ_SIDIS(N, M)
         G1gq = G1gq + asDIS * G11GQ_SIDIS(N, M)
         G1qg = G1qg + asDIS * G11QG_SIDIS(N, M)
      endif
*
*     Assemble PDF-FF combinations
*
      Lqq = (0d0, 0d0)
      Lgq = (0d0, 0d0)
      Lqg = (0d0, 0d0)
      do i = 1, 6
         ch = bq(i)
         if (i.eq.1) ch = bq(2)
         if (i.eq.2) ch = bq(1)
         Lqq = Lqq + ch * (xfphN(-i) * xdphN(-i) + xfphN(i) * xdphN(i))
         Lqg = Lqg + ch * ( xdphN(-i) + xdphN(i) ) * xfphN(0)
         Lgq = Lgq + ch * xdphN(0) * ( xfphN(-i) + xfphN(i) )
      enddo
*
*     Compute N-space structure functions
*
      SFN = G1qq * Lqq + G1gq * Lgq + G1qg * Lqg
*
      return
      end
