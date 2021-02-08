************************************************************************
*
*     SIDISNStructureFunctions.f:
*
*     Routine that returns the semi-inclusive structure functions in
*     N-space.
*
************************************************************************
      subroutine SIDISNStructureFunctions(N, M, Q0, Q, SFN)
*
      implicit none
*
      include "../commons/evol.h"
      include "../commons/renfacscales.h"
      include "../commons/alphas.h"
      include "../commons/ipt.h"
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
      double complex C21QQ_SIDIS, C21GQ_SIDIS, C21QG_SIDIS
      double complex CL1QQ_SIDIS, CL1GQ_SIDIS, CL1QG_SIDIS
      double complex C2qq, C2gq, C2qg
      double complex CLqq, CLgq, CLqg
      double complex Lqq, Lgq, Lqg
      double complex ch
      character*5 evolbkp
**
*     Output Variables
*
      double complex SFN(2)
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
      evol = "SPACE"
      call NDistributions(N, 2, Qv, xfevN)
      call evln2lhac(xfevN, xfphN)
*
*     Call FFs at the final scale and rotate them into the physical
*     basis
*
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
      C2qq = (1d0, 0d0)
      C2gq = (0d0, 0d0)
      C2qg = (0d0, 0d0)
      CLqq = (0d0, 0d0)
      CLgq = (0d0, 0d0)
      CLqg = (0d0, 0d0)
      if (ipt.ge.1) then
         C2qq = C2qq + asDIS * C21QQ_SIDIS(N, M)
         C2gq = C2gq + asDIS * C21GQ_SIDIS(N, M)
         C2qg = C2qg + asDIS * C21QG_SIDIS(N, M)
         CLqq = CLqq + asDIS * CL1QQ_SIDIS(N, M)
         CLgq = CLgq + asDIS * CL1GQ_SIDIS(N, M)
         CLqg = CLqg + asDIS * CL1QG_SIDIS(N, M)
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
      SFN(1) = C2qq * Lqq + C2gq * Lgq + C2qg * Lqg
      SFN(2) = CLqq * Lqq + CLgq * Lgq + CLqg * Lqg
*
      return
      end
