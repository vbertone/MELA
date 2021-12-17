************************************************************************
*
*     NDistributions.f:
*
*     Routine that returns the evolved N-space distributions.
*
************************************************************************
      subroutine NDistributions(N,nfi,nff,fn)
*
      implicit none
*
      include "../commons/alpha.h"
      include "../commons/massthrs.h"
      include "../commons/renscheme.h"
**
*     Input Variables
*
      integer nfi,nff
      double complex N
**
*     Internal Variables
*
      integer ifl
      double complex evf(19,19)
**
*     Output Variables
*
      double complex fn(19)
*
*     Call N-space distributions
*
      call electronPDFsn(N, fn) ! Electron PDFs
*
*     Run over sub-intervals
*
      do ifl = nfi, nff
*
*     Call evolution kernels
*
         if(renscheme.eq."MSBAR")then
            call path_ordering(N,ath(ifl),ath(ifl+1),ifl,evf)
         elseif(renscheme.eq."FIXED")then
            call alpha_fixed(N,q2th(ifl),q2th(ifl+1),ifl,evf)
         elseif(renscheme.eq."ALPMZ")then
            call alphaMZ_analytic(N,q2th(ifl),q2th(ifl+1),ifl,evf)
         endif
*
*     Convolute with vector of PDFs
*
         call mvmult(evf,19,19,fn,19)
      enddo
*
      return
      end
