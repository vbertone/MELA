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
      include "../commons/tecparam.h"
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
      double precision range,perc
      integer nintstep
      integer nminstep
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
*     Adjust the number of iterations in the path ordering
*     (if there are thresholds, we need more points where
*     the range of variation is larger).
*     No less than NMINSTEP points for step.
         nminstep = 50
*
         if(renscheme.eq."MSBAR")then
            range = ath(nff+1)-ath(nfi)
            perc  = dble(nint) * (ath(ifl+1)-ath(ifl))/range
         else
            range = q2th(nff+1)-q2th(nfi)
            perc  = dble(nint) * (q2th(ifl+1)-q2th(ifl))/range
         endif
         if (perc.lt.nminstep) then
            nintstep = nminstep
         else
            nintstep = int(perc)
         endif
c         write(6,*) "ifl",ifl,"nintstep",nintstep
*         
*     Call evolution kernels
*
         if(renscheme.eq."MSBAR")then
            call path_ordering(N,ath(ifl),ath(ifl+1),ifl,evf,nintstep)
         elseif(renscheme.eq."FIXED")then
            call alpha_fixed(N,q2th(ifl),q2th(ifl+1),ifl,evf)
         elseif(renscheme.eq."ALPMZ")then
            call alphaMZ_pathordering(N,q2th(ifl),q2th(ifl+1),
     .           ifl,evf,nintstep) 
c            call alphaMZ_analytic(N,q2th(ifl),q2th(ifl+1),ifl,evf)
         endif
*
*     Convolute with vector of PDFs
*
         call mvmult(evf,19,19,fn,19)
      enddo
*
      return
      end
