************************************************************************
*
*     xDistributions.f:
*
*     Routine that returns the evolved x-space distributions.
*
************************************************************************
      subroutine xDistributions(x,Q,xfph)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/alpha.h"
      include "../commons/massthrs.h"
      include "../commons/ns.h"
      include "../commons/nffn.h"
      include "../commons/nfmax.h"
**
*     Input Variables
*
      double precision x
      double precision Q
**
*     Internal Variables
*
      integer m,j,nfi,nff
      integer ip
      double precision r
      double precision theta,sigma,t
      double precision aQED
      double precision xfev(19)
      double complex s,tmp(19)
      double complex xfN(19)
**
*     Output Variables
*
      double precision xfph(-9:9)
*
*     Determine number of active flavours at the final scale
*
      if(ns.eq."FFNS")then
         nfi = nffn
         nff = nffn
      elseif(ns.eq."VFNS")then
         nfi = 1
         do nff = 1, 9
            if (Q**2.ge.q2th(nff)) exit
         enddo
         if(nff.gt.nfmax) nff = nfmax
      endif
*
*     Coupling at the final scale
*
      ath(nff+1) = aQED(Q**2d0)
*
      t = - dlog(x)
      r = 2d0 * 43 / 5d0 / t
*
      do ip=1,19
         tmp(ip) = (0d0,0d0)
      enddo
      m = 83                  ! Must be odd
      do j=1,m-1
         theta = - pi + dble(j) * ( 2d0 * pi / dble(m) )
         sigma = theta + (theta/tan(theta)-1d0)/tan(theta)
         s     = r * theta * dcmplx(1d0/tan(theta),1d0) + 1d0
*
         call NDistributions(s,nfi,nff,xfN)
         do ip=1,19
            tmp(ip) = tmp(ip)
     1              + exp( t * s ) * dcmplx(1d0,sigma) * xfN(ip)
         enddo
      enddo
*
      do ip=1,19
         xfev(ip) = x * r * dble(tmp(ip)) / m
      enddo
*
*     Rotate back to the physical basis
*
      call evln2lha(xfev,xfph)
*
      return
      end
