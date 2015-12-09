************************************************************************
*
*     xDistibutions.f:
*
*     Routine that returns the evolved x-space distributions.
*
************************************************************************
      subroutine xDistibutions(x,nQ,Q,xfph)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/alphas.h"
**
*     Input Variables
*
      integer nQ
      double precision x
      double complex Q(100)
**
*     Internal Variables
*
      integer m,j
      integer ip,iQ
      double precision r,rmax
      double precision theta,sigma,t
      double complex as
      double complex s,tmp(13)
      double complex xfN(13)
      double complex xfev(13)
**
*     Output Variables
*
      double complex xfph(-6:6)
*
*     Precompute alphas
*
      asEvolIni(1) = as(Q(1)**2d0)
      asEvolFin(1) = as(Q(2)**2d0)
      do iQ=2,nQ-1
         asEvolIni(iQ) = asEvolFin(iQ-1)
         asEvolFin(iQ) = as(Q(iQ+1)**2d0)
      enddo
*
      t = - dlog(x)
      m = 33                  ! Must be odd
      r = 2d0 * m / 5d0 / t
*
      rmax = 10d0
*
      if(r.gt.rmax) r = rmax
*
      do ip=1,13
         tmp(ip) = (0d0,0d0)
      enddo
      do j=1,m-1
         theta = - pi + dble(j) * ( 2d0 * pi / dble(m) )
         sigma = theta + (theta/tan(theta)-1d0)/tan(theta)
         s     = r * theta * dcmplx(1d0/tan(theta),1d0) + 1d0
*
         call NDistributions(s,nQ,Q,xfN)
         do ip=1,13
            tmp(ip) = tmp(ip)
     1              + exp( t * s ) * dcmplx(1d0,sigma) * xfN(ip)
         enddo
      enddo
*
      do ip=1,13
         xfev(ip) = x * r * tmp(ip) / m
         if(abs(xfev(ip)).lt.1d-12) xfev(ip) = (0d0,0d0)
      enddo
*
*     Rotate back to the physical basis
*
      call evln2lhac(xfev,xfph)
*
      return
      end