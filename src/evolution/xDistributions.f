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
**
*     Input Variables
*
      double precision x
      double precision Q
**
*     Internal Variables
*
      integer m,j
      integer ip,iQ
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
*     Coupling at the final scale
*
      aq = aQED(Q**2d0)
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
         call NDistributions(s,Q,xfN)
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
