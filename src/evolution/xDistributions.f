************************************************************************
*
*     xDistributions.f:
*
*     Routine that returns the evolved x-space distributions.
*
************************************************************************
      subroutine xDistributions(x,nQ,Q,xfph)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/alpha.h"
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
      double precision r
      double precision theta,sigma,t
      double complex aQED
      double complex s,tmp(7)
      double complex xfN(7)
      double complex xfev(7)
**
*     Output Variables
*
      double complex xfph(-3:3)
*
*     Precompute alpha
*
      aEvolIni(1) = aQED(Q(1)**2d0)
      aEvolFin(1) = aQED(Q(2)**2d0)
      do iQ=2,nQ-1
         aEvolIni(iQ) = aEvolFin(iQ-1)
         aEvolFin(iQ) = aQED(Q(iQ+1)**2d0)
      enddo
*
      t = - dlog(x)
      r = 2d0 * 43 / 5d0 / t
*
      do ip=1,7
         tmp(ip) = (0d0,0d0)
      enddo
      m = 83                  ! Must be odd
      do j=1,m-1
         theta = - pi + dble(j) * ( 2d0 * pi / dble(m) )
         sigma = theta + (theta/tan(theta)-1d0)/tan(theta)
         s     = r * theta * dcmplx(1d0/tan(theta),1d0) + 1d0
*
         call NDistributions(s,nQ,Q,xfN)
         do ip=1,7
            tmp(ip) = tmp(ip)
     1              + exp( t * s ) * dcmplx(1d0,sigma) * xfN(ip)
         enddo
      enddo
*
      do ip=1,7
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
*
************************************************************************
      subroutine xDistributionsReal(x,Q0,Q,xfph)
*
      implicit none
**
*     Input Variables
*
      double precision x
      double precision Q0, Q
**
*     Internal Variables
*
      integer i
      double complex Qv(100)
      double complex xfphc(-3:3)
**
*     Output Variables
*
      double precision xfph(-3:3)
*
*     Precompute alpha
*
      Qv(1) = dcmplx(Q0, 0d0)
      Qv(2) = dcmplx(Q, 0d0)
      call xDistributions(x, 2, Qv, xfphc)
*
      do i=-3,3
         xfph(i) = dble(xfphc(i))
      enddo
*
      return
      end
