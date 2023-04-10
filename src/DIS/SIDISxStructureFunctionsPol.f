************************************************************************
*
*     SIDISxStructureFunctionsPol.f:
*
*     Routine that returns the x-space semi-inclusive polarised
*     structure functions.
*
************************************************************************
      subroutine SIDISxStructureFunctionsPol(x, z, Q0, Q, SFx)
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/renfacscales.h"
      include "../commons/consts.h"
      include "../commons/ns.h"
      include "../commons/massthrs.h"
      include "../commons/ipt.h"
**
*     Input Variables
*
      double precision x, z
      double precision Q0, Q
**
*     Internal Variables
*
      integer m, j, k
      integer iQ
      double precision rmax
      double precision rx, thetax, sigmax(32), tx
      double precision rz, thetaz, sigmaz(32), tz
      double complex sx(32), sz(32)
      double complex aQCD
      double complex tmp1, tmp2
      double complex SFN
**
*     Output Variables
*
      double precision SFx
*
      asDIS = aQCD((kfacQ*Q)**2d0)
      asEvolIni(1) = aQCD(Q0**2d0)
      asEvolFin(1) = aQCD((kfacQ*Q)**2d0)
*
*     Couplings
*
      call ComputeChargesDIS_MELA(dcmplx(Q**2d0, 0d0), bq, dq)
*
      m = 33                  ! Must be odd
      rmax = 10d0
*
      tx = - dlog(x)
      rx = 2d0 * m / 5d0 / tx
      if(rx.gt.rmax) rx = rmax
*
      tz = - dlog(z)
      rz = 2d0 * m / 5d0 / tz
      if(rz.gt.rmax) rz = rmax
*
      do j=1,m-1
         thetax    = - pi + dble(j) * ( 2d0 * pi / dble(m) )
         sigmax(j) = thetax + (thetax/tan(thetax)-1d0)/tan(thetax)
         sx(j)     = rx * thetax * dcmplx(1d0/tan(thetax),1d0) + 1d0
      enddo
      do k=1,m-1
         thetaz    = - pi + dble(k) * ( 2d0 * pi / dble(m) )
         sigmaz(k) = thetaz + (thetaz/tan(thetaz)-1d0)/tan(thetaz)
         sz(k)     = rz * thetaz * dcmplx(1d0/tan(thetaz),1d0) + 1d0
      enddo
*
      tmp1 = (0d0, 0d0)
      do k=1,m-1
         tmp2 = (0d0,0d0)
         do j=1,m-1
            call SIDISNStructureFunctions(sx(j), sz(k), Q0, Q, SFN)
            tmp2 = tmp2
     1           + exp( tx * sx(j) ) * dcmplx(1d0,sigmax(j))
     2           * SFN
         enddo
         tmp2 = x * rx * tmp2 / m
         tmp1 = tmp1
     1        + exp( tz * sz(k) ) * dcmplx(1d0,sigmaz(k))
     2        * tmp2
      enddo
      SFx = z * rz * dreal(tmp1) / dble(m)
*
      if(abs(SFx).lt.1d-12) SFx = 0d0
*
      return
      end
