************************************************************************
*
*     SIDISxStructureFunctions.f:
*
*     Routine that returns the x-space semi-inclusive structure
*     functions.
*
************************************************************************
      subroutine SIDISxStructureFunctions(x, z, Q0, Q, SFx)
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
      integer isf,ip,iQ
      double precision rmax
      double precision rx, thetax, sigmax(32), tx
      double precision rz, thetaz, sigmaz(32), tz
      double complex sx(32), sz(32)
      double complex aQCD
      double complex tmp1(2), tmp2(2)
      double complex SFN(2)
**
*     Output Variables
*
      double precision SFx(2)
*
      asDIS  = aQCD((kfacQ*Q)**2d0)
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
      do ip=1,2
         tmp1(ip) = (0d0, 0d0)
      enddo
      do k=1,m-1
         do ip = 1, 2
            tmp2(ip) = (0d0,0d0)
         enddo
         do j=1,m-1
            call SIDISNStructureFunctions(sx(j), sz(k), Q0, Q, SFN)
            do ip=1,2
               tmp2(ip) = tmp2(ip)
     1              + exp( tx * sx(j) ) * dcmplx(1d0,sigmax(j))
     2              * SFN(ip)
            enddo
         enddo
         do ip=1,2
            tmp2(ip) = x * rx * tmp2(ip) / m
         enddo
         do ip=1,2
            tmp1(ip) = tmp1(ip)
     1           + exp( tz * sz(k) ) * dcmplx(1d0,sigmaz(k))
     2           * tmp2(ip)
         enddo
      enddo
      do ip=1,2
         SFx(ip) = z * rz * tmp1(ip) / m
      enddo
*
      do ip=1,2
         if(abs(SFx(ip)).lt.1d-12) SFx(ip) = (0d0,0d0)
      enddo
*
      return
      end
*
************************************************************************
      double precision function F2LO(x, z, Q0, Q)
*
      implicit none
*
      include "../commons/evol.h"
*
      integer k
      double precision x, z
      double complex Q0, Q, Qv(2)
      double complex xf(-6:6), xd(-6:6)
      double complex bq(6), dq(6), ch
*
      Qv(1) = Q0
      Qv(2) = Q
      evol = "SPACE"
      call xDistributions(x, 2, Qv ,xf)
      evol = "TIME"
      call xDistributions(z, 2, Qv ,xd)
      call ComputeChargesDIS_MELA(Q**2, bq, dq)
      F2LO = (0d0, 0d0)
      do k = 1, 6
         ch = bq(k)
         if (k.eq.1) ch = bq(2)
         if (k.eq.2) ch = bq(1)
         F2LO = F2LO + ch * ( xd(-k) * xf(-k) + xd(k) * xf(k) )
      enddo
*
      return
      end

