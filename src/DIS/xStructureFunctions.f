************************************************************************
*
*     xStructureFunctions.f:
*
*     Routine that returns the x-space structure functions.
*
************************************************************************
      subroutine xStructureFunctions(x,nQ,Q,SFx)
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
      integer nQ
      double precision x
      double complex Q(100)
**
*     Internal Variables
*
      integer m,j
      integer isf,ip,iQ
      integer iptbkp
      integer maxfl
      double precision r,rmax
      double precision theta,sigma,t
      double precision rho(6),W2,kappa
      double complex aQCD
      double complex s,tmp(3,6)
      double complex SFN(3,6)
**
*     Output Variables
*
      double complex SFx(3,0:6)
*
*     Precompute alphas
*
      maxfl = 3
      if(ns.eq."FFNS")then
         iptbkp = ipt
         ipt    = max(0,ipt-1)
         asDIS = aQCD((krenQ*krf*Q(nQ))**2d0)
         asEvolIni(1) = aQCD(krf**2 * Q(1)**2d0)
         asEvolFin(1) = aQCD(kfacQ**2 * krf**2 * Q(2)**2d0)
         do iQ=2,nQ-1
            asEvolIni(iQ) = asEvolFin(iQ-1)
            asEvolFin(iQ) = aQCD(kfacQ**2 * krf**2 * Q(iQ+1)**2d0)
         enddo
         ipt    = iptbkp
         W2     = abs(Q(nq)**2d0) * ( 1d0 - x ) / x
         kappa  = 4d0
         do ip=4,6
            if(W2.ge.kappa*abs(q2th(ip))) maxfl = maxfl + 1 
         enddo
      elseif(ns.eq."VFNS")then
         asDIS = aQCD((krenQ*krf*Q(nQ))**2d0)
         asEvolIni(1) = aQCD(krf**2 * Q(1)**2d0)
         asEvolFin(1) = aQCD(kfacQ**2 * krf**2 * Q(2)**2d0)
         do iQ=2,nQ-1
            asEvolIni(iQ) = asEvolFin(iQ-1)
            asEvolFin(iQ) = aQCD(kfacQ**2 * krf**2 * Q(iQ+1)**2d0)
         enddo
         do ip=4,6
            if(abs(Q(nq)**2d0).ge.abs(q2th(ip))) maxfl = maxfl + 1 
         enddo
      endif
*
*     Couplings
*
      call ComputeChargesDIS_MELA(Q(nQ)**2d0,bq,dq)
*
*     Threasholds
*
      do ip=1,6
         rho(ip) = 1d2
      enddo
*
      if(ns.eq."FFNS")then
         do ip=4,6
            rho(ip) = 1d0 / ( 1d0
     1              + 4d0 * abs( q2th(ip) / Q(nQ) / Q(nQ) ) )
         enddo
*     Initialize Adler coefficients
         call adlersr_coeff(Q(nQ))
      endif
*
      t = - dlog(x)
      m = 17                  ! Must be odd
      r = 2d0 * m / 5d0 / t
*
      rmax = 10d0
*
      if(r.gt.rmax) r = rmax
*
      do isf=1,3
         do ip=1,6
            tmp(isf,ip) = (0d0,0d0)
         enddo
      enddo
      do j=1,m-1
         theta = - pi + dble(j) * ( 2d0 * pi / dble(m) )
         sigma = theta + (theta/tan(theta)-1d0)/tan(theta)
         s     = r * theta * dcmplx(1d0/tan(theta),1d0) + 1d0
*
         call NStructureFunctions(s,nQ,Q,SFN)
         do isf=1,3
            do ip=1,maxfl
               tmp(isf,ip) = tmp(isf,ip)
     1              + exp( t * s ) * dcmplx(1d0,sigma) * SFN(isf,ip)
            enddo
         enddo
      enddo
*
*     F2, FL, xF3
*
      do isf=1,3
         SFx(isf,0) = (0d0,0d0)
         do ip=1,6
            SFx(isf,ip) = x * r * tmp(isf,ip) / m
            if(x.ge.rho(ip)) SFx(isf,ip) = (0d0,0d0)
            if(abs(SFx(isf,ip)).lt.1d-12) SFx(isf,ip) = (0d0,0d0)
            SFx(isf,0) = SFx(isf,0) + SFx(isf,ip)
         enddo
      enddo
*
      return
      end
*
************************************************************************
      subroutine xStructureFunctionsReal(x,Q0,Q,SFx)
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
      integer i,isf
      double complex Qv(100)
      double complex SFxc(3,0:6)
**
*     Output Variables
*
      double precision SFx(3)
*
*     Precompute alphas
*
      Qv(1) = dcmplx(Q0, 0d0)
      Qv(2) = dcmplx(Q, 0d0)
      call xStructureFunctions(x, 2, Qv, SFxc)
*
      do isf=1,3
         SFx(isf) = dreal(SFxc(isf, 0))
      enddo
*
      return
      end
