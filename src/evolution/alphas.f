****************************************************************
*
*     alphas.f
*     
*     Routine that returns alpha_s/(4*pi).
*
*     Note that the only argument of routine AQCD is the energy Q2
*     of the process which is convereted into the renormalization
*     scale inside the routine itself.
*
****************************************************************
      function aQCD(q2)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/alphas.h"
      include "../commons/massthrs.h"
      include "../commons/hqmass.h"
      include "../commons/nffn.h"
      include "../commons/nfmax.h"
      include "../commons/renfacscales.h"
      include "../commons/ns.h"
      include "../commons/modev.h"
      include "../commons/ipt.h"
**
*     Input Variables
*
      double complex q2
**
*     Internal Variables
*
      integer i
      integer nfi,nff
      integer dnf,snf
      double precision c1,c2,kappa,ln
      double complex mur2th(4:6)
      double complex mur2,mur20
      double complex asi,asr0
      double complex as_expanded_MELA,as_exact_MELA
      external as_expanded_MELA,as_exact_MELA
**
*     Output Variables
*
      double complex aQCD
*
      asr0 = asref / 4d0 / pi
*
      mur20 = q2ref
      mur2  = krf**2d0 * q2
      do i=4,6
         mur2th(i) = krf**2d0 * q2th(i)
c         mur2th(i) = q2th(i)
      enddo
*
      if(ns.eq."FFNS")then
         nfi = nffn
         nff = nffn
      elseif(ns.eq."VFNS")then
         if(abs(mur2).ge.abs(mur2th(6)))then
            nff = 6
         elseif(abs(mur2).ge.abs(mur2th(5)))then
            nff = 5
         elseif(abs(mur2).ge.abs(mur2th(4)))then
            nff = 4
         else
            nff = 3
         endif
         if(nff.gt.nfmax) nff = nfmax
*
         if(abs(mur20).gt.abs(mur2th(6)))then
            nfi = 6
         elseif(abs(mur20).gt.abs(mur2th(5)))then
            nfi = 5
         elseif(abs(mur20).gt.abs(mur2th(4)))then
            nfi = 4
         else
            nfi = 3
         endif
      endif
*
*     c1 and c2 are the same coefficients used in eq. (2.42) of hep-ph/0408244 and 
*     obtained in eq. (10) of hep-ph/9706430. In the following they are divided by 
*     (4*pi) and (4*pi)^2 respectively to match the notations. Note that in terms 
*     of the MSbar mass this coefficients change.
*
      kappa = krf**2d0     ! mu_R / mu_F
c      kappa = 1d0      ! mu_R / mu_F
      ln = dlog(kappa)
*     Pole Mass
      if(hqmass.eq.0)then
         if(nff.gt.nfi)then
            c1 = 2d0 / 3d0 * ln
            c2 = 4d0 / 9d0 * ln**2d0 + 38d0 / 3d0 * ln + 14d0 / 3d0
         elseif(nff.lt.nfi)then
            c1 = - 2d0 / 3d0 * ln
            c2 = 4d0 / 9d0 * ln**2d0 - 38d0 / 3d0 * ln - 14d0 / 3d0
         endif
*     MSbar mass
      elseif(hqmass.eq.1)then
         if(nff.gt.nfi)then
            c1 = 0d0
            c2 = - 22d0 / 9d0
         elseif(nff.lt.nfi)then
            c1 = 0d0
            c2 = 22d0 / 9d0
         endif
      endif
*
 10   if(nff.eq.nfi) then
         if(modev.eq."TRN".or.modev.eq."GFN")then
            aQCD = as_expanded_MELA(nfi,mur20,asr0,mur2,ipt)
         elseif(modev.eq."ITE".or.modev.eq."PTH")then
            aQCD = as_exact_MELA(nfi,mur20,asr0,mur2,ipt)
         endif
         return
      else
         if(nff.gt.nfi)then
            dnf = 1
            snf = 1
         else
            dnf = -1
            snf = 0
         endif
*
         if(modev.eq."TRN".or.modev.eq."GFN")then
            asi = as_expanded_MELA(nfi,mur20,asr0,mur2th(nfi+snf),ipt)
         elseif(modev.eq."ITE".or.modev.eq."PTH")then
            asi = as_exact_MELA(nfi,mur20,asr0,mur2th(nfi+snf),ipt)
         endif
*
*     NLO and NNLO threshold matchings
*
         if(ipt.eq.1)then
            asi = asi * ( 1d0 + c1 * asi )
         elseif(ipt.eq.2)then
            asi = asi * ( 1d0 + c1 * asi + c2 * asi**2d0 )
         endif
*
         asr0  = asi
         mur20 = mur2th(nfi+snf)
         nfi  = nfi + dnf
         goto 10
      endif
*
      end
*
****************************************************************
*
*     Routines for the computation of alpha_s with fixed number of 
*     flavours
* 
*     - as_expanded_MELA: computes alpha_s aQCD function of alpha_s at a given
*                    refernce scale.
*
*     - as_exact_MELA: Exact solution of the QCD beta function 
*                 equation using fourth order Runge-Kutta 
*                 algorithm.
*                 (Used in PEGASUS and Les Houches tables)
*      
****************************************************************
      function as_expanded_MELA(nf,mu20,as0,mu2,ipt)
*
      implicit none
*
      include "../commons/beta.h"
**
*     Input Variables
*
      integer nf,ipt
      double complex mu2
      double complex mu20,as0
**
*     Internal Variables
*
      double complex asi
      double complex alo,t,aQCD,den
**
*     Output Variables
*
      double complex as_expanded_MELA
*
      asi = as0
      t = zlog(mu2/mu20)
      den = 1d0 + beta0(nf) * asi * t
      alo = asi / den
*
*     LO
*
      aQCD = alo
*
*     NLO
*
      if(ipt.ge.1)then
         aQCD = alo * ( 1d0 - b1(nf) * alo * zlog(den) )
      endif
*
*     NNLO
*
      if(ipt.eq.2)then
         aQCD = alo * ( 1d0 
     1      + ( alo * ( alo - asi ) * ( b2(nf) - b1(nf)**2d0 )
     2      + aQCD * b1(nf) * zlog(aQCD/asi) ) )
      endif
*
      as_expanded_MELA = aQCD
*
      return
      end
*
****************************************************************************
      FUNCTION AS_EXACT_MELA(NF,MU20,AS0,MU2,IPT)
*
      IMPLICIT NONE
*
      include "../commons/beta.h" 
      include "../commons/consts.h"
**
*     Input Variables
*
      INTEGER NF,IPT
      DOUBLE COMPLEX MU2
      DOUBLE COMPLEX AS0,MU20
**
*     Inernal Variables
*
      INTEGER NSTEP,K1
      DOUBLE COMPLEX AQCD
      DOUBLE PRECISION SXTH
      DOUBLE COMPLEX FBETAMELA
      DOUBLE COMPLEX XK0,XK1,XK2,XK3
      DOUBLE COMPLEX DLR,LRRAT
      
      PARAMETER(NSTEP=10)
      PARAMETER(SXTH=0.166666666666666D0 )
**
*     Output Variables
*
      DOUBLE COMPLEX AS_EXACT_MELA
*
*     ..The beta functions FBETAMELAn at N^nLO for n = 1, 2, and 3
*
      AQCD    = AS0
      LRRAT = ZLOG (MU2/MU20)
      DLR   = LRRAT / NSTEP
*     
*     ..Solution of the evolution equation depending on  NAORD
*   (fourth-order Runge-Kutta beyond the leading order)
*     
      IF(IPT.EQ.0)THEN
         AQCD = AS0 / ( 1D0 + BETA0(NF) * AS0 * LRRAT )
      ELSEIF(IPT.EQ.1)THEN
         DO 2 K1=1,NSTEP
            XK0 = DLR * FBETAMELA(AQCD,NF,IPT)
            XK1 = DLR * FBETAMELA(AQCD + 0.5d0 * XK0,NF,IPT)
            XK2 = DLR * FBETAMELA(AQCD + 0.5d0 * XK1,NF,IPT)
            XK3 = DLR * FBETAMELA(AQCD + XK2,NF,IPT)
            AQCD = AQCD + SXTH * ( XK0 + 2D0 * XK1 + 2d0 * XK2 + XK3 )
 2        CONTINUE
      ELSEIF(IPT.EQ.2)THEN
         DO 3 K1 = 1, NSTEP
            XK0 = DLR * FBETAMELA(AQCD,NF,IPT)
            XK1 = DLR * FBETAMELA(AQCD + 0.5d0 * XK0,NF,IPT)
            XK2 = DLR * FBETAMELA(AQCD + 0.5d0 * XK1,NF,IPT)
            XK3 = DLR * FBETAMELA(AQCD + XK2,NF,IPT)
            AQCD = AQCD + SXTH * ( XK0 + 2D0 * XK1 + 2D0 * XK2 + XK3 )
 3       CONTINUE
      ENDIF
*
      AS_EXACT_MELA = AQCD
*
      RETURN
      END
*
****************************************************************************
      function fbetaMELA(a,nf,ipt)
*
      implicit none
*
      include "../commons/beta.h"
*
      double complex fbetaMELA,a
      integer nf,ipt
*
      if(ipt.eq.1)then
         fbetaMELA = - A**2d0 * ( BETA0(NF) + A * BETA1(NF) )
      elseif(ipt.eq.2)then
         fbetaMELA = - A**2d0 * ( BETA0(NF) + A * ( BETA1(NF)
     1             + A * BETA2(NF) ) )
      endif
*
      return
      end
