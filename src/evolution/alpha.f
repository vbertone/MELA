****************************************************************
*
*     alpha.f
*     
*     Routine that returns alpha_s/(4*pi).
*
*     Note that the only argument of routine AQED is the energy Q2
*     of the process which is convereted into the renormalization
*     scale inside the routine itself.
*
****************************************************************
      function aQED(q2)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/alpha.h"
      include "../commons/massthrs.h"
      include "../commons/nffn.h"
      include "../commons/nfmax.h"
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
      double complex mur2th(3)
      double complex mur2,mur20
      double complex ai,ar0
      double complex a_expanded_MELA,a_exact_MELA
      external a_expanded_MELA,a_exact_MELA
**
*     Output Variables
*
      double complex aQED
*
      ar0 = aref / 4d0 / pi
*
      mur20 = q2ref
      mur2  = q2
      do i=1,3
         mur2th(i) = q2th(i)
      enddo
*
      if(ns.eq."FFNS")then
         nfi = nffn
         nff = nffn
      elseif(ns.eq."VFNS")then
         if(abs(mur2).ge.abs(mur2th(3)))then
            nff = 3
         elseif(abs(mur2).ge.abs(mur2th(2)))then
            nff = 2
         elseif(abs(mur2).ge.abs(mur2th(1)))then
            nff = 1
         else
            nff = 3
         endif
         if(nff.gt.nfmax) nff = nfmax
*
         if(abs(mur20).gt.abs(mur2th(3)))then
            nfi = 3
         elseif(abs(mur20).gt.abs(mur2th(2)))then
            nfi = 2
         elseif(abs(mur20).gt.abs(mur2th(1)))then
            nfi = 1
         else
            nfi = 3
         endif
      endif
*
 10   if(nff.eq.nfi) then
         if(modev.eq."TRN".or.modev.eq."GFN")then
            aQED = a_expanded_MELA(nfi,mur20,ar0,mur2,ipt)
         elseif(modev.eq."ITE".or.modev.eq."PTH")then
            aQED = a_exact_MELA(nfi,mur20,ar0,mur2,ipt)
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
            ai = a_expanded_MELA(nfi,mur20,ar0,mur2th(nfi+snf),ipt)
         elseif(modev.eq."ITE".or.modev.eq."PTH")then
            ai = a_exact_MELA(nfi,mur20,ar0,mur2th(nfi+snf),ipt)
         endif
*
         ar0  = ai
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
*     - a_expanded_MELA: computes alpha_s aQED function of alpha_s at a given
*                    refernce scale.
*
*     - a_exact_MELA: Exact solution of the QED beta function 
*                 equation using fourth order Runge-Kutta 
*                 algorithm.
*                 (Used in PEGASUS and Les Houches tables)
*      
****************************************************************
      function a_expanded_MELA(nf,mu20,a0,mu2,ipt)
*
      implicit none
*
      include "../commons/beta.h"
**
*     Input Variables
*
      integer nf,ipt
      double complex mu2
      double complex mu20,a0
**
*     Internal Variables
*
      double complex ai
      double complex alo,t,aQED,den
**
*     Output Variables
*
      double complex a_expanded_MELA
*
      ai = a0
      t = zlog(mu2/mu20)
      den = 1d0 + beta0(nf) * ai * t
      alo = ai / den
*
*     LO
*
      aQED = alo
*
*     NLO
*
      if(ipt.ge.1)then
         aQED = alo * ( 1d0 - b1(nf) * alo * zlog(den) )
      endif
*
      a_expanded_MELA = aQED
*
      return
      end
*
****************************************************************************
      FUNCTION A_EXACT_MELA(NF,MU20,A0,MU2,IPT)
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
      DOUBLE COMPLEX A0,MU20
**
*     Inernal Variables
*
      INTEGER NSTEP,K1
      DOUBLE COMPLEX AQED
      DOUBLE PRECISION SXTH
      DOUBLE COMPLEX FBETAMELA
      DOUBLE COMPLEX XK0,XK1,XK2,XK3
      DOUBLE COMPLEX DLR,LRRAT
      
      PARAMETER(NSTEP=10)
      PARAMETER(SXTH=0.166666666666666D0 )
**
*     Output Variables
*
      DOUBLE COMPLEX A_EXACT_MELA
*
*     ..The beta functions FBETAMELAn at N^nLO for n = 1, 2, and 3
*
      AQED  = A0
      LRRAT = ZLOG (MU2/MU20)
      DLR   = LRRAT / NSTEP
*     
*     ..Solution of the evolution equation depending on  NAORD
*   (fourth-order Runge-Kutta beyond the leading order)
*     
      IF(IPT.EQ.0)THEN
         AQED = A0 / ( 1D0 + BETA0(NF) * A0 * LRRAT )
      ELSEIF(IPT.GE.1)THEN
         DO 2 K1=1,NSTEP
            XK0 = DLR * FBETAMELA(AQED,NF,IPT)
            XK1 = DLR * FBETAMELA(AQED + 0.5d0 * XK0,NF,IPT)
            XK2 = DLR * FBETAMELA(AQED + 0.5d0 * XK1,NF,IPT)
            XK3 = DLR * FBETAMELA(AQED + XK2,NF,IPT)
            AQED = AQED + SXTH * ( XK0 + 2D0 * XK1 + 2d0 * XK2 + XK3 )
 2        CONTINUE
      ENDIF
*
      A_EXACT_MELA = AQED
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
      if(ipt.eq.0)then
         fbetaMELA = - A**2d0 * BETA0(NF)
      elseif(ipt.eq.1)then
         fbetaMELA = - A**2d0 * ( BETA0(NF) + A * BETA1(NF) )
      endif
*
      return
      end
