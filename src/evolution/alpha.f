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
      include "../commons/ipt.h"
**
*     Input Variables
*
      double precision q2
**
*     Internal Variables
*
      integer nfi,nff
      integer dnf,snf
      double precision mur20
      double precision ai,ar0
      double precision a_exact_MELA
      external a_exact_MELA
**
*     Output Variables
*
      double precision aQED
*
      ar0 = aref / 4d0 / pi
*
      mur20 = q2ref
*
      if(ns.eq."FFNS")then
         nfi = nffn
         nff = nffn
      elseif(ns.eq."VFNS")then
         do nff=1,9
            if (q2.le.q2th(nff)) exit
         enddo
         if(nff.gt.nfmax) nff = nfmax

         do nfi=1,9
            if (mur20.le.q2th(nfi)) exit
         enddo
      endif
*
 10   if(nff.eq.nfi) then
         aQED = a_exact_MELA(nfi,mur20,ar0,q2,ipt)
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
         ai = a_exact_MELA(nfi,mur20,ar0,q2th(nfi+snf),ipt)
*
         ar0   = ai
         mur20 = q2th(nfi+snf)
         nfi   = nfi + dnf
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
*     - a_exact_MELA: Exact solution of the QED beta function equation
*                 using fourth order Runge-Kutta algorithm. (Used in
*                 PEGASUS and Les Houches tables)
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
      DOUBLE PRECISION MU2
      DOUBLE PRECISION A0,MU20
**
*     Inernal Variables
*
      INTEGER NSTEP,K1
      DOUBLE PRECISION AQED
      DOUBLE PRECISION SXTH
      DOUBLE PRECISION FBETAMELA
      DOUBLE PRECISION XK0,XK1,XK2,XK3
      DOUBLE PRECISION DLR,LRRAT
      
      PARAMETER(NSTEP=10)
      PARAMETER(SXTH=0.166666666666666D0)
**
*     Output Variables
*
      DOUBLE PRECISION A_EXACT_MELA
*
*     ..The beta functions FBETAMELAn at N^nLO for n = 1, 2, and 3
*
      AQED  = A0
      LRRAT = DLOG (MU2/MU20)
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
      double precision fbetaMELA,a
      integer nf,ipt
*
      if(ipt.eq.0)then
         fbetaMELA = - A**2 * BETA0(NF)
      elseif(ipt.eq.1)then
         fbetaMELA = - A**2 * ( BETA0(NF) + A * BETA1(NF) )
      endif
*
      return
      end
