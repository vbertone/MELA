***********************************************************************
*
*     gfuncs.f:
*
*     This subroutine returns the QCD evolution factors, in N space, from
*     Q2I to Q2F with NF active flavours, for Non Singlet and Singlet-Gluon
*     using the g-function approach.
*
***********************************************************************
      SUBROUTINE GFUNCS(ZN,ASI,QI,QF,NF,EFNNS,EFNSG)
*
      IMPLICIT NONE
*
      include "../commons/ipt.h"
      include "../commons/renfacscales.h"
      include "../commons/beta.h"
      include "../commons/evol.h"
      include "../commons/pol.h"
**
*     Input Variables
*
      INTEGER NF
      DOUBLE COMPLEX ASI
      DOUBLE COMPLEX QI,QF
      DOUBLE COMPLEX ZN
**
*     Internal Variables
*
      INTEGER I,J,K
      DOUBLE PRECISION BT0
      DOUBLE COMPLEX LAMBDAB,LN1ML,AQCD
      DOUBLE COMPLEX GAMMA0(2,2),GAMMA1(2,2),GAMMA2(2,2)
      DOUBLE COMPLEX GAMMA0NS,GAMMA1NS(3),GAMMA2NS(3)
      DOUBLE COMPLEX G0(2,2),G1(2,2),G2(2,2),G3(2,2)
      DOUBLE COMPLEX G0NS(3),G1NS(3),G2NS(3),G3NS(3)
      DOUBLE COMPLEX LNSUD(2,2),LNSUDNS(3)
      DOUBLE COMPLEX SUD(2,2),SUDNS(3)
      DOUBLE COMPLEX DELTA,SH,CH
**
*     Output Variables
*
      DOUBLE COMPLEX EFNNS(3),EFNSG(2,2)
*
*     LO
*
      IF(EVOL.EQ."SPACE")THEN
         IF(POL.EQ."OFF")THEN
            CALL ANDIM_LO(ZN,NF,GAMMA0NS,GAMMA0)
         ELSE
            CALL ANDIM_LO_POL(ZN,NF,GAMMA0NS,GAMMA0)
         ENDIF
      ELSEIF(EVOL.EQ."TIME")THEN
         IF(POL.EQ."OFF")THEN
            CALL ANDIM_LO_TL(ZN,NF,GAMMA0NS,GAMMA0)
         ELSE
C            CALL ANDIM_LO_TL_POL(ZN,NF,GAMMA0NS,GAMMA0)
         ENDIF
      ENDIF
*
*     NLO
*
      IF(IPT.GE.1)THEN
         IF(EVOL.EQ."SPACE")THEN
            IF(POL.EQ."OFF")THEN
               CALL ANDIM_NLO(ZN,NF,GAMMA1NS,GAMMA1)
            ELSE
               CALL ANDIM_NLO_POL(ZN,NF,GAMMA1NS,GAMMA1)
            ENDIF
         ELSEIF(EVOL.EQ."TIME")THEN
            IF(POL.EQ."OFF")THEN
               CALL ANDIM_NLO_TL(ZN,NF,GAMMA1NS,GAMMA1)
            ELSE
C               CALL ANDIM_NLO_TL_POL(ZN,NF,GAMMA1NS,GAMMA1)
            ENDIF
         ENDIF
      ENDIF
*
*     NNLO
*
      IF(IPT.GE.2)THEN
         IF(EVOL.EQ."SPACE")THEN
            IF(POL.EQ."OFF")THEN
               CALL ANDIM_NNLO(ZN,NF,GAMMA2NS,GAMMA2)
            ELSE
C               CALL ANDIM_NNLO_POL(ZN,NF,GAMMA2NS,GAMMA2)
            ENDIF
         ELSEIF(EVOL.EQ."TIME")THEN
            IF(POL.EQ."OFF")THEN
               CALL ANDIM_NNLO_TL(ZN,NF,GAMMA2NS,GAMMA2)
            ELSE
C               CALL ANDIM_NNLO_TL_POL(ZN,NF,GAMMA2NS,GAMMA2)
            ENDIF
         ENDIF
      ENDIF
*
      BT0 = - 2D0 * BETA0(NF)
      LAMBDAB = ASI * BT0 * ZLOG( KRF * QF / QI )
      LN1ML = ZLOG( 1D0 - LAMBDAB )
*
*     LO g-functions
*
*     Singlet
*
      G0(1,1) = (1D0, 0D0)
      G0(1,2) = (0D0, 0D0)
      G0(2,1) = (0D0, 0D0)
      G0(2,2) = (1D0, 0D0)
      DO I=1,2
         DO J=1,2
            G1(I,J) = - 2D0 * GAMMA0(I,J) * LN1ML / BT0
            G2(I,J) = (0D0, 0D0)
            G3(I,J) = (0D0, 0D0)
         ENDDO
      ENDDO
*
*     Non-singlet
*
      DO I=1,3
         G0NS(I) = (1D0, 0D0)
         G1NS(I) = - 2D0 * GAMMA0NS * LN1ML / BT0
         G2NS(I) = (0D0, 0D0)
         G3NS(I) = (0D0, 0D0)
      ENDDO
*
*     NLO g-functions
*
      IF(IPT.GE.1)THEN
*
*     Singlet
*
         DO I=1,2
            DO J=1,2
               G0(I,J) = G0(I,J)
     1              + ASI * 2D0 * ( GAMMA1(I,J) - B1(NF) * GAMMA0(I,J) )
     2              * LAMBDAB / ( 1D0 - LAMBDAB ) / BT0
               G2(I,J) = - 2D0 * GAMMA0(I,J) * ( B1(NF) / BT0 * LN1ML
     1              + DLOG(KRF) ) / ( 1D0 - LAMBDAB )
            ENDDO
         ENDDO
*
*     Non-singlet
*
         DO I=1,3
            G0NS(I) = G0NS(I)
     1           + ASI * 2D0 * ( GAMMA1NS(I) - B1(NF) * GAMMA0NS )
     2           * LAMBDAB / ( 1D0 - LAMBDAB ) / BT0
            G2NS(I) = - 2D0 * GAMMA0NS
     1           * ( B1(NF) / BT0 * LN1ML + DLOG(KRF) )
     2           / ( 1D0 - LAMBDAB )
         ENDDO
      ENDIF
*
*     NNLO g-functions
*
      IF(IPT.GE.2)THEN
*
*     Singlet
*
         DO I=1,2
            DO J=1,2
C               G0(I,J) = G0(I,J)
C     1              + ASI**2 * ( 0D0 ) / 2D0 / BT0**2 / ( 1 - LAMBDAB )**2
               G3(I,J) = 2D0 * GAMMA0(I,J) / BT0 * ( B2(NF) * LAMBDAB
     1              - B1(NF)**2 * ( LAMBDAB + LN1ML - LN1ML**2 / 2D0 )
     2              + BT0 * B1(NF) * ( LN1ML - 1D0 ) * DLOG(KRF)
     3              + BT0**2 * DLOG(KRF)**2 / 2D0 )
     4              / ( 1D0 - LAMBDAB )**2
            ENDDO
         ENDDO
*
*     Non-singlet
*
         DO I=1,3
            G0NS(I) = G0NS(I)
     1   + ASI**2 * ( B1(NF) * ( 2D0 * GAMMA1NS(I)
     2   * ( BT0 * ( ( LAMBDAB - 2D0 ) * LAMBDAB - 2D0 * LN1ML )
     3   - 2D0 * LAMBDAB**2 * 2D0 * GAMMA0NS )
     4   + 2D0 * BT0**2 * 2D0 * GAMMA0NS * DLOG(KRF) )
     5   + B1(NF)**2 * 2D0 * GAMMA0NS * ( LAMBDAB**2 * 2D0 * GAMMA0NS
     6   - BT0 * ( ( LAMBDAB - 2D0 ) * LAMBDAB - 2D0 * LN1ML ) )
     7   + B2(NF)*BT0 * (LAMBDAB - 2D0) * LAMBDAB * 2D0 * GAMMA0NS
     8   + 2D0 * GAMMA1NS(I) * ( LAMBDAB**2 * 2D0 * GAMMA1NS(I)
     9   - 2D0 * BT0**2 * DLOG(KRF) ) - BT0 * ( LAMBDAB - 2D0 )
     1   * LAMBDAB * 2D0 * GAMMA2NS(I) )
     2   / 2D0 / BT0**2 / ( 1 - LAMBDAB )**2
            G3NS(I) = 2D0 * GAMMA0NS / BT0 * ( B2(NF) * LAMBDAB
     1           - B1(NF)**2 * ( LAMBDAB + LN1ML - LN1ML**2 / 2D0 )
     2           + BT0 * B1(NF) * ( LN1ML - 1D0 ) * DLOG(KRF)
     3           + BT0**2 * DLOG(KRF)**2 / 2D0 )
     4           / ( 1D0 - LAMBDAB )**2
         ENDDO
      ENDIF
*
*     Argument of the exponential
*
*     Singlet
*
      DO I=1,2
         DO J=1,2
            LNSUD(I,J) = G1(I,J) + ASI * G2(I,J) + ASI**2 * G3(I,J)
         ENDDO
      ENDDO
*
*     Non-singlet
*
      DO I=1,3
         LNSUDNS(I) = G1NS(I) + ASI * G2NS(I) + ASI**2 * G3NS(I)
      ENDDO
*
*     Now take the exponential
*
*     Singlet
*
*     Formulae (15)-(19) on http://mathworld.wolfram.com/MatrixExponential.html
*
      DELTA = ZSQRT( ( LNSUD(1,1) - LNSUD(2,2) )**2D0 
     1     + 4D0 * LNSUD(1,2) * LNSUD(2,1) )
      SH = ZEXP( ( LNSUD(1,1) + LNSUD(2,2) ) / 2D0 ) 
     1     * ( ZEXP( DELTA / 2D0 ) - ZEXP( - DELTA / 2D0 ) ) / 2D0
      CH = ZEXP( ( LNSUD(1,1) + LNSUD(2,2) ) / 2D0 )
     1     * ( ZEXP( DELTA / 2D0 ) + ZEXP( - DELTA / 2D0 ) ) / 2D0
*
      SUD(1,1) = CH + SH * ( LNSUD(1,1) - LNSUD(2,2) ) / DELTA
      SUD(1,2) = 2D0 * LNSUD(1,2) * SH / DELTA
      SUD(2,1) = 2D0 * LNSUD(2,1) * SH / DELTA
      SUD(2,2) = CH - SH * ( LNSUD(1,1) - LNSUD(2,2) ) / DELTA
*
*     Non-singlet
*
      DO I=1,3
         SUDNS(I) = ZEXP(LNSUDNS(I))
      ENDDO
*
*     Finally multiply by g0 to obtain the evolution kernels
*
*     Singlet
*
      CALL MMULT(SUD,2,2,G0,2,2,EFNSG)
*
*     Non-singlet
*
      DO I=1,3
         EFNNS(I) = G0NS(I) * SUDNS(I)
      ENDDO
*
      RETURN
      END
 
