***********************************************************************
*
*     path_ordering.f:
*
*     This subroutine returns the QCD evolution factors, in N space, from
*     Q2I to Q2F with NF active flavours, for Non Singlet and Singlet-Gluon
*     using the path ordering method.
*
***********************************************************************
      SUBROUTINE PATH_ORDERING(ZN,ASI,ASF,NF,EFNNS,EFNSG)
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
      DOUBLE COMPLEX ASI,ASF
      DOUBLE COMPLEX ZN
**
*     Internal Variables
*
      INTEGER NINT
      INTEGER I,J,K,L
      DOUBLE PRECISION AS_CHECK,PREC
      DOUBLE PRECISION BT0,BT1,BT2
      DOUBLE COMPLEX AF,AI,DA,AK
      DOUBLE COMPLEX LOGR
      DOUBLE COMPLEX G0(2,2),G1(2,2),G2(2,2)
      DOUBLE COMPLEX G0NS,G1NS(3),G2NS(3)
      DOUBLE COMPLEX PSG(2,2),PNS(3)
      DOUBLE COMPLEX SPSG(2,2),SPNS(3)
      DOUBLE COMPLEX DEFNSG(2,2),EFNSG_TMP(2,2)
      DOUBLE COMPLEX DELTA,SH,CH
      PARAMETER(NINT=20)
      PARAMETER(PREC=1D-7)
**
*     Output Variables
*
      DOUBLE COMPLEX EFNNS(3),EFNSG(2,2)
*
      DO I=1,2
         DO J=1,2
            G0(I,J) = (0D0,0D0)
            G1(I,J) = (0D0,0D0)
            G2(I,J) = (0D0,0D0)
*
            SPSG(I,J) = (0D0,0D0)
         ENDDO
      ENDDO
*
      G0NS = (0D0,0D0)
      DO I=1,3
         G1NS(I) = (0D0,0D0)
         G2NS(I) = (0D0,0D0)
*
         SPNS(I) = (0D0,0D0)
      ENDDO
*
      BT0 = BETA0(NF)
      BT1 = 0D0
      BT2 = 0D0
*
      LOGR = - 2D0 * DLOG(KRF)
*
*     LO
*
      IF(EVOL.EQ."SPACE")THEN
         IF(POL.EQ."OFF")THEN
            CALL ANDIM_LO(ZN,NF,G0NS,G0)
         ELSE
            CALL ANDIM_LO_POL(ZN,NF,G0NS,G0)
         ENDIF
      ELSEIF(EVOL.EQ."TIME")THEN
         IF(POL.EQ."OFF")THEN
            CALL ANDIM_LO_TL(ZN,NF,G0NS,G0)
         ELSE
C            CALL ANDIM_LO_TL_POL(ZN,NF,G0NS,G0)
         ENDIF
      ENDIF
*
*     NLO
*
      IF(IPT.GE.1)THEN
         IF(EVOL.EQ."SPACE")THEN
            IF(POL.EQ."OFF")THEN
               CALL ANDIM_NLO(ZN,NF,G1NS,G1)
            ELSE
               CALL ANDIM_NLO_POL(ZN,NF,G1NS,G1)
            ENDIF
         ELSEIF(EVOL.EQ."TIME")THEN
            IF(POL.EQ."OFF")THEN
               CALL ANDIM_NLO_TL(ZN,NF,G1NS,G1)
            ELSE
C               CALL ANDIM_NLO_TL_POL(ZN,NF,G1NS,G1)
            ENDIF
         ENDIF
         BT1 = BETA1(NF)
*        Anomalous dimensions at NLO for muR.ne.muF
         IF(KRF.NE.1D0)THEN
            DO I=1,2
               DO J=1,2
                  G1(I,J) = G1(I,J) - BETA0(NF) * LOGR * G0(I,J) !Second line of eq. (2.8) of hep-ph/0408244)
               ENDDO
            ENDDO
*
            DO I=1,3
               G1NS(I) = G1NS(I) - BETA0(NF) * LOGR * G0NS       !Second line of eq. (2.8) of hep-ph/0408244)
            ENDDO
         ENDIF
      ENDIF
*
*     NNLO
*
      IF(IPT.GE.2)THEN
         IF(EVOL.EQ."SPACE")THEN
            IF(POL.EQ."OFF")THEN
               CALL ANDIM_NNLO(ZN,NF,G2NS,G2)
            ELSE
C               CALL ANDIM_NNLO_POL(ZN,NF,G2NS,G2)
            ENDIF
         ELSEIF(EVOL.EQ."TIME")THEN
            IF(POL.EQ."OFF")THEN
               CALL ANDIM_NNLO_TL(ZN,NF,G2NS,G2)
            ELSE
C               CALL ANDIM_NNLO_TL_POL(ZN,NF,G2NS,G2)
            ENDIF
         ENDIF
         BT2 = BETA2(NF)
         IF(KRF.NE.1D0)THEN
            DO I=1,2
               DO J=1,2
                  G2(I,J) = G2(I,J) - 2D0 * BETA0(NF) * LOGR * G1(I,J)
     1                    - ( BETA1(NF) * LOGR
     2                    +   BETA0(NF)**2D0 * LOGR**2D0 ) * G0(I,J)       !Third line of eq. (2.8) hep-ph/0408244 
               ENDDO
            ENDDO
*
            DO I=1,3
               G2NS(I) = G2NS(I) - 2D0 * BETA0(NF) * LOGR * G1NS(I)
     1                 - ( BETA1(NF) * LOGR
     2                 +   BETA0(NF)**2D0 * LOGR**2D0 ) * G0NS             !Third line of eq. (2.8) hep-ph/0408244
            ENDDO
         ENDIF
      ENDIF
*
*     Singlet
*
*     Solution at leading order (easy)
*
      DA = ( ASF - ASI ) / DBLE(NINT)
      IF(IPT.EQ.0)THEN
         AI = ASI
         AF = ASF
         DO I=1,2
            DO J=1,2
               PSG(I,J) = G0(I,J) * ( 1D0 / AF + 1D0 / AI ) / BT0
*
               SPSG(I,J) = - DA * PSG(I,J) / 2D0
            ENDDO
         ENDDO
*
         AK = AI
         DO K=1,NINT-1
            AK = AK + DA
            DO I=1,2
               DO J=1,2
                  PSG(I,J)  = G0(I,J) / BT0 / AK
                  SPSG(I,J) = SPSG(I,J) - DA * PSG(I,J)
               ENDDO
            ENDDO
         ENDDO
*
*     Check
*
         AS_CHECK = ABS( AK + DA - AF )
         IF(AS_CHECK.GT.PREC)THEN
            WRITE(6,*) "In path_ordering.f:"
            WRITE(6,*) "Mismatch, ASF =",AF,", AK =",AK + DA
            CALL EXIT(-10)
         ENDIF
*
*     Formulae (15)-(19) on http://mathworld.wolfram.com/MatrixExponential.html
*
         DELTA = ZSQRT( ( SPSG(1,1) - SPSG(2,2) )**2D0 
     1         + 4D0 * SPSG(1,2) * SPSG(2,1) )
         SH = ZEXP( ( SPSG(1,1) + SPSG(2,2) ) / 2D0 ) 
     1      * ( ZEXP( DELTA / 2D0 ) - ZEXP( - DELTA / 2D0 ) ) / 2D0
         CH = ZEXP( ( SPSG(1,1) + SPSG(2,2) ) / 2D0 )
     1      * ( ZEXP( DELTA / 2D0 ) + ZEXP( - DELTA / 2D0 ) ) / 2D0
*
         EFNSG(1,1) = CH + SH * ( SPSG(1,1) - SPSG(2,2) ) / DELTA
         EFNSG(1,2) = 2D0 * SPSG(1,2) * SH / DELTA
         EFNSG(2,1) = 2D0 * SPSG(2,1) * SH / DELTA
         EFNSG(2,2) = CH - SH * ( SPSG(1,1) - SPSG(2,2) ) / DELTA
*
*     Solution at NLO and NNLO
*
      ELSE
         EFNSG(1,1) = (1D0,0D0)
         EFNSG(1,2) = (0D0,0D0)
         EFNSG(2,1) = (0D0,0D0)
         EFNSG(2,2) = (1D0,0D0)
*
         AI = ASI
         DO K=1,NINT
            AF = AI + DA
            DO I=1,2
               DO J=1,2
                  PSG(I,J) = ( G0(I,J) + AI*G1(I,J) + AI**2D0*G2(I,J) ) 
     1                     / ( BT0 + AI * BT1 + AI**2D0 * BT2 ) / AI
*
                  PSG(I,J) = PSG(I,J)
     1                     + ( G0(I,J) + AF*G1(I,J) + AF**2D0*G2(I,J) )
     2                     / ( BT0 + AF * BT1 + AF**2D0 * BT2 ) / AF
*
                  SPSG(I,J) = - DA * PSG(I,J) / 2D0
               ENDDO
            ENDDO
*
*     Formulae (15)-(19) on http://mathworld.wolfram.com/MatrixExponential.html
*
            DELTA = ZSQRT( ( SPSG(1,1) - SPSG(2,2) )**2D0 
     1            + 4D0 * SPSG(1,2) * SPSG(2,1) )
            SH = ZEXP( ( SPSG(1,1) + SPSG(2,2) ) / 2D0 ) 
     1         * ( ZEXP( DELTA / 2D0 ) - ZEXP( - DELTA / 2D0 ) ) / 2D0
            CH = ZEXP( ( SPSG(1,1) + SPSG(2,2) ) / 2D0 )
     1         * ( ZEXP( DELTA / 2D0 ) + ZEXP( - DELTA / 2D0 ) ) / 2D0
*
            DEFNSG(1,1) = CH + SH * ( SPSG(1,1) - SPSG(2,2) ) / DELTA
            DEFNSG(1,2) = 2D0 * SPSG(1,2) * SH / DELTA
            DEFNSG(2,1) = 2D0 * SPSG(2,1) * SH / DELTA
            DEFNSG(2,2) = CH - SH * ( SPSG(1,1) - SPSG(2,2) ) / DELTA
*
            EFNSG_TMP(1,1) = (0D0,0D0)
            EFNSG_TMP(1,2) = (0D0,0D0)
            EFNSG_TMP(2,1) = (0D0,0D0)
            EFNSG_TMP(2,2) = (0D0,0D0)
*
            DO I=1,2
               DO J=1,2
                  DO L=1,2
                     EFNSG_TMP(I,J) = EFNSG_TMP(I,J)
     1                              + DEFNSG(I,L) * EFNSG(L,J)
                  ENDDO
               ENDDO
            ENDDO
*
            DO I=1,2
               DO J=1,2
                  EFNSG(I,J) = EFNSG_TMP(I,J)
               ENDDO
            ENDDO
*
            AI = AI + DA
         ENDDO
      ENDIF
*
*     Non Singlet
*
      AI = ASI
      AF = ASF
      DO I=1,3
         PNS(I) = ( G0NS + AI * G1NS(I) + AI**2D0 * G2NS(I) )
     1          / ( BT0 + AI * BT1 + AI**2D0 * BT2 ) / AI
*
         PNS(I) = PNS(I)
     1          + ( G0NS + AF * G1NS(I) + AF**2D0 * G2NS(I) )
     2          / ( BT0 + AF * BT1 + AF**2D0 * BT2 ) / AF
*
         SPNS(I) = - DA * PNS(I) / 2D0
      ENDDO
*
      AK = AI
      DO K=1,NINT-1
         AK = AK + DA
*
         DO I=1,3
            PNS(I) = ( G0NS + AK * G1NS(I) + AK**2D0 * G2NS(I) )
     2             / ( BT0 + AK * BT1 + AK**2D0 * BT2 ) / AK
*
            SPNS(I) = SPNS(I) - DA * PNS(I)
         ENDDO
      ENDDO
*
*     Check
*
      AS_CHECK = ABS( AK + DA - AF )
      IF(AS_CHECK.GT.PREC)THEN
         WRITE(6,*) "In path_ordering.f:"
         WRITE(6,*) "Mismatch, AF =",AF,", AK =",AK + DA
         CALL EXIT(-10)
      ENDIF
*
      DO I=1,3
         EFNNS(I) = ZEXP( SPNS(I) )
      ENDDO
*
      RETURN
      END
 
