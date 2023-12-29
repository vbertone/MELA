************************************************************************
*
*     umatrix.f:
*
*     It computes the matrices used for the singlet evolution. 
*     Equation numbers refer to the "Notes on PT evolution".
*
*     Input variables: IPT : perturbative order
*                      ZN  : Mellin variable
*                      NF  : number of flavours
*
*     Functions:
*     Output: LP,LM     : LO matrix eigenvalues
*             EP,EM     : LO matrix eigenvectors
*             U(1..20)  : NLO or NNLO evolution matrices appearing into the 
*                         iterative solution for the singlet
*             R(0..20)  : NLO or NNLO evolution kernels appearing in the DGLAG
*                         written in terms of variation of alpha_s for singlet
*                         (eq. (2.20) hep-ph/0408244)                          
*             UNS(1..20): NLO or NNLO evolution coefficients appearing into the 
*                         iterative solution for the non-singlet
*             RNS(0..20): NLO or NNLO evolution kernels appearing in the DGLAG
*                         written in terms of variation of alpha_s for non-singlet
*
*     Modified to allow different renormalization and factorization scales 
*     in the anomalous dimensions.
*
*     KRF = MUR / MUF
*
*     Notation follows EQ. (2.8) of "NNLO evolution of DIS structure functions"
*     Van Neerven and Vogt hep-ph/0006154.
*
*     See also Eqns. (18-20) of the "Notes on PT evolution"
*
************************************************************************
      SUBROUTINE UMATRIX(ZN,NF,LP,LM,EP,EM,U,R,UNS,RNS)
*
      IMPLICIT NONE
*
      include "../commons/beta.h"
      include "../commons/colfact.h"
      include "../commons/ipt.h"
      include "../commons/renfacscales.h"
      include "../commons/evol.h"
      include "../commons/pol.h"
*
      INTEGER NF
      DOUBLE COMPLEX ZN
      DOUBLE COMPLEX U(20,2,2),UNS(20,3)
      DOUBLE COMPLEX R(0:20,2,2),RNS(0:20,3)
      DOUBLE COMPLEX LP,LM
      DOUBLE COMPLEX EP(2,2),EM(2,2)
*
      INTEGER I,J,K,II,JJ,KK,L
      DOUBLE PRECISION LOGR,B(20)
      DOUBLE COMPLEX G0(2,2),G1(2,2),G2(2,2)
      DOUBLE COMPLEX G0NS,G1NS(3),G2NS(3)
      DOUBLE COMPLEX P(0:20,2,2),PNS(0:20,3)
      DOUBLE COMPLEX RT(20,2,2)
      DOUBLE COMPLEX QQ0,QG0,GQ0,GG0
      DOUBLE COMPLEX SQ,LDIFF
*
*     Anomalous dimensions calculation
*
      LOGR = - 2D0 * DLOG(KRF)
*
*******************
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
*        Anomalous dimensions at NLO for muR.ne.muF
         IF(KRF.NE.1D0)THEN
            DO I=1,2
               DO J=1,2
                  G1(I,J) = G1(I,J) - BETA0(NF) * LOGR * G0(I,J) !Second line of eq. (2.8) of hep-ph/0408244)
               ENDDO
            ENDDO
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
         IF(KRF.NE.1D0)THEN
            DO I=1,2
               DO J=1,2
                  G2(I,J) = G2(I,J) - 2D0 * BETA0(NF) * LOGR * G1(I,J)
     1                    - ( BETA1(NF) * LOGR
     2                    +   BETA0(NF)**2D0 * LOGR**2D0 ) * G0(I,J)       !Third line of eq. (2.8) hep-ph/0408244 
               ENDDO
            ENDDO
            DO I=1,3
               G2NS(I) = G2NS(I) - 2D0 * BETA0(NF) * LOGR * G1NS(I)
     1                 - ( BETA1(NF) * LOGR
     2                 +   BETA0(NF)**2D0 * LOGR**2D0 ) * G0NS             !Third line of eq. (2.8) hep-ph/0408244
            ENDDO
         ENDIF
      ENDIF
*
*******************
*     b_i coefficients and DGLAP kernels assignment according to the perturbation order
*
      B(1) = B1(NF)
      B(2) = B2(NF)
      DO I=1,2
         DO J=1,2
            P(0,I,J) = G0(I,J)
            P(1,I,J) = G1(I,J)
            P(2,I,J) = G2(I,J)
         ENDDO
      ENDDO
      DO I=1,3
         PNS(0,I) = G0NS
         PNS(1,I) = G1NS(I)
         PNS(2,I) = G2NS(I)
      ENDDO

      DO K=IPT+1,20
         B(K) = 0.D0
         DO I=1,2
            DO J=1,2
               P(K,I,J) = (0.D0,0.D0)
            ENDDO
         ENDDO
         DO I=1,3
            PNS(K,I) = (0.D0,0.D0)
         ENDDO
      ENDDO
*
*******************************
*          SINGLET            *
*******************************
      DO K=1,20
         DO I=1,2
            DO J=1,2
               R(K,I,J)  = (0.D0,0.D0)
               RT(K,I,J) = (0.D0,0.D0)
               U(K,I,J)  = (0.D0,0.D0)
            ENDDO
         ENDDO
      ENDDO
*
*     Computation of R_0 (eq. (2.21) of hep-ph/0408244)
*    
      DO I=1,2
         DO J=1,2
            R(0,I,J) = P(0,I,J) / BETA0(NF)
         ENDDO
      ENDDO
*
      QQ0 = R(0,1,1)
      QG0 = R(0,1,2)
      GQ0 = R(0,2,1)
      GG0 = R(0,2,2)
*
*     Computation of the eigenvalues of R_0 (eq. (2.27) of hep-ph/0408244)
*
      SQ = SQRT( ( QQ0 - GG0 )**2D0 + 4D0 * QG0 * GQ0 )
      LP = 0.5d0 * ( QQ0 + GG0 + SQ )
      LM = 0.5d0 * ( QQ0 + GG0 - SQ )
*
*     Computation of the projectors of R_0 (eq. (2.28) of hep-ph/0408244)
*
      EM(1,1) = ( QQ0 - LP ) / ( LM - LP )
      EM(1,2) = QG0 / ( LM - LP )
      EM(2,1) = GQ0 / ( LM - LP )
      EM(2,2) = ( GG0 - LP ) / ( LM - LP )
*
      EP(1,1) = ( QQ0 - LM ) / ( LP - LM )
      EP(1,2) = QG0 / ( LP - LM )
      EP(2,1) = GQ0 / ( LP - LM )
      EP(2,2) = ( GG0 - LM ) / ( LP - LM )
*
*     Computation of R_{1...20} matrices (eq. (2.21) of hep-ph/0408244)
*
      DO K=1,20
         DO I=1,2
            DO J=1,2
               R(K,I,J) = P(K,I,J) / BETA0(NF)
               DO L=1,K
                  R(K,I,J) = R(K,I,J) - B(L) * R(K-L,I,J)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
*
      DO K=1,20
*
*     Computation of \widetilde{R}_k matrix (eq. (2.25) of hep-ph/0408244)
*
         DO I=1,2
            DO J=1,2
               RT(K,I,J) = R(K,I,J)
               IF(K.GT.1)THEN
                  DO L=1,K-1
                     DO KK=1,2
                        RT(K,I,J) = RT(K,I,J) + R(L,I,KK) * U(K-L,KK,J) 
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
*
*     Computation of U_k matrix from \widetilde{R}_k (eq. (2.31) of hep-ph/0408244)
*
         LDIFF = LM - LP
         DO I=1,2
            DO J=1,2
               DO II=1,2
                  DO JJ=1,2                       
                     U(K,I,J) = U(K,I,J)
     1               - ( EM(I,II) * RT(K,II,JJ) * EM(JJ,J) / DBLE(K) )
     2               - ( EP(I,II) * RT(K,II,JJ) * EP(JJ,J) / DBLE(K) )
     3               - ( EM(I,II) * RT(K,II,JJ) * EP(JJ,J) /
     4                 ( DBLE(K) + LDIFF ) )
     5               - ( EP(I,II) * RT(K,II,JJ) * EM(JJ,J) /
     6                 ( DBLE(K) - LDIFF ) )
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
*
*******************************
*        NON-SINGLET          *
*******************************
*
*     Computation of R_0 (eq. (2.21) of hep-ph/0408244)
*    
      DO I=1,3
         RNS(0,I) = PNS(0,I) / BETA0(NF)
      ENDDO
*
*     Computation of R_{1...20} matrices
*
      DO K=1,20
         DO I=1,3
            RNS(K,I) = PNS(K,I) / BETA0(NF)
            DO L=1,K
               RNS(K,I) = RNS(K,I) - B(L) * RNS(K-L,I)
            ENDDO
         ENDDO
      ENDDO
*
*     Computation of U_{1...20} matrices
*
      DO K=1,20
         DO I=1,3
            UNS(K,I) = - RNS(K,I) / DBLE(K)
            IF(K.GT.1)THEN
               DO L=1,K-1
                  UNS(K,I) = UNS(K,I) - RNS(L,I) * UNS(K-L,I) / DBLE(K)
               ENDDO
            ENDIF
         ENDDO
      ENDDO
*
      RETURN
      END
