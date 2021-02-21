***********************************************************************
*
*     evolfactn.f:
*
*     This subroutine returns the QCD evolution factors, in N space, from
*     Q2I to Q2F with NF active flavours, for Non Singlet and Singlet-Gluon.
*
***********************************************************************
      SUBROUTINE EVOLFACTN(ZN,Q2I,Q2F,NF,EFNNS,EFNSG)
*
      IMPLICIT NONE
*
      include "../commons/evscale.h"
      include "../commons/beta.h"
      include "../commons/consts.h"
      include "../commons/alphas.h"
      include "../commons/massthrs.h"
      include "../commons/ipt.h"
      include "../commons/modev.h"
**
*     Input Variables
*
      INTEGER NF
      DOUBLE COMPLEX Q2I,Q2F
      DOUBLE COMPLEX ZN
**
*     Internal Variables
*
      INTEGER I,J,K
      DOUBLE PRECISION TMP3
      
      DOUBLE COMPLEX TMP
      DOUBLE COMPLEX TMP2,TMP5,TMP6
      DOUBLE COMPLEX ASI,ASF,T
      DOUBLE COMPLEX TMP1(3),TMP4(3)
      DOUBLE COMPLEX EFNSGTMP(2,2)
      DOUBLE COMPLEX EXPM,EXPP,LP,LM
      DOUBLE COMPLEX U(20,2,2),UNS(20,3)
      DOUBLE COMPLEX R(0:20,2,2),RNS(0:20,3)
      DOUBLE COMPLEX EM(2,2),EP(2,2)
      DOUBLE COMPLEX U1(2,2)
      DOUBLE COMPLEX L(2,2),LU1(2,2),U1L(2,2)
      DOUBLE COMPLEX U2(2,2)
      DOUBLE COMPLEX U2L(2,2),U1LU1(2,2),LU1U1(2,2),LU2(2,2)
      DOUBLE COMPLEX LNS(3)
      DOUBLE COMPLEX USUM(2,2),USUMTMP(2,2),UINV(2,2),USUML(2,2)
      DOUBLE COMPLEX DETINV
**
*     Output Variables
*
      DOUBLE COMPLEX EFNNS(3),EFNSG(2,2)
*
*     If initial and final scales are equal
*
      IF(ABS(Q2I).EQ.ABS(Q2F))THEN
         DO I=1,3
            EFNNS(I) = (1D0,0D0)
         ENDDO
         EFNSG(1,1) = (1D0,0D0)
         EFNSG(1,2) = (0D0,0D0)
         EFNSG(2,1) = (0D0,0D0)
         EFNSG(2,2) = (1D0,0D0)
         RETURN
      ENDIF
*
*     Compute alpha_s
*
      IF(ABS(Q2I).EQ.ABS(Q20)) THEN
         ASI = AS0
      ELSEIF(ABS(Q2I).EQ.ABS(Q2TH(4))) THEN
         ASI = ASC
      ELSEIF(ABS(Q2I).EQ.ABS(Q2TH(5))) THEN
         ASI = ASB
      ELSEIF(ABS(Q2I).EQ.ABS(Q2TH(6))) THEN
         ASI = AST
      ENDIF
*
      IF(ABS(Q2F).EQ.ABS(Q20)) THEN
         ASF = AS0
      ELSEIF(ABS(Q2F).EQ.ABS(Q2TH(4))) THEN
         ASF = ASCM
      ELSEIF(ABS(Q2F).EQ.ABS(Q2TH(5))) THEN
         ASF = ASBM
      ELSEIF(ABS(Q2F).EQ.ABS(Q2TH(6))) THEN
         ASF = ASTM
      ELSE
         ASF = ASQ
      ENDIF
*
*     Path-ordering solution
*
      IF(MODEV.EQ."PTH")THEN
         CALL PATH_ORDERING(ZN,ASI,ASF,NF,EFNNS,EFNSG)
         RETURN
      ENDIF
*
*     g-functions solution
*
      !write(6,*) ZN, Q2I, asi
      IF(MODEV.EQ."GFN")THEN
         CALL GFUNCS(ZN,ASI,ZSQRT(Q2I),ZSQRT(Q2F),NF,EFNNS,EFNSG)
         RETURN
      ENDIF
*
*     U matrices evaluation for singlet and non-singlet
*
      CALL UMATRIX(ZN,NF,LP,LM,EP,EM,U,R,UNS,RNS)
*
*     LO evolution factors 
*
      T = ZLOG(ASF/ASI)
*
*     Non singlet
*
      DO I=1,3
         LNS(I) = EXP ( - RNS(0,I) * T )    !Note anyway that RNS(0,1) = RNS(0,2) = RNS(0,3)
      ENDDO
*
*     Singlet
*
      EXPM = EXP ( - LM * T )
      EXPP = EXP ( - LP * T )
      DO I=1,2
         DO J=1,2 
            L(I,J) = EXPM * EM(I,J) + EXPP * EP(I,J) !Eq. (2.29) of hep-ph/0408244
         ENDDO
      ENDDO
*
      DO I=1,3
         EFNNS(I) = (0D0,0D0)
      ENDDO
*
      DETINV=0D0
      DO I=1,2
         DO J=1,2
            UINV(I,J)     = (0D0,0D0)
            EFNSG(I,J)    = (0D0,0D0)
            EFNSGTMP(I,J) = (0D0,0D0)
         ENDDO
      ENDDO 

      IF(IPT.EQ.0)THEN
*
*******************
*     LO solution
*
         DO I= 1,3
            EFNNS(I) = LNS(I)
         ENDDO
         DO I=1,2
            DO J=1,2
               EFNSG(I,J) = L(I,J)
            ENDDO
         ENDDO
      ELSEIF(IPT.EQ.1)THEN
*
*******************
*     NLO solution
*
         IF(MODEV.EQ."TRN")THEN
*
*     Truncated solution (Eq. (2.24) of hep-ph/0408244 for the non-singlet)
*
*     Non singlet
*
            DO I=1,3
               EFNNS(I) = LNS(I) * ( (1D0,0D0) + UNS(1,I) * (ASF-ASI) )
            ENDDO
*
*     Singlet
*
            DO I=1,2
               DO J=1,2
                  U1(I,J) = U(1,I,J)
               ENDDO
            ENDDO
*
            CALL MMULT(U1,2,2,L,2,2,U1L)
            CALL MMULT(L,2,2,U1,2,2,LU1)
            DO I=1,2
               DO J=1,2
                  EFNSG(I,J) = L(I,J) + ASF * U1L(I,J) - ASI * LU1(I,J) !Eq. (2.24) of hep-ph/0408244 for the singlet 
               ENDDO
            ENDDO
         ELSEIF(MODEV.EQ."ITE")THEN
*
*     Iterated solution
*
*     Non singlet
*
            TMP = ZLOG( ( 1D0 + B1(NF) * ASF ) 
     1                / ( 1D0 + B1(NF) * ASI ) ) / B1(NF)
            DO I=1,3
               EFNNS(I) = LNS(I) * EXP( TMP * UNS(1,I) )   !Eq. (2.34) of hep-ph/0408244 (analytical solution)
            ENDDO
*
*     Singlet
*
            USUM(1,1)    = (1D0,0D0)
            USUM(1,2)    = (0D0,0D0)
            USUM(2,1)    = (0D0,0D0)
            USUM(2,2)    = (1D0,0D0)
            USUMTMP(1,1) = (1D0,0D0)
            USUMTMP(1,2) = (0D0,0D0)
            USUMTMP(2,1) = (0D0,0D0)
            USUMTMP(2,2) = (1D0,0D0)
*
*     Calculation of U(a_s) and U(a_0) up to a^20 (which is an approximation of the
*     complete U matrix fot the iterated solution)
*
            DO K=1,20
               DO I=1,2
                  DO J=1,2
                     USUM(I,J) = USUM(I,J)
     1                         + ( ASF**(DBLE(K)) * U(K,I,J) )
                     USUMTMP(I,J) = USUMTMP(I,J)
     1                            + ( ASI**(DBLE(K)) * U(K,I,J) )
                   ENDDO
               ENDDO
            ENDDO
*
*     Inversion of U(a_0)
*
            DETINV    = 1D0 / ( USUMTMP(1,1) * USUMTMP(2,2)
     1                - USUMTMP(1,2) * USUMTMP(2,1) )
            UINV(1,1) = DETINV * USUMTMP(2,2)
            UINV(1,2) = - DETINV * USUMTMP(1,2)
            UINV(2,1) = - DETINV * USUMTMP(2,1)
            UINV(2,2) = DETINV * USUMTMP(1,1)
*
            CALL MMULT(USUM,2,2,L,2,2,USUML)
            CALL MMULT(USUML,2,2,UINV,2,2,EFNSGTMP)
*
            DO I=1,2
               DO J=1,2
                  EFNSG(I,J) = EFNSGTMP(I,J)
               ENDDO
            ENDDO
         ENDIF
      ELSEIF(IPT.EQ.2)THEN
*
*******************
*     NNLO solution
*
         IF(MODEV.EQ."TRN")THEN
*
*     Truncated solution
*
*     Non-singlet (Eq. (2.24) of hep-ph/0408244 with commuting matrices)
*
            DO I=1,3
               EFNNS(I) = LNS(I) * ( (1D0,0D0)                            !LO   
     1                             + UNS(1,I) * ( ASF - ASI )             !NLO
     2                             + UNS(1,I)**2 * ( ASI**2 - ASF*ASI )   !NNLO
     3                             + UNS(2,I) * ( ASF**2 - ASI**2 ) )
            ENDDO
*
*     Singlet (Eq. (2.24) of hep-ph/0408244)
*
            DO I=1,2
               DO J=1,2
                  U1(I,J) = U(1,I,J)
                  U2(I,J) = U(2,I,J)
               ENDDO
            ENDDO
*
*     NLO terms
*
            CALL MMULT(U1,2,2,L,2,2,U1L)
            CALL MMULT(L,2,2,U1,2,2,LU1)
*
*     NNLO terms
*
            CALL MMULT(U2,2,2,L,2,2,U2L)
            CALL MMULT(U1L,2,2,U1,2,2,U1LU1)
            CALL MMULT(LU1,2,2,U1,2,2,LU1U1)
            CALL MMULT(L,2,2,U2,2,2,LU2)
            DO I=1,2
               DO J=1,2
                  EFNSG(I,J) = L(I,J)                                       !LO
     1                       + ASF * U1L(I,J) - ASI * LU1(I,J)              !NLO  
     2                       + ASF**2D0 * U2L(I,J) - ASI**2D0 * LU2(I,J)    !NNLO
     3                       + ASI * ( ASI*LU1U1(I,J) - ASF*U1LU1(I,J) )
               ENDDO
            ENDDO
         ELSEIF(MODEV.EQ."ITE")THEN
*
*     Iterated solution
*
*     Non-singlet 
*
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                     c
c               Analytical NNLO solution              c
c                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*
            DO I=1,3
               TMP1(I) = - ( RNS(1,I) * B1(NF) + RNS(2,I) )
     1                 / ( 2D0 * B2(NF) )
               TMP4(I) = ( RNS(1,I) * B1(NF)**2. - 2D0 * RNS(1,I)
     1                 * B2(NF) + RNS(2,I) * B1(NF) ) / B2(NF)
            ENDDO

            TMP2 = ZLOG( (1D0 + B1(NF) * ASF + B2(NF) * ASF**2.)
     1                 / (1D0 + B1(NF) * ASI + B2(NF) * ASI**2.) )   
            TMP3 = DSQRT( - B1(NF)**2. + 4D0 * B2(NF))
            TMP5 = B1(NF) + 2D0 * ASF * B2(NF)
            TMP6 = B1(NF) + 2D0 * ASI * B2(NF)
*
            DO I=1,3
               EFNNS(I) = LNS(I) * EXP( TMP1(I) * TMP2 )
     1                           * EXP( ( TMP4(I) / TMP3 ) 
     2                  * ( ATAN(TMP5 / TMP3) - ATAN(TMP6 / TMP3) ) )
            ENDDO
*
*     Iterated solution non-singlet-like
*
C            DO I=1,3
C               UNSSUMF(I) = (1D0,0D0)
C               UNSSUMI(I) = (1D0,0D0)
C            ENDDO
C*
C            DO K=1,20
C               DO I=1,3
C                  UNSSUMF(I) = UNSSUMF(I) + UNS(K,I) * ASF**DBLE(K)
C                  UNSSUMI(I) = UNSSUMI(I) + UNS(K,I) * ASI**DBLE(K)
C               ENDDO
C            ENDDO
C*
C            DO I=1,3
C               EFNNS(I) = ( UNSSUMF(I) / UNSSUMI(I) ) * LNS(I)
C            ENDDO
*
*     Singlet (equal to NLO above)
*
            USUM(1,1)    = (1d0,0d0)
            USUM(1,2)    = (0d0,0d0)
            USUM(2,1)    = (0d0,0d0)
            USUM(2,2)    = (1d0,0d0)
            USUMTMP(1,1) = (1d0,0d0)
            USUMTMP(1,2) = (0d0,0d0)
            USUMTMP(2,1) = (0d0,0d0)
            USUMTMP(2,2) = (1d0,0d0)
*
*     Calculation of U(a_s) and U(a_0) up to a^20 (which is an approximation of the 
*     complete U matrix fot the iterated solution)
*
            DO K=1,20
               DO I=1,2
                  DO J=1,2
                     USUM(I,J) = USUM(I,J)
     1                         + ( U(K,I,J) * ASF**(DBLE(K)) )
                     USUMTMP(I,J) = USUMTMP(I,J)
     1                            + ( U(K,I,J) * ASI**(DBLE(K)) )
                   ENDDO
               ENDDO
            ENDDO
*
*     Inversion of U(a_0)
*
            DETINV    = 1D0 / ( USUMTMP(1,1) * USUMTMP(2,2)
     1                - USUMTMP(1,2) * USUMTMP(2,1) )
            UINV(1,1) = DETINV * USUMTMP(2,2)
            UINV(1,2) = - DETINV * USUMTMP(1,2)
            UINV(2,1) = - DETINV * USUMTMP(2,1)
            UINV(2,2) = DETINV * USUMTMP(1,1)
*
            CALL MMULT(USUM,2,2,L,2,2,USUML)
            CALL MMULT(USUML,2,2,UINV,2,2,EFNSGTMP)
*
            DO I=1,2
               DO J=1,2
                  EFNSG(I,J) = EFNSGTMP(I,J)
               ENDDO
            ENDDO
         ENDIF
      ENDIF
*
      RETURN
      END
