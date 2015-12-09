************************************************************************
*
*     The following routines return the netral current coefficient
*     functions in the Mellin space in the FFNS
*
************************************************************************
      FUNCTION C2GH_FFNS_NC(JPT,N,Q2,M2,RUN)
*
      IMPLICIT NONE
*
      include "../commons/colfact.h"
      include "../commons/renfacscales.h"
      include "../commons/hqmass.h"
**
*     Input Variables
*
      INTEGER JPT,RUN
      DOUBLE COMPLEX Q2,M2
      DOUBLE COMPLEX N
**
*     Internal Variables
*
      DOUBLE PRECISION H1
      DOUBLE COMPLEX EPS,M,A,LN
      DOUBLE COMPLEX IQ(0:2),JQ(0:2)
      DOUBLE COMPLEX DIQ(0:2),DJQ(0:2)
      DOUBLE COMPLEX C2GH,C2GHBAR
      DOUBLE COMPLEX CH_FFNS_NC_AB_LIN
      DOUBLE COMPLEX DC2GH_FFNS_NC
**
*     Output Variables
*
      DOUBLE COMPLEX C2GH_FFNS_NC
*
      EPS = M2 / Q2
      IF(JPT.EQ.1)THEN
         CALL IQN(N,EPS,IQ)
         CALL JQN(N,EPS,JQ)
*
         C2GH_FFNS_NC = 4D0 * TR * ( 2D0 * ( 1D0 - 6D0 * EPS 
     1                - 4D0 * EPS**2D0 ) * IQ(2) - 2D0 * ( 1D0 
     2                - 2D0 * EPS ) * IQ(1) + IQ(0) - 4D0 * (2D0 - EPS) 
     3                * JQ(2) + 4D0 * ( 2D0 - EPS ) * JQ(1) - JQ(0) )
      ELSEIF(JPT.EQ.2)THEN
         C2GH    = CH_FFNS_NC_AB_LIN(3,N,Q2,M2)
         C2GHBAR = CH_FFNS_NC_AB_LIN(9,N,Q2,M2)
*
         C2GH_FFNS_NC = C2GH + C2GHBAR * ZLOG( KFACQ**2D0 * Q2 / M2 )
*
         IF(HQMASS.EQ.1)THEN
            IF(RUN.EQ.1)THEN
               LN = ZLOG( KRENQ**2D0 * Q2 / M2 )
            ELSE
               LN = 0D0
            ENDIF
            H1 = CF * ( 4D0 + 3D0 * LN )
            M  = SQRT(M2)
            A  = 1D0 / ( 1D0 + 4D0 * EPS )
*
            CALL IQN(N,EPS,IQ)
            CALL JQN(N,EPS,JQ)
            CALL DIQN(N,EPS,DIQ)
            CALL DJQN(N,EPS,DJQ)
*
            DC2GH_FFNS_NC = 4D0 * TR * 2D0 * EPS * ( 2D0 * ( - 6D0 
     1                    - 8D0 * EPS ) * IQ(2) + 4D0 * IQ(1) 
     2                    + 4D0 * JQ(2) - 4D0 * JQ(1)
     3                    - 4D0 * A**2D0 * ( 2D0 * ( 1D0 - 6D0 * EPS 
     4                    - 4D0 * EPS**2D0 ) * DIQ(2) - 2D0 * ( 1D0 
     5                    - 2D0 * EPS ) * DIQ(1) + DIQ(0) 
     6                    - 4D0 * ( 2D0-EPS ) * DJQ(2) + 4D0 
     7                    * ( 2D0 - EPS ) * DJQ(1) - DJQ(0) ) ) / M
*
            C2GH_FFNS_NC = C2GH_FFNS_NC + M * H1 * DC2GH_FFNS_NC
         ENDIF
      ENDIF
*
      RETURN
      END
*
************************************************************************
      FUNCTION C2QH_FFNS_NC(JPT,N,Q2,M2)
*
      IMPLICIT NONE
*
      include "../commons/renfacscales.h"
**
*     Input Variables
*
      INTEGER JPT
      DOUBLE COMPLEX Q2,M2
      DOUBLE COMPLEX N
**
*     Internal Variables
*
      DOUBLE COMPLEX C2QH,C2QHBAR
      DOUBLE COMPLEX CH_FFNS_NC_AB_LIN
**
*     Output Variables
*
      DOUBLE COMPLEX C2QH_FFNS_NC
*
      IF(JPT.EQ.1)THEN
         C2QH_FFNS_NC = (0D0,0D0)
      ELSEIF(JPT.EQ.2)THEN
         C2QH    = CH_FFNS_NC_AB_LIN(4,N,Q2,M2)
         C2QHBAR = CH_FFNS_NC_AB_LIN(10,N,Q2,M2)
*
         C2QH_FFNS_NC = C2QH + C2QHBAR * ZLOG( KFACQ**2D0 * Q2 / M2 )
      ENDIF
*
      RETURN
      END
*
************************************************************************
      FUNCTION C2QGRH_FFNS_NC(JPT,N,Q2,M2)
*
      IMPLICIT NONE
*
      include "../commons/adlersr.h"
      include "../commons/massthrs.h"
**
*     Input Variables
*
      INTEGER JPT
      DOUBLE COMPLEX Q2,M2
      DOUBLE COMPLEX N
**
*     Internal Variables
*
      INTEGER K
      DOUBLE COMPLEX C2QGRH
      DOUBLE COMPLEX CH_FFNS_NC_AB_LIN
**
*     Output Variables
*
      DOUBLE COMPLEX C2QGRH_FFNS_NC
*
      IF(JPT.EQ.1)THEN
         C2QGRH_FFNS_NC = (0D0,0D0)
      ELSEIF(JPT.EQ.2)THEN
         C2QGRH = CH_FFNS_NC_AB_LIN(7,N,Q2,M2)
*
         IF(M2.EQ.Q2TH(4))THEN
            K = 4
         ELSEIF(M2.EQ.Q2TH(5))THEN
            K = 5
         ELSEIF(M2.EQ.Q2TH(6))THEN
            K = 6
         ELSE
            WRITE(6,*) "In modules/evolutionQCD/wc_FFNS_NC.f:"
            WRITE(6,*) "M2 does not match with any HQ mass"
            WRITE(6,*) "M2 =",M2,"GeV^2"
            WRITE(6,*) "mc2(mu) =",Q2TH(4),"GeV^2,  mb2(mu) =",
     1                 Q2TH(5),"GeV^2,  mt2(mu) =",Q2TH(6)
            CALL EXIT(-10)
         ENDIF
*     Costant addition needed to fulfil the Adler sum rule (see eq. (97) of ArXiv.1001.2312)
         C2QGRH_FFNS_NC = C2QGRH - ADLERSR(K)
      ENDIF
*
      RETURN
      END
*
************************************************************************
      FUNCTION CLGH_FFNS_NC(JPT,N,Q2,M2,RUN)
*
      IMPLICIT NONE
*
      include "../commons/colfact.h"
      include "../commons/renfacscales.h"
      include "../commons/hqmass.h"
**
*     Input Variables
*
      INTEGER JPT,RUN
      DOUBLE COMPLEX Q2,M2
      DOUBLE COMPLEX N
**
*     Internal Variables
*
      DOUBLE PRECISION H1
      DOUBLE COMPLEX EPS,M,A,LN
      DOUBLE COMPLEX IQ(0:2),JQ(0:2)
      DOUBLE COMPLEX DIQ(0:2),DJQ(0:2)
      DOUBLE COMPLEX CLGH,CLGHBAR
      DOUBLE COMPLEX CH_FFNS_NC_AB_LIN
      DOUBLE COMPLEX DCLGH_FFNS_NC
**
*     Output Variables
*
      DOUBLE COMPLEX CLGH_FFNS_NC
*
      EPS = M2 / Q2
      IF(JPT.EQ.1)THEN
         CALL IQN(N,EPS,IQ)
         CALL JQN(N,EPS,JQ)
*
         CLGH_FFNS_NC = 4D0 * TR * ( - 8D0 * EPS * IQ(2) - 4D0 * JQ(2) 
     1                               + 4D0 * JQ(1) )
      ELSEIF(JPT.EQ.2)THEN
         CLGH    = CH_FFNS_NC_AB_LIN(5,N,Q2,M2)
         CLGHBAR = CH_FFNS_NC_AB_LIN(11,N,Q2,M2)
*
         CLGH_FFNS_NC = CLGH + CLGHBAR * ZLOG( KFACQ**2D0 * Q2 / M2 )
*
         IF(HQMASS.EQ.1)THEN
            IF(RUN.EQ.1)THEN
               LN = ZLOG( KRENQ**2D0 * Q2 / M2 )
            ELSE
               LN = 0D0
            ENDIF
            H1 = CF * ( 4D0 + 3D0 * LN )
            M  = SQRT(M2)
            A  = 1D0 / ( 1D0 + 4D0 * EPS )
*
            CALL IQN(N,EPS,IQ)
            CALL JQN(N,EPS,JQ)
            CALL DIQN(N,EPS,DIQ)
            CALL DJQN(N,EPS,DJQ)
*
            DCLGH_FFNS_NC = 4D0 * TR * 8D0 * EPS * ( - 2D0 * IQ(2) 
     1                    - A**2D0 * ( - 8D0 * EPS * DIQ(2) 
     2                    - 4D0 * DJQ(2) + 4D0 * DJQ(1) ) ) / M
*
            CLGH_FFNS_NC = CLGH_FFNS_NC + M * H1 * DCLGH_FFNS_NC
         ENDIF
      ENDIF
*
      RETURN
      END
*
************************************************************************
      FUNCTION CLQH_FFNS_NC(JPT,N,Q2,M2)
*
      IMPLICIT NONE
*
      include "../commons/renfacscales.h"
**
*     Input Variables
*
      INTEGER JPT
      DOUBLE COMPLEX Q2,M2
      DOUBLE COMPLEX N
**
*     Internal Variables
*
      DOUBLE COMPLEX CLQH,CLQHBAR
      DOUBLE COMPLEX CH_FFNS_NC_AB_LIN
**
*     Output Variables
*
      DOUBLE COMPLEX CLQH_FFNS_NC
*
      IF(JPT.EQ.1)THEN
         CLQH_FFNS_NC = (0D0,0D0)
      ELSEIF(JPT.EQ.2)THEN
         CLQH    = CH_FFNS_NC_AB_LIN(6,N,Q2,M2)
         CLQHBAR = CH_FFNS_NC_AB_LIN(12,N,Q2,M2)
*
         CLQH_FFNS_NC = CLQH + CLQHBAR * ZLOG( KFACQ**2D0 * Q2 / M2 )
      ENDIF
*
      RETURN
      END
*
************************************************************************
      FUNCTION CLQGRH_FFNS_NC(JPT,N,Q2,M2)
*
      IMPLICIT NONE
**
*     Input Variables
*
      INTEGER JPT
      DOUBLE COMPLEX Q2,M2
      DOUBLE COMPLEX N
**
*     Internal Variables
*
      DOUBLE COMPLEX CLQGRH
      DOUBLE COMPLEX CH_FFNS_NC_AB_LIN
**
*     Output Variables
*
      DOUBLE COMPLEX CLQGRH_FFNS_NC
*
      IF(JPT.EQ.1)THEN
         CLQGRH_FFNS_NC = (0D0,0D0)
      ELSEIF(JPT.EQ.2)THEN
         CLQGRH = CH_FFNS_NC_AB_LIN(8,N,Q2,M2)
*
         CLQGRH_FFNS_NC = CLQGRH
      ENDIF
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE IQN(N,EPS,IQ)
*
      IMPLICIT NONE
*
      include "../commons/consts.h"
**
*     Input Variable
*
      DOUBLE COMPLEX EPS
      DOUBLE COMPLEX N
**
*     Internal Variables
*
      INTEGER Q
      DOUBLE COMPLEX HYPGEO
      DOUBLE COMPLEX CONST,CONST1,CONST2,CONST3
      DOUBLE COMPLEX SUM
      DOUBLE COMPLEX A,B,C,Z
      DOUBLE COMPLEX G1,G2,G3,Z1,Z2,Z3
**
*     Output Variables
*
      DOUBLE COMPLEX IQ(0:2)
*
      Z = 1D0 / ( 1D0 + 4D0 * EPS )
*
      DO Q=0,2
         CONST1 = EXP( - N * ZLOG( 1D0 + 4D0 * EPS ) )
         CONST2 = EXP( - Q * ZLOG( 1D0 + 4D0 * EPS ) )
         CONST3 = N + Q
         CONST = CONST1 * CONST2 / CONST3
*
         A = 0.5D0
         B = N + Q
         C = N + Q + 0.5D0
*
         Z1 = B
         Z2 = C - B
         Z3 = C
*
         CALL GAMMAL(Z1,G1)
C         CALL GAMMAL(Z2,G2)
         G2 = DCMPLX(DLOG(DSQRT(PI)),0D0)
         CALL GAMMAL(Z3,G3)
*
         SUM = EXP( G1 + G2 - G3 ) * HYPGEO(A,B,C,Z)
*
         IQ(Q) = CONST * SUM
      ENDDO
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE JQN(N,EPS,JQ)
*
      IMPLICIT NONE
*
      include "../commons/consts.h"
**
*     Input Variable
*
      DOUBLE COMPLEX EPS
      DOUBLE COMPLEX N
**
*     Internal Variables
*
      INTEGER Q
      DOUBLE COMPLEX HYPGEO
      DOUBLE COMPLEX CONST,CONST1,CONST2
      DOUBLE COMPLEX A1,B1,C1,A2,B2,C2,Z
      DOUBLE COMPLEX G1,G2,G3,Z1,Z2,Z3
**
*     Output Variables
*
      DOUBLE COMPLEX JQ(0:2)
*
      Z = 1D0 / ( 1D0 + 4D0 * EPS )
*
      DO Q=0,2
         CONST1 = EXP( - N * ZLOG( 1D0 + 4D0 * EPS ) ) 
         CONST2 = EXP( - Q * ZLOG( 1D0 + 4D0 * EPS ) )
         CONST = CONST1 * CONST2
*
         A1 = 0.5D0
         B1 = N + Q
         C1 = N + Q + 0.5D0
*
         A2 = 0.5D0
         B2 = N + Q + 1D0
         C2 = N + Q + 1.5D0
*
         Z1 = B1
         Z2 = C1 - B1
         Z3 = C1
*
         CALL GAMMAL(Z1,G1)
C         CALL GAMMAL(Z2,G2)
         G2 = DCMPLX(DLOG(DSQRT(PI)),0D0)
         CALL GAMMAL(Z3,G3)
*
         JQ(Q) = CONST * EXP( G1 + G2 - G3 ) * ( HYPGEO(A1,B1,C1,Z)
     1         - ( N + Q ) / ( N + Q + 0.5D0 ) * HYPGEO(A2,B2,C2,Z) )
      ENDDO
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE DIQN(N,EPS,DIQ)
*
      IMPLICIT NONE
*
      include "../commons/consts.h"
**
*     Input Variable
*
      DOUBLE COMPLEX EPS
      DOUBLE COMPLEX N
**
*     Internal Variables
*
      INTEGER Q
      DOUBLE COMPLEX HYPGEO
      DOUBLE COMPLEX CONST,CONST1,CONST2
      DOUBLE COMPLEX SUM
      DOUBLE COMPLEX A,B,C,Z
      DOUBLE COMPLEX G1,G2,G3,Z1,Z2,Z3
**
*     Output Variables
*
      DOUBLE COMPLEX DIQ(0:2)
*
      Z = 1D0 / ( 1D0 + 4D0 * EPS )
*
      DO Q=0,2
         CONST1 = EXP( - N * ZLOG( 1D0 + 4D0 * EPS ) )
         CONST2 = EXP( - Q * ZLOG( 1D0 + 4D0 * EPS ) )
         CONST = CONST1 * CONST2 / Z
*
         A = 0.5D0
         B = N + Q + 1D0
         C = N + Q + 0.5D0
*
         Z1 = N + Q
         Z2 = 1D0 / 2D0
         Z3 = N + Q + 0.5D0
*
         CALL GAMMAL(Z1,G1)
C         CALL GAMMAL(Z2,G2)
         G2 = DCMPLX(DLOG(DSQRT(PI)),0D0)
         CALL GAMMAL(Z3,G3)
*
         SUM = EXP( G1 + G2 - G3 ) * HYPGEO(A,B,C,Z)
*
         DIQ(Q) = CONST * SUM
      ENDDO
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE DJQN(N,EPS,DJQ)
*
      IMPLICIT NONE
*
      include "../commons/consts.h"
**
*     Input Variable
*
      DOUBLE COMPLEX EPS
      DOUBLE COMPLEX N
**
*     Internal Variables
*
      INTEGER Q
      DOUBLE COMPLEX HYPGEO
      DOUBLE COMPLEX CONST,CONST1,CONST2
      DOUBLE COMPLEX A1,B1,C1,A2,B2,C2,A3,B3,C3,Z
      DOUBLE COMPLEX G1,G2,G3,Z1,Z2,Z3
**
*     Output Variables
*
      DOUBLE COMPLEX DJQ(0:2)
*
      Z = 1D0 / ( 1D0 + 4D0 * EPS )
*
      DO Q=0,2
         CONST1 = EXP( - N * ZLOG( 1D0 + 4D0 * EPS ) ) 
         CONST2 = EXP( - Q * ZLOG( 1D0 + 4D0 * EPS ) )
         CONST = CONST1 * CONST2 / Z
*
         A1 = 0.5D0
         B1 = N + Q + 1D0
         C1 = N + Q + 0.5D0
*
         A2 = 0.5D0
         B2 = N + Q + 2D0
         C2 = N + Q + 1.5D0
*
         A3 = 0.5D0
         B3 = N + Q + 1D0
         C3 = N + Q + 1.5D0
*
         Z1 = N + Q + 1D0
         Z2 = 1D0 / 2D0
         Z3 = N + Q + 0.5D0
*
         CALL GAMMAL(Z1,G1)
C         CALL GAMMAL(Z2,G2)
         G2 = DCMPLX(DLOG(DSQRT(PI)),0D0)
         CALL GAMMAL(Z3,G3)
*
         DJQ(Q) = CONST * EXP( G1 + G2 - G3 ) * ( HYPGEO(A1,B1,C1,Z)
     1          - ( N + Q + 1D0 ) * HYPGEO(A2,B2,C2,Z) / (N + Q + 0.5D0)
     2          + HYPGEO(A3,B3,C3,Z) / ( N + Q + 0.5D0 ) )
      ENDDO
*
      RETURN
      END
*
************************************************************************
*
*     Leading Order O(as):
*
*        - ico = 1  F_2,g (LO, gluon)
*        - ico = 2  F_L,g (LO, gluon)
*
*     Next to Leading Order O(as^2):
*
*        - ico = 3  F_2,g (NLO, gluon)
*        - ico = 4  F_2,q (NLO, pure singlet)
*        - ico = 5  F_L,g (NLO, gluon)
*        - ico = 6  F_L,q (NLO, pure singlet)
*
*        - ico = 7  F_2,q (NLO, gluon-radiation)
*        - ico = 8  F_L,q (NLO, gluon-radiation)
*
*     * Coefficients of the terms proportional to ln(q2/m2):
*
*        - ico = 9  F_2,g (NLO, gluon, bar term)
*        - ico = 10 F_2,q (NLO, pure singlet, bar term)
*        - ico = 11 F_L,g (NLO, gluon, bar term)
*        - ico = 12 F_L,q (NLO, pure singlet, bar term)
*
*        - ico = 13 F_2,q (NLO, gluon-radiation, bar term)
*        - ico = xx F_L,q (NLO, gluon-radiation, bar term) = 0
*
************************************************************************
*
*     Linear Iterpolation
*
************************************************************************
      FUNCTION CH_FFNS_NC_AB_LIN(ICO,N,Q2,M2)
*
      IMPLICIT NONE
*
      include "../commons/coeffhqmellin.h"
      include "../commons/consts.h"
**
*     Input Variables
*
      INTEGER ICO,IXI
      DOUBLE COMPLEX Q2,M2
      DOUBLE COMPLEX N
**
*     Internal Variables
*
      INTEGER I
      DOUBLE PRECISION XI,EPS,LOGXI
      DOUBLE PRECISION DIFF(NXI),SGN
      DOUBLE PRECISION Y,YY(0:1)
      DOUBLE PRECISION TOLL
      DOUBLE COMPLEX TRANS,TRANS_TMP(0:1)
      PARAMETER(TOLL=1D-10)
**
*     Output Variables
*
      DOUBLE COMPLEX CH_FFNS_NC_AB_LIN
*
      IF(DABS(DIMAG(Q2)).GT.TOLL.OR.DABS(DIMAG(M2)).GT.TOLL)THEN
         WRITE(6,*) "In CH_FFNS_NC_AB_LIN: Q2 and M2 must be real"
         WRITE(6,*) "Q2 = ",Q2,", M2 = ",M2
         CALL EXIT(-10)
      ENDIF
*
      XI = ABS( Q2 / M2 )
      LOGXI = DLOG(XI)
      EPS = 1D0 / XI
*
*     Find IXI such that XIGRID(IXI) < XI < XIGRID(IXI+1)
*
      DIFF(1) = XI - XIGRID(1)
      DO I=2,NXI
         DIFF(I) = XI - XIGRID(I)
         SGN = DIFF(I-1) * DIFF(I)
         IF (SGN.LT.0D0) THEN
             IXI = I - 1
         ENDIF
      ENDDO
*
*     Check that the value of XI is within the range of the AB CFs or on the borders
*
      IF(XI.LT.XIMIN.OR.XI.GT.XIMAX)THEN
*
*     Set to zero the O(as2) coefficient function if XI
*     is outside the parametrization range
*
         CH_FFNS_NC_AB_LIN = (0D0,0D0)
         RETURN
      ELSEIF(XI.EQ.XIMAX)THEN
         CALL MELLIN(ICO,NXI,N-1,TRANS)
         IF(ICO.EQ.1.OR.ICO.EQ.2)THEN
            CH_FFNS_NC_AB_LIN = TRANS / ( PI * EPS )        !NLO (LO)
         ELSE
            CH_FFNS_NC_AB_LIN = 16D0 * PI * TRANS / EPS     !NNLO (NLO)
         ENDIF
         RETURN
      ELSEIF(XI.EQ.XIMIN)THEN
         CALL MELLIN(ICO,1,N-1,TRANS)
         IF(ICO.EQ.1.OR.ICO.EQ.2)THEN
            CH_FFNS_NC_AB_LIN = TRANS / ( PI * EPS )        !NLO (LO)
         ELSE
            CH_FFNS_NC_AB_LIN = 16D0 * PI * TRANS / EPS     !NNLO (NLO)
         ENDIF
         RETURN
      ENDIF
*
*     Linear Interpolation
*
      DO I=0,1
         CALL MELLIN(ICO,IXI+I,N-1,TRANS_TMP(I))
         YY(I) = DLOG(XIGRID(IXI+I))
      ENDDO
      Y = DLOG(XI)
*
      TRANS = ( TRANS_TMP(1) - TRANS_TMP(0) ) * ( Y - YY(0) ) 
     1      / ( YY(1) - YY(0) ) + TRANS_TMP(0)
*
*     Different normalization according to the perturbative order
*
      IF(ICO.EQ.1.OR.ICO.EQ.2)THEN
         CH_FFNS_NC_AB_LIN = TRANS / ( PI * EPS )        !NLO  (LO)
      ELSE
         CH_FFNS_NC_AB_LIN = 16D0 * PI * TRANS / EPS     !NNLO (NLO)
      ENDIF
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE MELLIN(ICOEF,IXI,CZ,TRANS)
*
      IMPLICIT NONE
*
      include "../commons/coeffhqmellin.h"
**
*     Input Variables
*
      INTEGER ICOEF,IXI
      DOUBLE COMPLEX CZ
**
*     Internal Variables
*
      INTEGER I
      DOUBLE COMPLEX XTH
      DOUBLE COMPLEX A,B
      DOUBLE COMPLEX BETA,BETAC
**
*     Output Variables
*
      DOUBLE COMPLEX TRANS
*
      XTH = 1D0 / ( 4D0 / XIGRID(IXI) + 1D0 )
      A = CZ - COEF_P1(ICOEF)
      B = 1D0 - COEF_P2(ICOEF)
      BETA = BETAC(A,B)
*
      TRANS = (0D0,0D0)
*
      DO I=0,M_COEF(ICOEF)-1
         TRANS = TRANS + COEF(ICOEF,IXI,I+1) * XTH**(A + B - 1D0) * BETA
         BETA = BETA * A / ( A + B )
         A = A + 1D0
      ENDDO
*
      RETURN
      END
*
c$$$***********************************************************************
c$$$*
c$$$*     Hermite Cubic Iterpolation
c$$$*
c$$$***********************************************************************
c$$$      FUNCTION CH_FFNS_NC_AB_CUB(ICO,N,Q2,M2)
c$$$*
c$$$      IMPLICIT NONE
c$$$*
c$$$      include "../commons/coeffhqmellin.h"
c$$$      include "../commons/consts.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      INTEGER ICO,IXI
c$$$      DOUBLE PRECISION Q2,M2
c$$$      DOUBLE COMPLEX N
c$$$**
c$$$*     Internal Variables
c$$$*
c$$$      INTEGER I
c$$$      DOUBLE PRECISION XI,EPS,LOGXI
c$$$      DOUBLE PRECISION DIFF(NXI),SGN
c$$$      DOUBLE PRECISION Y,T,YY(-1:2),H(-1:2)
c$$$      DOUBLE PRECISION H00,H10,H01,H11
c$$$      DOUBLE PRECISION A,B,C,D
c$$$      DOUBLE COMPLEX TRANS,TRANS_TMP(-1:2)
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      DOUBLE COMPLEX CH_FFNS_NC_AB_CUB
c$$$*
c$$$      XI = Q2 / M2
c$$$      LOGXI = DLOG(XI)
c$$$      EPS = 1D0 / XI
c$$$*
c$$$*     Find IXI such that XIGRID(IXI) < XI < XIGRID(IXI+1)
c$$$*
c$$$      DIFF(1) = XI - XIGRID(1)
c$$$      DO I=2,NXI
c$$$         DIFF(I) = XI - XIGRID(I)
c$$$         SGN = DIFF(I-1) * DIFF(I)
c$$$         IF (SGN.LT.0D0) THEN
c$$$             IXI = I
c$$$         ENDIF
c$$$      ENDDO
c$$$*
c$$$*     Check that the value of XI is within the range of the AB CFs or on the borders
c$$$*
c$$$      IF(XI.LT.XIMIN.OR.XI.GT.XIMAX)THEN
c$$$C         WRITE(6,*) "In modules/evolutionQCD/wc_FFNS_NC.f:"
c$$$C         WRITE(6,*) "Alekhin and Bluemlein CFs out of range:"
c$$$C         WRITE(6,*) "Q2, M2 = ",Q2,M2
c$$$C         WRITE(6,*) XIMIN,"< XI=Q2/M2 <",XIMAX,", while XI =",XI
c$$$C         CALL EXIT(-10)
c$$$*
c$$$*     Set to zero the O(as2) coefficient function if the range
c$$$*     is outside that of the parametrization
c$$$*
c$$$         CH_FFNS_NC_AB_CUB = DCMPLX(0D0,0D0)
c$$$         RETURN
c$$$      ELSEIF(XI.EQ.XIMAX)THEN
c$$$         CALL MELLIN(ICO,NXI,N-1,TRANS)
c$$$         IF(ICO.EQ.1.OR.ICO.EQ.2)THEN
c$$$            CH_FFNS_NC_AB_CUB = TRANS / ( PI * EPS )        !NLO (LO)
c$$$         ELSE
c$$$            CH_FFNS_NC_AB_CUB = 16D0 * PI * TRANS / EPS     !NNLO (NLO)
c$$$         ENDIF
c$$$         RETURN
c$$$      ELSEIF(XI.EQ.XIMIN)THEN
c$$$         CALL MELLIN(ICO,1,N-1,TRANS)
c$$$         IF(ICO.EQ.1.OR.ICO.EQ.2)THEN
c$$$            CH_FFNS_NC_AB_CUB = TRANS / ( PI * EPS )        !NLO (LO)
c$$$         ELSE
c$$$            CH_FFNS_NC_AB_CUB = 16D0 * PI * TRANS / EPS     !NNLO (NLO)
c$$$         ENDIF
c$$$         RETURN
c$$$      ENDIF
c$$$*
c$$$*     Cubic Hermite Interpolation
c$$$*
c$$$      IF(IXI.GT.1.AND.IXI.LT.NXI-2)THEN
c$$$         DO I=-1,2
c$$$            CALL MELLIN(ICO,IXI+I,N-1,TRANS_TMP(I))
c$$$            YY(I) = DLOG(XIGRID(IXI+I))
c$$$            H(I) = DLOG(XIGRID(IXI+I+1)) - DLOG(XIGRID(IXI+I))
c$$$         ENDDO
c$$$      ELSEIF(IXI.EQ.1)THEN
c$$$         DO I=0,2
c$$$            CALL MELLIN(ICO,IXI+I,N-1,TRANS_TMP(I))
c$$$            YY(I) = DLOG(XIGRID(IXI+I))
c$$$            H(I) = DLOG(XIGRID(IXI+I+1)) - DLOG(XIGRID(IXI+I))
c$$$         ENDDO
c$$$      ELSEIF(IXI.EQ.NXI-2)THEN
c$$$         DO I=-1,1
c$$$            CALL MELLIN(ICO,IXI+I,N-1,TRANS_TMP(I))
c$$$            YY(I) = DLOG(XIGRID(IXI+I))
c$$$            H(I) = DLOG(XIGRID(IXI+I+1)) - DLOG(XIGRID(IXI+I))
c$$$         ENDDO
c$$$      ELSEIF(IXI.EQ.NXI-1)THEN
c$$$         DO I=-1,0
c$$$            CALL MELLIN(ICO,IXI+I,N-1,TRANS_TMP(I))
c$$$            YY(I) = DLOG(XIGRID(IXI+I))
c$$$            H(I) = DLOG(XIGRID(IXI+I+1)) - DLOG(XIGRID(IXI+I))
c$$$         ENDDO
c$$$      ENDIF
c$$$      Y = DLOG(XI)
c$$$      T = ( Y - YY(0) ) / H(0)
c$$$*
c$$$      H00 = 2D0 * T**3D0 - 3D0 * T**2D0 + 1D0
c$$$      H10 = T**3D0 -2D0 * T**2D0 + T
c$$$      H01 = - 2D0 * T**3D0 + 3d0* T**2D0
c$$$      H11 = T**3D0 -T**2D0
c$$$*
c$$$      A = 0D0 
c$$$      B = 0D0
c$$$      C = 0D0
c$$$      D = 0D0
c$$$*
c$$$*     A(Y)
c$$$*
c$$$      IF(IXI.GE.2)THEN
c$$$         A = - 0.5D0 * H10 * ( H(0) / H(-1) ) 
c$$$      ENDIF
c$$$*
c$$$*     B(Y)
c$$$*
c$$$      IF(IXI.EQ.1)THEN
c$$$         B = H00 - H10 - 0.5D0 * H11
c$$$      ELSEIF(IXI.EQ.NXI-1)THEN
c$$$         B = H00 - 0.5D0 * H10 * ( 1D0 - H(0) / H(-1) ) - H11
c$$$      ELSE
c$$$         B = H00 - 0.5D0 * H10 * ( 1D0 - H(0) / H(-1) ) - 0.5D0 * H11
c$$$      ENDIF
c$$$*
c$$$*     C(Y)
c$$$*
c$$$      IF(IXI.EQ.1)THEN
c$$$         C = H01 + 0.5D0 * H11 * ( 1D0 - H(0) / H(1) ) + H10
c$$$      ELSEIF(IXI.EQ.NXI-1)THEN
c$$$         C = H01 + H11 + 0.5D0 *  H10
c$$$      ELSE
c$$$         C = H01 + 0.5D0 * H11 * ( 1D0 - H(0) / H(1) ) + 0.5D0 * H10
c$$$      ENDIF
c$$$*
c$$$*     D(Y)
c$$$*
c$$$      IF(IXI.LE.NXI-2)THEN
c$$$         D = 0.5D0 * H11 * ( H(0) / H(1) )
c$$$      ENDIF 
c$$$*      
c$$$      TRANS = ( YY(-1) * TRANS_TMP(-1) * A 
c$$$     1      + YY(0) * TRANS_TMP(0) * B 
c$$$     2      + YY(1) * TRANS_TMP(1) * C 
c$$$     3      + YY(2) * TRANS_TMP(2) * D ) / Y
c$$$*
c$$$*     Different normalization according to the perturbative order
c$$$*
c$$$      IF(ICO.EQ.1.OR.ICO.EQ.2)THEN
c$$$         CH_FFNS_NC_AB_CUB = TRANS / ( PI * EPS )        !NLO  (LO)
c$$$      ELSE
c$$$         CH_FFNS_NC_AB_CUB = 16D0 * PI * TRANS / EPS     !NNLO (NLO)
c$$$      ENDIF
c$$$*
c$$$      RETURN
c$$$      END
c$$$*
