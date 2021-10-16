************************************************************************
*
*     No crossing:
*
*     - CROSS01 --> Q2I,Q2F < M_E^2           (or FFNS with NF = 0)
*     - CROSS02 --> M_E^2 < Q2I, Q2F < M_M^2  (or FFNS with NF = 1)
*     - CROSS03 --> M_M^2 < Q2I, Q2F < M_T^2  (or FFNS with NF = 2)
*     - CROSS04 --> Q2I,Q2F > M_T^2           (or FFNS with NF = 3)
*
*     One crossing:
*
*     - CROSS11 --> Q2I < M_E^2 < Q2F
*     - CROSS12 --> Q2I < M_M^2 < Q2F
*     - CROSS13 --> Q2I < M_T^2 < Q2F
*
*     Two crossings:
*
*     - CROSS21 --> Q2I < M_E^2 < M_M^2 < Q2F
*     - CROSS22 --> Q2I < M_M^2 < M_T^2 < Q2F
*
*     Three crossings:
*
*     - CROSS31 --> Q2I < M_E^2 < M_M^2 < M_T^2 < Q2F
*
************************************************************************
*
*     NO CROSSING
*
************************************************************************
      SUBROUTINE CROSS01(ZN,Q2I,Q2F,NFI,EFNNS,EFNSG,
     1                   ZFUNCNS3,ZFUNCNS8,
     2                   ZFUNCNSV3,ZFUNCNSV8)
*
      IMPLICIT NONE
**
*     Input Variables 
*     
      INTEGER NFI
      DOUBLE COMPLEX ZN
      DOUBLE COMPLEX Q2I,Q2F
**
*     Output Variables
* 
      DOUBLE COMPLEX EFNNS(3),EFNSG(2,2)
      DOUBLE COMPLEX ZFUNCNS3(2),ZFUNCNS8(2)
      DOUBLE COMPLEX ZFUNCNSV3,ZFUNCNSV8
*
      EFNNS(1) = (1D0, 0D0)
      EFNNS(2) = (1D0, 0D0)
      EFNNS(3) = (1D0, 0D0)
*
      EFNSG(1,1) = (1D0, 0D0)
      EFNSG(1,2) = (0D0, 0D0)
      EFNSG(2,1) = (0D0, 0D0)
      EFNSG(2,2) = (1D0, 0D0)
*
      ZFUNCNS3(1) = EFNSG(1,1)
      ZFUNCNS3(2) = EFNSG(1,2)
      ZFUNCNS8(1) = EFNSG(1,1)
      ZFUNCNS8(2) = EFNSG(1,2)
      ZFUNCNSV3   = EFNNS(3)
      ZFUNCNSV8   = EFNNS(3)
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE CROSS02(ZN,Q2I,Q2F,NFI,EFNNS,EFNSG,
     1                   ZFUNCNS3,ZFUNCNS8,
     2                   ZFUNCNSV3,ZFUNCNSV8)
*
      IMPLICIT NONE
**
*     Input Variables 
*     
      INTEGER NFI
      DOUBLE COMPLEX ZN
      DOUBLE COMPLEX Q2I,Q2F
**
*     Output Variables
* 
      DOUBLE COMPLEX EFNNS(3),EFNSG(2,2)
      DOUBLE COMPLEX ZFUNCNS3(2),ZFUNCNS8(2)
      DOUBLE COMPLEX ZFUNCNSV3,ZFUNCNSV8
*
      CALL EVOLFACTN(ZN,Q2I,Q2F,NFI,EFNNS,EFNSG)
*
      ZFUNCNS3(1) = EFNSG(1,1)
      ZFUNCNS3(2) = EFNSG(1,2)
      ZFUNCNS8(1) = EFNSG(1,1)
      ZFUNCNS8(2) = EFNSG(1,2)
      ZFUNCNSV3   = EFNNS(3)
      ZFUNCNSV8   = EFNNS(3)
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE CROSS03(ZN,Q2I,Q2F,NFI,EFNNS,EFNSG,
     1                   ZFUNCNS3,ZFUNCNS8,
     2                   ZFUNCNSV3,ZFUNCNSV8)
*
      IMPLICIT NONE
**
*     Input Variables 
*     
      INTEGER NFI
      DOUBLE COMPLEX ZN
      DOUBLE COMPLEX Q2I,Q2F
**
*     Output Variables
* 
      DOUBLE COMPLEX EFNNS(3),EFNSG(2,2)
      DOUBLE COMPLEX ZFUNCNS3(2),ZFUNCNS8(2)
      DOUBLE COMPLEX ZFUNCNSV3,ZFUNCNSV8
*
      CALL EVOLFACTN(ZN,Q2I,Q2F,NFI,EFNNS,EFNSG)
*
      ZFUNCNS3(1) = EFNNS(1)
      ZFUNCNS3(2) = (0D0, 0D0)
      ZFUNCNS8(1) = EFNSG(1,1)
      ZFUNCNS8(2) = EFNSG(1,2)
      ZFUNCNSV3   = EFNNS(2)
      ZFUNCNSV8   = EFNNS(3)
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE CROSS04(ZN,Q2I,Q2F,NFI,EFNNS,EFNSG,
     1                   ZFUNCNS3,ZFUNCNS8,
     2                   ZFUNCNSV3,ZFUNCNSV8)
*
      IMPLICIT NONE
**
*     Input Variables 
*     
      INTEGER NFI
      DOUBLE COMPLEX ZN
      DOUBLE COMPLEX Q2I,Q2F
**
*     Output Variables
* 
      DOUBLE COMPLEX EFNNS(3),EFNSG(2,2)
      DOUBLE COMPLEX ZFUNCNS3(2),ZFUNCNS8(2)
      DOUBLE COMPLEX ZFUNCNSV3,ZFUNCNSV8
*
      CALL EVOLFACTN(ZN,Q2I,Q2F,NFI,EFNNS,EFNSG)
*
      ZFUNCNS3(1) = EFNNS(1)
      ZFUNCNS3(2) = (0D0, 0D0)
      ZFUNCNS8(1) = EFNNS(1)
      ZFUNCNS8(2) = (0D0, 0D0)
      ZFUNCNSV3   = EFNNS(2)
      ZFUNCNSV8   = EFNNS(2)
*
      RETURN
      END
*
************************************************************************
*
*     1 CROSSING
*
************************************************************************
      SUBROUTINE CROSS11(ZN,Q2I,Q2F,NFI,EFNNS,EFNSG,
     1                   ZFUNCNS3,ZFUNCNS8,
     2                   ZFUNCNSV3,ZFUNCNSV8)
*
      IMPLICIT NONE
*
      include "../commons/alpha.h"
      include "../commons/massthrs.h"
      include "../commons/ipt.h"
      include "../commons/evol.h"
**
*     Input Variables 
*     
      INTEGER NFI
      DOUBLE COMPLEX ZN
      DOUBLE COMPLEX Q2I,Q2F
**
*     Internal Variables
* 
      INTEGER I,J,K,L
      DOUBLE COMPLEX ATH(3)
      DOUBLE COMPLEX MCSG(2,2),MCNS,VCNS(2)
      DOUBLE COMPLEX EFNNS1(3),EFNSG1(2,2),EFNNS2(3),EFNSG2(2,2)
      DOUBLE COMPLEX A2SG(2,2),A1SG,A2NS
**
*     Output Variables
* 
      DOUBLE COMPLEX EFNNS(3),EFNSG(2,2)
      DOUBLE COMPLEX ZFUNCNS3(2),ZFUNCNS8(2)
      DOUBLE COMPLEX ZFUNCNSV3,ZFUNCNSV8
*
      ATH(1) = AE
      ATH(2) = AM
      ATH(3) = AT
*
      CALL EVOLFACTN(ZN,Q2I,Q2TH(1),NFI,EFNNS1,EFNSG1)
      IF(IPT.EQ.1.AND.EVOL.EQ."TIME")THEN
         CALL MATCHCOEF1SGT(ZN,A1SG)
*
         MCSG(1,1) = (1D0,0D0)
         MCSG(1,2) = ATH(NFI+1) * A1SG
         MCSG(2,1) = (0D0,0D0)
         MCSG(2,2) = (1D0,0D0)
*
         MCNS = (1D0,0D0)
*
         VCNS(1) = (1D0,0D0)
         VCNS(2) = (0D0,0D0)
      ELSE
         MCSG(1,1) = (1D0,0D0)
         MCSG(1,2) = (0D0,0D0)
         MCSG(2,1) = (0D0,0D0)
         MCSG(2,2) = (1D0,0D0)
*
         MCNS = (1D0,0D0)
*
         VCNS(1) = (1D0,0D0)
         VCNS(2) = (0D0,0D0)
      ENDIF
      CALL EVOLFACTN(ZN,Q2TH(1),Q2F,NFI+1,EFNNS2,EFNSG2)
*
      DO I=1,2
         DO J=1,2
            EFNSG(I,J) = (0D0,0D0)
         ENDDO
         ZFUNCNS3(I) = (0D0,0D0)
         ZFUNCNS8(I) = (0D0,0D0)
      ENDDO
*     Singlet and gluon
      DO I=1,2
         DO J=1,2
            DO K=1,2
               DO L=1,2
                  EFNSG(I,J) = EFNSG(I,J) + EFNSG2(I,K) * MCSG(K,L) * 
     1                                      EFNSG1(L,J)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
*     Plus, Minus, Valence
      DO I=1,3
         EFNNS(I) = EFNNS2(I) * MCNS * EFNNS1(I)
      ENDDO
*     T_{3, 8}
      DO I=1,2
         DO J=1,2
            DO K=1,2
               ZFUNCNS3(I) = ZFUNCNS3(I) + EFNSG2(1,J) * MCSG(J,K) * 
     1                                     EFNSG1(K,I)
               ZFUNCNS8(I) = ZFUNCNS8(I) + EFNSG2(1,J) * MCSG(J,K) * 
     1                                     EFNSG1(K,I)
            ENDDO
         ENDDO
      ENDDO
*     V_{3, 8}
      ZFUNCNSV3 = EFNNS2(3) * MCNS * EFNNS1(3)
      ZFUNCNSV8 = EFNNS2(3) * MCNS * EFNNS1(3)
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE CROSS12(ZN,Q2I,Q2F,NFI,EFNNS,EFNSG,
     1                   ZFUNCNS3,ZFUNCNS8,
     2                   ZFUNCNSV3,ZFUNCNSV8)
*
      IMPLICIT NONE
*
      include "../commons/alpha.h"
      include "../commons/massthrs.h"
      include "../commons/ipt.h"
      include "../commons/evol.h"
**
*     Input Variables 
*     
      INTEGER NFI
      DOUBLE COMPLEX ZN
      DOUBLE COMPLEX Q2I,Q2F
**
*     Internal Variables
* 
      INTEGER I,J,K,L
      DOUBLE COMPLEX ATH(3)
      DOUBLE COMPLEX MCSG(2,2),MCNS,VCNS(2)
      DOUBLE COMPLEX EFNNS1(3),EFNSG1(2,2),EFNNS2(3),EFNSG2(2,2)
      DOUBLE COMPLEX A2SG(2,2),A1SG,A2NS
**
*     Output Variables
* 
      DOUBLE COMPLEX EFNNS(3),EFNSG(2,2)
      DOUBLE COMPLEX ZFUNCNS3(2),ZFUNCNS8(2)
      DOUBLE COMPLEX ZFUNCNSV3,ZFUNCNSV8
*
      ATH(1) = AE
      ATH(2) = AM
      ATH(3) = AT
*
      CALL EVOLFACTN(ZN,Q2I,Q2TH(2),NFI,EFNNS1,EFNSG1)
      IF(IPT.EQ.1.AND.EVOL.EQ."TIME")THEN
         CALL MATCHCOEF1SGT(ZN,A1SG)
*
         MCSG(1,1) = (1D0,0D0)
         MCSG(1,2) = ATH(NFI+1) * A1SG
         MCSG(2,1) = (0D0,0D0)
         MCSG(2,2) = (1D0,0D0)
*
         MCNS = (1D0,0D0)
*
         VCNS(1) = (1D0,0D0)
         VCNS(2) = - ATH(NFI+1) * A1SG
      ELSE
         MCSG(1,1) = (1D0,0D0)
         MCSG(1,2) = (0D0,0D0)
         MCSG(2,1) = (0D0,0D0)
         MCSG(2,2) = (1D0,0D0)
*
         MCNS = (1D0,0D0)
*
         VCNS(1) = (1D0,0D0)
         VCNS(2) = (0D0,0D0)
      ENDIF
      CALL EVOLFACTN(ZN,Q2TH(2),Q2F,NFI+1,EFNNS2,EFNSG2)
*
      DO I=1,2
         DO J=1,2
            EFNSG(I,J) = (0D0,0D0)
         ENDDO
         ZFUNCNS3(I) = (0D0,0D0)
         ZFUNCNS8(I) = (0D0,0D0)
      ENDDO
*     Singlet and gluon
      DO I=1,2
         DO J=1,2
            DO K=1,2
               DO L=1,2
                  EFNSG(I,J) = EFNSG(I,J) + EFNSG2(I,K) * MCSG(K,L) * 
     1                                      EFNSG1(L,J)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
*     Plus, Minus, Valence
      DO I=1,3
         EFNNS(I) = EFNNS2(I) * MCNS * EFNNS1(I)
      ENDDO
*     T_3
      DO I=1,2
         DO J=1,2
            ZFUNCNS3(I) = ZFUNCNS3(I) + EFNNS2(1) * VCNS(J) * 
     1                                  EFNSG1(J,I)
         ENDDO
      ENDDO
*     T_8
      DO I=1,2
         DO J=1,2
            DO K=1,2
               ZFUNCNS8(I) = ZFUNCNS8(I) + EFNSG2(1,J) * MCSG(J,K) * 
     1                                     EFNSG1(K,I)
            ENDDO
         ENDDO
      ENDDO
*     V_{3, 8}
      ZFUNCNSV3 = EFNNS2(2) * MCNS * EFNNS1(3)
      ZFUNCNSV8 = EFNNS2(3) * MCNS * EFNNS1(3)
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE CROSS13(ZN,Q2I,Q2F,NFI,EFNNS,EFNSG,
     1                   ZFUNCNS3,ZFUNCNS8,
     2                   ZFUNCNSV3,ZFUNCNSV8)
*
      IMPLICIT NONE
*
      include "../commons/alpha.h"
      include "../commons/massthrs.h"
      include "../commons/ipt.h"
      include "../commons/evol.h"
**
*     Input Variables 
*     
      INTEGER NFI
      DOUBLE COMPLEX ZN
      DOUBLE COMPLEX Q2I,Q2F
**
*     Internal Variables
* 
      INTEGER I,J,K,L
      DOUBLE COMPLEX ATH(3)
      DOUBLE COMPLEX MBSG(2,2),MBNS,VBNS(2)
      DOUBLE COMPLEX EFNNS1(3),EFNSG1(2,2),EFNNS2(3),EFNSG2(2,2)
      DOUBLE COMPLEX A2SG(2,2),A1SG,A2NS
**
*     Output Variables
* 
      DOUBLE COMPLEX EFNNS(3),EFNSG(2,2)
      DOUBLE COMPLEX ZFUNCNS3(2),ZFUNCNS8(2)
      DOUBLE COMPLEX ZFUNCNSV3,ZFUNCNSV8
*
      ATH(1) = AE
      ATH(2) = AM
      ATH(3) = AT
*
      CALL EVOLFACTN(ZN,Q2I,Q2TH(3),NFI,EFNNS1,EFNSG1)
      IF(IPT.EQ.1.AND.EVOL.EQ."TIME")THEN
         CALL MATCHCOEF1SGT(ZN,A1SG)
*
         MBSG(1,1) = (1D0,0D0)
         MBSG(1,2) = ATH(NFI+1) * A1SG
         MBSG(2,1) = (0D0,0D0)
         MBSG(2,2) = (1D0,0D0)
*
         MBNS = (1D0,0D0)
*
         VBNS(1) = (1D0,0D0)
         VBNS(2) = - 2D0 * ATH(NFI+1) * A1SG
      ELSE
         MBSG(1,1) = (1D0,0D0)
         MBSG(1,2) = (0D0,0D0)
         MBSG(2,1) = (0D0,0D0)
         MBSG(2,2) = (1D0,0D0)
*
         MBNS = (1D0,0D0)
*
         VBNS(1) = (1D0,0D0)
         VBNS(2) = (0D0,0D0)
      ENDIF
      CALL EVOLFACTN(ZN,Q2TH(3),Q2F,NFI+1,EFNNS2,EFNSG2)
*
      DO I=1,2
         DO J=1,2
            EFNSG(I,J) = (0D0,0D0)
         ENDDO
         ZFUNCNS3(I) = (0D0,0D0)
         ZFUNCNS8(I) = (0D0,0D0)
      ENDDO
*     Singlet and gluon
      DO I=1,2
         DO J=1,2
            DO K=1,2
               DO L=1,2
                  EFNSG(I,J) = EFNSG(I,J) + EFNSG2(I,K) * MBSG(K,L) * 
     1                                      EFNSG1(L,J)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
*     Plus, Minus, Valence
      DO I=1,3
         EFNNS(I) = EFNNS2(I) * MBNS * EFNNS1(I)
      ENDDO
*     T_3
      ZFUNCNS3(1) = EFNNS2(1) * MBNS * EFNNS1(1)
      ZFUNCNS3(2) = (0D0,0D0)
*     T_8
      DO I=1,2
         DO J=1,2
            ZFUNCNS8(I) = ZFUNCNS8(I) + EFNNS2(1) * VBNS(J) * 
     1                                  EFNSG1(J,I)
         ENDDO
      ENDDO
*     V_{3, 8}
      ZFUNCNSV3 = EFNNS2(2) * MBNS * EFNNS1(2)
      ZFUNCNSV8 = EFNNS2(2) * MBNS * EFNNS1(3)
*
      RETURN
      END
*
************************************************************************
*
*     2 CROSSINGS
*
************************************************************************
      SUBROUTINE CROSS21(ZN,Q2I,Q2F,NFI,EFNNS,EFNSG,
     1                   ZFUNCNS3,ZFUNCNS8,
     2                   ZFUNCNSV3,ZFUNCNSV8)
*
      IMPLICIT NONE
*
      include "../commons/alpha.h"
      include "../commons/massthrs.h"
      include "../commons/ipt.h"
      include "../commons/evol.h"
**
*     Input Variables 
*     
      INTEGER NFI
      DOUBLE COMPLEX ZN
      DOUBLE COMPLEX Q2I,Q2F
**
*     Internal Variables
* 
      INTEGER I,J,K,L,G,H
      DOUBLE COMPLEX ATH(3)
      DOUBLE COMPLEX MCSG(2,2),MCNS,VCNS(2)
      DOUBLE COMPLEX MBSG(2,2),MBNS,VBNS(2)
      DOUBLE COMPLEX EFNNS1(3),EFNSG1(2,2),EFNNS2(3),EFNSG2(2,2)
      DOUBLE COMPLEX EFNNS3(3),EFNSG3(2,2)
      DOUBLE COMPLEX A2SG(2,2),A1SG,A2NS
**
*     Output Variables
* 
      DOUBLE COMPLEX EFNNS(3),EFNSG(2,2)
      DOUBLE COMPLEX ZFUNCNS3(2),ZFUNCNS8(2)
      DOUBLE COMPLEX ZFUNCNSV3,ZFUNCNSV8
*
      ATH(1) = AE
      ATH(2) = AM
      ATH(3) = AT
*
      CALL EVOLFACTN(ZN,Q2I,Q2TH(1),NFI,EFNNS1,EFNSG1)
      IF(IPT.EQ.1.AND.EVOL.EQ."TIME")THEN
         CALL MATCHCOEF1SGT(ZN,A1SG)
*
         MCSG(1,1) = (1D0,0D0)
         MCSG(1,2) = ATH(NFI+1) * A1SG
         MCSG(2,1) = (0D0,0D0)
         MCSG(2,2) = (1D0,0D0)
*
         MCNS = (1D0,0D0)
*
         VCNS(1) = (1D0,0D0)
         VCNS(2) = (0D0,0D0)
      ELSE
         MCSG(1,1) = (1D0,0D0)
         MCSG(1,2) = (0D0,0D0)
         MCSG(2,1) = (0D0,0D0)
         MCSG(2,2) = (1D0,0D0)
*
         MCNS = (1D0,0D0)
*
         VCNS(1) = (1D0,0D0)
         VCNS(2) = (0D0,0D0)
      ENDIF
      CALL EVOLFACTN(ZN,Q2TH(1),Q2TH(2),NFI+1,EFNNS2,EFNSG2)
      IF(IPT.EQ.1.AND.EVOL.EQ."TIME")THEN
         CALL MATCHCOEF1SGT(ZN,A1SG)
*
         MBSG(1,1) = (1D0,0D0)
         MBSG(1,2) = ATH(NFI+2) * A1SG
         MBSG(2,1) = (0D0,0D0)
         MBSG(2,2) = (1D0,0D0)
*
         MBNS = (1D0,0D0)
*
         VBNS(1) = (1D0,0D0)
         VBNS(2) = - ATH(NFI+2) * A1SG
      ELSE
         MBSG(1,1) = (1D0,0D0)
         MBSG(1,2) = (0D0,0D0)
         MBSG(2,1) = (0D0,0D0)
         MBSG(2,2) = (1D0,0D0)
*
         MBNS = (1D0,0D0)
*
         VBNS(1) = (1D0,0D0)
         VBNS(2) = (0D0,0D0)
      ENDIF
      CALL EVOLFACTN(ZN,Q2TH(2),Q2F,NFI+2,EFNNS3,EFNSG3)
*
      DO I=1,2
         DO J=1,2
            EFNSG(I,J) = (0D0,0D0)
         ENDDO
         ZFUNCNS3(I) = (0D0,0D0)
         ZFUNCNS8(I) = (0D0,0D0)
      ENDDO
*     Singlet and gluon
      DO I=1,2
         DO J=1,2
            DO K=1,2
               DO L=1,2
                  DO H=1,2
                     DO G=1,2
                        EFNSG(I,J) = EFNSG(I,J) +
     1                               EFNSG3(I,K) * MBSG(K,L) * 
     2                               EFNSG2(L,H) * MCSG(H,G) * 
     3                               EFNSG1(G,J)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
*     Plus, Minus, Valence
      DO I=1,3
         EFNNS(I) = EFNNS3(I) * MBNS * EFNNS2(I) * MCNS * EFNNS1(I)
      ENDDO
*     T_3
      DO I=1,2
         DO J=1,2
            DO K=1,2
               DO L=1,2
                  ZFUNCNS3(I) = ZFUNCNS3(I) + EFNNS3(1) * VBNS(J) * 
     1                                        EFNSG2(J,K) * MCSG(K,L)*
     2                                        EFNSG1(L,I)
               ENDDO
            ENDDO 
         ENDDO
      ENDDO
*     T_8
      DO I=1,2
         DO J=1,2
            DO K=1,2
               DO L=1,2
                  DO H=1,2
                     ZFUNCNS8(I) = ZFUNCNS8(I) + 
     1                             EFNSG3(1,J) * MBSG(J,K) * 
     2                             EFNSG2(K,L) * MCSG(L,H) *
     3                             EFNSG1(H,I)
                  ENDDO
               ENDDO
            ENDDO 
         ENDDO
      ENDDO
*     V_{3, 8}
      ZFUNCNSV3 = EFNNS3(2) * MBNS * EFNNS2(3) * MCNS * EFNNS1(3)
      ZFUNCNSV8 = EFNNS3(3) * MBNS * EFNNS2(3) * MCNS * EFNNS1(3)
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE CROSS22(ZN,Q2I,Q2F,NFI,EFNNS,EFNSG,
     1                   ZFUNCNS3,ZFUNCNS8,
     2                   ZFUNCNSV3,ZFUNCNSV8)
*
      IMPLICIT NONE
*
      include "../commons/alpha.h"
      include "../commons/massthrs.h"
      include "../commons/ipt.h"
      include "../commons/evol.h"
**
*     Input Variables 
*     
      INTEGER NFI
      DOUBLE COMPLEX ZN
      DOUBLE COMPLEX Q2I,Q2F
**
*     Internal Variables
* 
      INTEGER I,J,K,L,G,H
      DOUBLE COMPLEX ATH(3)
      DOUBLE COMPLEX MBSG(2,2),MBNS,VBNS(2)
      DOUBLE COMPLEX MTSG(2,2),MTNS,VTNS(2)
      DOUBLE COMPLEX EFNNS1(3),EFNSG1(2,2),EFNNS2(3),EFNSG2(2,2)
      DOUBLE COMPLEX EFNNS3(3),EFNSG3(2,2)
      DOUBLE COMPLEX A2SG(2,2),A1SG,A2NS
**
*     Output Variables
* 
      DOUBLE COMPLEX EFNNS(3),EFNSG(2,2)
      DOUBLE COMPLEX ZFUNCNS3(2),ZFUNCNS8(2)
      DOUBLE COMPLEX ZFUNCNSV3,ZFUNCNSV8
*
      ATH(1) = AE
      ATH(2) = AM
      ATH(3) = AT
*
      CALL EVOLFACTN(ZN,Q2I,Q2TH(2),NFI,EFNNS1,EFNSG1)
      IF(IPT.EQ.1.AND.EVOL.EQ."TIME")THEN
         CALL MATCHCOEF1SGT(ZN,A1SG)
*
         MBSG(1,1) = (1D0,0D0)
         MBSG(1,2) = ATH(NFI+1) * A1SG
         MBSG(2,1) = (0D0,0D0)
         MBSG(2,2) = (1D0,0D0)
*
         MBNS = (1D0,0D0)
*
         VBNS(1) = (1D0,0D0)
         VBNS(2) = - ATH(NFI+1) * A1SG
      ELSE
         MBSG(1,1) = (1D0,0D0)
         MBSG(1,2) = (0D0,0D0)
         MBSG(2,1) = (0D0,0D0)
         MBSG(2,2) = (1D0,0D0)
*
         MBNS = (1D0,0D0)
*
         VBNS(1) = (1D0,0D0)
         VBNS(2) = (0D0,0D0)
      ENDIF
      CALL EVOLFACTN(ZN,Q2TH(2),Q2TH(3),NFI+1,EFNNS2,EFNSG2)
      IF(IPT.EQ.1.AND.EVOL.EQ."TIME")THEN
         CALL MATCHCOEF1SGT(ZN,A1SG)
*
         MTSG(1,1) = (1D0,0D0)
         MTSG(1,2) = ATH(NFI+2) * A1SG
         MTSG(2,1) = (0D0,0D0)
         MTSG(2,2) = (1D0,0D0)
*
         MTNS = (1D0,0D0)
*
         VTNS(1) = (1D0,0D0)
         VTNS(2) = - 2D0 * ATH(NFI+2) * A1SG
      ELSE
         MTSG(1,1) = (1D0,0D0)
         MTSG(1,2) = (0D0,0D0)
         MTSG(2,1) = (0D0,0D0)
         MTSG(2,2) = (1D0,0D0)
*
         MTNS = (1D0,0D0)
*
         VTNS(1) = (1D0,0D0)
         VTNS(2) = (0D0,0D0)
      ENDIF
      CALL EVOLFACTN(ZN,Q2TH(3),Q2F,NFI+2,EFNNS3,EFNSG3)
*
      DO I=1,2
         DO J=1,2
            EFNSG(I,J) = (0D0,0D0)
         ENDDO
         ZFUNCNS3(I) = (0D0,0D0)
         ZFUNCNS8(I) = (0D0,0D0)
      ENDDO
*     Singlet and gluon
      DO I=1,2
         DO J=1,2
            DO K=1,2
               DO L=1,2
                  DO H=1,2
                     DO G=1,2
                        EFNSG(I,J) = EFNSG(I,J) +
     1                               EFNSG3(I,K) * MTSG(K,L) * 
     2                               EFNSG2(L,H) * MBSG(H,G) * 
     3                               EFNSG1(G,J)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
*     Plus, Minus, Valence
      DO I=1,3
         EFNNS(I) = EFNNS3(I) * MTNS * EFNNS2(I) * MBNS * EFNNS1(I)
      ENDDO
*     T_3
      DO I=1,2
         DO J=1,2
            ZFUNCNS3(I) = ZFUNCNS3(I) + EFNNS3(1) * MTNS * 
     1                                  EFNNS2(1) * VBNS(J) * 
     2                                  EFNSG1(J,I)
         ENDDO
      ENDDO
*     T_8
      DO I=1,2
         DO J=1,2
            DO K=1,2
               DO L=1,2
                  ZFUNCNS8(I) = ZFUNCNS8(I) + EFNNS3(1) * VTNS(J) * 
     1                                        EFNSG2(J,K) * MBSG(K,L)*
     2                                        EFNSG1(L,I)
               ENDDO
            ENDDO 
         ENDDO
      ENDDO
*     V_{3, 8}
      ZFUNCNSV3 = EFNNS3(2) * MTNS * EFNNS2(2) * MBNS * EFNNS1(3)
      ZFUNCNSV8 = EFNNS3(2) * MTNS * EFNNS2(3) * MBNS * EFNNS1(3)
*
      RETURN
      END
*
************************************************************************
*
*     3 CROSSINGS
*
************************************************************************
      SUBROUTINE CROSS31(ZN,Q2I,Q2F,NFI,EFNNS,EFNSG,
     1                   ZFUNCNS3,ZFUNCNS8,
     2                   ZFUNCNSV3,ZFUNCNSV8)
*
      IMPLICIT NONE
*
      include "../commons/alpha.h"
      include "../commons/massthrs.h"
      include "../commons/ipt.h"
      include "../commons/evol.h"
**
*     Input Variables 
*     
      INTEGER NFI
      DOUBLE COMPLEX ZN
      DOUBLE COMPLEX Q2I,Q2F
**
*     Internal Variables
* 
      INTEGER I,J,K,L,H,G,S,T
      DOUBLE COMPLEX ATH(3)
      DOUBLE COMPLEX MCSG(2,2),MCNS,VCNS(2)
      DOUBLE COMPLEX MBSG(2,2),MBNS,VBNS(2)
      DOUBLE COMPLEX MTSG(2,2),MTNS,VTNS(2)
      DOUBLE COMPLEX EFNNS1(3),EFNSG1(2,2),EFNNS2(3),EFNSG2(2,2)
      DOUBLE COMPLEX EFNNS3(3),EFNSG3(2,2),EFNNS4(3),EFNSG4(2,2)
      DOUBLE COMPLEX A2SG(2,2),A1SG,A2NS
**
*     Output Variables
* 
      DOUBLE COMPLEX EFNNS(3),EFNSG(2,2)
      DOUBLE COMPLEX ZFUNCNS3(2),ZFUNCNS8(2)
      DOUBLE COMPLEX ZFUNCNSV3,ZFUNCNSV8
*
      ATH(1) = AE
      ATH(2) = AM
      ATH(3) = AT
*
      CALL EVOLFACTN(ZN,Q2I,Q2TH(1),NFI,EFNNS1,EFNSG1)
      IF(IPT.EQ.1.AND.EVOL.EQ."TIME")THEN
         CALL MATCHCOEF1SGT(ZN,A1SG)
*
         MCSG(1,1) = (1D0,0D0)
         MCSG(1,2) = ATH(NFI+1) * A1SG
         MCSG(2,1) = (0D0,0D0)
         MCSG(2,2) = (1D0,0D0)
*
         MCNS = (1D0,0D0)
*
         VCNS(1) = (1D0,0D0)
         VCNS(2) = (0D0,0D0)
      ELSE
         MCSG(1,1) = (1D0,0D0)
         MCSG(1,2) = (0D0,0D0)
         MCSG(2,1) = (0D0,0D0)
         MCSG(2,2) = (1D0,0D0)
*
         MCNS = (1D0,0D0)
*
         VCNS(1) = (1D0,0D0)
         VCNS(2) = (0D0,0D0)
      ENDIF
      CALL EVOLFACTN(ZN,Q2TH(1),Q2TH(2),NFI+1,EFNNS2,EFNSG2)
      IF(IPT.EQ.1.AND.EVOL.EQ."TIME")THEN
         CALL MATCHCOEF1SGT(ZN,A1SG)
*
         MBSG(1,1) = (1D0,0D0)
         MBSG(1,2) = ATH(NFI+2) * A1SG
         MBSG(2,1) = (0D0,0D0)
         MBSG(2,2) = (1D0,0D0)
*
         MBNS = (1D0,0D0)
*
         VBNS(1) = (1D0,0D0)
         VBNS(2) = - ATH(NFI+2) * A1SG
      ELSE
         MBSG(1,1) = (1D0,0D0)
         MBSG(1,2) = (0D0,0D0)
         MBSG(2,1) = (0D0,0D0)
         MBSG(2,2) = (1D0,0D0)
*
         MBNS = (1D0,0D0)
*
         VBNS(1) = (1D0,0D0)
         VBNS(2) = (0D0,0D0)
      ENDIF
      CALL EVOLFACTN(ZN,Q2TH(2),Q2TH(3),NFI+2,EFNNS3,EFNSG3)
      IF(IPT.EQ.1.AND.EVOL.EQ."TIME")THEN
         CALL MATCHCOEF1SGT(ZN,A1SG)
*
         MTSG(1,1) = (1D0,0D0)
         MTSG(1,2) = ATH(NFI+3) * A1SG
         MTSG(2,1) = (0D0,0D0)
         MTSG(2,2) = (1D0,0D0)
*
         MTNS = (1D0,0D0)
*
         VTNS(1) = (1D0,0D0)
         VTNS(2) = - 2D0 * ATH(NFI+3) * A1SG
      ELSE
         MTSG(1,1) = (1D0,0D0)
         MTSG(1,2) = (0D0,0D0)
         MTSG(2,1) = (0D0,0D0)
         MTSG(2,2) = (1D0,0D0)
*
         MTNS = (1D0,0D0)
*
         VTNS(1) = (1D0,0D0)
         VTNS(2) = (0D0,0D0)
      ENDIF
      CALL EVOLFACTN(ZN,Q2TH(3),Q2F,NFI+3,EFNNS4,EFNSG4)
*
      DO I=1,2
         DO J=1,2
            EFNSG(I,J) = (0D0,0D0)
         ENDDO
         ZFUNCNS3(I) = (0D0,0D0)
         ZFUNCNS8(I) = (0D0,0D0)
      ENDDO
*     Singlet and gluon
      DO I=1,2
         DO J=1,2
            DO K=1,2
               DO L=1,2
                  DO H=1,2
                     DO G=1,2
                        DO S=1,2
                           DO T=1,2
                              EFNSG(I,J) = EFNSG(I,J) +
     1                                     EFNSG4(I,K) * MTSG(K,L) *
     2                                     EFNSG3(L,H) * MBSG(H,G) *
     3                                     EFNSG2(G,S) * MCSG(S,T) *
     4                                     EFNSG1(T,J)
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
*     Plus, Minus, Valence
      DO I=1,3
         EFNNS(I) = EFNNS4(I) * MTNS * 
     1              EFNNS3(I) * MBNS * 
     2              EFNNS2(I) * MCNS * 
     3              EFNNS1(I)
      ENDDO
*     T_3
      DO I=1,2
         DO J=1,2
            DO K=1,2
               DO L=1,2
                  ZFUNCNS3(I) = ZFUNCNS3(I) + EFNNS4(1) * MTNS *
     1                                        EFNNS3(1) * VBNS(J) * 
     2                                        EFNSG2(J,K) * MCSG(K,L)*
     3                                        EFNSG1(L,I)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
*     T_8
      DO I=1,2
         DO J=1,2
            DO K=1,2
               DO L=1,2
                  DO H=1,2
                     DO S=1,2
                        ZFUNCNS8(I) = ZFUNCNS8(I) + 
     1                                EFNNS4(1) * VTNS(J) *
     2                                EFNSG3(J,K) * MBSG(K,L) * 
     3                                EFNSG2(L,H) * MCSG(H,S) *
     4                                EFNSG1(S,I)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
*     V_{3, 8}
      ZFUNCNSV3 = EFNNS4(2) * MTNS * 
     1            EFNNS3(2) * MBNS * 
     2            EFNNS2(3) * MCNS * 
     3            EFNNS1(3)
      ZFUNCNSV8 = EFNNS4(2) * MTNS * 
     1            EFNNS3(3) * MBNS * 
     2            EFNNS2(3) * MCNS * 
     3            EFNNS1(3)
*
      RETURN
      END
