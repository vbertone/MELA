************************************************************************
*
*     No crossing:
*
*     - CROSS01 --> Q2I,Q2F < M_C^2          (or FFNS with NF = 3)
*     - CROSS02 --> M_C^2 < Q2I,Q2F < M_B^2  (or FFNS with NF = 4)
*     - CROSS03 --> M_B^2 < Q2I,Q2F < M_T^2  (or FFNS with NF = 5)
*     - CROSS04 --> Q2I,Q2F > M_T^2          (or FFNS with NF = 6)
*
*     One crossing:
*
*     - CROSS11 --> Q2I < M_C^2 < Q2F
*     - CROSS12 --> Q2I < M_B^2 < Q2F
*     - CROSS13 --> Q2I < M_T^2 < Q2F
*
*     Two crossings:
*
*     - CROSS21 --> Q2I < M_C^2 < M_B^2 < Q2F
*     - CROSS22 --> Q2I < M_B^2 < M_T^2 < Q2F
*
*     Three crossings:
*
*     - CROSS31 --> Q2I < M_C^2 < M_B^2 < M_T^2 < Q2F
*
************************************************************************
*
*     NO CROSSING
*
************************************************************************
      SUBROUTINE CROSS01(ZN,Q2I,Q2F,NFI,EFNNS,EFNSG,
     1                   ZFUNCNS15,ZFUNCNS24,ZFUNCNS35,
     2                   ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35)
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
      DOUBLE COMPLEX ZFUNCNS15(2),ZFUNCNS24(2),ZFUNCNS35(2)
      DOUBLE COMPLEX ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35
*
      CALL EVOLFACTN(ZN,Q2I,Q2F,NFI,EFNNS,EFNSG)
*
      ZFUNCNS15(1) = EFNSG(1,1)
      ZFUNCNS15(2) = EFNSG(1,2)
      ZFUNCNS24(1) = EFNSG(1,1)
      ZFUNCNS24(2) = EFNSG(1,2)
      ZFUNCNS35(1) = EFNSG(1,1)
      ZFUNCNS35(2) = EFNSG(1,2)
      ZFUNCNSV15   = EFNNS(3)
      ZFUNCNSV24   = EFNNS(3)
      ZFUNCNSV35   = EFNNS(3)      
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE CROSS02(ZN,Q2I,Q2F,NFI,EFNNS,EFNSG,
     1                   ZFUNCNS15,ZFUNCNS24,ZFUNCNS35,
     2                   ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35)
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
      DOUBLE COMPLEX ZFUNCNS15(2),ZFUNCNS24(2),ZFUNCNS35(2)
      DOUBLE COMPLEX ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35
*
      CALL EVOLFACTN(ZN,Q2I,Q2F,NFI,EFNNS,EFNSG)
*
      ZFUNCNS15(1) = EFNNS(1)
      ZFUNCNS15(2) = 0D0
      ZFUNCNS24(1) = EFNSG(1,1)
      ZFUNCNS24(2) = EFNSG(1,2)
      ZFUNCNS35(1) = EFNSG(1,1)
      ZFUNCNS35(2) = EFNSG(1,2)
      ZFUNCNSV15   = EFNNS(2)
      ZFUNCNSV24   = EFNNS(3)
      ZFUNCNSV35   = EFNNS(3)      
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE CROSS03(ZN,Q2I,Q2F,NFI,EFNNS,EFNSG,
     1                   ZFUNCNS15,ZFUNCNS24,ZFUNCNS35,
     2                   ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35)
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
      DOUBLE COMPLEX ZFUNCNS15(2),ZFUNCNS24(2),ZFUNCNS35(2)
      DOUBLE COMPLEX ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35
*
      CALL EVOLFACTN(ZN,Q2I,Q2F,NFI,EFNNS,EFNSG)
*
      ZFUNCNS15(1) = EFNNS(1)
      ZFUNCNS15(2) = 0D0
      ZFUNCNS24(1) = EFNNS(1)
      ZFUNCNS24(2) = 0D0
      ZFUNCNS35(1) = EFNSG(1,1)
      ZFUNCNS35(2) = EFNSG(1,2)
      ZFUNCNSV15   = EFNNS(2)
      ZFUNCNSV24   = EFNNS(2)
      ZFUNCNSV35   = EFNNS(3)      
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE CROSS04(ZN,Q2I,Q2F,NFI,EFNNS,EFNSG,
     1                   ZFUNCNS15,ZFUNCNS24,ZFUNCNS35,
     2                   ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35)
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
      DOUBLE COMPLEX ZFUNCNS15(2),ZFUNCNS24(2),ZFUNCNS35(2)
      DOUBLE COMPLEX ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35
*
      CALL EVOLFACTN(ZN,Q2I,Q2F,NFI,EFNNS,EFNSG)
*
      ZFUNCNS15(1) = EFNNS(1)
      ZFUNCNS15(2) = 0D0
      ZFUNCNS24(1) = EFNNS(1)
      ZFUNCNS24(2) = 0D0
      ZFUNCNS35(1) = EFNNS(1)
      ZFUNCNS35(2) = 0D0
      ZFUNCNSV15   = EFNNS(2)
      ZFUNCNSV24   = EFNNS(2)
      ZFUNCNSV35   = EFNNS(2)      
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
     1                   ZFUNCNS15,ZFUNCNS24,ZFUNCNS35,
     2                   ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35)
*
      IMPLICIT NONE
*
      include "../commons/alphas.h"
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
      DOUBLE COMPLEX ASTH(4:6)
      DOUBLE COMPLEX MCSG(2,2),MCNS,VCNS(2)
      DOUBLE COMPLEX EFNNS1(3),EFNSG1(2,2),EFNNS2(3),EFNSG2(2,2)
      DOUBLE COMPLEX A2SG(2,2),A1SG,A2NS
**
*     Output Variables
* 
      DOUBLE COMPLEX EFNNS(3),EFNSG(2,2)
      DOUBLE COMPLEX ZFUNCNS15(2),ZFUNCNS24(2),ZFUNCNS35(2)
      DOUBLE COMPLEX ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35
*
      ASTH(4) = ASC
      ASTH(5) = ASB
      ASTH(6) = AST
*
      CALL EVOLFACTN(ZN,Q2I,Q2TH(4),NFI,EFNNS1,EFNSG1)
      IF(IPT.EQ.2)THEN
         CALL MATCHCOEF2SG(ZN,A2SG)
         CALL MATCHCOEF2NS(ZN,A2NS)
*
         MCSG(1,1) = (1D0,0D0) + ASTH(NFI+1)**2D0 * (A2NS + A2SG(1,1))
         MCSG(1,2) = ASTH(NFI+1)**2D0 * A2SG(1,2)
         MCSG(2,1) = ASTH(NFI+1)**2D0 * A2SG(2,1)
         MCSG(2,2) = (1D0,0D0) + ASTH(NFI+1)**2D0 * A2SG(2,2)
*
         MCNS = (1D0,0D0) + ASTH(NFI+1)**2D0 * A2NS
*
         VCNS(1) = (1D0,0D0)+ASTH(NFI+1)**2D0*(A2NS - 3D0 * A2SG(1,1))
         VCNS(2) = - 3D0 * ASTH(NFI+1)**2D0 * A2SG(1,2)
      ELSEIF(IPT.EQ.1.AND.EVOL.EQ."TIME")THEN
         CALL MATCHCOEF1SGT(ZN,A1SG)
*
         MCSG(1,1) = (1D0,0D0)
         MCSG(1,2) = ASTH(NFI+1) * A1SG
         MCSG(2,1) = (0D0,0D0)
         MCSG(2,2) = (1D0,0D0)
*
         MCNS = (1D0,0D0)
*
         VCNS(1) = (1D0,0D0)
         VCNS(2) = - 3D0 * ASTH(NFI+1) * A1SG
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
      CALL EVOLFACTN(ZN,Q2TH(4),Q2F,NFI+1,EFNNS2,EFNSG2)
*
      DO I=1,2
         DO J=1,2
            EFNSG(I,J) = (0D0,0D0)
         ENDDO
         ZFUNCNS15(I) = (0D0,0D0)
         ZFUNCNS24(I) = (0D0,0D0)
         ZFUNCNS35(I) = (0D0,0D0)
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
*     T_{3,8}(I=1),V_{3,8}(I=2),V(I=3)
      DO I=1,3
         EFNNS(I) = EFNNS2(I) * MCNS * EFNNS1(I)
      ENDDO
*     T_15
      DO I=1,2
         DO J=1,2
            ZFUNCNS15(I) = ZFUNCNS15(I) + EFNNS2(1) * VCNS(J) * 
     1                                    EFNSG1(J,I)
         ENDDO
      ENDDO
*     T_{24,35}
      DO I=1,2
         DO J=1,2
            DO K=1,2
               ZFUNCNS24(I) = ZFUNCNS24(I) + EFNSG2(1,J) * MCSG(J,K) * 
     1                                       EFNSG1(K,I)
               ZFUNCNS35(I) = ZFUNCNS35(I) + EFNSG2(1,J) * MCSG(J,K) * 
     1                                       EFNSG1(K,I)
            ENDDO
         ENDDO
      ENDDO
*     V_{15,24,35}
      ZFUNCNSV15 = EFNNS2(2) * MCNS * EFNNS1(3)
      ZFUNCNSV24 = EFNNS2(3) * MCNS * EFNNS1(3)
      ZFUNCNSV35 = EFNNS2(3) * MCNS * EFNNS1(3)
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE CROSS12(ZN,Q2I,Q2F,NFI,EFNNS,EFNSG,
     1                   ZFUNCNS15,ZFUNCNS24,ZFUNCNS35,
     2                   ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35)
*
      IMPLICIT NONE
*
      include "../commons/alphas.h"
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
      DOUBLE COMPLEX ASTH(4:6)
      DOUBLE COMPLEX MBSG(2,2),MBNS,VBNS(2)
      DOUBLE COMPLEX EFNNS1(3),EFNSG1(2,2),EFNNS2(3),EFNSG2(2,2)
      DOUBLE COMPLEX A2SG(2,2),A1SG,A2NS
**
*     Output Variables
* 
      DOUBLE COMPLEX EFNNS(3),EFNSG(2,2)
      DOUBLE COMPLEX ZFUNCNS15(2),ZFUNCNS24(2),ZFUNCNS35(2)
      DOUBLE COMPLEX ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35
*
      ASTH(4) = ASC
      ASTH(5) = ASB
      ASTH(6) = AST
*
      CALL EVOLFACTN(ZN,Q2I,Q2TH(5),NFI,EFNNS1,EFNSG1)
      IF(IPT.EQ.2)THEN
         CALL MATCHCOEF2SG(ZN,A2SG)
         CALL MATCHCOEF2NS(ZN,A2NS)
*
         MBSG(1,1) = (1D0,0D0) + ASTH(NFI+1)**2D0 * (A2NS + A2SG(1,1))
         MBSG(1,2) = ASTH(NFI+1)**2D0 * A2SG(1,2)
         MBSG(2,1) = ASTH(NFI+1)**2D0 * A2SG(2,1)
         MBSG(2,2) = (1D0,0D0) + ASTH(NFI+1)**2D0 * A2SG(2,2)
*
         MBNS = (1D0,0D0) + ASTH(NFI+1)**2D0 * A2NS
*
         VBNS(1) = (1D0,0D0)+ASTH(NFI+1)**2D0*(A2NS - 4D0 * A2SG(1,1))
         VBNS(2) = - 4D0 * ASTH(NFI+1)**2D0 * A2SG(1,2)
      ELSEIF(IPT.EQ.1.AND.EVOL.EQ."TIME")THEN
         CALL MATCHCOEF1SGT(ZN,A1SG)
*
         MBSG(1,1) = (1D0,0D0)
         MBSG(1,2) = ASTH(NFI+1) * A1SG
         MBSG(2,1) = (0D0,0D0)
         MBSG(2,2) = (1D0,0D0)
*
         MBNS = (1D0,0D0)
*
         VBNS(1) = (1D0,0D0)
         VBNS(2) = - 4D0 * ASTH(NFI+1) * A1SG
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
      CALL EVOLFACTN(ZN,Q2TH(5),Q2F,NFI+1,EFNNS2,EFNSG2)
*
      DO I=1,2
         DO J=1,2
            EFNSG(I,J) = (0D0,0D0)
         ENDDO
         ZFUNCNS24(I) = (0D0,0D0)
         ZFUNCNS35(I) = (0D0,0D0)
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
*     T_{3,8}(I=1),V_{3,8}(I=2),V(I=3)
      DO I=1,3
         EFNNS(I) = EFNNS2(I) * MBNS * EFNNS1(I)
      ENDDO
*     T_15
      ZFUNCNS15(1) = EFNNS2(1) * MBNS * EFNNS1(1)
      ZFUNCNS15(2) = (0D0,0D0)
*     T_24
      DO I=1,2
         DO J=1,2
            ZFUNCNS24(I) = ZFUNCNS24(I) + EFNNS2(1) * VBNS(J) * 
     1                                    EFNSG1(J,I)
         ENDDO
      ENDDO
*     T_35
      DO I=1,2
         DO J=1,2
            DO K=1,2
               ZFUNCNS35(I) = ZFUNCNS35(I) + EFNSG2(1,J) * MBSG(J,K) *
     1                                       EFNSG1(K,I)
            ENDDO 
         ENDDO
      ENDDO
*     V_{15,24,35}
      ZFUNCNSV15 = EFNNS2(2) * MBNS * EFNNS1(2)
      ZFUNCNSV24 = EFNNS2(2) * MBNS * EFNNS1(3)
      ZFUNCNSV35 = EFNNS2(3) * MBNS * EFNNS1(3)
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE CROSS13(ZN,Q2I,Q2F,NFI,EFNNS,EFNSG,
     1                   ZFUNCNS15,ZFUNCNS24,ZFUNCNS35,
     2                   ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35)
*
      IMPLICIT NONE
*
      include "../commons/alphas.h"
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
      DOUBLE COMPLEX ASTH(4:6)
      DOUBLE COMPLEX MTSG(2,2),MTNS,VTNS(2)
      DOUBLE COMPLEX EFNNS1(3),EFNSG1(2,2),EFNNS2(3),EFNSG2(2,2)
      DOUBLE COMPLEX A2SG(2,2),A1SG,A2NS
**
*     Output Variables
* 
      DOUBLE COMPLEX EFNNS(3),EFNSG(2,2)
      DOUBLE COMPLEX ZFUNCNS15(2),ZFUNCNS24(2),ZFUNCNS35(2)
      DOUBLE COMPLEX ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35
*
      ASTH(4) = ASC
      ASTH(5) = ASB
      ASTH(6) = AST
*
      CALL EVOLFACTN(ZN,Q2I,Q2TH(6),NFI,EFNNS1,EFNSG1)
      IF(IPT.EQ.2)THEN
         CALL MATCHCOEF2SG(ZN,A2SG)
         CALL MATCHCOEF2NS(ZN,A2NS)
*
         MTSG(1,1) = (1D0,0D0) + ASTH(NFI+1)**2D0 * (A2NS + A2SG(1,1))
         MTSG(1,2) = ASTH(NFI+1)**2D0 * A2SG(1,2)
         MTSG(2,1) = ASTH(NFI+1)**2D0 * A2SG(2,1)
         MTSG(2,2) = (1D0,0D0) + ASTH(NFI+1)**2D0 * A2SG(2,2)
*
         MTNS = (1D0,0D0) + ASTH(NFI+1)**2D0 * A2NS
*
         VTNS(1) = (1D0,0D0)+ASTH(NFI+1)**2D0*(A2NS - 5D0 * A2SG(1,1))
         VTNS(2) = - 5D0 * ASTH(NFI+1)**2D0 * A2SG(1,2)
      ELSEIF(IPT.EQ.1.AND.EVOL.EQ."TIME")THEN
         CALL MATCHCOEF1SGT(ZN,A1SG)
*
         MTSG(1,1) = (1D0,0D0)
         MTSG(1,2) = ASTH(NFI+1) * A1SG
         MTSG(2,1) = (0D0,0D0)
         MTSG(2,2) = (1D0,0D0)
*
         MTNS = (1D0,0D0)
*
         VTNS(1) = (1D0,0D0)
         VTNS(2) = - 5D0 * ASTH(NFI+1) * A1SG
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
      CALL EVOLFACTN(ZN,Q2TH(6),Q2F,NFI+1,EFNNS2,EFNSG2)
*
      DO I=1,2
         DO J=1,2
            EFNSG(I,J) = (0D0,0D0)
         ENDDO
         ZFUNCNS35(I) = (0D0,0D0)
      ENDDO
*     Singlet and gluon
      DO I=1,2
         DO J=1,2
            DO K=1,2
               DO L=1,2
                  EFNSG(I,J) = EFNSG(I,J) + EFNSG2(I,K) * MTSG(K,L) * 
     1                                      EFNSG1(L,J)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
*     T_{3,8}(I=1),V_{3,8}(I=2),V(I=3)
      DO I=1,3
         EFNNS(I) = EFNNS2(I) * MTNS * EFNNS1(I)
      ENDDO
*     T_15
      ZFUNCNS15(1) = EFNNS2(1) * MTNS * EFNNS1(1)
      ZFUNCNS15(2) = (0D0,0D0)
*     T_24
      ZFUNCNS24(1) = EFNNS2(1) * MTNS * EFNNS1(1)
      ZFUNCNS24(2) = (0D0,0D0)
*     T_35
      DO I=1,2
         DO J=1,2
            ZFUNCNS35(I) = ZFUNCNS35(I) + EFNNS2(1) * VTNS(J) * 
     1                                    EFNSG1(J,I)
         ENDDO
      ENDDO
*     V_{15,24,35}
      ZFUNCNSV15 = EFNNS2(2) * MTNS * EFNNS1(2)
      ZFUNCNSV24 = EFNNS2(2) * MTNS * EFNNS1(2)
      ZFUNCNSV35 = EFNNS2(2) * MTNS * EFNNS1(3)
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
     1                   ZFUNCNS15,ZFUNCNS24,ZFUNCNS35,
     2                   ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35)
*
      IMPLICIT NONE
*
      include "../commons/alphas.h"
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
      DOUBLE COMPLEX ASTH(4:6)
      DOUBLE COMPLEX MCSG(2,2),MCNS,VCNS(2)
      DOUBLE COMPLEX MBSG(2,2),MBNS,VBNS(2)
      DOUBLE COMPLEX EFNNS1(3),EFNSG1(2,2),EFNNS2(3),EFNSG2(2,2)
      DOUBLE COMPLEX EFNNS3(3),EFNSG3(2,2)
      DOUBLE COMPLEX A2SG(2,2),A1SG,A2NS
**
*     Output Variables
* 
      DOUBLE COMPLEX EFNNS(3),EFNSG(2,2)
      DOUBLE COMPLEX ZFUNCNS15(2),ZFUNCNS24(2),ZFUNCNS35(2)
      DOUBLE COMPLEX ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35
*
      ASTH(4) = ASC
      ASTH(5) = ASB
      ASTH(6) = AST
*
      CALL EVOLFACTN(ZN,Q2I,Q2TH(4),NFI,EFNNS1,EFNSG1)
      IF(IPT.EQ.2)THEN
         CALL MATCHCOEF2SG(ZN,A2SG)
         CALL MATCHCOEF2NS(ZN,A2NS)
*
         MCSG(1,1) = (1D0,0D0) + ASTH(NFI+1)**2D0 * (A2NS + A2SG(1,1))
         MCSG(1,2) = ASTH(NFI+1)**2D0 * A2SG(1,2)
         MCSG(2,1) = ASTH(NFI+1)**2D0 * A2SG(2,1)
         MCSG(2,2) = (1D0,0D0) + ASTH(NFI+1)**2D0 * A2SG(2,2)
*
         MCNS = (1D0,0D0) + ASTH(NFI+1)**2D0 * A2NS
*
         VCNS(1) = (1D0,0D0)+ASTH(NFI+1)**2D0*(A2NS - 3D0 * A2SG(1,1))
         VCNS(2) = - 3D0 * ASTH(NFI+1)**2D0 * A2SG(1,2)
      ELSEIF(IPT.EQ.1.AND.EVOL.EQ."TIME")THEN
         CALL MATCHCOEF1SGT(ZN,A1SG)
*
         MCSG(1,1) = (1D0,0D0)
         MCSG(1,2) = ASTH(NFI+1) * A1SG
         MCSG(2,1) = (0D0,0D0)
         MCSG(2,2) = (1D0,0D0)
*
         MCNS = (1D0,0D0)
*
         VCNS(1) = (1D0,0D0)
         VCNS(2) = - 3D0 * ASTH(NFI+1) * A1SG
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
      CALL EVOLFACTN(ZN,Q2TH(4),Q2TH(5),NFI+1,EFNNS2,EFNSG2)
      IF(IPT.EQ.2)THEN
         MBSG(1,1) = (1D0,0D0) + ASTH(NFI+2)**2D0 * ( A2NS + A2SG(1,1))
         MBSG(1,2) = ASTH(NFI+2)**2D0 * A2SG(1,2)
         MBSG(2,1) = ASTH(NFI+2)**2D0 * A2SG(2,1)
         MBSG(2,2) = (1D0,0D0) + ASTH(NFI+2)**2D0 * A2SG(2,2)
*
         MBNS = (1D0,0D0) + ASTH(NFI+2)**2D0 * A2NS
*
         VBNS(1) = (1D0,0D0)+ASTH(NFI+2)**2D0*(A2NS - 4D0 * A2SG(1,1))
         VBNS(2) = - 4D0 * ASTH(NFI+2)**2D0 * A2SG(1,2)
      ELSEIF(IPT.EQ.1.AND.EVOL.EQ."TIME")THEN
         CALL MATCHCOEF1SGT(ZN,A1SG)
*
         MBSG(1,1) = (1D0,0D0)
         MBSG(1,2) = ASTH(NFI+2) * A1SG
         MBSG(2,1) = (0D0,0D0)
         MBSG(2,2) = (1D0,0D0)
*
         MBNS = (1D0,0D0)
*
         VBNS(1) = (1D0,0D0)
         VBNS(2) = - 4D0 * ASTH(NFI+2) * A1SG
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
      CALL EVOLFACTN(ZN,Q2TH(5),Q2F,NFI+2,EFNNS3,EFNSG3)
*
      DO I=1,2
         DO J=1,2
            EFNSG(I,J) = (0D0,0D0)
         ENDDO
         ZFUNCNS15(I) = (0D0,0D0)
         ZFUNCNS24(I) = (0D0,0D0)
         ZFUNCNS35(I) = (0D0,0D0)
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
*     T_{3,8}(I=1),V_{3,8}(I=2),V(I=3)
      DO I=1,3
         EFNNS(I) = EFNNS3(I) * MBNS * EFNNS2(I) * MCNS * EFNNS1(I)
      ENDDO
*     T_15
      DO I=1,2
         DO J=1,2
            ZFUNCNS15(I) = ZFUNCNS15(I) + EFNNS3(1) * MBNS * 
     1                                    EFNNS2(1) * VCNS(J) * 
     2                                    EFNSG1(J,I)
         ENDDO
      ENDDO
*     T_24
      DO I=1,2
         DO J=1,2
            DO K=1,2
               DO L=1,2
                  ZFUNCNS24(I) = ZFUNCNS24(I) + EFNNS3(1) * VBNS(J) * 
     1                                          EFNSG2(J,K) * MCSG(K,L)*
     2                                          EFNSG1(L,I)
               ENDDO
            ENDDO 
         ENDDO
      ENDDO
*     T_35
      DO I=1,2
         DO J=1,2
            DO K=1,2
               DO L=1,2
                  DO H=1,2
                     ZFUNCNS35(I) = ZFUNCNS35(I) + 
     1                              EFNSG3(1,J) * MBSG(J,K) * 
     2                              EFNSG2(K,L) * MCSG(L,H) *
     3                              EFNSG1(H,I)
                  ENDDO
               ENDDO
            ENDDO 
         ENDDO
      ENDDO
*     V_{15,24,35}
      ZFUNCNSV15 = EFNNS3(2) * MBNS * EFNNS2(2) * MCNS * EFNNS1(3)
      ZFUNCNSV24 = EFNNS3(2) * MBNS * EFNNS2(3) * MCNS * EFNNS1(3)
      ZFUNCNSV35 = EFNNS3(3) * MBNS * EFNNS2(3) * MCNS * EFNNS1(3)
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE CROSS22(ZN,Q2I,Q2F,NFI,EFNNS,EFNSG,
     1                   ZFUNCNS15,ZFUNCNS24,ZFUNCNS35,
     2                   ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35)
*
      IMPLICIT NONE
*
      include "../commons/alphas.h"
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
      DOUBLE COMPLEX ASTH(4:6)
      DOUBLE COMPLEX MBSG(2,2),MBNS,VBNS(2)
      DOUBLE COMPLEX MTSG(2,2),MTNS,VTNS(2)
      DOUBLE COMPLEX EFNNS1(3),EFNSG1(2,2),EFNNS2(3),EFNSG2(2,2)
      DOUBLE COMPLEX EFNNS3(3),EFNSG3(2,2)
      DOUBLE COMPLEX A2SG(2,2),A1SG,A2NS
**
*     Output Variables
* 
      DOUBLE COMPLEX EFNNS(3),EFNSG(2,2)
      DOUBLE COMPLEX ZFUNCNS15(2),ZFUNCNS24(2),ZFUNCNS35(2)
      DOUBLE COMPLEX ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35
*
      ASTH(4) = ASC
      ASTH(5) = ASB
      ASTH(6) = AST
*
      CALL EVOLFACTN(ZN,Q2I,Q2TH(5),NFI,EFNNS1,EFNSG1)
      IF(IPT.EQ.2)THEN
         CALL MATCHCOEF2SG(ZN,A2SG)
         CALL MATCHCOEF2NS(ZN,A2NS)
*
         MBSG(1,1) = (1D0,0D0) + ASTH(NFI+1)**2D0 * (A2NS + A2SG(1,1))
         MBSG(1,2) = ASTH(NFI+1)**2D0 * A2SG(1,2)
         MBSG(2,1) = ASTH(NFI+1)**2D0 * A2SG(2,1)
         MBSG(2,2) = (1D0,0D0) + ASTH(NFI+1)**2D0 * A2SG(2,2)
*
         MBNS = (1D0,0D0) + ASTH(NFI+1)**2D0 * A2NS
*
         VBNS(1) = (1D0,0D0)+ASTH(NFI+1)**2D0*(A2NS - 4D0 * A2SG(1,1))
         VBNS(2) = - 4D0 * ASTH(NFI+1)**2D0 * A2SG(1,2)
      ELSEIF(IPT.EQ.1.AND.EVOL.EQ."TIME")THEN
         CALL MATCHCOEF1SGT(ZN,A1SG)
*
         MBSG(1,1) = (1D0,0D0)
         MBSG(1,2) = ASTH(NFI+1) * A1SG
         MBSG(2,1) = (0D0,0D0)
         MBSG(2,2) = (1D0,0D0)
*
         MBNS = (1D0,0D0)
*
         VBNS(1) = (1D0,0D0)
         VBNS(2) = - 4D0 * ASTH(NFI+1) * A1SG
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
      CALL EVOLFACTN(ZN,Q2TH(5),Q2TH(6),NFI+1,EFNNS2,EFNSG2)
      IF(IPT.EQ.2)THEN
         CALL MATCHCOEF2SG(ZN,A2SG)
         CALL MATCHCOEF2NS(ZN,A2NS)
*
         MTSG(1,1) = (1D0,0D0) + ASTH(NFI+2)**2D0 * ( A2NS + A2SG(1,1))
         MTSG(1,2) = ASTH(NFI+2)**2D0 * A2SG(1,2)
         MTSG(2,1) = ASTH(NFI+2)**2D0 * A2SG(2,1)
         MTSG(2,2) = (1D0,0D0) + ASTH(NFI+2)**2D0 * A2SG(2,2)
*
         MTNS = (1D0,0D0) + ASTH(NFI+2)**2D0 * A2NS
*
         VTNS(1) = (1D0,0D0)+ASTH(NFI+2)**2D0*(A2NS - 5D0 * A2SG(1,1))
         VTNS(2) = - 5D0 * ASTH(NFI+2)**2D0 * A2SG(1,2)
      ELSEIF(IPT.EQ.1.AND.EVOL.EQ."TIME")THEN
         CALL MATCHCOEF1SGT(ZN,A1SG)
*
         MTSG(1,1) = (1D0,0D0)
         MTSG(1,2) = ASTH(NFI+2) * A1SG
         MTSG(2,1) = (0D0,0D0)
         MTSG(2,2) = (1D0,0D0)
*
         MTNS = (1D0,0D0)
*
         VTNS(1) = (1D0,0D0)
         VTNS(2) = - 5D0 * ASTH(NFI+2) * A1SG
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
      CALL EVOLFACTN(ZN,Q2TH(6),Q2F,NFI+2,EFNNS3,EFNSG3)
*
      DO I=1,2
         DO J=1,2
            EFNSG(I,J) = (0D0,0D0)
         ENDDO
         ZFUNCNS24(I) = (0D0,0D0)
         ZFUNCNS35(I) = (0D0,0D0)
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
*     T_{3,8}(I=1),V_{3,8}(I=2),V(I=3)
      DO I=1,3
         EFNNS(I) = EFNNS3(I) * MTNS * EFNNS2(I) * MBNS * EFNNS1(I)
      ENDDO
*     T_15
      ZFUNCNS15(1) = EFNNS3(1) * MTNS * EFNNS2(1) * MBNS * EFNNS1(1)
      ZFUNCNS15(2) = (0D0,0D0)
*     T_24
      DO I=1,2
         DO J=1,2
            ZFUNCNS24(I) = ZFUNCNS24(I) + EFNNS3(1) * MTNS * 
     1                                    EFNNS2(1) * VBNS(J) * 
     2                                    EFNSG1(J,I)
         ENDDO
      ENDDO
*     T_35
      DO I=1,2
         DO J=1,2
            DO K=1,2
               DO L=1,2
                  ZFUNCNS35(I) = ZFUNCNS35(I) + EFNNS3(1) * VTNS(J) * 
     1                                          EFNSG2(J,K) * MBSG(K,L)*
     2                                          EFNSG1(L,I)
               ENDDO
            ENDDO 
         ENDDO
      ENDDO
*     V_{15,24,35}
      ZFUNCNSV15 = EFNNS3(2) * MTNS * EFNNS2(2) * MBNS * EFNNS1(2)
      ZFUNCNSV24 = EFNNS3(2) * MTNS * EFNNS2(2) * MBNS * EFNNS1(3)
      ZFUNCNSV35 = EFNNS3(2) * MTNS * EFNNS2(3) * MBNS * EFNNS1(3)
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
     1                   ZFUNCNS15,ZFUNCNS24,ZFUNCNS35,
     2                   ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35)
*
      IMPLICIT NONE
*
      include "../commons/alphas.h"
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
      DOUBLE COMPLEX ASTH(4:6)
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
      DOUBLE COMPLEX ZFUNCNS15(2),ZFUNCNS24(2),ZFUNCNS35(2)
      DOUBLE COMPLEX ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35
*
      ASTH(4) = ASC
      ASTH(5) = ASB
      ASTH(6) = AST
*
      CALL EVOLFACTN(ZN,Q2I,Q2TH(4),NFI,EFNNS1,EFNSG1)
      IF(IPT.EQ.2)THEN
         CALL MATCHCOEF2SG(ZN,A2SG)
         CALL MATCHCOEF2NS(ZN,A2NS)
*
         MCSG(1,1) = (1D0,0D0) + ASTH(NFI+1)**2D0 * (A2NS + A2SG(1,1))
         MCSG(1,2) = ASTH(NFI+1)**2D0 * A2SG(1,2)
         MCSG(2,1) = ASTH(NFI+1)**2D0 * A2SG(2,1)
         MCSG(2,2) = (1D0,0D0) + ASTH(NFI+1)**2D0 * A2SG(2,2)
*
         MCNS = (1D0,0D0) + ASTH(NFI+1)**2D0 * A2NS
*
         VCNS(1) = (1D0,0D0)+ASTH(NFI+1)**2D0*(A2NS - 3D0 * A2SG(1,1))
         VCNS(2) = - 3D0 * ASTH(NFI+1)**2D0 * A2SG(1,2)
      ELSEIF(IPT.EQ.1.AND.EVOL.EQ."TIME")THEN
         CALL MATCHCOEF1SGT(ZN,A1SG)
*
         MCSG(1,1) = (1D0,0D0)
         MCSG(1,2) = ASTH(NFI+1) * A1SG
         MCSG(2,1) = (0D0,0D0)
         MCSG(2,2) = (1D0,0D0)
*
         MCNS = (1D0,0D0)
*
         VCNS(1) = (1D0,0D0)
         VCNS(2) = - 3D0 * ASTH(NFI+1) * A1SG
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
      CALL EVOLFACTN(ZN,Q2TH(4),Q2TH(5),NFI+1,EFNNS2,EFNSG2)
      IF(IPT.EQ.2)THEN
         CALL MATCHCOEF2SG(ZN,A2SG)
         CALL MATCHCOEF2NS(ZN,A2NS)
*
         MBSG(1,1) = (1D0,0D0) + ASTH(NFI+2)**2D0 * ( A2NS + A2SG(1,1))
         MBSG(1,2) = ASTH(NFI+2)**2D0 * A2SG(1,2)
         MBSG(2,1) = ASTH(NFI+2)**2D0 * A2SG(2,1)
         MBSG(2,2) = (1D0,0D0) + ASTH(NFI+2)**2D0 * A2SG(2,2)
*
         MBNS = (1D0,0D0) + ASTH(NFI+2)**2D0 * A2NS
*
         VBNS(1) = (1D0,0D0)+ASTH(NFI+2)**2D0*(A2NS - 4D0 * A2SG(1,1))
         VBNS(2) = - 4D0 * ASTH(NFI+2)**2D0 * A2SG(1,2)
      ELSEIF(IPT.EQ.1.AND.EVOL.EQ."TIME")THEN
         CALL MATCHCOEF1SGT(ZN,A1SG)
*
         MBSG(1,1) = (1D0,0D0)
         MBSG(1,2) = ASTH(NFI+2) * A1SG
         MBSG(2,1) = (0D0,0D0)
         MBSG(2,2) = (1D0,0D0)
*
         MBNS = (1D0,0D0)
*
         VBNS(1) = (1D0,0D0)
         VBNS(2) = - 4D0 * ASTH(NFI+2) * A1SG
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
      CALL EVOLFACTN(ZN,Q2TH(5),Q2TH(6),NFI+2,EFNNS3,EFNSG3)
      IF(IPT.EQ.2)THEN
         CALL MATCHCOEF2SG(ZN,A2SG)
         CALL MATCHCOEF2NS(ZN,A2NS)
*
         MTSG(1,1) = (1D0,0D0) + ASTH(NFI+3)**2D0 * ( A2NS + A2SG(1,1))
         MTSG(1,2) = ASTH(NFI+3)**2D0 * A2SG(1,2)
         MTSG(2,1) = ASTH(NFI+3)**2D0 * A2SG(2,1)
         MTSG(2,2) = (1D0,0D0) + ASTH(NFI+3)**2D0 * A2SG(2,2)
*
         MTNS = (1D0,0D0) + ASTH(NFI+3)**2D0 * A2NS
*
         VTNS(1) = (1D0,0D0)+ASTH(NFI+3)**2D0*(A2NS - 5D0 * A2SG(1,1))
         VTNS(2) = - 5D0 * ASTH(NFI+3)**2D0 * A2SG(1,2)
      ELSEIF(IPT.EQ.1.AND.EVOL.EQ."TIME")THEN
         CALL MATCHCOEF1SGT(ZN,A1SG)
*
         MTSG(1,1) = (1D0,0D0)
         MTSG(1,2) = ASTH(NFI+3) * A1SG
         MTSG(2,1) = (0D0,0D0)
         MTSG(2,2) = (1D0,0D0)
*
         MTNS = (1D0,0D0)
*
         VTNS(1) = (1D0,0D0)
         VTNS(2) = - 5D0 * ASTH(NFI+3) * A1SG
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
      CALL EVOLFACTN(ZN,Q2TH(6),Q2F,NFI+3,EFNNS4,EFNSG4)
*
      DO I=1,2
         DO J=1,2
            EFNSG(I,J) = (0D0,0D0)
         ENDDO
         ZFUNCNS15(I) = (0D0,0D0)
         ZFUNCNS24(I) = (0D0,0D0)
         ZFUNCNS35(I) = (0D0,0D0)
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
*     T_{3,8}(I=1),V_{3,8}(I=2),V(I=3)
      DO I=1,3
         EFNNS(I) = EFNNS4(I) * MTNS * 
     1              EFNNS3(I) * MBNS * 
     2              EFNNS2(I) * MCNS * 
     3              EFNNS1(I)
      ENDDO
*     T_15
      DO I=1,2
         DO J=1,2
            ZFUNCNS15(I) = ZFUNCNS15(I) + EFNNS4(1) * MTNS *
     1                                    EFNNS3(1) * MBNS * 
     2                                    EFNNS2(1) * VCNS(J) * 
     3                                    EFNSG1(J,I)
         ENDDO
      ENDDO
*     T_24
      DO I=1,2
         DO J=1,2
            DO K=1,2
               DO L=1,2
                  ZFUNCNS24(I) = ZFUNCNS24(I) + EFNNS4(1) * MTNS *
     1                                          EFNNS3(1) * VBNS(J) * 
     2                                          EFNSG2(J,K) * MCSG(K,L)*
     3                                          EFNSG1(L,I)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
*     T_35
      DO I=1,2
         DO J=1,2
            DO K=1,2
               DO L=1,2
                  DO H=1,2
                     DO S=1,2
                        ZFUNCNS35(I) = ZFUNCNS35(I) + 
     1                                 EFNNS4(1) * VTNS(J) *
     2                                 EFNSG3(J,K) * MBSG(K,L) * 
     3                                 EFNSG2(L,H) * MCSG(H,S) *
     4                                 EFNSG1(S,I)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
*     V_{15,24,35}
      ZFUNCNSV15 = EFNNS4(2) * MTNS * 
     1             EFNNS3(2) * MBNS * 
     2             EFNNS2(2) * MCNS * 
     3             EFNNS1(3)
      ZFUNCNSV24 = EFNNS4(2) * MTNS * 
     1             EFNNS3(2) * MBNS * 
     2             EFNNS2(3) * MCNS * 
     3             EFNNS1(3)
      ZFUNCNSV35 = EFNNS4(2) * MTNS * 
     1             EFNNS3(3) * MBNS * 
     2             EFNNS2(3) * MCNS * 
     3             EFNNS1(3)
*
      RETURN
      END
