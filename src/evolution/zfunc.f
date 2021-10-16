************************************************************************
*
*     zfunc.f:
*     It returns the evolution kernels in N space from Q20 to Q2
*
************************************************************************
      SUBROUTINE ZFUNC(ZN,Q20,Q2,
     1                 ZFUNCNS,ZFUNCNS3,ZFUNCNS8,
     2                 ZFUNCSG,ZFUNCNSV3,ZFUNCNSV8)
*
      IMPLICIT NONE
*
      include "../commons/massthrs.h"
      include "../commons/nffn.h"
      include "../commons/nfmax.h"
      include "../commons/ns.h"
      include "../commons/nf.h"
**
*     Input Variables
*
      DOUBLE COMPLEX Q20,Q2
      DOUBLE COMPLEX ZN
**
*     Output Variables
*
      DOUBLE COMPLEX ZFUNCNS(3),ZFUNCSG(2,2)
      DOUBLE COMPLEX ZFUNCNS3(2),ZFUNCNS8(2)
      DOUBLE COMPLEX ZFUNCNSV3,ZFUNCNSV8
*
*     Evolution Kernels
*
*     Fixed Flavour Number Schemes
*
      IF(NS.EQ."FFNS")THEN
         NFI = NFFN
         NFF = NFFN
         IF(NFFN.EQ.0)THEN
            CALL CROSS01(ZN,Q20,Q2,NFFN,ZFUNCNS,ZFUNCSG,
     1                   ZFUNCNS3,ZFUNCNS8,
     2                   ZFUNCNSV3,ZFUNCNSV8)
         ELSEIF(NFFN.EQ.1)THEN
            CALL CROSS02(ZN,Q20,Q2,NFFN,ZFUNCNS,ZFUNCSG,
     1                   ZFUNCNS3,ZFUNCNS8,
     2                   ZFUNCNSV3,ZFUNCNSV8)
         ELSEIF(NFFN.EQ.2)THEN
            CALL CROSS03(ZN,Q20,Q2,NFFN,ZFUNCNS,ZFUNCSG,
     1                   ZFUNCNS3,ZFUNCNS8,
     2                   ZFUNCNSV3,ZFUNCNSV8)
         ELSEIF(NFFN.EQ.3)THEN
            CALL CROSS04(ZN,Q20,Q2,NFFN,ZFUNCNS,ZFUNCSG,
     1                   ZFUNCNS3,ZFUNCNS8,
     2                   ZFUNCNSV3,ZFUNCNSV8)
         ENDIF
*
*     Variable Flavour Number Schemes
*
      ELSEIF(NS.EQ."VFNS")THEN
*
*     Determine number of active flavours at the initial (NFI) and
*     final (NFF) scale.
*
         IF(ABS(Q20).GT.ABS(Q2TH(3)))THEN
            NFI = 3
         ELSEIF(ABS(Q20).GT.ABS(Q2TH(2)))THEN
            NFI = 2
         ELSEIF(ABS(Q20).GT.ABS(Q2TH(1)))THEN
            NFI = 1
         ELSE
            NFI = 0
         ENDIF
         IF(NFI.GT.NFMAX) NFI = NFMAX
*     
         IF(ABS(Q2).GT.ABS(Q2TH(3)))THEN
            NFF = 3
         ELSEIF(ABS(Q2).GT.ABS(Q2TH(2)))THEN
            NFF = 2
         ELSEIF(ABS(Q2).GT.ABS(Q2TH(1)))THEN
            NFF = 1
         ELSE
            NFF = 0
         ENDIF
         IF(NFF.GT.NFMAX) NFF = NFMAX
*
         IF(NFI.EQ.0)THEN
            IF(NFF.EQ.0)THEN
*     NO CROSSING
               CALL CROSS01(ZN,Q20,Q2,NFI,ZFUNCNS,ZFUNCSG,
     1                      ZFUNCNS3,ZFUNCNS8,
     2                      ZFUNCNSV3,ZFUNCNSV8)
            ELSEIF(NFF.EQ.1)THEN
*     ME2 CROSSING
               CALL CROSS11(ZN,Q20,Q2,NFI,ZFUNCNS,ZFUNCSG,
     1                      ZFUNCNS3,ZFUNCNS8,
     2                      ZFUNCNSV3,ZFUNCNSV8)
            ELSEIF(NFF.EQ.2)THEN
*     ME2 AND MM2 CROSSING
               CALL CROSS21(ZN,Q20,Q2,NFI,ZFUNCNS,ZFUNCSG,
     1                      ZFUNCNS3,ZFUNCNS8,
     2                      ZFUNCNSV3,ZFUNCNSV8)
            ELSEIF(NFF.EQ.3)THEN
*     ME2, MM2 AND MT2 CROSSING
               CALL CROSS31(ZN,Q20,Q2,NFI,ZFUNCNS,ZFUNCSG,
     1                      ZFUNCNS3,ZFUNCNS8,
     2                      ZFUNCNSV3,ZFUNCNSV8)
            ELSE
               WRITE(6,*) "In src/evolution/zfunc.f:"
               WRITE(6,*) 'Undefined Final Flavour Number, NFF =',NFF
               CALL EXIT(-10)
            ENDIF
         ELSEIF(NFI.EQ.1)THEN
            IF(NFF.EQ.1)THEN
*     NO CROSSING
               CALL CROSS02(ZN,Q20,Q2,NFI,ZFUNCNS,ZFUNCSG,
     1                      ZFUNCNS3,ZFUNCNS8,
     2                      ZFUNCNSV3,ZFUNCNSV8)
            ELSEIF(NFF.EQ.2)THEN
*     MM2 CROSSING
               CALL CROSS12(ZN,Q20,Q2,NFI,ZFUNCNS,ZFUNCSG,
     1                      ZFUNCNS3,ZFUNCNS8,
     2                      ZFUNCNSV3,ZFUNCNSV8)
            ELSEIF(NFF.EQ.3)THEN
*     MM2 AND MT2 CROSSING
               CALL CROSS22(ZN,Q20,Q2,NFI,ZFUNCNS,ZFUNCSG,
     1                      ZFUNCNS3,ZFUNCNS8,
     2                      ZFUNCNSV3,ZFUNCNSV8)
            ELSE
               WRITE(6,*) "In src/evolution/zfunc.f:"
               WRITE(6,*) 'Undefined Final Flavour Number, NFF =',NFF
               CALL EXIT(-10)
            ENDIF
         ELSEIF(NFI.EQ.2)THEN
            IF(NFF.EQ.2)THEN
*     NO CROSSING
               CALL CROSS03(ZN,Q20,Q2,NFI,ZFUNCNS,ZFUNCSG,
     1                      ZFUNCNS3,ZFUNCNS8,
     2                      ZFUNCNSV3,ZFUNCNSV8) 
            ELSEIF(NFF.EQ.3)THEN
*     MT2 CROSSING
               CALL CROSS13(ZN,Q20,Q2,NFI,ZFUNCNS,ZFUNCSG,
     1                      ZFUNCNS3,ZFUNCNS8,
     2                      ZFUNCNSV3,ZFUNCNSV8)
            ELSE
               WRITE(6,*) "In src/evolution/zfunc.f:"
               WRITE(6,*) 'Undefined Final Flavour Number, NFF =',NFF
               CALL EXIT(-10)
            ENDIF
         ELSEIF(NFI.EQ.3)THEN
*     NO CROSSING
            CALL CROSS04(ZN,Q20,Q2,NFI,ZFUNCNS,ZFUNCSG,
     1                   ZFUNCNS3,ZFUNCNS8,
     2                   ZFUNCNSV3,ZFUNCNSV8)
         ELSE
            WRITE(6,*) "In src/evolution/zfunc.f:"
            WRITE(6,*) 'Undefined Initial Flavour Number, NFI =',NFI
            CALL EXIT(-10)
         ENDIF
      ENDIF
*
      RETURN
      END
