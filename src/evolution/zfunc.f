************************************************************************
*
*     zfunc.f:
*     It returns the evolution kernels in N space from Q20 to Q2, 
*     according to eq. (20) of the "Notes on Perturbative Evolution"
*     (only forward evolution allowed).
*
*      - Fixed Flavour Number Scheme ("FFNS")
*      - Variable Flavour Number Scheme ("VFNS")
*
************************************************************************
      SUBROUTINE ZFUNC(ZN,Q20,Q2,
     1                 ZFUNCNS,ZFUNCNS15,ZFUNCNS24,ZFUNCNS35,
     2                 ZFUNCSG,ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35)
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
      DOUBLE COMPLEX ZFUNCNS15(2),ZFUNCNS24(2),ZFUNCNS35(2)
      DOUBLE COMPLEX ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35
*
*     Evolution Kernels
*
*     Fixed Flavour Number Schemes
*
      IF(NS.EQ."FFNS")THEN
         NFI = NFFN
         NFF = NFFN
         IF(NFFN.EQ.3)THEN
            CALL CROSS01(ZN,Q20,Q2,NFFN,ZFUNCNS,ZFUNCSG,
     1                   ZFUNCNS15,ZFUNCNS24,ZFUNCNS35,
     2                   ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35)
         ELSEIF(NFFN.EQ.4)THEN
            CALL CROSS02(ZN,Q20,Q2,NFFN,ZFUNCNS,ZFUNCSG,
     1                   ZFUNCNS15,ZFUNCNS24,ZFUNCNS35,
     2                   ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35)
         ELSEIF(NFFN.EQ.5)THEN
            CALL CROSS03(ZN,Q20,Q2,NFFN,ZFUNCNS,ZFUNCSG,
     1                   ZFUNCNS15,ZFUNCNS24,ZFUNCNS35,
     2                   ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35)
         ELSEIF(NFFN.EQ.6)THEN
            CALL CROSS04(ZN,Q20,Q2,NFFN,ZFUNCNS,ZFUNCSG,
     1                   ZFUNCNS15,ZFUNCNS24,ZFUNCNS35,
     2                   ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35)
         ENDIF
*
*     Variable Flavour Number Schemes
*
      ELSEIF(NS.EQ."VFNS")THEN
*
*     Determine number of active flavours at the initial (NFI) and
*     final (NFF) scale.
*
         IF(ABS(Q20).GT.ABS(Q2TH(6)))THEN
            NFI = 6
         ELSEIF(ABS(Q20).GT.ABS(Q2TH(5)))THEN
            NFI = 5
         ELSEIF(ABS(Q20).GT.ABS(Q2TH(4)))THEN
            NFI = 4
         ELSE
            NFI = 3
         ENDIF
         IF(NFI.GT.NFMAX) NFI = NFMAX
*     
         IF(ABS(Q2).GT.ABS(Q2TH(6)))THEN
            NFF = 6
         ELSEIF(ABS(Q2).GT.ABS(Q2TH(5)))THEN
            NFF = 5
         ELSEIF(ABS(Q2).GT.ABS(Q2TH(4)))THEN
            NFF = 4
         ELSE
            NFF = 3
         ENDIF
         IF(NFF.GT.NFMAX) NFF = NFMAX
*
         IF(NFI.EQ.3)THEN
            IF(NFF.EQ.3)THEN
*     NO CROSSING
               CALL CROSS01(ZN,Q20,Q2,NFI,ZFUNCNS,ZFUNCSG,
     1                      ZFUNCNS15,ZFUNCNS24,ZFUNCNS35,
     2                      ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35)
            ELSEIF(NFF.EQ.4)THEN
*     MC2 CROSSING
               CALL CROSS11(ZN,Q20,Q2,NFI,ZFUNCNS,ZFUNCSG,
     1                      ZFUNCNS15,ZFUNCNS24,ZFUNCNS35,
     2                      ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35)
            ELSEIF(NFF.EQ.5)THEN
*     MC2 AND MB2 CROSSING
               CALL CROSS21(ZN,Q20,Q2,NFI,ZFUNCNS,ZFUNCSG,
     1                      ZFUNCNS15,ZFUNCNS24,ZFUNCNS35,
     2                      ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35)
            ELSEIF(NFF.EQ.6)THEN
*     MC2, MB2 AND MT2 CROSSING
               CALL CROSS31(ZN,Q20,Q2,NFI,ZFUNCNS,ZFUNCSG,
     1                      ZFUNCNS15,ZFUNCNS24,ZFUNCNS35,
     2                      ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35)
            ELSE
               WRITE(6,*) "In src/evolution/zfunc.f:"
               WRITE(6,*) 'Undefined Final Flavour Number, NFF =',NFF
               CALL EXIT(-10)
            ENDIF
         ELSEIF(NFI.EQ.4)THEN
            IF(NFF.EQ.4)THEN
*     NO CROSSING
               CALL CROSS02(ZN,Q20,Q2,NFI,ZFUNCNS,ZFUNCSG,
     1                      ZFUNCNS15,ZFUNCNS24,ZFUNCNS35,
     2                      ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35)
            ELSEIF(NFF.EQ.5)THEN
*     MB2 CROSSING
               CALL CROSS12(ZN,Q20,Q2,NFI,ZFUNCNS,ZFUNCSG,
     1                      ZFUNCNS15,ZFUNCNS24,ZFUNCNS35,
     2                      ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35)
            ELSEIF(NFF.EQ.6)THEN
*     MB2 AND MT2 CROSSING
               CALL CROSS22(ZN,Q20,Q2,NFI,ZFUNCNS,ZFUNCSG,
     1                      ZFUNCNS15,ZFUNCNS24,ZFUNCNS35,
     2                      ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35)
            ELSE
               WRITE(6,*) "In src/evolution/zfunc.f:"
               WRITE(6,*) 'Undefined Final Flavour Number, NFF =',NFF
               CALL EXIT(-10)
            ENDIF
         ELSEIF(NFI.EQ.5)THEN
            IF(NFF.EQ.5)THEN
*     NO CROSSING
               CALL CROSS03(ZN,Q20,Q2,NFI,ZFUNCNS,ZFUNCSG,
     1                      ZFUNCNS15,ZFUNCNS24,ZFUNCNS35,
     2                      ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35) 
            ELSEIF(NFF.EQ.6)THEN
*     MT2 CROSSING
               CALL CROSS13(ZN,Q20,Q2,NFI,ZFUNCNS,ZFUNCSG,
     1                      ZFUNCNS15,ZFUNCNS24,ZFUNCNS35,
     2                      ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35)
            ELSE
               WRITE(6,*) "In src/evolution/zfunc.f:"
               WRITE(6,*) 'Undefined Final Flavour Number, NFF =',NFF
               CALL EXIT(-10)
            ENDIF
         ELSEIF(NFI.EQ.6)THEN
*     NO CROSSING
            CALL CROSS04(ZN,Q20,Q2,NFI,ZFUNCNS,ZFUNCSG,
     1                   ZFUNCNS15,ZFUNCNS24,ZFUNCNS35,
     2                   ZFUNCNSV15,ZFUNCNSV24,ZFUNCNSV35)
         ELSE
            WRITE(6,*) "In src/evolution/zfunc.f:"
            WRITE(6,*) 'Undefined Initial Flavour Number, NFI =',NFI
            CALL EXIT(-10)
         ENDIF
      ENDIF
*
      RETURN
      END
