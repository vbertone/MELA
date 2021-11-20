************************************************************************
*
*     zfunc.f:
*     It returns the evolution kernels in N space from Q20 to Q2
*
************************************************************************
      SUBROUTINE ZFUNC(ZN,Q2,EFSG,EFV,EFT1,EFT2)
*
      IMPLICIT NONE
*
      include "../commons/massthrs.h"
      include "../commons/nffn.h"
      include "../commons/nfmax.h"
      include "../commons/ns.h"
      include "../commons/alpha.h"
**
*     Input Variables
*
      DOUBLE PRECISION Q2
      DOUBLE COMPLEX ZN
**
*     Output Variables
*
      INTEGER I,J
      DOUBLE COMPLEX EFSG(4,4),EFNS(2,3)
      DOUBLE COMPLEX EFV(3)
      DOUBLE COMPLEX EFT1(3,4),EFT2(3,4)
*
*     Evolution Kernels
*
*     1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
*     g Sgl T1l T2l  Vl V1l V2l Sgu T1u T2u  Vu V1u V2u Sgd T1d T2d  Vd V1d V2d
*
*     Fixed Flavour Number Schemes
*
      IF(NS.EQ."FFNS")THEN
*
*     Singlet and photon
*
         CALL PATH_ORDERING(ZN,AE,AQ,NFFN,EFNS,EFSG)
*
*     Total and non-singlet valence
*
         DO I=1,3
            EFV(I) = EFNS(2,I)
         ENDDO
*
*     Non singlet non-valence
*
         IF (NFFN.GE.1.AND.NFFN.LE.3) THEN
            DO I=1,3
               DO J=1,4
                  EFT1(I,J) = EFSG(I+1,J)
                  EFT2(I,J) = EFSG(I+1,J)
               ENDDO
            ENDDO

         ELSEIF (NFFN.EQ.4) THEN
            DO I=1,2
               DO J=1,4
                  EFT1(I,J) = EFSG(I+1,J)
               ENDDO
            ENDDO
            EFT1(3,1) = (0D0,0D0)
            EFT1(3,2) = (0D0,0D0)
            EFT1(3,3) = (0D0,0D0)
            EFT1(3,4) = EFNS(1,3)

            DO I=1,3
               DO J=1,4
                  EFT2(I,J) = EFSG(I+1,J)
               ENDDO
            ENDDO

         ELSEIF (NFFN.EQ.5) THEN
            EFT1(1,1) = (0D0,0D0)
            EFT1(1,2) = EFNS(1,1)
            EFT1(1,3) = (0D0,0D0)
            EFT1(1,4) = (0D0,0D0)
            DO J=1,4
               EFT1(2,J) = EFSG(3,J)
            ENDDO
            EFT1(3,1) = (0D0,0D0)
            EFT1(3,2) = (0D0,0D0)
            EFT1(3,3) = (0D0,0D0)
            EFT1(3,4) = EFNS(1,3)
            DO I=1,3
               DO J=1,4
                  EFT2(I,J) = EFSG(I+1,J)
               ENDDO
            ENDDO

         ELSEIF (NFFN.EQ.6) THEN
            EFT1(1,1) = (0D0,0D0)
            EFT1(1,2) = EFNS(1,1)
            EFT1(1,3) = (0D0,0D0)
            EFT1(1,4) = (0D0,0D0)
            EFT1(2,1) = (0D0,0D0)
            EFT1(2,2) = (0D0,0D0)
            EFT1(2,3) = EFNS(1,2)
            EFT1(2,4) = (0D0,0D0)
            EFT1(3,1) = (0D0,0D0)
            EFT1(3,2) = (0D0,0D0)
            EFT1(3,3) = (0D0,0D0)
            EFT1(3,4) = EFNS(1,3)
            DO I=1,3
               DO J=1,4
                  EFT2(I,J) = EFSG(I+1,J)
               ENDDO
            ENDDO

         ELSEIF (NFFN.EQ.7) THEN
            EFT1(1,1) = (0D0,0D0)
            EFT1(1,2) = EFNS(1,1)
            EFT1(1,3) = (0D0,0D0)
            EFT1(1,4) = (0D0,0D0)
            EFT1(2,1) = (0D0,0D0)
            EFT1(2,2) = (0D0,0D0)
            EFT1(2,3) = EFNS(1,2)
            EFT1(2,4) = (0D0,0D0)
            EFT1(3,1) = (0D0,0D0)
            EFT1(3,2) = (0D0,0D0)
            EFT1(3,3) = (0D0,0D0)
            EFT1(3,4) = EFNS(1,3)
            EFT2(1,1) = (0D0,0D0)
            EFT2(1,1) = EFNS(1,1)
            EFT2(1,1) = (0D0,0D0)
            EFT2(1,1) = (0D0,0D0)
            DO I=2,3
               DO J=1,4
                  EFT2(I,J) = EFSG(I+1,J)
               ENDDO
            ENDDO

         ELSEIF (NFFN.EQ.8) THEN
            EFT1(1,1) = (0D0,0D0)
            EFT1(1,2) = EFNS(1,1)
            EFT1(1,3) = (0D0,0D0)
            EFT1(1,4) = (0D0,0D0)
            EFT1(2,1) = (0D0,0D0)
            EFT1(2,2) = (0D0,0D0)
            EFT1(2,3) = EFNS(1,2)
            EFT1(2,4) = (0D0,0D0)
            EFT1(3,1) = (0D0,0D0)
            EFT1(3,2) = (0D0,0D0)
            EFT1(3,3) = (0D0,0D0)
            EFT1(3,4) = EFNS(1,3)
            EFT2(1,1) = (0D0,0D0)
            EFT2(1,1) = EFNS(1,1)
            EFT2(1,1) = (0D0,0D0)
            EFT2(1,1) = (0D0,0D0)
            DO J=1,4
               EFT2(2,J) = EFSG(3,J)
            ENDDO
            EFT2(3,1) = (0D0,0D0)
            EFT2(3,1) = (0D0,0D0)
            EFT2(3,1) = (0D0,0D0)
            EFT2(3,1) = EFNS(1,3)

         ELSEIF (NFFN.EQ.9) THEN
            EFT1(1,1) = (0D0,0D0)
            EFT1(1,2) = EFNS(1,1)
            EFT1(1,3) = (0D0,0D0)
            EFT1(1,4) = (0D0,0D0)
            EFT1(2,1) = (0D0,0D0)
            EFT1(2,2) = (0D0,0D0)
            EFT1(2,3) = EFNS(1,2)
            EFT1(2,4) = (0D0,0D0)
            EFT1(3,1) = (0D0,0D0)
            EFT1(3,2) = (0D0,0D0)
            EFT1(3,3) = (0D0,0D0)
            EFT1(3,4) = EFNS(1,3)
            EFT2(1,1) = (0D0,0D0)
            EFT2(1,1) = EFNS(1,1)
            EFT2(1,1) = (0D0,0D0)
            EFT2(1,1) = (0D0,0D0)
            EFT2(2,1) = (0D0,0D0)
            EFT2(2,1) = (0D0,0D0)
            EFT2(2,1) = EFNS(1,2)
            EFT2(2,1) = (0D0,0D0)
            EFT2(3,1) = (0D0,0D0)
            EFT2(3,1) = (0D0,0D0)
            EFT2(3,1) = (0D0,0D0)
            EFT2(3,1) = EFNS(1,3)

         ELSE
            WRITE(6,*) "In src/evolution/zfunc.f:"
            WRITE(6,*) 'Undefined NFFN =',NFFN
            CALL EXIT(-10)
         ENDIF
*
*     Variable Flavour Number Schemes
*
      ELSEIF(NS.EQ."VFNS")THEN

*=========================================================
*=========================================================
*=========================================================

      ENDIF
*
      RETURN
      END
