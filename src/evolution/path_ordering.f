***********************************************************************
*
*     path_ordering.f:
*
*     This subroutine returns the QED evolution factors, in N space, from
*     Q2I to Q2F with NF active flavours, for Non Singlet and Singlet-Gluon
*     using the path ordering method.
*
***********************************************************************
      SUBROUTINE PATH_ORDERING(ZN,AII,AFF,NF,EVF)
*
      IMPLICIT NONE
*
      include "../commons/ipt.h"
      include "../commons/beta.h"
      include "../commons/activeflavours.h"
**
*     Input Variables
*
      INTEGER NF
      DOUBLE COMPLEX AII,AFF
      DOUBLE COMPLEX ZN
**
*     Internal Variables
*
      INTEGER NINT,NEXP
      INTEGER I,J,K,NA
      INTEGER MP(0:3)
      DOUBLE PRECISION PREC
      DOUBLE PRECISION BT0,BT1
      DOUBLE COMPLEX AF,AI,DA,AK
      DOUBLE COMPLEX G0(4,4),G1(4,4)
      DOUBLE COMPLEX G0NS(3),G1NS(2,3)
      DOUBLE COMPLEX SPSG(4,4),SPNS(2,3)
      DOUBLE COMPLEX DEFSG(4,4)
      DOUBLE COMPLEX EFNS(2,3),EFSG(4,4)
      PARAMETER(NINT=100)
      PARAMETER(PREC=1D-7)
      PARAMETER(NEXP=5)
      DATA MP / 1, 2, 8, 14 /
**
*     Output Variables
*
      DOUBLE COMPLEX EVF(19,19)
*
      BT0 = BETA0(NF)
      BT1 = 0D0
*
*     LO
*
      CALL ANDIM_LO(ZN,NF,G0NS,G0)
*
*     NLO
*
      IF(IPT.GE.1)THEN
         CALL ANDIM_NLO(ZN,NF,G1NS,G1)
         BT1 = BETA1(NF)
      ENDIF
*
*     Singlet
*
*     Solution at LO
*
      DA = ( AFF - AII ) / DBLE(NINT)
      IF(IPT.EQ.0)THEN
         AI = AII
         AF = AFF
         DO I=1,4
            DO J=1,4
               SPSG(I,J) = - DA * G0(I,J) * ( 1D0 / AF + 1D0 / AI )
     1              / BT0 / 2D0
            ENDDO
         ENDDO
*
         AK = AI
         DO K=1,NINT-1
            AK = AK + DA
            DO I=1,4
               DO J=1,4
                  SPSG(I,J) = SPSG(I,J) - DA * G0(I,J) / BT0 / AK
               ENDDO
            ENDDO
         ENDDO
*
*     Now we need to exponentiate SPSG that is a 4x4 matrix. In this
*     case, analytical formulas are too involved, therefore we adopt an
*     expanded approach truncating the series to the first NEXP+1 terms.
*
         CALL MATRIXEXP(NEXP,4,SPSG,EFSG)
*
*     Solution at NLO
*
      ELSE
         EFSG(1,1) = (1D0,0D0)
         EFSG(1,2) = (0D0,0D0)
         EFSG(2,1) = (0D0,0D0)
         EFSG(2,2) = (1D0,0D0)
*
         AI = AII
         DO K=1,NINT
            AF = AI + DA
            DO I=1,2
               DO J=1,2
                  SPSG(I,J) = - DA * ( ( G0(I,J) + AI * ( G1(I,J)
     1                 - G0(I,J) * BT1 / BT0 ) ) / BT0 / AI
     2                 + ( G0(I,J) + AF * ( G1(I,J)
     3                 - G0(I,J) * BT1 / BT0 ) ) / BT0 / AF ) / 2D0
               ENDDO
            ENDDO
*
*     Now we need to exponentiate SPSG that is a 4x4 matrix. In this
*     case, analytical formulas are too involved, therefore we adopt an
*     expanded approach truncating the series to the first NEXP+1 terms.
*
            CALL MATRIXEXP(NEXP,4,SPSG,DEFSG)
            CALL MMULT(DEFSG,4,4,EFSG,4,4,EFSG)
*
            AI = AI + DA
         ENDDO
      ENDIF
*
*     Non Singlet
*
      AI = AII
      AF = AFF
      DO I=1,2
         DO J=1,3
            SPNS(I,J) = - DA * ( ( G0NS(J)
     1           + AI * ( G1NS(I,J) - G0NS(J) * BT1/BT0 ) ) / BT0 / AI         
     2           + ( G0NS(J) + AF * ( G1NS(I,J) - G0NS(J) * BT1/BT0 ) )
     3           / BT0 / AF ) / 2D0
         ENDDO
      ENDDO
*
      AK = AI
      DO K=1,NINT-1
         AK = AK + DA
*
         DO I=1,2
            DO J=1,3
               SPNS(I,J) = SPNS(I,J)
     1              - DA * ( G0NS(J)
     2              + AK * ( G1NS(I,J) - G0NS(J) * BT1/BT0 ) )/BT0/AK
            ENDDO
         ENDDO
      ENDDO
*
      DO I=1,2
         DO J=1,3
            EFNS(I,J) = ZEXP( SPNS(I,J) )
         ENDDO
      ENDDO
*
*     Contruct evolution matrix according to nf
*
*     1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
*     g Sgl T1l T2l  Vl V1l V2l Sgu T1u T2u  Vu V1u V2u Sgd T1d T2d  Vd V1d V2d
*
*     Initialise to zero
*     
      DO I=1,19
         DO J=1,19
            EVF(I,J) = (0D0, 0D0)
         ENDDO
      ENDDO
*
*     Singlet matrix
*
      DO I=1,4
         DO J=1,4
            EVF(MP(I-1),MP(J-1)) = EFSG(I,J)
         ENDDO
      ENDDO
*
*     Total valence
*
      DO I=1,3
         EVF(MP(I)+3,MP(I)+3) = EFSG(I,J)
      ENDDO
*
      DO K=1,3
         IF (K.EQ.1) NA = NL(NF)
         IF (K.EQ.2) NA = NU(NF)
         IF (K.EQ.3) NA = ND(NF)
         DO I=1,2
            IF (NA.GT.I) THEN
               EVF(MP(K)+I,  MP(K)+I)   = EFNS(1,K)
               EVF(MP(K)+I+3,MP(K)+I+3) = EFNS(2,K)
            ELSE
               DO J=1,4
                  EVF(MP(K)+I,MP(J-1)) = EFSG(2,J)
               ENDDO
               EVF(MP(K)+I+3,MP(K)+3)  = EFNS(2,K)
            ENDIF
         ENDDO
      ENDDO
*
      RETURN
      END
 
