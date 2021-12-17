***********************************************************************
*
*     alphaMZ_analytic.f:
*
*     This subroutine returns the QED evolution factors, in N space,
*     from Q2II to Q2FF with NF active flavours, for Non Singlet and
*     Singlet-Gluon in the alphaMZ scheme using a (semi)analytical
*     approach.
*
***********************************************************************
      SUBROUTINE ALPHAMZ_ANALYTIC(ZN,Q2II,Q2FF,NF,EVF)
*
      IMPLICIT NONE
*
      include "../commons/ipt.h"
      include "../commons/activeflavours.h"
      include "../commons/facscheme.h"
      include "../commons/tecparam.h"
      include "../commons/alpha.h"
      include "../commons/beta.h"
      include "../commons/consts.h"
      include "../commons/nfsum.h"
**
*     Input Variables
*
      INTEGER NF
      DOUBLE PRECISION Q2II,Q2FF
      DOUBLE COMPLEX ZN
**
*     Internal Variables
*
      INTEGER I,J,K,IQ,NA
      INTEGER MP(0:3)
      DOUBLE PRECISION LN2
      DOUBLE PRECISION AFPI
      DOUBLE PRECISION MZ2
      DOUBLE PRECISION LNQ2I,LNQ2F,DLNQ2
      DOUBLE COMPLEX G0(4,4),G1(4,4)
      DOUBLE COMPLEX G0NS(3),G1NS(2,3)
      DOUBLE COMPLEX JLL,JGL,FDEN,FDGMN
      DOUBLE COMPLEX JM(4,4),JS(4,4),JSINV(4,4)
      DOUBLE COMPLEX SPSG(4,4),SPNS(2,3),EXPSPSG(4,4)
      DOUBLE COMPLEX SGTMP1(4,4)
      DOUBLE COMPLEX SGTMP2(4,4)
      DOUBLE COMPLEX SG(4,4)
      DOUBLE COMPLEX EFNS(2,3),EFSG(4,4),EFSGTMP(4,4)
      DATA MP / 1, 2, 8, 14 /
      PARAMETER(MZ2 = 91.1876D0**2)
**
*     Output Variables
*
      DOUBLE COMPLEX EVF(19,19)
*
*     LO splitting functions
*
      CALL ANDIM_LO(ZN,NF,G0NS,G0)
*
*     NLO splitting functions
*
      IF(IPT.GE.1)THEN
         CALL ANDIM_NLO(ZN,NF,G1NS,G1)
      ENDIF
*
*     Log of the scales
      LN2 = DLOG(Q2FF/Q2II)
*
*     Our splitting functions require alpha/(4*pi)
      AFPI = AREF / 4D0 / PI
*      
*     Singlet
*
*     Solution at LO
*     
      IF(IPT.EQ.0)THEN
         DO I=1,4
            DO J=1,4
               SPSG(I,J) = AFPI * G0(I,J) * LN2
            ENDDO
         ENDDO
*
*     Now we need to exponentiate SPSG that is a 4x4 matrix. In this
*     case, analytical formulas are too involved, therefore we adopt an
*     expanded approach truncating the series to the first NEXP+1 terms.
*
         CALL MATRIXEXP(NEXP,4,SPSG,EFSG)
*     
*     Non Singlet
*     
         DO I=1,2
            DO J=1,3
               SPNS(I,J) = AFPI * G0NS(J) * LN2
            ENDDO
         ENDDO
*     
*     Solution at NLO
*
      ELSE
*
*     Initialize evolution operator to unity
*
         DO I=1,4
            DO J=1,4
               EFSG(I,J) = (0D0, 0D0)
            ENDDO
            EFSG(I,I) = (1D0, 0D0)
         ENDDO
*
*     Delta-scheme corrections
*
         DO I=1,4
            DO J=1,4
               JM(I,J) = (0D0,0D0)
            ENDDO
         ENDDO
*
         JLL = (0D0, 0D0)
         JGL = (0D0, 0D0)
         IF (FACSCHEME.EQ."DELTA") THEN
            JLL = FDEN(ZN)
            JGL = FDGMN(ZN)
            JM(1,2) = JGL
            JM(2,2) = JLL
         ENDIF
*
         DLNQ2 = LN2 / DBLE(NINT)
         LNQ2I = DLOG(Q2II)
         DO IQ=1,NINT
            LNQ2F = LNQ2I + DLNQ2
            DO I=1,4
               DO J=1,4
                  SGTMP1(I,J) = G0(I,J) + AFPI * (
     1                 G1(I,J)
     2                 + ( 20D0 / 9D0 * NFSUM2(NF)
     3                 - BETA0(NF) * ( ( LNQ2I + LNQ2F ) / 2D0
     4                 - DLOG(MZ2) ) )
     5                 * G0(I,J)
     6                 )
                  JS(I,J) = AFPI * JM(I,J)
               ENDDO
               JS(I,I) = 1D0 + JS(I,I)
            ENDDO
*     
            DO I=1,4
               DO J=1,4
                  JSINV(I,J) = (0D0, 0D0)
               ENDDO
               JSINV(I,I) = (1D0, 0D0)
            ENDDO
            JSINV(1,2) = JSINV(1,2) - AFPI * JGL / ( 1D0 + AFPI * JLL )
            JSINV(2,2) = JSINV(2,2) / ( 1D0 + AFPI * JLL )
*     
            CALL MMULT(JS,4,4,SGTMP1,4,4,SGTMP2)
            CALL MMULT(SGTMP2,4,4,JSINV,4,4,SG)
*     
            DO I=1,4
               DO J=1,4
                  SPSG(I,J) = AFPI * SG(I,J) * DLNQ2
               ENDDO
            ENDDO
*
*     Exponentiate
*
            CALL MATRIXEXP(NEXP,4,SPSG,EXPSPSG)
*
*     Multiply evolutiion operator
*
            DO I=1,4
               DO J=1,4
                  EFSGTMP(I,J) = EFSG(I,J)
               ENDDO
            ENDDO
            CALL MMULT(EFSGTMP,4,4,EXPSPSG,4,4,EFSG)
*
            LNQ2I = LNQ2F
         ENDDO
*
*     Non Singlet (equal in MSbar and Delta scheme for alpha fixed)
*     
         DO I=1,2
            DO J=1,3
               SPNS(I,J) = AFPI * ( G0NS(J) + AFPI * ( G1NS(I,J) 
     1              + ( 20D0 / 9D0 * NFSUM2(NF)
     2              - BETA0(NF) * DLOG(DSQRT(Q2II * Q2FF) / MZ2) )
     3              * G0NS(J)
     4              ) ) * LN2
            ENDDO
         ENDDO
*         
      ENDIF
*
*     Exponentiate the non-singlet (here because same at LO and NLO)
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
         EVF(MP(I)+3,MP(I)+3) = EFNS(2,I)
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
                  EVF(MP(K)+I,MP(J-1)) = EFSG(K+1,J)
               ENDDO
               EVF(MP(K)+I+3,MP(K)+3)  = EFNS(2,K)
            ENDIF
         ENDDO
      ENDDO
*
      RETURN
      END
 
