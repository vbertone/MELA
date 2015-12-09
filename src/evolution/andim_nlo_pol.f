*************************************************************
*
*     File: andim_nlo_pol.f
*     
*     Returns the polarized NLO anomalous dimensions in the 
*     a_s=alpha_s/4pi notation
*
*     Input variables: N  : Mellin variable
*                      NF : Number of flavours
*
*
*     Output: P1NS : Complex valued Non-singlet LO anomalous
*                    dimension in N space
*
*             P1SG(2,2) : matrix of the complex valued LO
*                         singlet anomalous dimensions in 
*                         N space
*
*************************************************************
      SUBROUTINE ANDIM_NLO_POL(N,NF,P1NS,P1SG) 
*
      IMPLICIT NONE
*
      include "../commons/colfact.h"
      include "../commons/consts.h"
*
* ---------------------------------------------------------------------
*
*     Internal variables
*
      INTEGER I,J
      CHARACTER*2 ISCH
      DOUBLE COMPLEX N2,N3,NP1,NM1,NP2,NNP13,NNP1,NNP12
      DOUBLE COMPLEX S1,S2,S3,S1T,S2T,S3T,S12S21,S12T,S1N,S2N
      DOUBLE COMPLEX CSU1P,CSU2P,CSU3P,CSU1M,CSU2M,CSU3M
      DOUBLE COMPLEX CSUM1,CSUM2,CSUM3,CSUM1P,CSUM2P,CSUM3P,
     1               CSUM1M,CSUM2M,CSUM3M,CSUMTP,CSUMTM
      COMMON/SOMM/CSUM1,CSUM2,CSUM3,CSUM1P,CSUM2P,CSUM3P,
     1            CSUM1M,CSUM2M,CSUM3M,CSUMTP,CSUMTM
      DOUBLE COMPLEX FQ,FG,FQMS,FGMS,ZETA(4),P1QQ,P1GG,P1QG,P1GQ
      DOUBLE COMPLEX P1QQNSP,P1QQNSM,P1QQPS,CZERO,P0(2,2),XB0,XB1
      REAL*8 K1,K2,ETAP,ETAM
*
* ---------------------------------------------------------------------
*
*     Input variables
*
      DOUBLE COMPLEX N
      INTEGER NF       
*
*     Output variables  
*
      DOUBLE COMPLEX P1NS(3)
      DOUBLE COMPLEX P1SG(2,2)
*
* ---------------------------------------------------------------------
*     
      CZERO=CMPLX(0D0,0D0)

      DO I=1,3
         P1NS(I)=CZERO
      ENDDO
*
      DO I=1,2
         DO J=1,2
            P1SG(I,J)=CZERO
         ENDDO
      ENDDO

* factorization and renormalization scale factors
      K1=1
      K2=1
* factorization scheme
      ISCH='MS'
*
      N2=N*N
      N3=N*N2
      NM1=N-1
      NP1=N+1
      NP2=N+2
      NNP1=N*NP1
      NNP12=NNP1*NNP1
      NNP13=NNP12*NNP1
      CALL SOMME(NM1)
      S1=CSUM1
      S2=CSUM2
      S3=CSUM3
      S1T=CSUM1P-S1
      S2T=0.50D0*CSUM2P-S2
      S3T=0.25D0*CSUM3P-S3
      S12S21=S1*S2+S3
      S12T=S1*S2T+S3T-CSUMTP
      S1N=CSUM1+1/N
      S2N=CSUM2+1/N2
C PQQ-PQQBAR
      CSU1M=CSUM1M+2D0/N
      CSU2M=CSUM2M+4D0/N2
      CSU3M=CSUM3M+8D0/N3
      CSUMTM=CSUMTM+S1N/N2
      ETAM=1D0
C PQQ+PQQBAR
      CSU1P=CSUM1P
      CSU2P=CSUM2P
      CSU3P=CSUM3P
      CSUMTP=CSUMTP-S1N/N2
      ETAP=-1D0

      XB0=(33-2*NF)/6D0
      XB1=(153-19*NF)/(33D0-2D0*NF)
c splitting functions: p0(i,j) + alpha/(2*pi)*p1(i,j)
c one-loop splitting functions: 
      P0(1,1)=CF*(1.5D0-1/N-1/NP1-2*S1)
      P0(1,2)=2*NF*TR*NM1/NNP1
      P0(2,1)=CF*NP2/NNP1
      P0(2,2)=XB0+CA*(-2*NM1/NNP1-2*S1)
c two-loop non-singlet (Floratos et al divided by -8)
      P1QQNSM=-CF**2*(S1N*(4*N+2)/NNP12+2*(2*S1N-1/NNP1)*(S2N-CSU2M)
     #         +3*S2N+8*CSUMTM-CSU3M-3D0/8D0-((3*N+1)*N2-1)/NNP13
     #         -ETAM*((4*N+4)*N+2)/NNP13)
     #       -CF*CA*(67/9D0*S1N-(2*S1N-1/NNP1)*(2*S2N-CSU2M)
     #         -11/3D0*S2N-4*CSUMTM+0.5D0*CSU3M-17/24D0
     #         -1/18D0*((((151*N+236)*N+88)*N+3)*N+18)/NNP13
     #         +ETAM*((2*N+2)*N+1)/NNP13)
     #       -CF*TR*NF*(-20*S1N/9D0+4D0/3D0*S2N+1D0/6D0
     #         +((22*N+10)*N-6)/(9*NNP12))
      P1QQNSP=-CF**2*(S1N*(4*N+2)/NNP12+2*(2*S1N-1/NNP1)*(S2N-CSU2P)
     #         +3*S2N+8*CSUMTP-CSU3P-3D0/8D0-((3*N+1)*N2-1)/NNP13
     #         -ETAP*((4*N+4)*N+2)/NNP13)
     #       -CF*CA*(67/9D0*S1N-(2*S1N-1/NNP1)*(2*S2N-CSU2P)
     #         -11/3D0*S2N-4*CSUMTP+0.5D0*CSU3P-17/24D0
     #         -1/18D0*((((151*N+236)*N+88)*N+3)*N+18)/NNP13
     #         +ETAP*((2*N+2)*N+1)/NNP13)
     #       -CF*TR*NF*(-20*S1N/9D0+4D0/3D0*S2N+1D0/6D0
     #         +((22*N+10)*N-6)/(9*NNP12))
C two-loop singlet (Mertig-Van Neerven divided by -8, erratum included)
      P1QQPS=-2*CF*TR*NF*NP2*(N3+2*N+1)/NNP13
      P1QQ = P1QQNSP+P1QQPS
      P1QG=-2*CA*TR*NF*(-S1**2/N+2*S1**2/NP1-2*S1/N2+4*S1/NP1**2-S2/N
     #        +2*S2/NP1-2*S2T/N+4*S2T/NP1-4/N+3/NP1-3/N2+8/NP1**2
     #        +2/N3+12/NP1**3 )
     #     -CF*TR*NF*(2*S1**2/N-4*S1**2/NP1-2*S2/N+4*S2/NP1
     #        +14/N-19/NP1-1/N2-8/NP1**2-2/N3+4/NP1**3 ) 
      P1GQ=-CF*CA*(-2*S1**2/N+S1**2/NP1+16*S1/(3*N)-5*S1/(3*NP1)+2*S2/N
     #        -S2/NP1+4*S2T/N-2*S2T/NP1-56/(9*N)-20/(9*NP1)+28/(3*N2)
     #        -38/(3*NP1**2)-4/N3-6/NP1**3 )
     #     -0.5D0*CF**2*(4*S1**2/N-2*S1**2/NP1-8*S1/N+2*S1/NP1+8*S1/N2
     #        -4*S1/NP1**2+4*S2/N-2*S2/NP1+15/N-6/NP1-12/N2+3/NP1**2
     #        +4/N3-2/NP1**3 )
     #     -4*CF*TR*NF*(-2*S1/(3*N)+S1/(3*NP1)+7/(9*N)-2/(9*NP1)
     #        -2/(3*N2)+1/(3*NP1**2) )
      P1GG=-0.5D0*CA**2*(134/9D0*S1+8*S1/N2-16*S1/NP1**2+8*S2/N
     #        -16*S2/NP1+4*S3-8*S12S21+8*S2T/N-16*S2T/NP1+4*S3T-8*S12T
     #        -107/(9*N)+241/(9*NP1)+58/(3*N2)-86/(3*NP1**2)-8/N3
     #        -48/NP1**3-16/3D0   )
     #     -4*CA*TR*NF*(-5*S1/9D0+14/(9*N)-19/(9*NP1)-1/(3*N2)
     #        -1/(3*NP1**2)+1/3D0)
     #     -CF*TR*NF*(-10/NP1+2/NP1**2+4/NP1**3+1+10/N-10/N2+4/N3)
c The change of scheme from MSbar to any other is defined as follows:
c *** in the ABFR notation, powers of as/(2 pi) ***
c   [Delta Sigma]   [ 1+a*zeta(1)    a*zeta(2) ][Delta Sigma(MS)]
c   [           ] = [                          ][               ]
c   [Delta   g  ]   [   a*zeta(3)  1+a*zeta(4) ][  Delta g(MS)  ]
c and
c   Delta q_NS = ( 1+a*zeta(1) ) Delta q_NS(MS)
c where a=alfa/(2*pi).
c Coefficient functions are transformed according to
c                                [ 1-a*zeta(1)   -a*zeta(2) ]
c   [ cq   cg  ] = [ cqms  cgms ][                          ] + O(a^**2)
c                                [  -a*zeta(3)  1-a*zeta(4) ]
c                = [ cqms-a*zeta(1)  cgms-a*zeta(2) ] + O(a^**2),
c so that cq*Delta Sigma + cg*Delta g is unchanged to order a.
c We define fq and fg by
c   cq = 1 + a fq  -> zeta(1) = fqms - fq
c   cg = a fg      -> zeta(2) = fgms - fg
      FQMS=CF*(S1*(1.5D0+S1+1/N+1/NP1)-S2-9D0/2D0+3/N)
      FGMS=-2*TR*NF*NM1/NNP1*(S1+1)
      IF(ISCH.EQ.'MS')THEN
C MSbar:
         FQ=FQMS
         FG=FGMS
      ELSEIF(ISCH.EQ.'DI')THEN
C DIS scheme:
         FQ=CZERO
         FG=CZERO
      ELSEIF(ISCH.EQ.'AB')THEN
c Adler-Bardeen scheme: MSbar minimally modified
         FQ=FQMS
         FG=FGMS-2*NF*TR/N
      ELSEIF(ISCH.EQ.'OS')THEN
c off-shell initial parton (Kodaira)
         FQ=CF*(-1-3/(2*N)+4/NP1+2/N2-2/NP1**2+3/2D0*S1N-4*S2N)
         FG=-4*NF*TR*(N3-N2+N+1)/NNP12
      ELSEIF(ISCH.EQ.'AR')THEN
C massive quark (Altarelli-Ross)
         FQ=CF*(-5/2D0-1/(2*N)+2/NP1+1/N2-2/NP1**2
     #     +(3.5D0+1/(N*NP1)-S1N)*S1N-3*S2N)
         FG=-2*NF*TR*(N2+1+N*NM1*S1N)/(N2*NP1)
      ELSE
         WRITE(*,*)'You cannot choose isch=',ISCH
         STOP
      ENDIF
      ZETA(1)=FQMS-FQ
      ZETA(2)=FGMS-FG
      ZETA(3)=0
      ZETA(4)=0

      P1QQNSP=P1QQNSP-XB0*ZETA(1)
      P1QQNSM=P1QQNSM-XB0*ZETA(1)

* switch to NNPDF notation, powers of as/(4 pi)
      P1NS(1)=4*P1QQNSP
      P1NS(2)=4*P1QQNSM
      P1NS(3)=P1NS(2)


      P1SG(1,1)=P1QQ-XB0*ZETA(1)+ZETA(2)*P0(2,1)-ZETA(3)*P0(1,2)
     #     +XB0*P0(1,1)*LOG(K1*K2)
      P1SG(1,2)=P1QG+ZETA(2)*(P0(2,2)-P0(1,1)-XB0)
     #     +(ZETA(1)-ZETA(4))*P0(1,2)
     #     +XB0*P0(1,2)*LOG(K1*K2)
      P1SG(2,1)=P1GQ+ZETA(3)*(P0(1,1)-P0(2,2)-XB0)
     #     +(ZETA(4)-ZETA(1))*P0(2,1)
     #     +XB0*P0(2,1)*LOG(K1*K2)
      P1SG(2,2)=P1GG-XB0*ZETA(4)+ZETA(3)*P0(1,2)-ZETA(2)*P0(2,1)
     #     +XB0*P0(2,2)*LOG(K1*K2)



* switch to NNPDF notation, powers of as/(4 pi)
      DO I=1,2
         DO J=1,2
            P1SG(I,J)=4*P1SG(I,J)
         ENDDO
      ENDDO


****************************************************************************
      
*
       RETURN
       END
*
c------------------------------------------------------------------
      SUBROUTINE SOMME(XN)
c---calculation of all sums
      IMPLICIT NONE
      DOUBLE COMPLEX CSUM1,CSUM2,CSUM3,CSUM1P,CSUM2P,CSUM3P,
     .  CSUM1M,CSUM2M,CSUM3M,CSUMTP,CSUMTM,
     .  XN,XNON2,XNON2MH
      REAL*8 HALF
      PARAMETER(HALF=0.5d0)
      COMMON/SOMM/CSUM1,CSUM2,CSUM3,CSUM1P,CSUM2P,CSUM3P,
     .  CSUM1M,CSUM2M,CSUM3M,CSUMTP,CSUMTM
      CALL SUMEV(XN,CSUM1,CSUM2,CSUM3)
      XNON2=HALF*XN
      CALL SUMEV(XNON2,CSUM1P,CSUM2P,CSUM3P)
      XNON2MH=XNON2-HALF
      CALL SUMEV(XNON2MH,CSUM1M,CSUM2M,CSUM3M)
      CALL STILDE(XN,CSUMTP,CSUMTM)
      RETURN    
      END   
c-----------------------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION SUM1(N)
      IMPLICIT NONE
      external psi
      REAL*8 EULER,B2,B4,B6,B8,B10
      DOUBLE COMPLEX N,N1,N1I,N1SQI,PSI
      PARAMETER(EULER=.5772156649d0)
      PARAMETER(B2=1./6.,B4=-1./30.,B6=1./42.,B8=-1./30.,B10=5./66.)

      SUM1=CMPLX(0d0,0d0)
      N1=N+1
      SUM1 = EULER + PSI(N1)
c 8    IF (DBLE(N1) .LT. 10.) THEN
c         S1=S1-1./N1
c         N1=N1+1
c         GO TO 8
c      ENDIF
c      N1I=1./N1
c      N1SQI=N1I*N1I
c      S1=S1+EULER+LOG(N1)-(60.+(60.*B2+(30.*B4+(20.*B6+
c     . (15.*B8+12.*B10*N1SQI)*N1SQI)*N1SQI)*N1SQI)*N1I)/120.*N1I
      RETURN
      END
c-----------------------------------------------------------------------------
      SUBROUTINE SUMEV(N,S1,S2,S3)
      IMPLICIT NONE
      DOUBLE COMPLEX N,N1,N1I,N1SQI,S1,S2,S3
      REAL*8 EULER,ZETA2,ZETA3,B2,B4,B6,B8,B10
      PARAMETER(EULER=.5772156649)
      PARAMETER (ZETA2=1.644934067,ZETA3=1.2020569031)
      PARAMETER(B2=1./6.,B4=-1./30.,B6=1./42.,B8=-1./30.,B10=5./66.)
      S1=CMPLX(0d0,0d0)
      S2=CMPLX(0d0,0d0)
      S3=CMPLX(0d0,0d0)
      N1=N+1.
 10   IF(DBLE(N1) .LT. 10.) THEN
         S1=S1-1./N1
         S2=S2-1./N1**2
         S3=S3-1./N1**3
         N1=N1+1.      
         GO TO 10
      ENDIF
      N1I=1./N1
      N1SQI=N1I*N1I
      S1=S1+EULER+LOG(N1)-(60.+(60.*B2+(30.*B4+(20.*B6+
     . (15.*B8+12.*B10*N1SQI)*N1SQI)*N1SQI)*N1SQI)*N1I)/120.*N1I
      S2=S2+ZETA2-(1.+(0.5+(B2+(B4+
     . (B6+(B8+B10*N1SQI)*N1SQI)*N1SQI)*N1SQI)*N1I)*N1I)*N1I
      S3=S3+ZETA3-0.5*(1.+N1I+(3.*B2+(5.*B4+(7.*B6+
     . (9.*B8+11.*B10*N1SQI)*N1SQI)*N1SQI)*N1SQI)*N1SQI)*N1SQI
      RETURN
      END
c---------------------------------------------------------------------------
      SUBROUTINE STILDE(N,SP,SM)
      IMPLICIT NONE
      DOUBLE COMPLEX N,N1,N1I,N1SQI,SP,SM,TEMP,SUM1
      REAL*8 ETA,ZETA3,EULER,A1,A4,A5,A6,A7
      PARAMETER (ZETA3=1.2020569031,EULER=.5772156649)
      PARAMETER(A1=5./8.,A4=2./3.,A5=11./12.,A6=149./120.,A7=469./120.)
      SP=CMPLX(0d0,0d0)
      SM=CMPLX(0d0,0d0)
      N1=N+1.
      ETA=1.
C---eta=(-1)**(N+1)
   40 IF (DBLE(N1) .LT. 10.) THEN
         SP=SP+ETA*SUM1(N1)/N1**2
         SM=SM-ETA*SUM1(N1)/N1**2
         N1=N1+1.      
         ETA=-ETA
         GO TO 40
      ENDIF
      N1I=1./N1
      N1SQI=N1I*N1I
      TEMP=0.5*(((LOG(N1)+EULER)*(1.+(1.-(1.-3.*N1SQI)*N1SQI)*N1I)
     . +(A4+(A5-(A6+A7*N1I)*N1I)*N1I)*N1SQI)*N1SQI)
      SP=SP-A1*ZETA3+ETA*TEMP
      SM=SM-A1*ZETA3-ETA*TEMP
      RETURN
      END
