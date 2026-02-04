************************************************************************
*
*     KL series implementation of MT of Log(2-z)
*     z^m Log(2-z) = Log(2) - sum_1^inf (1/n)(z/2)^n
*     log[2]/(N+m) - sum_(k=1)^40 (1/k)(1/2)^k 1/(N+k+m)
*
************************************************************************
      FUNCTION LOG2KL(N,m)
*
      IMPLICIT NONE
*
*     Input Variables
      DOUBLE COMPLEX N
      DOUBLE PRECISION m
*
*     Output Variable
      DOUBLE COMPLEX LOG2KL
*
*     Internal Variables
      INTEGER k,kcheck
      DOUBLE COMPLEX TERM,SUM
      DOUBLE PRECISION THRESH
*
*     Initialization
      SUM = (0.0D0, 0.0D0)
      THRESH = 1.0D-6
*
*     Series Summation
      DO k = 1, 40
         TERM = 1d0/k*0.5d0**k/(N+m+k)
         SUM = SUM + TERM
      END DO
      LOG2KL=0.693147d0/(N+m)-SUM
*
*     Check the N=41 term
      kcheck=41
      TERM = 1/kcheck*0.5d0**kcheck/(N+m+kcheck)
C
      IF (ABS(TERM) .GT. THRESH) THEN
         WRITE(6,*) 'Warning: Term 101 =', TERM, ' exceeds threshold.'
      END IF
C
      RETURN
      END     
************************************************************************
*
*     KL series implementation 
*     sum_n=0^infty (9+k)!/k! (1/2)^k/(N+m+k) * 1/2^10 * 1/9!
*     Mellin Transform of (2-z)^-10
*
************************************************************************
      FUNCTION SMKL10(N,m)
*
      IMPLICIT NONE
*
*     Input Variables
      DOUBLE COMPLEX N
      DOUBLE PRECISION m
*
*     Output Variable
      DOUBLE COMPLEX SMKL10
*
*     Internal Variables
      INTEGER k,kcheck
      DOUBLE COMPLEX TERM
      DOUBLE PRECISION THRESH
*
*     Initialization
      SMKL10 = (0.0D0, 0.0D0)
      THRESH = 1.0D-6
*
*     Series Summation
      DO k = 0, 100
         TERM = (9d0+k)*(8d0+k)*(7d0+k)*(6d0+k)*(5d0+k)
     &       *(4d0+k)*(3d0+k)*(2d0+k)
     &       *(1d0+k)*0.5d0**k/(N+m+k)
         SMKL10 = SMKL10 + TERM
      END DO
      SMKL10=SMKL10*1d0/362880d0/1024d0
*
*     Check the N=101 term
      kcheck=101
      TERM = (9d0+kcheck)*(8d0+kcheck)
     &       *(7d0+kcheck)*(6d0+kcheck)
     &       *(5d0+kcheck)*(4d0+kcheck)*(3d0+kcheck)
     &       *(2d0+kcheck)*(1d0+kcheck)*0.5d0**kcheck/(N+m+kcheck)
C
      IF (ABS(TERM) .GT. THRESH) THEN
         WRITE(6,*) 'Warning: Term 101 =', TERM, ' exceeds threshold.'
      END IF
C
      RETURN
      END       
************************************************************************
*
*     KL series implementation 
*     sum_n=0^infty (7+k)!/k! (1/2)^k/(N+m+k) * 1/2^8 * 1/7!
*     Mellin Transform of (2-z)^-8
*
************************************************************************
      FUNCTION SMKL8(N,m)
*
      IMPLICIT NONE
*
*     Input Variables
      DOUBLE COMPLEX N
      DOUBLE PRECISION m
*
*     Output Variable
      DOUBLE COMPLEX SMKL8
*
*     Internal Variables
      INTEGER k,kcheck
      DOUBLE COMPLEX TERM
      DOUBLE PRECISION THRESH
*
*     Initialization
      SMKL8 = (0.0D0, 0.0D0)
      THRESH = 1.0D-6
*
*     Series Summation
      DO k = 0, 100
         TERM = (7d0+k)*(6d0+k)*(5d0+k)
     &       *(4d0+k)*(3d0+k)*(2d0+k)
     &       *(1d0+k)*0.5d0**k/(N+m+k)
         SMKL8 = SMKL8 + TERM
      END DO
      SMKL8=SMKL8*1d0/5040d0/256d0
*
*     Check the N=101 term
      kcheck=101
      TERM = (7d0+kcheck)*(6d0+kcheck)
     &       *(5d0+kcheck)*(4d0+kcheck)*(3d0+kcheck)
     &       *(2d0+kcheck)*(1d0+kcheck)*0.5d0**kcheck/(N+m+kcheck)
C
      IF (ABS(TERM) .GT. THRESH) THEN
         WRITE(6,*) 'Warning: Term 101 =', TERM, ' exceeds threshold.'
      END IF
C
      RETURN
      END     
************************************************************************
*
*     KL series implementation 
*     sum_n=0^infty (4+k)!/k! (1/2)^k/(N+m+k) * 1/2^5 * 1/4!
*     Mellin Transform of (2-z)^-5
*
************************************************************************
      FUNCTION SMKL5(N,m)
*
      IMPLICIT NONE
*
*     Input Variables
      DOUBLE COMPLEX N
      DOUBLE PRECISION m
*
*     Output Variable
      DOUBLE COMPLEX SMKL5
*
*     Internal Variables
      INTEGER k,kcheck
      DOUBLE COMPLEX TERM
      DOUBLE PRECISION THRESH
*
*     Initialization
      SMKL5 = (0.0D0, 0.0D0)
      THRESH = 1.0D-6
*
*     Series Summation
      DO k = 0, 100
         TERM = (4d0+k)*(3d0+k)*(2d0+k)
     &       *(1d0+k)*0.5d0**k/(N+m+k)
         SMKL5 = SMKL5 + TERM
      END DO
      SMKL5=SMKL5*1d0/24d0/32d0
*
*     Check the N=101 term
      kcheck=101
      TERM = (4d0+kcheck)*(3d0+kcheck)
     &       *(2d0+kcheck)*(1d0+kcheck)*0.5d0**kcheck/(N+m+kcheck)
C
      IF (ABS(TERM) .GT. THRESH) THEN
         WRITE(6,*) 'Warning: Term 101 =', TERM, ' exceeds threshold.'
      END IF
C
      RETURN
      END     
************************************************************************
*
*     KL series implementation 
*     sum_n=0^infty (5+k)!/k! (1/2)^k/(N+m+k) * 1/2^6 * 1/5!
*     Mellin Transform of (2-z)^-6
*
************************************************************************
      FUNCTION SMKL6(N,m)
*
      IMPLICIT NONE
*
*     Input Variables
      DOUBLE COMPLEX N
      DOUBLE PRECISION m
*
*     Output Variable
      DOUBLE COMPLEX SMKL6
*
*     Internal Variables
      INTEGER k,kcheck
      DOUBLE COMPLEX TERM
      DOUBLE PRECISION THRESH
*
*     Initialization
      SMKL6 = (0.0D0, 0.0D0)
      THRESH = 1.0D-6
*
*     Series Summation
      DO k = 0, 100
         TERM = (5d0+k)*(4d0+k)*(3d0+k)*(2d0+k)
     &    	*(1d0+k)*0.5d0**k/(N+m+k)
         SMKL6 = SMKL6 + TERM
      END DO
      SMKL6=SMKL6*1d0/120d0/64d0
*
*     Check the N=101 term
      kcheck=101
      TERM = (5d0+kcheck)*(4d0+kcheck)*(3d0+kcheck)
     &       *(2d0+kcheck)*(1d0+kcheck)*0.5d0**kcheck/(N+m+kcheck)
*      Part that should be here to check poles (N+m+k)=0
*      IF (A.EQ.INT(REAL(A)).AND.REAL(A).LE.0D0.AND.AIMAG(A).EQ.0D0)THEN
*      	WRITE(6,*) 'Warning: A =', A, ' A!=0,-1,-2,...'
*      END IF
C
      IF (ABS(TERM) .GT. THRESH) THEN
         WRITE(6,*) 'Warning: Term 101 =', TERM, ' exceeds threshold.'
      END IF
C
      RETURN
      END     
************************************************************************
*
*     LerchTranscendent.f:
*
*     It returns the Lerch Transcendent function.
*
*     The function is written in terms of an infinate series
*     we truncate at N=100
*
************************************************************************
      FUNCTION LERCH(Z,S,A)
*
      IMPLICIT NONE
*
*     Input Variables
      DOUBLE COMPLEX Z, S, A
*
*     Output Variable
      DOUBLE COMPLEX LERCH
*
*     Internal Variables
      INTEGER N
      DOUBLE COMPLEX TERM
      DOUBLE PRECISION THRESH
*
*     Initialization
      LERCH = (0.0D0, 0.0D0)
      THRESH = 1.0D-6
*
*     Series Summation
      DO N = 0, 100
         TERM = Z**N / (DBLE(N) + A)**S
         LERCH = LERCH + TERM
      END DO
*
*     Check the N=101 term
      TERM = Z**101 / (101.0D0 + A)**S
      IF (A.EQ.INT(REAL(A)).AND.REAL(A).LE.0D0.AND.AIMAG(A).EQ.0D0)THEN
      	WRITE(6,*) 'Warning: A =', A, ' A!=0,-1,-2,...'
      END IF
C
      IF (ABS(TERM) .GT. THRESH) THEN
         WRITE(6,*) 'Warning: Term 101 =', TERM, ' exceeds threshold.'
      END IF
C
      RETURN
      END      
* =================================================================av==
* =====================================================================
*
* KL implementation S_3 = zeta(3) - (-1)^l/(l-1)! Gamma^(l-1)(N-1)
*
* =====================================================================
*
      FUNCTION STHREE (Z)
*
      IMPLICIT NONE
*
*     Input Variables
      DOUBLE COMPLEX Z
*
*     Output Variable
      DOUBLE COMPLEX STHREE
*
*     Internal Variables
      DOUBLE COMPLEX LERCH,L1,L2,L3
      DOUBLE PRECISION THRESH
       
       
* ..Shift of the argument using the functional equation
*
*
      L1 =(1d0,0d0)
      L2 =(3d0,0d0)
      L3= Z-L1
      STHREE = 1.2020569031d0 - LERCH(L1,L2,L3)
*
      RETURN
      END
************************************************************************
*
*     hypgeo.f:
*
*     It returns the hypergeometric function/series.
*
*     It takes the complex arguments A, B, C and Z with the constrain:
*     C = A + B or C = A + B - 1.
*
*     It merges the hypergeometric series around z = 0 and an expansion 
*     around z = 1.
*
************************************************************************
      FUNCTION HYPGEO(A,B,C,Z)
*
      IMPLICIT NONE
**
*     Input Variables
*
      DOUBLE COMPLEX A,B,C,Z
**
*     Internal Variables
*
      INTEGER SEL
      DOUBLE PRECISION EPS
      DOUBLE PRECISION PREC,PREC1
      DOUBLE PRECISION MODSUM1,MODSUM2,MODC

      PARAMETER(PREC=1D-8)
      PARAMETER(PREC1=1D-8)
**
*     Output Variables
*
      DOUBLE COMPLEX HYPGEO
*
      MODSUM1 = DSQRT( DREAL(A+B-C)**2D0 + DIMAG(A+B-C)**2D0 )
      MODSUM2 = DSQRT( DREAL(A+B-C-1D0)**2D0 + DIMAG(A+B-C-1D0)**2D0 )
      MODC = DSQRT( DREAL(C)**2D0 + DIMAG(C)**2D0 )
      IF(MODSUM1.LT.PREC1*MODC)THEN
         SEL = 0
      ELSEIF(MODSUM2.LT.PREC1*MODC)THEN
         SEL = 1
      ELSE
         WRITE(6,*) "In modules/utilities/hypgeo.f:"
         WRITE(6,*) "HYPGEO(A,B,C,Z) can be used only for:"
         WRITE(6,*) " 1) C = A + B"
         WRITE(6,*) " 2) C = A + B - 1"
         WRITE(6,*) "   "
         WRITE(6,*) "A =",A,",   B =",B
         WRITE(6,*) "But C =",C
         CALL EXIT(-10)
      ENDIF
*
      EPS = ( 1D0 - DREAL(Z) ) / ( 4D0 * DREAL(Z) )
      IF(EPS.GT.0.1)THEN                            !Selection of the transition point between z=0 and z=1
         CALL HGZ0(A,B,C,Z,PREC,HYPGEO)
      ELSE
         IF(SEL.EQ.0)THEN
            CALL HGZ1M0(A,B,Z,PREC,HYPGEO)
         ELSEIF(SEL.EQ.1)THEN
            CALL HGZ1MM1(A,B,Z,PREC,HYPGEO)
         ELSE
            WRITE(6,*) "In modules/utilities/hypgeo.f:"
            WRITE(6,*) "Unknow value of SEL =",SEL
            CALL EXIT(-10)
         ENDIF
      ENDIF
*
      RETURN
      END
*
************************************************************************
*
*     Hypergeometric Series around z = 0.
*
************************************************************************
      SUBROUTINE HGZ0(A,B,C,Z,PREC,HYPGEOMZ0)
*
      IMPLICIT NONE
**
*     Input Variables
*
      DOUBLE PRECISION PREC
      DOUBLE COMPLEX A,B,C,Z
**
*     Internal Variables
*
      INTEGER I,IMAX
      DOUBLE PRECISION PRECZ0
      DOUBLE COMPLEX PART
      DOUBLE COMPLEX AN,BN,CN,IC

      PARAMETER(IMAX=80)
**
*     Output Variables
*
      DOUBLE COMPLEX HYPGEOMZ0
*
      HYPGEOMZ0 = (1D0,0D0)
      PART = (1D0,0D0)
*
      DO I=1,IMAX
         AN = A + I - 1D0
         BN = B + I - 1D0
         CN = C + I - 1D0
         IC = I
*
         PART = PART * ( AN * BN / CN / IC ) * Z
         HYPGEOMZ0 = HYPGEOMZ0 + PART
*
         PRECZ0 = SQRT( DREAL(PART)**2D0 + AIMAG(PART)**2D0 )
         IF(PRECZ0.LT.PREC)THEN
            RETURN
         ELSEIF(PRECZ0.GE.1D30)THEN
            WRITE(6,*) "In modules/utilities/hypgeo.f:"
            WRITE(6,*) "Hypergeometric series around z=0 not converging"
            WRITE(6,*) "Precision =",PRECZ0
            CALL EXIT(-10)
         ENDIF
      ENDDO
*
      RETURN
      END
*
************************************************************************
*
*     Hypergeometric Series around z = 1 for c = a + b.
*     See eq. 15.3.10 pg. 559 of M. Abramowitz and I. A. Stegun. 
*     "Handbook of mathematical functions" .
*     http://people.math.sfu.ca/~cbm/aands/page_559.htm
*
************************************************************************
      SUBROUTINE HGZ1M0(A,B,Z,PREC,HYPGEOMZ1)
*
      IMPLICIT NONE
**
*     Input Variables
*
      DOUBLE PRECISION PREC
      DOUBLE COMPLEX A,B,Z
**
*     Internal Variables
*
      INTEGER I,IMAX
      DOUBLE PRECISION PRECZ1
      DOUBLE COMPLEX FACT
      DOUBLE COMPLEX PART
      DOUBLE COMPLEX GA,GB,GAB
      DOUBLE COMPLEX PSI
      DOUBLE COMPLEX PSI1,PSI2,PSI3
      DOUBLE COMPLEX AN,BN,IC

      PARAMETER(IMAX=80)
**
*     Output Variables
*
      DOUBLE COMPLEX HYPGEOMZ1
*
      CALL GAMMAL(A,GA)
      CALL GAMMAL(B,GB)
      CALL GAMMAL(A+B,GAB)
*
      FACT = EXP( GAB - GA - GB )
*
      PSI1 = PSI((1D0,0D0))
      PSI2 = PSI(A)
      PSI3 = PSI(B)
*
      HYPGEOMZ1 = FACT * ( 2D0 * PSI1 - PSI2 - PSI3 - ZLOG(1D0-Z) )
      PART = (1D0,0.0D0)
*
      DO I=1,IMAX
         AN = A + I - 1D0
         BN = B + I - 1D0
         IC = I
*
         PSI1 = PSI1 + 1D0 / DBLE(I)
         PSI2 = PSI2 + 1D0 / AN
         PSI3 = PSI3 + 1D0 / BN
*
         PART = PART * ( AN * BN / IC / IC ) * ( 1D0 - Z )
*
         HYPGEOMZ1 = HYPGEOMZ1 + FACT * PART 
     1             * ( 2D0 * PSI1 - PSI2 - PSI3 - ZLOG(1D0-Z) )
*
         PRECZ1 = SQRT( DREAL(PART)**2D0 + AIMAG(PART)**2D0 )
         IF(PRECZ1.LT.PREC)THEN
            RETURN
         ELSEIF(PRECZ1.GE.1D30)THEN
            WRITE(6,*) "In modules/utilities/hypgeo.f:"
            WRITE(6,*) "Hypergeometric series around z=1 not converging"
            WRITE(6,*) "Precision =",PRECZ1
            CALL EXIT(-10)
         ENDIF
      ENDDO
*
      RETURN
      END
*
************************************************************************
*
*     Hypergeometric Series around z = 1 for c = a + b - 1.
*     See eq. 15.3.12 pg. 560 of M. Abramowitz and I. A. Stegun. 
*     "Handbook of mathematical functions".
*     http://people.math.sfu.ca/~cbm/aands/page_560.htm
*
************************************************************************
      SUBROUTINE HGZ1MM1(A,B,Z,PREC,HYPGEOMZ1)
*
      IMPLICIT NONE
**
*     Input Variables
*
      DOUBLE PRECISION PREC
      DOUBLE COMPLEX A,B,Z
**
*     Internal Variables
*
      INTEGER I,IMAX
      DOUBLE PRECISION PRECZ1
      DOUBLE COMPLEX FACT
      DOUBLE COMPLEX PART
      DOUBLE COMPLEX GA,GB,GAB
      DOUBLE COMPLEX PSI
      DOUBLE COMPLEX PSI1,PSI2,PSI3
      DOUBLE COMPLEX AN,BN,IC

      PARAMETER(IMAX=80)
**
*     Output Variables
*
      DOUBLE COMPLEX HYPGEOMZ1
*
      CALL GAMMAL(A-1D0,GA)
      CALL GAMMAL(B-1D0,GB)
      CALL GAMMAL(A+B-1D0,GAB)
*
      FACT = EXP( GAB - GA - GB )
*
      PSI1 = PSI((1D0,0D0))
      PSI2 = PSI(A)
      PSI3 = PSI(B)
*
      HYPGEOMZ1 = FACT * ( 1D0 / ( A - 1D0) / ( B - 1D0 ) / ( 1D0 - Z )
     1          + ZLOG(1D0-Z) - 2D0 * PSI1 + PSI2 + PSI3 - 1D0 )
      PART = (1D0,0D0)
*
      DO I=1,IMAX
         AN = A + I - 1D0
         BN = B + I - 1D0
         IC = I
*
         PSI1 = PSI1 + 1D0 / DBLE(I)
         PSI2 = PSI2 + 1D0 / AN
         PSI3 = PSI3 + 1D0 / BN
*
         PART = PART * ( AN * BN / IC / IC ) * ( 1D0 - Z )
*
         HYPGEOMZ1 = HYPGEOMZ1 + FACT * PART
     1             * ( ZLOG(1D0-Z) - 2D0 * PSI1 + PSI2 + PSI3 
     2             - 1D0 / ( IC + 1D0 ) ) / ( IC + 1D0 )
*
         PRECZ1 = SQRT( DREAL(PART)**2D0 + AIMAG(PART)**2D0 )
         IF(PRECZ1.LT.PREC)THEN
            RETURN
         ELSEIF(PRECZ1.GE.1D30)THEN
            WRITE(6,*) "In modules/utilities/hypgeo.f:"
            WRITE(6,*) "Precision =",PRECZ1
            CALL EXIT(-10)
         ENDIF
      ENDDO
*
      RETURN
      END
