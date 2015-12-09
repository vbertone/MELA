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
