************************************************************************
*
*      mmult.f:
*
*      Multply a matrix by another matrix.
*
************************************************************************
      SUBROUTINE MMULT(A,ROWSA,COLSA,B,ROWSB,COLSB,C)
*
      IMPLICIT NONE
**
*     Input Variables
*
      INTEGER ROWSA,COLSA,ROWSB,COLSB
      DOUBLE COMPLEX A(ROWSA,COLSA),B(ROWSB,COLSB)
**
*     Internal Variables
*
      INTEGER I,J,K
**
*     Output Variables
*
      DOUBLE COMPLEX C(ROWSA,COLSB)
*
*     Check that the matrices to multiply have the  correct dimensions.
*
      IF(COLSA.NE.ROWSB)THEN
          WRITE(6,*) "In mmult.f:"
          WRITE(6,*) "Input matrices do not match"
          CALL EXIT(-10)
      ENDIF
*
*     Initialize the output matrix to zero.
*
      DO I=1,ROWSA
         DO J=1,COLSB
            C(I,J) = (0d0,0d0)
         ENDDO
      ENDDO
*
      DO I=1,ROWSA
         DO J=1,COLSB
            DO K=1,COLSA
	       C(I,J) = C(I,J) + A(I,K) * B(K,J)
	    ENDDO
	 ENDDO
      ENDDO
*
      RETURN
      END
