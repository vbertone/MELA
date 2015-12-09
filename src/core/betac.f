************************************************************************
*
*     BETA(A,B) FOR COMPLEX ARGUMENTS
*
************************************************************************
      FUNCTION BETAC(A,B)
*
      IMPLICIT NONE
*
      DOUBLE COMPLEX T1,T2,T3,T,BETAC,A,B
*
      CALL GAMMAL(A,T1)
      CALL GAMMAL(B,T2)
      CALL GAMMAL(A+B,T3)
*
      T = T1 + T2 - T3
*
      BETAC = ZEXP(T)
*
      RETURN
      END
