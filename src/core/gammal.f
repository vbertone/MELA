************************************************************************
*
*     Log of Gamma function with a complex argument 
*
************************************************************************
      SUBROUTINE GAMMAL(ZZ,RES)
*
      implicit double precision (a-h, o-z)
*
      DOUBLE COMPLEX T1,T2,ZZ,Z,X,X2,RES,t
*
      Z=ZZ
      PI = 3.14159 26535 89793 23846 2643 D+0
*
      ONE=DCMPLX(1.0D0,0.0D0)
      T=DCMPLX(0.0D0,0.0D0)
2     R=SQRT(DREAL(Z)**2+DIMAG(Z)**2)
      IF(R.GT.10.0D0) GOTO 1
      T=T-LOG(Z)
      Z=Z+ONE
      GOTO 2
1     CONTINUE
*
      T1=Z*(LOG(Z)-1.0D0)+LOG(2.0D0*PI/Z)/2.0D0
*
      X=ONE/Z
      X2=X*X
      T2 = (1.D0/12.D0+(-1.D0/360.D0+(1.D0/1260.D0+(-1.D0/1680.D0+(1.D0/
     #1188.D0+(-691.D0/360360.D0+(1.D0/156.D0-3617.D0/122400.D0*X2)*X2
     #)*X2)*X2)*X2)*X2)*X2)*X
*
      RES=T1+T2+T
*
      RETURN
      END
