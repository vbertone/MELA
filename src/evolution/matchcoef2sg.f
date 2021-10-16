************************************************************************
*
*     NLO Time-like matching conditions (FFs)
*
*     Reference: hep-ph/0504192 (Eq. 27)
*
************************************************************************
      SUBROUTINE MATCHCOEF1SGT(N,A1SG)
*
      IMPLICIT NONE
*
      include "../commons/colfact.h"
**
*     Input Variables
*
      DOUBLE COMPLEX N
**
*     Internal Variables
*
      DOUBLE COMPLEX NM,N1,NS,NMS,N1S
**
*     Output Variables  
*
      DOUBLE COMPLEX A1SG
*
      NM  = N - 1D0
      N1  = N + 1D0
      NS  = N * N
      NMS = NM * NM
      N1S = N1 * N1
*
c      A1SG = 2D0 * CF * ( - 2D0 / NM +  2D0 / N - 1D0 / N1 ! Mine
c     1     + 4D0 / NMS - 4D0 / NS + 2D0 / N1S )
      A1SG = 2D0 * CF * ( - ( 2D0 + N + NS ) / N / ( NS - 1D0 ) ! Paper
     1      + 4D0 / NMS - 4D0 / NS + 2D0 / N1S )

*
      RETURN
      END
