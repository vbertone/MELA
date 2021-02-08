************************************************************************
*
*     The following routines return the Zero-Mass coefficient functions
*     for semi-nclusive DIS structure functions in Mellin
*     space. Expressions are taken from
*     https://arxiv.org/pdf/hep-ph/0107064.pdf.
*
************************************************************************
      DOUBLE COMPLEX FUNCTION CL1QQ_SIDIS(N, M)
*
      IMPLICIT NONE
*
      include "../commons/colfact.h"
*
*     Input variables
*
      DOUBLE COMPLEX N, M
*
      CL1QQ_SIDIS = 8D0 * CF / ( M + 1D0 ) / ( N + 1D0 ) 
*
      RETURN
      END
*
************************************************************************
      DOUBLE COMPLEX FUNCTION CL1GQ_SIDIS(N, M)
*
      IMPLICIT NONE
*
      include "../commons/colfact.h"
*
*     Input variables
*
      DOUBLE COMPLEX N, M
*
      CL1GQ_SIDIS = 8D0 * CF / M / ( M + 1D0 ) / ( N + 1D0 ) 
*
      RETURN
      END
*
************************************************************************
      DOUBLE COMPLEX FUNCTION CL1QG_SIDIS(N)
*
      IMPLICIT NONE
*
      include "../commons/colfact.h"
*
*     Input variables
*
      DOUBLE COMPLEX N
*
      CL1QG_SIDIS = 2D0 * TR / ( N + 1D0 ) / ( N + 2D0 ) 
*
      RETURN
      END
*
************************************************************************
      DOUBLE COMPLEX FUNCTION C21QQ_SIDIS(N, M)
*
      IMPLICIT NONE
*
      include "../commons/consts.h"
      include "../commons/colfact.h"
*
*     Input variables
*
      DOUBLE COMPLEX N, M
*
*     Internal variables
*
      DOUBLE COMPLEX N1, M1, NS, MS
      DOUBLE COMPLEX S1N,S1M,S2N,S2M
      DOUBLE COMPLEX PSI,DPSI
      DOUBLE COMPLEX CL1QQ_SIDIS
*
      N1 = N + 1D0
      M1 = M + 1D0
      NS = N * N
      MS = M * M
*
      S1N = EMC + PSI(N1)
      S2N = ZETA2 - DPSI(N1,1)
      S1M = EMC + PSI(M1)
      S2M = ZETA2 - DPSI(M1,1)
*
      C21QQ_SIDIS = 2 * CF * ( - 8D0 - 1D0 / MS + 2D0 / M1**2 + 1D0 / NS
     1     + ( ( 1D0 + M + N )**2 - 1D0 ) / M / M1 / N / N1 + 3D0 * S2M
     2     - S2N + ( S1M + S1N ) * ( S1M + S1N - 1D0 / M / M1
     2     - 1D0 / N / N1 ) + 2D0 / M / M1 / N / N1 )
     3     + CL1QQ_SIDIS(N, M)
*
      RETURN
      END
*
************************************************************************
      DOUBLE COMPLEX FUNCTION C21GQ_SIDIS(N, M)
*
      IMPLICIT NONE
*
      include "../commons/consts.h"
      include "../commons/colfact.h"
*
*     Input variables
*
      DOUBLE COMPLEX N, M
*
*     Internal variables
*
      DOUBLE COMPLEX N1, M1, NS, MS, MM
      DOUBLE COMPLEX S1N, S1M
      DOUBLE COMPLEX PSI
      DOUBLE COMPLEX CL1GQ_SIDIS
*
      N1 = N + 1D0
      M1 = M + 1D0
      MM = M - 1D0
      NS = N * N
      MS = M * M
*
      S1N = EMC + PSI(N1)
      S1M = EMC + PSI(M1)
*
      C21GQ_SIDIS = 2 * CF * ( ( 2D0 - 2d0 * M -9D0 * MS + M**3 - M**4
     1     + M**5 ) / MS / MM**2 / M1**2 + 2D0 * M / N / M1 / MM
     2     - ( 2D0 - M - MS ) / M / M1 / MM / N1
     3     - ( 2D0 + M + MS ) * ( S1N + S1M ) / M / M1 / MM
     4     - 2D0 / M1 / N / N1 + 2D0 / M1 / N / N1 )
     5     + CL1GQ_SIDIS(N, M)
*
      RETURN
      END
*
************************************************************************
      DOUBLE COMPLEX FUNCTION C21QG_SIDIS(N, M)
*
      IMPLICIT NONE
*
      include "../commons/consts.h"
      include "../commons/colfact.h"
*
*     Input variables
*
      DOUBLE COMPLEX N, M
*
*     Internal variables
*
      DOUBLE COMPLEX N1, M1, NM, MM
      DOUBLE COMPLEX S1N, S1M
      DOUBLE COMPLEX PSI
      DOUBLE COMPLEX CL1QG_SIDIS
*
      N1 = N + 1D0
      M1 = M + 1D0
      NM = N - 1D0
      MM = M - 1D0
*
      S1N = EMC + PSI(N1)
      S1M = EMC + PSI(M1)
*
      C21QG_SIDIS = 2 * TR * ( NM * ( 1D0 / MM - 1D0 / M + 1D0 / N
     1     - S1M - S1N ) / N / N1 + ( 2D0 + N + N**2 )
     2     * ( 1D0 / MM - 1D0 / M - S1M - S1N ) / N / N1 / ( N + 2D0 )
     3     + 1D0 / N / N )
     4     + CL1QG_SIDIS(N)
*
      RETURN
      END
