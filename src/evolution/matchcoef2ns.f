*
* ..File: matchcoef2ns.f     
*
*
* ..The subroutine  MATCHCOEF2NS  for the NNLO (alpha_s^2) heavy quark 
*    contribution  A2NS  to the non-singlet operator matrix element 
*    (OME) in N-space in the MS(bar) scheme for mu_f^2 = m_H^2.
*    The coupling constant is normalized as  a_s = alpha_s/(4*pi).
*
* ..This quantity, presented in Appendix B of Buza, Matiounine, Smith 
*    and van Neerven, Eur. Phys. J. C1 (1998) 301 (BSMN), is required 
*    for the N_f matching of the NNLO parton densities.
*
* =====================================================================
*
*
      SUBROUTINE MATCHCOEF2NS(N,A2NS)
      IMPLICIT NONE 
*
      include "../commons/colfact.h"
      include "../commons/consts.h"
*
      DOUBLE COMPLEX DPSI,PSI,S1,S2,S3
*
* ---------------------------------------------------------------------
*
*     Internal variables
*
      DOUBLE COMPLEX NI,N1,N1I
      DOUBLE COMPLEX S1M,S2M,S3M,S21,S31
      DOUBLE COMPLEX A0
      DOUBLE COMPLEX C0,C1
      DOUBLE COMPLEX D1,D11
      DOUBLE COMPLEX G1,G12,G2,G22
      DOUBLE COMPLEX A2QQ
*
* ---------------------------------------------------------------------
*
*     Input variables
*
      DOUBLE COMPLEX N      
*
*     Output variables  
*
      DOUBLE COMPLEX A2NS
*
* ---------------------------------------------------------------------
*     
*     Some useful definition

      S1 = EMC + PSI(N+1.d0)
      S2 = ZETA2 - DPSI(N+1.d0,1)
      S3 = ZETA3 + (1.d0/2.d0)*DPSI(N+1.d0,2)
*
      N1 = N + 1.D0
      NI = 1.D0/N
      N1I = 1.D0/N1
*
      S1M = S1 - NI
      S2M = S2 - NI*NI
      S3M = S3 - NI**3.
      S21 = S2 + N1I*N1I
      S31 = S3 + N1I**3.
*
* ---------------------------------------------------------------------
*
*  ..Moments of the basic x-space functions 
*
      A0 = - S1M
*
      C0 = NI
      C1 = N1I
*
      D1  = - NI*NI
      D11 = - N1I*N1I
*
      G1  = S2M - ZETA2
      G12 = S21 - ZETA2
      G2  = - 2.D0* ( S3M - ZETA3 )  
      G22 = - 2.D0* ( S31 - ZETA3 )  
*
* ---------------------------------------------------------------------
*
* ..The moments of the OME A_{qq,H}^{NS,(2)} given in Eq. (B.4) of BMSN 
*
      A2QQ = 224.D0/27.D0 * A0 - 8.D0/3.D0 * ZETA3 + 40.D0/9.D0 * ZETA2 
     1       + 73.D0/18.D0 + 44.D0/27.D0 * C0 - 268.D0/27.D0 * C1 
     2       + 8.D0/3.D0 * (D1 - D11) + 20.D0/9.D0 * (G1 + G12) 
     3       + 2.D0/3.D0 * (G2 + G22)
*
* ..Output to the array 
*
      A2NS = CF * TR * A2QQ
*
* ---------------------------------------------------------------------
*
      RETURN
      END
*
* =================================================================av==
*
*     The following function returns the Mellin transform of the 
*     complete matching condition reported in eq. (B.4) of hep-ph/9612398 
*     including the dependence on the scale.
*
*     This is needed for the implementation of FONLL.
*
* =====================================================================
*
      FUNCTION A2NSL(N,M2,Q2)
*
      IMPLICIT NONE
*
      include "../commons/colfact.h"
      include "../commons/consts.h"
**
*     Input Variables
*
      DOUBLE PRECISION M2,Q2
      DOUBLE COMPLEX N
**
*     Internal Variables
*
      DOUBLE PRECISION L
      DOUBLE COMPLEX N1,NS,N1S
      DOUBLE COMPLEX PSI,DPSI
      DOUBLE COMPLEX S1,S2
      DOUBLE COMPLEX CL2,CL
      DOUBLE COMPLEX A2NS
**
*     Output Variables
*
      DOUBLE COMPLEX A2NSL
*
      L = DLOG(Q2/M2)
*
      N1  = N + 1D0
      NS  = N * N
      N1S = N1 * N1
*
      S1 = EMC + PSI(N+1D0)
      S2 = ZETA2 - DPSI(N+1D0,1)
*
      CL2 = 8D0 / 3D0 * ( - S1 + 1D0 / N ) - 4D0 / 3D0 / N 
     1    - 4D0 / 3D0 / N1 + 2D0
*
      CL  = 80D0 / 9D0 * ( - S1 + 1D0 / N ) + 8D0 / 3D0 
     1    * ( 2D0 * ( S2 - ZETA2 ) - 1D0 / NS + 1D0 / N1S )
     2    + 8D0 / 9D0 / N - 88D0 / 9D0 / N1 + 16D0 * ZETA2 / 3D0
     3    + 2D0 / 3D0 
*
      CALL MATCHCOEF2NS(N,A2NS)
*
      A2NSL = CF * TR * ( CL2 * L**2D0 - CL * L ) + A2NS
*
      RETURN
      END
