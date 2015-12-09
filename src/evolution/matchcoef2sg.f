*
* ..File: matchcoef2sg.f     
*
*
* ..The subroutine  MATCHCOEF2SG  for the NNLO (alpha_s^2) heavy quark 
*    contribution  A2SG  to the singlet operator matrix element 
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
      SUBROUTINE MATCHCOEF2SG(N,A2SG)
      IMPLICIT NONE 
*
      include "../commons/colfact.h"
      include "../commons/consts.h"
      include "../commons/hqmass.h"
*
      DOUBLE COMPLEX DPSI,PSI,S1,S2,S3
*
* ---------------------------------------------------------------------
*
*     Internal variables
*
      DOUBLE PRECISION H1
      DOUBLE COMPLEX NM,N1,N2,NI,NMI,N1I,N2I
      DOUBLE COMPLEX S1M,S2M,S3M,S11,S21,S31,S22
      DOUBLE COMPLEX A0
      DOUBLE COMPLEX B1,B1M,B11,B2,B2M,B21,B3
      DOUBLE COMPLEX C0,CM,C1,C2
      DOUBLE COMPLEX D1,D11,D12,D2,D21,D22,d3,D31
      DOUBLE COMPLEX E2
      DOUBLE COMPLEX F1,F1M,F11,F12,F2,F21
      DOUBLE COMPLEX A2HQ,A2HG,A2GQ,A2GGF,A2GGA
      DOUBLE COMPLEX FN1,FN2
*
* ---------------------------------------------------------------------
*
*     Input variables
*
      DOUBLE COMPLEX N    
*
*     Output variables  
*
      DOUBLE COMPLEX A2SG(2,2)
*
* ---------------------------------------------------------------------
*     
*     Some useful definition

      S1 = EMC + PSI(N+1.d0)
      S2 = ZETA2 - DPSI(N+1.d0,1)
      S3 = ZETA3 + (1.d0/2.d0)*DPSI(N+1.d0,2)
*
      NM = N - 1.D0
      N1 = N + 1.D0
      N2 = N + 2.D0
      NI = 1.D0/N
      NMI = 1.D0/NM
      N1I = 1.D0/N1
      N2I = 1.D0/N2
*
      S1M = S1 - NI
      S2M = S2 - NI*NI
      S3M = S3 - NI**3.
      S11 = S1 + N1I
      S21 = S2 + N1I*N1I
      S31 = S3 + N1I**3.
      S22 = S21 + N2I*N2I
*
* ---------------------------------------------------------------------
*
*  ..Moments of the basic x-space functions 
*
      A0  = - S1M
*
      B1  = - S1 * NI
      B1M = - S1M * NMI
      B11 = - S11 * N1I
      B2  = (S1**2. + S2) * NI
      B2M = (S1M**2. + S2M) * NMI
      B21 = (S11**2. + S21) * N1I
      B3  = - (S1**3. + 3.D0*S1*S2 + 2.D0*S3) * NI
*
      C0 = NI
      CM = NMI
      C1 = N1I
      C2 = N2I
*
      D1  = - NI*NI
      D11 = - N1I*N1I
      D12 = - N2I*N2I
      D2  = 2.D0* NI**3.
      D21 = 2.D0* N1I**3.
      D22 = 2.D0* N2I**3.
      D3  = - 6.D0* NI**4.
      D31 = - 6.D0* N1I**4.
*
      E2 = 2.D0* NI * ( ZETA3 - S3 + NI * (ZETA2 - S2 - NI * S1) )
*
      F1  = NI  * ( ZETA2 - S2 )
      F1M = NMI * ( ZETA2 - S2M )
      F11 = N1I * ( ZETA2 - S21 )
      F12 = N2I * ( ZETA2 - S22 )
      F2  = - NI  * F1
      F21 = - N1I * F11
*
* ---------------------------------------------------------------------
*             
* ..The moments of the OME's A_Hq^{PS,(2)} and A_Hg^{S,(2)} given in 
*    Eqs. (B.1) and (B.3) of BMSN. For the latter quantity an accurate
*    x-space parametrization is used instead of the full expression.
*
      A2HQ = - (32.D0/3.D0*CM + 8.D0*(C0-C1) - 32.D0/3.D0 * C2) * ZETA2
     1       - 448.D0/27.D0 * CM - 4.D0/3.D0 * C0 - 124.D0/3.D0 * C1 
     2       + 1600.D0/27.D0 * C2 - 4.D0/3.D0 * (D3 + D31) + 2.D0* D2 
     3       + 10.D0* D21 + 16.D0/3.D0* D22 - 16.D0* ZETA2 * (D1 + D11) 
     4       - 56.D0/3.D0 * D1 - 88.D0/3.D0 * D11 - 448.D0/9.D0 * D12 
     5       + 32.D0/3.D0 * F1M + 8.D0* (F1 - F11) - 32.D0/3.D0 * F12 
     6       + 16.D0* (F2 + F21)
*
      A2HG = - 0.006D0 - 1.111D0 * B3 - 0.400D0 * B2 - 2.770D0 * B1
     1       - 24.89D0 * CM - 187.8D0 * C0 + 249.6D0 * C1
     2       - 1.556D0 * D3 - 3.292D0 * D2 - 93.68D0 * D1 - 146.8D0* E2
*
* ..The moments of the OME's A_{gq,H}^{S,(2)} and A_{gg,H}^{S,(2)} 
*    given in Eqs. (B.5) and (B.7) of BMSN.
*
      A2GQ =   4.D0/3.D0 * (2.D0* B2M - 2.D0* B2 + B21)
     1       + 8.D0/9.D0 * (10.D0* B1M - 10.D0* B1 + 8.D0* B11)
     2       + 1.D0/27.D0 * (448.D0* (CM - C0) + 344.D0* C1)  
*
      A2GGF = - 15.D0 - 8.D0* CM + 80.D0* C0 - 48.D0* C1 - 24.D0* C2 
     1        + 4.D0/3.D0 * (D3 + D31) + 6.D0*D2 + 10.D0*D21 + 32.D0*D1 
     2        + 48.D0* D11
      A2GGA =   224.D0/27.D0 * A0 + 10.D0/9.D0 - 4.D0/3.D0 * B11 
     1        + 1.D0/27.D0 * (556.D0* CM - 628.D0* C0 
     2        + 548.D0* C1 - 700.D0* C2) + 4.D0/3.D0 * (D2 + D21)
     3        + 1.D0/9.D0 * (52.D0* D1 + 88.D0* D11)
*
* ---------------------------------------------------------------------
*
* ..Output to the array 
*
      A2SG(1,1) = CF * TR * A2HQ
      A2SG(1,2) = A2HG
      A2SG(2,1) = CF * TR * A2GQ
      A2SG(2,2) = TR * ( CF * A2GGF + CA * A2GGA )
*
      IF(HQMASS.EQ.1)THEN
         H1 = 4D0 * CF
         FN1 = 4D0 * TR * ( 2D0 / N2 - 2D0 / N1 + 1D0 / N )
         FN2 = - 4D0 * TR / 3D0
*
         A2SG(1,2) = A2SG(1,2) - 2D0 * H1 * FN1
         A2SG(2,2) = A2SG(2,2) - 2D0 * H1 * FN2
      ENDIF
*
* ---------------------------------------------------------------------
*
*
       RETURN
       END
*
* =================================================================av==
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
