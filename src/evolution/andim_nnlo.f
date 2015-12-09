********************************************************************************
*                                                                              *
*     andim_nnlo.f                                                             *
*                                                                              *
********************************************************************************
*                                                                              *
*     Returns the NNLO space-like anomalous dimensions in the a_s=alpha_s/4pi  *
*     notation                                                                 *
*                                                                              *
*     Input variables:                                                         *
*     N  : Mellin variable                                                     *
*     NF : Number of flavours                                                  *
*                                                                              *
*     Output:                                                                  *
*     P2NS(3) : Complex valued non-singlet NNLO anomalous dimension in N space *
*     P2SG(2,2) : matrix of the complex-valued singlet                         *
*                 NNLO anomalous dimension in N space                          *
*                                                                              *
********************************************************************************
      SUBROUTINE ANDIM_NNLO(N,NF,P2NS,P2SG) 
      IMPLICIT NONE
*
      include "../commons/colfact.h"
      include "../commons/consts.h"
**
*     Input variables
*
      DOUBLE COMPLEX N
      INTEGER NF       
**
*     Internal variables
*
      DOUBLE COMPLEX NI,NI2,NI3,NM,NMI,NMI2,N1,N1I,N1I2,N1I3,N2,N2I
      DOUBLE COMPLEX S1M,S11,S2M,S21,S31,S12
      DOUBLE COMPLEX A0,B1,B1MNS,B1MSG,B11,B2,B2M,B21,B3,B31,B4,B12
      DOUBLE COMPLEX C0,CM,C1,C2,C3,C4
      DOUBLE COMPLEX D1,D1M,D11,D2,D21,D3,D31,D4,D41
      DOUBLE COMPLEX E1,E11,E2
      DOUBLE COMPLEX QG1,QG2,GQ0,GQ1,GQ2,GG0,GG1,GG2
      DOUBLE COMPLEX PP20,PP21,PM20,PM21,PS1SG,PS2NS,PS2SG,PF2
      DOUBLE COMPLEX DPSI,PSI,S1,S2,S3,S4
**
*     Output variables  
*
      DOUBLE COMPLEX P2NS(3)
      DOUBLE COMPLEX P2SG(2,2)
**     
*     Some useful definitions
*
      S1 = EMC + PSI(N+1.d0)
      S2 = ZETA2 - DPSI(N+1.d0,1)
      S3 = ZETA3 + (1.d0/2.d0)*DPSI(N+1.d0,2)
      S4 = ZETA4 - (1.d0/6.d0)*DPSI(N+1.d0,3)
*
      NI = 1.D0/N
      NI2 = NI*NI
      NI3 = NI*NI2
      NM = N - 1.D0
      NMI = 1./NM
      NMI2 = NMI*NMI
*
      N1 = N + 1.D0
      N1I = 1.D0/N1
      N1I2 = N1I*N1I
      N1I3 = N1I*N1I2
      N2 = N + 2.D0
      N2I = 1.D0/N2
*
      S1M = S1 - NI
      S11 = S1 + N1I
      S12 = S11 + N2I
      S2M = S2 - NI2
      S21 = S2 + N1I2
      S31 = S3 + N1I3
**
* ..The moments of the functions employed in the parametrizations:
*
* ...1/(1-x)_+ [A0]  and  x^a ln^b (1-x) [B`b(a)' with M for a = -1]
*
      A0  = - S1M
      B1  = - S1 * NI
*
* ...with special care for the first moment of x^-1 ln(1-x)
*
      IF ( ( DABS (DIMAG(N)) .LT. 1.D-5 ) .AND.
     ,      ( DABS ( DBLE(N) - 1.D0 ) .LT. 1.D-5 ) ) THEN
        B1MNS = - ZETA2
      ELSE
        B1MNS = - S1M * NMI
      ENDIF
      B1MSG = - S1M * NMI
      B11 = - S11 * N1I
      B12 = - S12 * N2I
      B2  = (S1**2. + S2) * NI
      B2M = (S1M**2. + S2M) * NMI
      B21 = (S11**2. + S21) * N1I
      B3  = - (S1**3. + 3.D0*S1*S2 + 2.D0*S3) * NI
      B31 = - (S11**3. + 3.D0*S11*S21 + 2.D0*S31) * N1I
      B4  = (S1**4. + 6.D0*S1**2.*S2 + 8.D0*S1*S3 + 3.D0*S2**2. + 
     1       6.D0*S4) * NI
*
* ...x^a [C`a']
*
      C0 = NI
      C1 = N1I
      CM = NMI
      C2 = N2I
      C3 = 1.D0/(N+3.D0)
      C4 = 1.D0/(N+4.D0)
*
* ...x^a ln^b x [D`b(a)']
*
      D1  = - NI2
      D1M = - NMI2
      D11 = - N1I2
      D2  = 2.D0* NI3
      D21 = 2.D0* N1I3
      D3  = - 6.D0* NI2*NI2
      D31 = - 6.D0* N1I2*N1I2
      D4  = 24.D0* NI2*NI3
      D41 = 24.D0* N1I2*N1I3
*
* ...x^a ln^b x ln(1-x) [E`b(a)']
*
      E1  = S1*NI2 + (S2-ZETA2)*NI
      E11 = S11*N1I2 + (S21-ZETA2)*N1I
      E2  = 2.D0* ( - S1*NI3 + (ZETA2-S2)*NI2 - (S3-ZETA3)*NI )
*
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                            NON-SINGLET
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*
* ..The parametrized n_f^{0,1} components of P_ns^(2)i, i = +, -, v
*    as given in S. Moch, J. Vermaseren and A. Vogt, hep-ph/0403192
*
      PP20 = + 1174.898D0*A0 + 1295.384D0 + 714.1D0*B1 - 522.1D0*C3
     1       + 243.6D0*C2 - 3135.D0*C1 + 1641.1D0*C0 + 1258.D0*D1
     2       + 294.9D0 * D2 + 800.D0/27.D0 * D3 + 128.D0/81.D0 * D4
     3       + 563.9D0 * E1 + 256.8D0 * E2
      PP21 = - 183.187D0 * A0 - 173.927D0 - 5120.D0/81.D0 * B1
     1       + 44.79D0 * C3 + 72.94D0 * C2 + 381.1D0 * C1 - 197.D0 * C0
     2       - 152.6D0 * D1 - 2608.D0/81.D0 * D2 - 192.D0/81.D0 * D3
     3       - 56.66D0 * E1 - 1.497D0 * D31 
*
      PM20 = + 1174.898D0 * A0 + 1295.470D0 + 714.1D0 * B1 - 433.2D0*C3
     1       + 297.D0 * C2 - 3505.D0*C1 + 1860.2D0 * C0 + 1465.2D0 * D1
     2       + 399.2D0 * D2 + 320.D0/9.D0 * D3 + 116.D0/81.D0 * D4
     3       + 684.D0 * E1 + 251.2D0 * E2
      PM21 = - 183.187D0 * A0 - 173.933D0 - 5120.D0/81.D0 * B1
     1       + 34.76D0 * C3 + 77.89D0 * C2 + 406.5D0 * C1 - 216.62D0*C0
     2       - 172.69D0 * D1 - 3216.D0/81.D0 * D2 - 256.D0/81.D0 * D3
     3       - 65.43D0 * E1 - 1.136D0 * D31
*
      PS2NS  = - 163.9D0 * (B1MNS-B1)-7.208D0*(B11-B12) + 4.82D0*(C3-C4)
     1       - 43.12D0 * (C2-C3) + 44.51D0 * (C1-C2) + 151.49D0*(C0-C1)
     2       + 178.04D0 * D1 + 6.892D0*D2 - 40.D0/27.D0*(2.D0*D3 - D4)
     2       - 173.1D0 * E1 + 46.18D0 * E2
*
* ..The exact n_f^2 contribution first determined by J.A. Gracey in
*    hep-ph/9401214
*
      PF2  = - ( 17.D0/72.D0 - 2.D0/27.D0 * S1 - 10.D0/27.D0 * S2
     1       + 2.D0/9.D0*S3 - (12.D0* N**4. + 2.D0*N**3. - 12.D0*N**2.
     2       - 2.D0* N + 3.D0)/(27.D0* N**3. * N1**3.) ) * 32.D0/3.D0
**
* ..Flavour-number loop and output to the array 
*
      P2NS(1) = PP20 + NF * (PP21 + NF * PF2)
      P2NS(2) = PM20 + NF * (PM21 + NF * PF2)
      P2NS(3) = P2NS(2) + NF * PS2NS 
*
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                                 SINGLET
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*
* ..The parametrized n_f-components of P_ps^(2) and P_qg^(2) ...
*   [ P_qq^(2) is obtained below by adding the non-singlet quantity
*     P_ns^(2)+ provided by the subroutine  P2NS2MOM  to P_ps^(2) ]
*
      PS1SG = 
     ,- 3584./27.D0* (D1M-D1) - 506.D0* (CM-C0) + 160./27.D0* (D4-D41)
     ,- 400.D0/9.D0 * (D3-D31) + 131.4D0 * (D2-D21) - 661.6D0 * (D1-D11)
     ,- 5.926D0 * (B3-B31) - 9.751D0 * (B2-B21) - 72.11D0 * (B1-B11)
     ,+ 177.4D0 * (C0-C1) + 392.9D0 * (C1-C2) - 101.4D0 * (C2-C3)
     ,- 57.04D0 * (E1-E11)
      PS2SG = 
     ,256.D0/81.D0* (CM-C0) + 32.D0/27.D0* (D3-D31) + 17.89D0 * (D2-D21)
     ,+ 61.75D0 * (D1-D11) + 1.778D0 * (B2-B21) + 5.944D0 * (B1-B11)
     ,+ 100.1D0 * (C0-C1) - 125.2D0 * (C1-C2) + 49.26D0 * (C2-C3)
     ,- 12.59D0 * (C3-C4) - 1.889D0 * (E1-E11)
*
      QG1 =
     ,- 896./3.D0* D1M - 1268.3D0 * CM + 536./27.D0* D4 - 44./3.D0* D3
     ,+ 881.5D0 * D2 + 424.9D0 * D1 + 100.D0/27.D0* B4 - 70.D0/9.D0* B3
     ,- 120.5D0 * B2 + 104.42D0 * B1 + 2522.D0* C0 - 3316.D0 * C1
     ,+ 2126.D0 * C2 + 1823.D0 * E1 - 25.22D0 * E2 - 252.5D0 * D31
      QG2 =
     ,1112.D0/243.D0* CM - 16.D0/9.D0* D4 - 376./27.D0* D3 - 90.8D0 * D2
     ,- 254.D0 * D1 + 20.D0/27.D0* B3 + 200.D0/27.D0* B2 - 5.496D0 * B1
     ,- 252.D0 * C0 + 158.D0 * C1 + 145.4D0 * C2 - 139.28D0 * C3
     ,- 53.09D0 * E1 - 80.616D0 * E2 - 98.07D0 * D21 + 11.7D0 * D31
*
* ...and of P^(2)_gq and P^(2)_gg  [GQ2 is exact], all as given by 
*     A. Vogt, S. Moch and J. Vermaseren in hep-ph/0404111
*
      GQ0 =
     ,1189.3D0 * D1M + 6163.1D0 * CM - 4288./81.D0 * D4 + 1568./9.D0* D3
     ,- 1794.D0* D2 + 4033.D0* D1 + 400.D0/81.D0* B4 + 2200.D0/27.D0* B3
     ,+ 606.3D0 * B2 + 2193.D0* B1 - 4307.D0* C0 + 489.3D0 * C1
     ,+ 1452.D0* C2 + 146.D0* C3 - 447.3D0 * E2 - 972.9D0 * D21
      GQ1 =
     ,71.082D0 * D1M - 46.41D0 * CM + 128.D0/27.D0* D4 + 704./81.D0* D3
     ,+ 20.39D0 * D2 + 174.8D0 * D1 - 400.D0/81.D0* B3 - 68.069D0 * B2
     ,- 296.7D0 * B1 - 183.8D0 * C0 + 33.35D0 * C1 - 277.9D0 * C2
     ,+ 108.6D0 * D21 - 49.68D0 * E1
      GQ2 =
     ,(64.D0* (- CM + C0 + 2.D0* C1)+ 320.D0* (B1MSG - B1 + 0.8D0 * B11)
     ,+ 96.D0* (B2M - B2 + 0.5D0 * B21) ) / 27.D0
*
      GG0 =
     ,2675.8D0 * D1M + 14214.D0* CM - 144.D0 * D4 + 72.D0 * D3
     ,- 7471.D0* D2 + 274.4D0 * D1 - 20852.D0* C0 + 3968.D0* C1
     ,- 3363.D0* C2 + 4848.D0* C3 + 7305.D0* E1 + 8757.D0* E2
     ,+ 3589.D0* B1 + 4425.894D0 + 2643.521D0 * A0
      GG1 =
     ,157.27D0 * D1M + 182.96D0 * CM + 512.D0/27.D0 * D4
     ,+ 832.D0/9.D0 * D3 + 491.3D0 * D2 + 1541.D0* D1 - 350.2D0 * C0
     ,+ 755.7D0 * C1 - 713.8D0 * C2 + 559.3D0 * C3 + 26.15D0 * E1
     ,- 808.7D0 * E2 - 320.D0 * B1 - 528.723D0 - 412.172D0 * A0
      GG2 =
     ,- 680.D0/243.D0 * CM - 32.D0/27.D0 * D3 + 9.68D0 * D2
     ,- 3.422D0 * D1 - 13.878D0 * C0 + 153.4D0 * C1 - 187.7D0 * C2
     ,+ 52.75D0 * C3 - 115.6D0 * E1 + 85.25D0 * E11 - 63.23D0 * E2
     ,+ 6.463D0 - 16.D0/9.D0 * A0
*
* ..Flavour-number loop and output to the array 
*
      P2SG(1,1) = P2NS(1) + NF * ( PS1SG + NF * PS2SG )
      P2SG(1,2) =       NF * ( QG1 + NF * QG2 )
      P2SG(2,1) = GQ0 + NF * ( GQ1 + NF * GQ2 )
      P2SG(2,2) = GG0 + NF * ( GG1 + NF * GG2 )
*
       RETURN
       END
*
* =================================================================av===
