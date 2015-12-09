********************************************************************************
*                                                                              *
*     andim_nnlo_tl.f                                                          *
*                                                                              *
********************************************************************************
*                                                                              *
*     Returns the NNLO time-like anomalous dimensions in the a_s=alpha_s/4pi   *
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
      SUBROUTINE ANDIM_NNLO_TL(N,NF,P2NS,P2SG) 
*
      IMPLICIT NONE
**
*     Input Variables
*
      INTEGER NF
      DOUBLE COMPLEX N
**
*     Output Variables
*
      DOUBLE COMPLEX P2NS(3)
      DOUBLE COMPLEX P2SG(2,2)
*
      CALL PT2NSP(P2NS(1),P2NS(2),P2NS(3),N,NF)
      CALL PT2IJP(P2SG(1,1),P2SG(2,1),P2SG(1,2),P2SG(2,2),N,NF)
*
      P2SG(1,1) = P2SG(1,1) + P2NS(1)
*
      RETURN
      END
*
************************************************************************
*
* ..File: pt2mom.f
*
*
* ..The subroutines  PT2NSP & PT2IJP for the (complex) moments of the
*    parametrized non-singlet and singlet MSbar splitting functions for 
*    the evolution of unpolarized fragmentation densities, mu_r = mu_f.
*
* ..The QCD colour factors have been hard-wired in the parametrizations.
*    The coupling constant is normalized as  a_s = alpha_s/(4*pi). 
*    
* ..The required routines for Euler's psi-function and its derivatives 
*    are appended at the end of the file.
*
* ..References:
*
*   A. Mitov, S. Moch and A. Vogt, PL B638 (2006) 61, hep-ph/0604053;        
*   S. Moch and A. Vogt, PL B659 (2008) 290, arXiv:0708.3899; 
*   A. Almasy, S. Moch and A. Vogt, NP B854 (2012) 133, arXiv:1107.2263
*
* ..for the ns; ps and gg; qg, gq parts and all parametrizations, resp.
*
* =====================================================================
*
*    
       SUBROUTINE PT2NSP ( PT2PLSN, PT2MINN, PT2VALN, N, NF )
*
       IMPLICIT DOUBLE COMPLEX (A - Z)
       DOUBLE PRECISION EMC, ZETA2, ZETA3, ZETA4
       INTEGER NF
       PARAMETER ( EMC  = 0.57721 56649D0, ZETA2 = 1.64493 40668D0,
     ,            ZETA3 = 1.20205 69032D0, ZETA4 = 1.08232 32337D0 )
*
* ..Analytic continuations of simple harmonic sums
*
       HS1(Z) =  EMC  + PSI (Z+1.) 
       HS2(Z) = ZETA2 - DPSI (Z+1.,1)
       HS3(Z) = ZETA3 + 0.5 * DPSI (Z+1.,2)
*
* ---------------------------------------------------------------------
*
* ..Some abbreviations
*
       NI = 1./N
       NI2 = NI*NI
       NI3 = NI*NI2
       NM = N - 1.
       NMI = 1./NM
*
       N1 = N + 1.
       N1I = 1./N1
       N1I2 = N1I*N1I
       N1I3 = N1I*N1I2
       N2 = N + 2.
       N2I = 1./N2
*
       S1 = HS1(N)
       S2 = HS2(N)
       S3 = HS3(N)
       S1M = S1 - NI
       S11 = S1 + N1I
       S12 = S11 + N2I
       S21 = S2 + N1I2
       S31 = S3 + N1I3
*
* ---------------------------------------------------------------------
*
* ..The moments of the functions employed in the parametrizations:
*
* ...1/(1-x)_+ [A0]  and  x^a ln^b (1-x) [B`b(a)' with M for a = -1]
*
       A0  = - S1M
       B1  = - S1 * NI
* ...with special care for the first moment of x^-1 ln(1-x)
       IF ( ( DABS (DIMAG(N)) .LT. 1.D-5 ) .AND. 
     ,      ( DABS ( DBLE(N) - 1.D0 ) .LT. 1.D-5 ) ) THEN
         B1M = - ZETA2
       ELSE
         B1M = - S1M * NMI
       ENDIF
       B11 = - S11 * N1I
       B12 = - S12 * N2I
*
* ...x^a [C`a'] 
*
       C0 = NI
       C1 = N1I
       C2 = N2I
       C3 = 1./(N+3.)
       C4 = 1./(N+4.)
*
* ...x^a ln^b x [D`b(a)'] 
*
       D1  = - NI2
       D11 = - N1I2
       D2  = 2.* NI3
       D3  = - 6.* NI2*NI2
       D31 = - 6.* N1I2*N1I2
       D4  = 24.* NI2*NI3
       D41 = 24.* N1I2*N1I3
*
* ...x^a ln^b x ln(1-x) [E`b(a)']
*
       E1  = S1*NI2 + (S2-ZETA2)*NI
       E11 = S11*N1I2 + (S21-ZETA2)*N1I
       E2  = 2.* ( - S1*NI3 + (ZETA2-S2)*NI2 - (S3-ZETA3)*NI )
       E21 = 2.* ( - S11*N1I3 + (ZETA2-S21)*N1I2 - (S31-ZETA3)*N1I )
*
* ---------------------------------------------------------------------
*
* ..The parametrized n_f^{0,1} components of Pt_ns^(2)i, i = +, -, v
*   (the latter is identical to the spacelike result)
*
       PP2 = + 1174.898 * A0 + 1295.625 - 707.67 * B1 + 593.9 * C3
     ,       - 1075.3 * C2 - 4249.4 * C1 + 1658.7 * C0 + 1327.5 * D1
     ,       - 189.37 * D2 - 352./9.D0 * D3 + 128./81.D0 * D4 
     ,       - 56.907 * E1 - 559.1 * E11 - 519.37 * E2 
     ,     + NF * ( - 183.187 * A0 - 173.935 + 5120/81.D0 * B1
     ,       - 31.84 * C3 + 181.18 * C2 + 466.29 * C1 - 198.10 * C0
     ,       - 168.89 * D1 - 176./81.D0 * D2 + 64./27.D0 * D3
     ,       - 50.758 * E1 + 85.72 * E11 + 28.551 * E2 - 23.102 * E21
     ,       - 39.113 * D11 )
*
       PM2 = + 1174.898 * A0 + 1295.622 - 707.94 * B1 + 407.89 * C3
     ,       - 577.42 * C2 - 4885.7 *C1 + 1981.3 * C0 + 1625.5 * D1
     ,       - 38.298 * D2 - 3072./81.D0 * D3 - 140./81.D0 * D4
     ,       + 4563.2 * E1 - 5140.6 * E11 + 1905.4 * E2 + 1969.5 * E21
     ,       - 437.03 * D31 - 34.683 * D41 
     ,     + NF * ( - 183.187 * A0 - 173.9376 + 5120./81.D0 * B1
     ,       - 85.786 * C3 + 209.19 * C2 + 511.92 * C1 - 217.84 * C0
     ,       - 188.99 * D1 - 784./81.D0 * D2 + 128./81.D0 * D3
     ,       + 71.428 * E1 - 23.722 * E11 + 30.554 * E2 - 18.975 * E21
     ,       + 92.453* D11 )
*
       PS2 = - 163.9 * (B1M-B1) -7.208 * (B11-B12) + 4.82 * (C3-C4)
     ,       - 43.12 * (C2-C3) + 44.51 * (C1-C2) + 151.49 * (C0-C1)
     ,       + 178.04 * D1 + 6.892 * D2 - 40./27.D0 * (2.*D3 - D4)
     ,       - 173.1 * E1 + 46.18 * E2
*
* ..The exact n_f^2 contribution, also identical to the spacelike case 
*   (first derived there by J.A. Gracey in Phys. Lett. B322 (1994) 141)
*
       PF2 = - ( 17./72.D0 - 2./27.D0 * S1 - 10./27.D0 * S2
     ,        + 2./9.D0 * S3 - (12.* N**4 + 2.* N**3 - 12.* N**2
     ,        - 2.* N + 3.)/(27.* N**3 * N1**3) ) * 32./3.D0
*
* ---------------------------------------------------------------------
*
* ..Assemble the pieces for the output
*
       PT2PLSN = PP2 + NF**2 * PF2
       PT2MINN = PM2 + NF**2 * PF2
       PT2VALN = PT2MINN + NF * PS2
*
       RETURN
       END
*
* =====================================================================
*
*    
       SUBROUTINE PT2IJP ( PT2PSN, PT2QGN, PT2GQN, PT2GGN, N, NF )
*
       IMPLICIT DOUBLE COMPLEX (A - Z)
       DOUBLE PRECISION EMC, ZETA2, ZETA3, ZETA4
       INTEGER NF
       PARAMETER ( EMC  = 0.57721 56649D0, ZETA2 = 1.64493 40668D0,
     ,            ZETA3 = 1.20205 69032D0, ZETA4 = 1.08232 32337D0 )
*
* ..Analytic continuations of simple harmonic sums
*
       HS1(Z) =  EMC  + PSI (Z+1.) 
       HS2(Z) = ZETA2 - DPSI (Z+1.,1)
       HS3(Z) = ZETA3 + 0.5 * DPSI (Z+1.,2)
       HS4(Z) = ZETA4 - 1./6.D0 * DPSI (Z+1.,3)
*
* ---------------------------------------------------------------------
*
* ..Some abbreviations
*
       NI = 1./N
       NI2 = NI*NI
       NI3 = NI*NI2
       NM = N - 1.
       NMI = 1./NM
       NMI2 = NMI*NMI
       NMI3 = NMI*NMI2
*
       N1 = N + 1.
       N1I = 1./N1
       N1I2 = N1I*N1I
       N1I3 = N1I*N1I2
       N2 = N + 2.
       N2I = 1./N2
       N2I2 = N2I*N2I
       N2I3 = N2I*N2I2
*
       S1 = HS1(N)
       S2 = HS2(N)
       S3 = HS3(N)
       S4 = HS4(N)
       S1M = S1 - NI
       S11 = S1 + N1I
       S12 = S11 + N2I
       S2M = S2 - NI2
       S21 = S2 + N1I2
       S22 = S21 + N2I2
       S31 = S3 + N1I3
*
* ---------------------------------------------------------------------
*
* ..The moments of the functions employed in the parametrizations:
*
* ...1/(1-x)_+ [A0]  and  x^a ln^b (1-x) [B`b(a)' with M for a = -1]
*
       A0  = - S1M
       B1  = - S1 * NI
       B1M = - S1M * NMI
       B11 = - S11 * N1I
       B12 = - S12 * N2I
       B2  = (S1**2 + S2) * NI
       B2M = (S1M**2 + S2M) * NMI
       B21 = (S11**2 + S21) * N1I
       B22 = (S12**2 + S22) * N2I
       B3  = - (S1**3 + 3.*S1*S2 + 2.*S3) * NI
       B31 = - (S11**3 + 3.*S11*S21 + 2.*S31) * N1I
       B4  = (S1**4 + 6.*S1**2*S2 + 8.*S1*S3 + 3.*S2**2 + 6.*S4) * NI
*
* ...x^a [C`a'] 
*
       C0 = NI
       CM = NMI
       C1 = N1I
       C2 = N2I
       C3 = 1./(N+3.)
       C4 = 1./(N+4.)
       C5 = 1./(N+5.)
*
* ...x^a ln^b x [D`b(a)'] 
*
       D1  = - NI2
       D1M = - NMI2
       D11 = - N1I2
       D12 = - N2I2
       D2  = 2.* NI3
       D2M = 2.* NMI3
       D21 = 2.* N1I3
       D22 = 2.* N2I3
       D3  = - 6.* NI2*NI2
       D3M = - 6.* NMI2*NMI2
       D31 = - 6.* N1I2*N1I2
       D32 = - 6.* N2I2*N2I2
       D4  = 24.* NI2*NI3
       D4M = 24.* NMI2*NMI3
       D41 = 24.* N1I2*N1I3
*
* ...x^a ln^b x ln(1-x) [E`b(a)'] and ln x ln^2(1-x) [F1]
*
       E1  = S1*NI2 + (S2-ZETA2)*NI
       E11 = S11*N1I2 + (S21-ZETA2)*N1I
       E12 = S12*N2I2 + (S22-ZETA2)*N2I
       E2  = 2.* ( - S1*NI3 + (ZETA2-S2)*NI2 - (S3-ZETA3)*NI )
       F1  = 2.*NI* ( ZETA3+ZETA2*S1 - 0.5*NI* (S1*S1+S2)- S1*S2-S3 )
*
* ---------------------------------------------------------------------
*
* ..The parametrized n_f-components of Pt_ps^(2) and Pt_qg^(2)
*   [ Pt_qq^(2) is obtained by adding the non-singlet quantity 
*     Pt_ns^(2) provided by the subroutine  PT2NSP  to Pt_ps^(2) ]
*
       PS1 = - 256./9.D0 * (D3M-D3) - 128./9.D0 * (D2M-D2) 
     ,       + 324.07 * (D1M-D1) + 479.87 * (CM-C0) + 9.072 * (D4-D41) 
     ,       + 47.322 * (D3-D31) + 425.14 * (D2-D21) + 656.49 * (D1-D11)
     ,       - 5.926 * (B3-B31) - 9.751 * (B2-B21) - 8.650 * (B1-B11) 
     ,       - 106.65 * (C0-C1) - 848.97 * (C1-C2) + 368.79 * (C2-C3) 
     ,       - 61.284 * (C3-C4) + 96.171 * (E1-E11)
       PS2 = - 128./81.D0 * (CM-C0) + 0.019122 * (D4-D41) 
     ,       - 1.900 * (D3-D31) + 9.1682 * (D2-D21) + 57.713 * (D1-D11) 
     ,       + 1.778 * (B2-B21) + 16.611 * (B1-B11) + 87.795 * (C0-C1) 
     ,       - 57.688 * (C1-C2) - 41.827 * (C2-C3) + 25.628 * (C3-C4) 
     ,       - 7.9934 * (C4-C5) - 2.1031 * (E1-E11) + 26.294 * (D11-D12)
     ,       - 7.8645 * (D31-D32)
*
       QG1 = - 64. * (D3M + D2M) + 675.83 * D1M + 1141.7 * CM 
     ,       + 42.328 * D4 + 361.28 * D3 + 1512.* D2 + 1864. * D1 
     ,       + 100./27. * B4 + 350./9.* B3 + 263.07 * B2 + 693.84 * B1
     ,       + 603.71 * C0 - 882.48 * C1 +  4723.2 * C2 - 4745.8 * C3
     ,       - 175.28 * C4 - 1809.4 * E1 - 107.59 * E11 - 885.5 * D41
       QG2 = - 32./27.* D2M - 3.1752 * D1M - 2.8986 * CM + 21.569 * D3
     ,       + 255.62 * D2 + 619.75 * D1 - 100./27.* B3 - 35.446 * B2
     ,       - 103.609 * B1 - 113.81 * C0 + 341.26 * C1 - 853.35 * C2 
     ,       + 492.1 * C3 + 14.803 * C4 + 966.96 * E1 - 709.1 * E11 
     ,       - 1.593 * F1 - 333.8 * D31
       QG3 = ( 4.* C0 + 6.*(D1 + B1) + 3.8696 * (C0 - 2.*C1 + 2.*C2)
     ,       + 4.* (D1 - 2.*D11 + 2.*D12 + B1 - 2.*B11 + 2.*B12)
     ,       + 3.* (D2 - 2.*D21 + 2.*D22 + B2 - 2.*B21 + 2.*B22)
     ,       + 6.* (E1 - 2.*E11 + 2.*E12) ) * 4./9.D0
*
* ...and of Pt^(2)_gq and Pt^(2)_gg  
*
       GQ0 =   256. * D4M + 3712./3.D0 * D3M + 1001.89 * D2M 
     ,       + 4776.5 * D1M + 5803.7 * CM - 30.062 * D4 - 126.38 * D3
     ,       - 0.71252 * D2 + 4.4136 * D1 + 400./81.D0 * B4 
     ,       + 520./27.D0 * B3 - 220.13 * B2 - 152.6 * B1 
     ,       + 272.85 * C0 - 7188.7 * C1 + 5693.2 * C2 + 146.98 * C3 
     ,       + 128.19 * C4 - 1300.6 * E1 - 71.23 * F1 + 543.8 * D31 
       GQ1 =   1280./81.D0 * D3M + 2912./27.D0 * D2M + 141.93 * D1M   
     ,       + 6.0041 * CM - 48.60 * D3 - 343.1 * D2 - 492.0 * D1
     ,       + 80./81.D0 * B3 + 1040./81.D0 * B2 - 16.914 * B1
     ,       - 871.3 * C0 + 790.13 * C1 - 241.23 * C2 + 43.252 * C3
     ,       - 4.3465 * D31 + 55.048 * E1
*
       GG0 =   576. * D4M + 3168. * D3M + 3651.1 * D2M + 10233. * D1M
     ,       + 14214.4 * CM + 191.99 * D4 + 3281.7 * D3 + 13528. * D2 
     ,       + 12258. * D1 - 28489. * C0 + 7469. * C1 + 30421. * C2 
     ,       - 53017. * C3 + 19556. * C4 - 186.4 * E1 - 21328. * E2 
     ,       + 5685.8 * D31 - 3590.1 * B1 + 4425.451 + 2643.521 * A0 
       GG1 =   448./9.D0 * D3M + 2368./9.D0 * D2M - 5.470 * D1M 
     ,       - 804.13 * CM + 18.085 * D4 + 155.10 * D3 + 482.94 * D2
     ,       + 4.9934 * D1 + 248.95 * C0 + 260.6 * C1 + 272.79 * C2
     ,       + 2133.2 * C3 - 926.87 * C4 + 1266.5 * E1 - 29.709 * E2 
     ,       + 87.771 * F1 + 485.18 * D31 + 319.97 * B1 - 528.719 
     ,       - 412.172 * A0
       GG2 =   32./27.D0 * D2M + 368./81.D0 * D1M + 472./243.D0 * CM
     ,       - 5.0372 * D3 - 44.80 * D2 - 69.712 * D1 - 77.190 * C0 
     ,       + 153.27 * C1 - 106.03 * C2 + 11.995 * C3 - 115.01 * E1 
     ,       + 96.522 * E11 - 62.908 * E2 + 6.4628 - 16./9.D0 * A0
*
* ---------------------------------------------------------------------
*
* ..Assemble the pieces for the output
*
       PT2PSN =       NF * ( PS1 + NF * PS2 )
c       PT2QGN =       NF * ( QG1 + NF * (QG2 + NF * QG3) )
c       PT2GQN = GQ0 + NF *   GQ1 
       PT2QGN = ( QG1 + NF * (QG2 + NF * QG3) ) / 2D0
       PT2GQN = 2D0 * NF * ( GQ0 + NF * GQ1 )
       PT2GGN = GG0 + NF * ( GG1 + NF * GG2 )
*
       RETURN
       END
*
* =============================================================aaa/av==
