********************************************************************************
*                                                                              *
*     andim_nlo_tl.f                                                           *
*                                                                              *
********************************************************************************
*                                                                              *
*     Returns the NLO time-like anomalous dimensions in the a_s=alpha_s/4pi    *
*     notation                                                                 *
*                                                                              *
*     Input variables:                                                         *
*     N  : Mellin variable                                                     *
*     NF : Number of flavours                                                  *
*                                                                              *
*     Output:                                                                  *
*     P1NS(3) : Complex valued non-singlet NLO anomalous dimension in N space  *
*     P1SG(2,2) : matrix of the complex-valued singlet                         *
*                 NLO anomalous dimension in N space                           *
*                                                                              *
********************************************************************************
      SUBROUTINE ANDIM_NLO_TL(N,NF,P1NS,P1SG)
*
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
      INTEGER I,J
      DOUBLE COMPLEX NS,NT,NFO,NFI,NSI,NSE,NE,NN
      DOUBLE COMPLEX N1,N2,NM,NM2,NMS,NMT,N1S,N1T,N2S,N2T
      DOUBLE COMPLEX N3,N4,N5,N6
      DOUBLE COMPLEX S1M,S11,S12,S13,S14,S15,S16,S2M,S21,S22
      DOUBLE COMPLEX SPMOM,SLC,SLV,SSCHLM,SSTR2M,SSTR3M,SSCHLP
      DOUBLE COMPLEX SSTR2P,SSTR3P
      DOUBLE COMPLEX DS2NM,DS2N,DS2N1,DS2N2
      DOUBLE COMPLEX PPSA,PQGA,PGQA,PGGA,PQGB,PGQB,PGGB,PQGC,PGGC,PPSTL
      DOUBLE COMPLEX PNPA,PNMA,PNSB,PNSC
      DOUBLE COMPLEX PQQATL,PQQBTL,PGGATL,PGGBTL,PGGCTL,PNSTL
      DOUBLE COMPLEX DPSI,PSI,S1,S2
**
*     Output variables
*
      DOUBLE COMPLEX P1NS(3)
      DOUBLE COMPLEX P1SG(2,2)
*
      NS = N * N
      NT = NS * N
      NFO = NT * N
      NFI = NFO * N
      NSI = NFI * N
      NSE = NSI * N
      NE = NSE * N
      NN = NE * N
*
      NM2 = N - 2d0
      NM = N - 1d0
      N1 = N + 1d0
      N2 = N + 2d0
      NMS = NM * NM
      NMT = NMS * NM
      N1S = N1 * N1
      N1T = N1S * N1
      N2S = N2 * N2
      N2T = N2S * N2
*
      S1 = EMC + PSI(N1)
      S2 = ZETA2 - DPSI(N1,1)
*
*     Analytic continuations of the occuring sums as given in GRV (1990) 
*     (with an improved parametrization of the moments of  Sp(x)/(1+x).)
*
      N3 = N + 3D0
      N4 = N + 4D0
      N5 = N + 5D0
      N6 = N + 6D0
*
      S11 = S1  + 1D0/N1
      S12 = S11 + 1D0/N2
      S13 = S12 + 1D0/N3
      S14 = S13 + 1D0/N4
      S15 = S14 + 1D0/N5
      S16 = S15 + 1D0/N6
      S1M = S1 - 1D0 / N
      S21 = S2  + 1D0/N1S
      S22 = S21 + 1D0/N2S
      S2M = S2 - 1D0 / NS
*
      SPMOM = 1.0000D0 * (ZETA2 - S1 / N ) / N  -
     1     0.9992D0 * (ZETA2 - S11/ N1) / N1 +
     2     0.9851D0 * (ZETA2 - S12/ N2) / N2 -
     3     0.9005D0 * (ZETA2 - S13/ N3) / N3 +
     4     0.6621D0 * (ZETA2 - S14/ N4) / N4 -
     5     0.3174D0 * (ZETA2 - S15/ N5) / N5 +
     6     0.0699D0 * (ZETA2 - S16/ N6) / N6  
*
      SLC = - 5D0/8D0 * ZETA3
      SLV = - ZETA2/2D0* (PSI(N1/2D0) - PSI(N/2D0)) 
     1     + S1/NS + SPMOM
      SSCHLM = SLC - SLV
      SSTR2M = ZETA2 - DPSI (N1/2D0,1)
      SSTR3M = 0.5D0 * DPSI (N1/2D0,2) + ZETA3
*
      SSCHLP = SLC + SLV
      SSTR2P = ZETA2 - DPSI(N2/2D0,1)
      SSTR3P = 0.5D0 * DPSI(N2/2D0,2) + ZETA3
*
      DS2NM = - DPSI(NM2/2D0+1D0,1) + DPSI(NM/2D0+1D0,1)
      DS2N  = - DPSI(NM/2D0+1D0,1) + DPSI(N/2D0+1D0,1)
      DS2N1 = - DPSI(N/2D0+1D0,1)  + DPSI(N1/2D0+1D0,1)
      DS2N2 = - DPSI(N1/2D0+1D0,1) + DPSI(N2/2D0+1D0,1)
*
*     The contributions to P1NS as given in Gonzalez-Arroyo et al. (1979) 
*     (Note that the anomalous dimensions in the literature often differ 
*     from these moments of the splitting functions by factors -1 or -2,
*     in addition to possible different normalizations of the coupling)
*
      PNMA = ( 16D0* S1 * (2D0* N + 1D0) / (NS * N1S) +
     1     16D0* (2D0* S1 - 1D0/(N * N1)) * ( S2 - SSTR2M ) +
     2     64D0* SSCHLM + 24D0* S2 - 3D0 - 8D0* SSTR3M -
     3     8D0* (3D0* NT + NS -1D0) / (NT * N1T) +
     4     16D0* (2D0* NS + 2D0* N +1D0) / (NT * N1T) ) * (-0.5D0)
      PNPA = ( 16D0* S1 * (2D0* N + 1D0) / (NS * N1S) +
     1     16D0* (2D0* S1 - 1D0/(N * N1)) * ( S2 - SSTR2P ) +
     2     64D0* SSCHLP + 24D0* S2 - 3D0 - 8D0* SSTR3P -
     3     8D0* (3D0* NT + NS -1D0) / (NT * N1T) -
     4     16D0* (2D0* NS + 2D0* N +1D0)/(NT * N1T) ) * (-0.5D0)
*
      PNSB = ( S1 * (536D0/9D0 + 8D0* (2D0* N + 1D0) / (NS * N1S)) -
     1     (16D0* S1 + 52D0/3D0- 8D0/(N * N1)) * S2 - 43D0/6D0 -
     2     (151D0* NFO + 263D0* NT + 97D0* NS + 3D0* N + 9D0) *
     3     4D0/ (9D0* NT * N1T) ) * (-0.5D0)
      PNSC = ( -160D0/9D0* S1 + 32D0/3.* S2 + 4D0/3D0 +
     1     16D0*(11D0*NS+5D0*N-3D0)/(9D0* NS * N1S))*(-0.5D0)
*
*     The contributions to P1SG as given in Floratos et al. (1981) 
*     Pure singlet (PS) and PGG
*
      PPSA = (5d0* NFI + 32d0* NFO + 49d0* NT+38d0* NS + 28d0* N + 8d0) 
     1     / (NM * NT * N1T * N2S) * 2d0     
      PGGA = - (2d0* NFI + 5d0* NFO + 8d0* NT + 7d0* NS- 2d0* N - 2d0)
     1     * 8d0* S1 / (NMS * NS * N1S * N2S) -  67d0/9d0* S1 + 8d0/3d0
     1     - 4d0* SSTR2P * (NS + N + 1d0) / (NM * N * N1 * N2)
     1     + 2d0* S1 * SSTR2P - 4d0* SSCHLP + 0.5d0 * SSTR3P
     1     + (457d0* NN + 2742d0* NE + 6040d0* NSE + 6098d0* NSI
     1     + 1567d0* NFI - 2344d0* NFO - 1632d0* NT + 560d0* NS
     1     + 1488d0* N + 576d0) / (18d0* NMS * NT * N1T * N2T)
      PGGB = (38d0* NFO + 76d0* NT + 94d0* NS + 56d0* N + 12d0) *(-2d0)
     1     / (9d0* NM * NS * N1S * N2)  +  20d0/9d0* S1  -  4d0/3d0
      PGGC = (2d0* NSI + 4d0* NFI + NFO - 10d0* NT - 5d0* NS - 4d0* N
     1     - 4d0) * (-2d0) / (NM * NT * N1T * N2)  -  1d0
*
      PPSTL = -40d0/9d0 * 1d0/NM + 4d0/NT + 10d0/NS - 16d0/N
     1     + 8d0/N1 + 112d0/9d0 * 1d0/N2 + 18d0/N1S
     1     + 4d0/N1T + 16d0/3d0 * 1d0/N2S
*
*     The contributions to P1QQ as given in Gluck Reya Vogt (1993)
*
      PQQATL = ( -4d0 * S1 + 3d0 + 2d0/(N*N1) ) 
     1     * ( 2d0*S2 - 2d0 * ZETA2 - (2d0*N + 1d0)/(NS*N1S) )
      PQQBTL = -80d0/9d0 * 1d0/NM + 8d0/NT + 12d0/NS - 12d0/N 
     1     + 8d0/N1T + 28d0/N1S - 4d0/N1 + 32d0/3d0 * 1d0/N2S 
     1     + 224d0/9d0 * 1/N2
*
*     The contributions to P1QG as given in Gluck Reya Vogt (1993)
*
      PQGA = S11 * (NS + N + 2)/(N * N1 * N2) + 1d0/NS - 5d0/3d0 * 1d0/N
     1     - 1d0/(N * N1) - 2d0/N1S + 4d0/3d0 * 1d0/N1 + 4d0/N2S
     1     - 4d0/3d0 * 1d0/N2
      PQGB = ( - 2d0 * S11**2d0 + 2d0 * S11 + 10d0 * S21 )
     1     * ( NS + N + 2d0 ) / ( N * N1 * N2 )
     1     + 4d0 * S11 
     1     * ( -1d0/NS + 1d0/N + 1d0/(N*N1) + 2d0/N1S - 4d0/N2S )
     1     - 2d0/NT + 5d0/NS - 12d0/N + 4d0/(NS*N1) - 12d0/(N*N1S)
     1     - 6d0/(N*N1) + 4d0/N1T - 4d0/N1S + 23d0/N1 - 20d0/N2
      PQGC = ( 2d0 * S11**2d0 - 10d0/3d0 * S11 - 6d0 * S21
     1     + 1d0 * ( DPSI(N2/2D0,1) - DPSI(N1/2d0,1) ) - 6d0 * ZETA2 ) 
     1     * ( NS + N + 2D0 ) / ( N * N1 * N2 )
     1     - 4d0 * S11 * 
     1     ( -2d0/NS + 1d0/N + 1d0/(N*N1) + 4d0/N1S - 6d0/N2S )
     1     - 40d0/9d0 * 1d0/NM + 4d0/NT + 8d0/3d0 * 1d0/NS 
     1     + 26d0/9d0 * 1/N - 8d0/(NS*N1S) + 22d0/3d0 * 1d0/(N*N1)
     1     + 16d0/N1T + 68d0/3d0 * 1d0/N1S - 190d0/9d0 * 1d0/N1
     1     + 8d0/(N1S*N2) - 4d0/N2S + 356d0/9d0 * 1d0/N2
*
*     The contributions to P1GQ as given in Gluck Reya Vogt (1993)
*
      PGQA = ( S1**2d0 - 3d0*S2 - 4d0 * ZETA2 )
     1     * ( NS + N + 2d0 ) / ( NM * N * N1 )
     1     + 2d0 * S1
     1     * ( 4d0/NMS - 2d0/(NM*N) - 4d0/NS + 3d0/N1S - 1d0/N1 )
     1     - 8d0/(NMS*N) + 8d0/(NM*NS) + 2d0/NT + 8d0/NS - 1d0/(2d0*N)
     1     + 1d0/N1T - 5d0/2d0 * 1d0/N1S + 9d0/2d0 * 1d0/N1
*
      PGQB = ( -1d0 * S1**2d0 + 5d0 * S2
     1     - 0.5d0 * ( DPSI(N1/2d0,1) - DPSI(N/2d0,1) ) + ZETA2 )
     1     * ( NS + N + 2d0 ) / ( NM * N * N1 )
     1     + 2d0 * S1 
     1     * ( -2d0/NMS + 2d0/(NM*N) + 2d0/NS - 2d0/N1S + 1d0/N1 )
     1     - 8d0/NMT + 6d0/NMS + 17d0/9d0 * 1d0/NM + 4d0/(NMS*N) 
     1     - 12d0/(NM*NS) - 8d0/NS + 5d0/N - 2d0/(NS*N1) - 2d0/N1T
     1     - 7d0/N1S - 1d0/N1 - 8d0/3d0 * 1d0/N2S - 44d0/9d0 * 1d0/N2
*
*     The contributions to P1GG as given in Gluck Reya Vogt (1993)
*
      PGGATL = - 16d0/3d0 * 1d0/NMS + 80d0/9d0 * 1d0/NM + 8d0/NT
     1     - 16d0/NS + 12d0/N + 8d0/N1T - 24d0/N1S + 4d0/N1
     1     - 16d0/3d0 * 1d0/N2S - 224d0/9d0 * 1d0/N2
      PGGBTL = S2 - 1d0/NMS + 1d0/NS - 1d0/N1S + 1d0/N2S - ZETA2
*
      PGGCTL = - 8D0 * S1 * S2 + 8D0 * S1
     1     * ( 1D0 / NMS - 1D0 / NS + 1D0 / N1S - 1D0 / N2S + ZETA2 )
     1     + ( 8D0 * S2 - 8D0 * ZETA2 )
     1     * ( 1D0 / NM - 1D0 / N + 1D0 / N1 - 1D0 / N2 + 11D0 / 12D0 )
     1     - 8D0 / NMT + 22D0 / 3D0 * 1D0 / NMS - 8D0 / ( NMS * N )
     1     - 8D0 / ( NM * NS ) - 8D0 / NT - 14D0 / 3D0 * 1D0 / NS
     1     - 8D0 / N1T + 14D0 / 3D0 * 1D0 / N1S - 8D0 / ( N1S * N2 )
     1     - 8D0 / ( N1 * N2S ) - 8D0 / N2T - 22D0 / 3D0 * 1D0 / N2S
*
*     The contributions to P1NS as given in Gluck Reya Vogt (1993)
*
      PNSTL = ( - 4d0*S1 + 3d0 + 2d0/(N*N1) )
     1     * ( 2d0*S2 - 2d0*ZETA2 - (2d0*N + 1d0)/(NS*N1S) )
*
*     Output to the array
*
      DO I=1,2
         DO J=1,2
            P1SG(I,J) = 0D0
         ENDDO
      ENDDO
*
      do I=1,3
         P1NS(I) = 0D0
      enddo
*
*     NON SINGLET     
*
*     Plus  
      P1NS(1) = CF *((CF-CA/2d0)* PNPA + CA* PNSB + TR*dble(NF)* PNSC)
     1        + CF**2d0 * PNSTL * 4d0 
*     Minus=Valence
      P1NS(2) = CF *((CF-CA/2d0)* PNMA + CA* PNSB + TR*dble(NF)* PNSC)
     1        + CF**2d0 * PNSTL * 4d0 
      P1NS(3) = P1NS(2)
*     SINGLET (time-like)
      P1SG(1,1) = P1NS(1) + TR * dble(NF) * CF * PPSTL * 4d0
      P1SG(1,2) = ( CF**2d0 * PGQA + CF * CA * PGQB ) 
     1          * 4d0 * 2d0 * dble(NF)
      P1SG(2,1) = ( 8d0/3d0 * (TR * dble(NF) )**2d0 * PQGA 
     1          + CF * TR * dble(NF) * PQGB 
     1          + CA * TR * dble(NF) * PQGC )
     1          * 4d0 / (2d0 * dble(NF) ) 
      P1SG(2,2) = (CA*CA*PGGA + TR*dble(NF)*(CA*PGGB+CF*PGGC))*4d0
     1          + CF * TR * dble(NF) * 4d0 * PGGATL
     1          + CA * TR * dble(NF) * 4d0 * (-8d0/3d0) * PGGBTL
     1          + CA * CA * 4d0 * PGGCTL
*
      RETURN
      END
