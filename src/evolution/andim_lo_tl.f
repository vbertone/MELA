********************************************************************************
*                                                                              *
*     andim_lo_tl.f                                                            *
*                                                                              *
********************************************************************************
*                                                                              *
*     Returns the LO anomalous dimensions in the a_s=alpha_s/4pi notation      *
*                                                                              *
*     Input variables:                                                         *
*     N  : Mellin variable                                                     *
*     NF : Number of flavours                                                  *
*                                                                              *
*     Output:                                                                  *
*     P0NS : Complex valued non-singlet LO anomalous dimension in N space      *
*     P0SG(2,2) : matrix of the complex-valued singlet                         *
*                 LO anomalous dimension in N space                            *
*                                                                              *
*     The explicit expressions for the singlet and non-singlet N-space         *
*     splitting functions are taken from hep-ph/0604160, eqs.(B.3)-(B.6),      *
*     with an overall "-" sign coming from the definition of the anomalous     *
*     dimension (which differs from our by a "-" sign, see eq.25)              *
*     Notice that the time-like splitting functions coincide, at LO, with      *
*     their space-like counterparts, provided the off-diagonal Pqg and Pgq     *
*     are exchanged                                                            *
*                                                                              *
********************************************************************************
      SUBROUTINE ANDIM_LO_TL(N,NF,P0NS,P0SG)
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
      DOUBLE COMPLEX NS,N1,N2,NM
      DOUBLE COMPLEX PQQA,PQGA,PGQA,PGGA,PGGB
      DOUBLE COMPLEX PSI, S1
**
*     Output variables
*
      DOUBLE COMPLEX P0NS
      DOUBLE COMPLEX P0SG(2,2)
*
      NS = N * N
      N1 = N + cmplx(1d0,0d0)
      N2 = N + (2d0,0d0)
      NM = N - (1d0,0d0)
*
      S1 = cmplx(EMC,0d0) + PSI(N1)
*
      PQQA = (3d0,0d0) - 4d0 * S1 + 2d0/(N * N1)
      PQGA = 2d0 * (NS + N + 2d0) / (N * N1 * N2)
      PGQA = 4d0 * (NS + N + 2d0) / (N * N1 * NM)
      PGGA = 11d0/3D0 - 4d0* S1 + 4d0/(N * NM) + 4d0/(N1 * N2) 
      PGGB = - 4d0/3D0
*
*     Output to the array
*
      P0NS      = CF * PQQA
*
      P0SG(1,1) = CF * PQQA
      P0SG(1,2) = dble(NF) * CF * PGQA
      P0SG(2,1) = TR * PQGA
      P0SG(2,2) = CA * PGGA + TR * dble(NF) * PGGB
*
      RETURN
      END
