*************************************************************
*
*     File: andim_lo_pol.f
*     
*     Returns the polarized LO anomalous dimensions in the 
*     a_s=alpha_s/4pi notation
*
*     Input variables: N  : Mellin variable
*                      NF : Number of flavours
*
*     Functions: PSI from data/theo/harmsums/psi.f
*
*     Output: P0NS : Complex valued Non-singlet LO anomalous
*                    dimension in N space
*
*             P0SG(2,2) : matrix of the complex valued LO
*                         singlet anomalous dimensions in 
*                         N space
*
*************************************************************
       SUBROUTINE ANDIM_LO_POL(N,NF,P0NS,P0SG)
*
       IMPLICIT NONE
*
       include "../commons/colfact.h"
       include "../commons/consts.h"
*
       DOUBLE COMPLEX PSI, S1
*
* ---------------------------------------------------------------------
*
*     Internal variables
*
       DOUBLE COMPLEX NS,N1,N2,NM
       DOUBLE COMPLEX PQQA,PQGA,PGQA,PGGA,PGGB
*
* ---------------------------------------------------------------------
*
*     Input variables
*
       DOUBLE COMPLEX N
       INTEGER NF       
*
*     Output variables  
*
       DOUBLE COMPLEX P0NS 
       DOUBLE COMPLEX P0SG(2,2)
*
* ---------------------------------------------------------------------
*
       NS = N * N
       N1 = N + cmplx(1d0,0d0)
       N2 = N + (2d0,0d0)
       NM = N - (1d0,0d0)
*       
       S1 = cmplx(EMC,0d0) + PSI(N1)
*
*     These are the LO polarized splitting functions taken
*     from Vogt's code PEGASUS
*
       PQQA = (3d0,0d0) - 4d0 * S1 + 2d0/(N * N1)
       PQGA = 4d0 * NM / (N * N1)
       PGQA = 2d0 * N2 / (N * N1)
       PGGA = 11D0/3D0 - 4D0 * S1 + 8D0/(N * N1)
       PGGB = - 4D0/3D0
*          
*     Output to the array
*     
       P0NS      = CF * PQQA
*
       P0SG(1,1) = CF * PQQA
       P0SG(1,2) = TR * dble(NF) * PQGA
       P0SG(2,1) = CF * PGQA
       P0SG(2,2) = CA * PGGA + TR * dble(NF) * PGGB
*
       RETURN
       END
