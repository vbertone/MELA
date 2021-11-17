*************************************************************
*
*     File: andim_lo.f
*     
*     Returns the LO anomalous dimensions in the 
*     a=alpha/4pi notation
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
      SUBROUTINE ANDIM_LO(N,NF,P0NS,P0SG)
*     
      IMPLICIT NONE
*     
      include "../commons/colfact.h"
      include "../commons/consts.h"
*     
*     ---------------------------------------------------------------------
*     
*     Internal variables
*     
      DOUBLE COMPLEX NS,N1,N2,NM
      DOUBLE COMPLEX PSI,S1
      DOUBLE COMPLEX PQQA,PQGA,PGQA,PGGA,PGGB
*     
*     ---------------------------------------------------------------------
*     
*     Input variables
*     
      DOUBLE COMPLEX N
      INTEGER NF
      DOUBLE PRECISION NFSUM2
*     
*     Output variables  
*     
      DOUBLE COMPLEX P0NS
      DOUBLE COMPLEX P0SG(2,2)
*     
*     ---------------------------------------------------------------------
*     
      NS = N * N
      N1 = N + cmplx(1d0,0d0)
      N2 = N + (2d0,0d0)
      NM = N - (1d0,0d0)
*     
      S1 = cmplx(EMC,0d0) + PSI(N1)
*     
      PQQA = (3d0,0d0) - 4d0* S1 + 2d0/(N * N1)
      PQGA = 4d0* (NS + N + 2d0) / (N * N1 * N2)
      PGQA = 2d0 * (NS + N + 2d0) / (N * N1 * NM)
      PGGB = - 4d0/3D0
*     
*     Output to the array
*     
      P0NS      = PQQA
*
**     gstagn: hacking to add quark contributions 3D0 * CH2(5)       
      NFSUM2    = dble(NF)
*      NFSUM2    = dble(NF) + 3D0 * CH2(5)
*      
      P0SG(1,1) = PQQA
      P0SG(1,2) = dble(NF) * PQGA
      P0SG(2,1) = PGQA
      P0SG(2,2) = NFSUM2 * PGGB
*     
      RETURN
      END
*     
*----------------------------------------------------------------------

