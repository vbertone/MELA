************************************************************************
*
*     ComputeChargesDIS_MELA.f:
*
*     It sets the charges used for the computation of the DIS observables.
*
************************************************************************
      subroutine ComputeChargesDIS_MELA(Q2,bq,dq)
*
      implicit none
*
      include "../commons/process.h"
      include "../commons/SinThetaW.h"
      include "../commons/ZedMass.h"
**
*     Input Variables
*
      double complex Q2
**
*     Input Variables
*
      integer i
      double precision eq(6),eq2(6)
      double precision vq(6),aq(6)
      double precision ve,ae
      double complex pz,pz2
**
*     Double precision
*
      double complex bq(6),dq(6)
*
*     Initialize charges and couplings
*
*     Electric Charges
*
      eq(1) = 2d0 / 3d0
      eq(2) = - 1d0 / 3d0
      eq(3) = - 1d0 / 3d0
      eq(4) = 2d0 / 3d0
      eq(5) = - 1d0 / 3d0
      eq(6) = 2d0 / 3d0
*
*     Squared Charges
*
      eq2(1) = eq(1) * eq(1) ! 4d0 / 9d0
      eq2(2) = eq(2) * eq(2) ! 1d0 / 9d0
      eq2(3) = eq(3) * eq(3) ! 1d0 / 9d0
      eq2(4) = eq(4) * eq(4) ! 4d0 / 9d0
      eq2(5) = eq(5) * eq(5) ! 1d0 / 9d0
      eq2(6) = eq(6) * eq(6) ! 4d0 / 9d0
*
*     Vector Couplings
*
      vq(1) = + 0.5d0 - 4d0 / 3d0 * SinThetaW
      vq(2) = - 0.5d0 + 2d0 / 3d0 * SinThetaW
      vq(3) = - 0.5d0 + 2d0 / 3d0 * SinThetaW
      vq(4) = + 0.5d0 - 4d0 / 3d0 * SinThetaW
      vq(5) = - 0.5d0 + 2d0 / 3d0 * SinThetaW
      vq(6) = + 0.5d0 - 4d0 / 3d0 * SinThetaW
*
*     Axial Couplings
*
      aq(1) = + 0.5d0
      aq(2) = - 0.5d0
      aq(3) = - 0.5d0
      aq(4) = + 0.5d0
      aq(5) = - 0.5d0
      aq(6) = + 0.5d0
*
*     Vector and Axial Electron Couplings
*
      ve = - 0.5d0 + 2d0 * SinThetaW
      ae = - 0.5d0
*
      if(proc.eq."EM")then
         do i=1,6
            bq(i) = eq2(i)
            dq(i) = (0d0,0d0)
         enddo
      elseif(proc.eq."NC")then
         pz  = Q2 / ( Q2 + MZ**2d0 ) 
     1        / ( 4d0 * SinThetaW * ( 1d0 - SinThetaW ) )
         pz2 = pz * pz
         do i=1,6
            bq(i) = eq2(i) 
     1            - 2d0 * eq(i) * vq(i) * ve * pz
     2            + ( ve**2d0 + ae**2d0 ) * ( vq(i)**2d0 + aq(i)**2d0 )
     3            * pz2
            dq(i) = - 2d0 * eq(i) * aq(i) * ae * pz
     1            + 2d0 * vq(i) * aq(i) * ( 2d0 * ve * ae ) * pz2
         enddo
      endif
*
      return
      end
