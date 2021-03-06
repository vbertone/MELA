************************************************************************
*
*     initBeta.f:
*
*     Initialize the coefficients of the QCD beta function.
*
************************************************************************
      subroutine initBeta
*
      implicit none
*
      include "../commons/beta.h"
**
*     Internal Variables
*
      integer i
*
*     Beta function and mass anomalous dimension coefficients
*
      do i=3,6
         beta0(i) = ( 33d0 - 2d0 * i ) / 3d0
         beta1(i) = 102d0 - 38d0 / 3d0 * i
         beta2(i) = 2857d0/2d0 - 5033d0/18d0 * i + 325d0/54d0 * i**2d0
         b1(i) = beta1(i) / beta0(i)
         b2(i) = beta2(i) / beta0(i)
      enddo
*
      return
      end
