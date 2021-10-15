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
         beta0(i) = - 4d0 * i / 3d0
         beta1(i) = - 4d0 * i
         b1(i) = beta1(i) / beta0(i)
      enddo
*
      return
      end
