************************************************************************
*
*     initBeta.f:
*
*     Initialize the coefficients of the QED beta function.
*
************************************************************************
      subroutine initBeta
*
      implicit none
*
      include "../commons/activeflavours.h"
      include "../commons/charges.h"
      include "../commons/beta.h"
**
*     Internal Variables
*
      integer i
*
*     Initialise numver of active lepton, up-type, and down-type
*     flavours.
*
      nl(0) = 0
      nl(1) = 1
      nl(2) = 1
      nl(3) = 1
      nl(4) = 1
      nl(5) = 2
      nl(6) = 2
      nl(7) = 3
      nl(8) = 3
      nl(9) = 3
*
      nu(0) = 0
      nu(1) = 0
      nu(2) = 1
      nu(3) = 1
      nu(4) = 1
      nu(5) = 1
      nu(6) = 2
      nu(7) = 2
      nu(8) = 2
      nu(9) = 3
*
      nd(0) = 0
      nd(1) = 0
      nd(2) = 0
      nd(3) = 1
      nd(4) = 2
      nd(5) = 2
      nd(6) = 2
      nd(7) = 2
      nd(8) = 3
      nd(9) = 3
*
*     Number of colours
*
      nc = 3
*
*     electric charges
*
      el2 = 1d0
      eu2 = 4d0 / 9d0
      ed2 = 1d0 / 9d0
*
      el4 = el2**2
      eu4 = eu2**2
      ed4 = ed2**2
*
*     Beta function and mass anomalous dimension coefficients
*
      do i = 1, 9
         beta0(i) = - 4d0 * ( nl(i)
     1        + nc * ( ed2 * nd(i) + eu2 * nu(i) ) ) / 3d0
         beta1(i) = - 4d0 * ( nl(i)
     1        + nc * ( ed4 * nd(i) + eu4 * nu(i) ) )
         b1(i) = beta1(i) / beta0(i)
      enddo
*
      return
      end
