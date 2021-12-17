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
      include "../commons/nfsum.h"
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
*     Beta function coefficients
*     
      do i = 1,9
         if (quarksalpha) then
            beta0(i) = - 4d0 / 3d0 *
     1           ( el2 * nl(i) + nc * ( eu2 * nu(i) + ed2 * nd(i) ) )
            beta1(i) = - 4d0 *
     1           ( el4 * nl(i) + nc * ( eu4 * nu(i) + ed4 * nd(i) ) )
         else
            beta0(i) = - 4d0 * el2 * nl(i) / 3d0
            beta1(i) = - 4d0 * el4 * nl(i)    
         endif
      enddo
*     
*     Turn off quark charges if option quark on
      if (.not.quarks) then
         eu2 = 0d0
         ed2 = 0d0
         eu4 = 0d0
         ed4 = 0d0
      endif
*     
      do i = 1,9      
         nfsum2(i) = el2 * nl(i) + nc * ( eu2 * nu(i) + ed2 * nd(i) )
         nfsum4(i) = el4 * nl(i) + nc * ( eu4 * nu(i) + ed4 * nd(i) )
      enddo
*
      return
      end
************************************************************************
      subroutine Getb0(b0)
      implicit none
      include "../commons/beta.h"
      include "../commons/consts.h"
      double precision b0(9)
      b0(:) = -beta0(:)/4d0/pi
      return
      end
************************************************************************
      subroutine Getb1(b1)
      implicit none
      include "../commons/beta.h"
      include "../commons/consts.h"      
      double precision b1(9)
      b1(:) = -beta1(:)/4d0/4d0/pi/pi
      return
      end      
************************************************************************
      subroutine GetC2(C2)
      implicit none
      include "../commons/nfsum.h"
      double precision C2(9)
      C2(:) = nfsum2(:)
      return
      end
************************************************************************
      subroutine GetC4(C4)
      implicit none
      include "../commons/nfsum.h"
      double precision C4(9)
      C4(:) = nfsum4(:)
      return
      end
