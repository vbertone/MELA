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
      integer nlaem(0:9), nuaem(0:9), ndaem(0:9)
      double precision C2aem(1:9), C4aem(1:9)
      double precision C2(1:9), C4(1:9)      
*
*     fill number of colours
*
      nc = 3
*
*     fill electric charges
*
      el2 = 1d0
      eu2 = 4d0 / 9d0
      ed2 = 1d0 / 9d0
*
      el4 = el2**2
      eu4 = eu2**2
      ed4 = ed2**2
*     
*     fill b0 and b1 according to nlmaxaem, numaxaem, ndmaxaem
*
      call setnlnund(nlmaxaem,numaxaem,ndmaxaem,nlaem,nuaem,ndaem)
      call fillC2(nlaem,nuaem,ndaem,C2aem)
      call fillC4(nlaem,nuaem,ndaem,C4aem)            
      do i = 1,9         
         beta0(i) = - 4d0 / 3d0 * C2aem(i)
         beta1(i) = - 4d0 * C4aem(i)
      enddo
*     
*     fill nl,nu,nd,nfsum2,nfsum4 according to nlmax, numax, ndmax
*
      call setnlnund(nlmax,numax,ndmax,nl,nu,nd)
      call fillC2(nl,nu,nd,C2)
      call fillC4(nl,nu,nd,C4)            
      do i = 1,9     
         nfsum2(i) = C2(i)
         nfsum4(i) = C4(i)
      enddo
*
      return
      end
************************************************************************
*     Calculate C2, sum of charges squared
      subroutine fillC2(nl,nu,nd,C2)
      implicit none
      include "../commons/charges.h"
*     Input
      integer nl(0:9), nu(0:9), nd(0:9)      
*     Output
      double precision C2(1:9)
*     Internal
      integer i
*
      do i = 1,9
         C2(i) = el2 * nl(i) + nc * ( eu2 * nu(i) + ed2 * nd(i) )  
      enddo
*
      return
      end
************************************************************************
*     Calculate C4, sum of charges to the fourth power
      subroutine fillC4(nl,nu,nd,C4)
      implicit none
      include "../commons/charges.h"
*     Input
      integer nl(0:9), nu(0:9), nd(0:9)      
*     Output
      double precision C4(1:9)
*     Internal
      integer i
*      
      do i = 1,9
         C4(i) = el4 * nl(i) + nc * ( eu4 * nu(i) + ed4 * nd(i) )
      enddo
*
      return
      end      
************************************************************************      
*     Utility subroutine to fill the nl,nu,nd vectors according to
*     the maximum number of allowed leptons, down- and up-quarks.
      subroutine setnlnund(nlmax,numax,ndmax,nl,nu,nd)
      implicit none
*     Input
      integer nlmax, numax, ndmax
*     Output
      integer nl(0:9), nu(0:9), nd(0:9)      
*     Internal
      integer i
*           
      do i=0,9
         nl(i) = 0
      enddo      
      if (nlmax.ge.1) then
         do i=1,9
            nl(i) = 1
         enddo
      endif
      if (nlmax.ge.2) then
         do i=5,9
            nl(i) = 2
         enddo
      endif
      if (nlmax.ge.3) then
         do i=7,9
            nl(i) = 3
         enddo
      endif
*
      do i=0,9
         nu(i) = 0
      enddo      
      if (numax.ge.1) then
         do i=2,9
            nu(i) = 1
         enddo
      endif
      if (numax.ge.2) then
         do i=6,9
            nu(i) = 2
         enddo
      endif
      if (numax.ge.3) then
         do i=9,9
            nu(i) = 3
         enddo
      endif
*      
      do i=0,9
         nd(i) = 0
      enddo      
      if (ndmax.ge.1) then
         do i=3,9
            nd(i) = 1
         enddo
      endif
      if (ndmax.ge.2) then
         do i=4,9
            nd(i) = 2
         enddo
      endif
      if (ndmax.ge.3) then
         do i=8,9
            nd(i) = 3
         enddo
      endif
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
