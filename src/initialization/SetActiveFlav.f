************************************************************************
*
*     SetNFmax.f:
*
*     Sets the maximum number of active flavours to be used in the VFNS
*
************************************************************************
*     this set also nfmax for alpha
      subroutine SetActiveFlav(nlmaxin,numaxin,ndmaxin)
      implicit none
      include "../commons/activeflavours.h"
      integer nlmaxin,numaxin,ndmaxin
*
      nlmax = nlmaxin
      numax = numaxin
      ndmax = ndmaxin
*
      nlmaxaem = nlmax
      numaxaem = numax
      ndmaxaem = ndmax      
*      
      return
      end
************************************************************************
      subroutine SetActiveFlavAem(nlmaxaemin,numaxaemin,ndmaxaemin)
      implicit none
      include "../commons/activeflavours.h"
      integer nlmaxaemin,numaxaemin,ndmaxaemin
*
      nlmaxaem = nlmaxaemin
      numaxaem = numaxaemin
      ndmaxaem = ndmaxaemin
*
      return
      end
************************************************************************
      subroutine GetActiveFlav(nlmaxout,numaxout,ndmaxout)
      implicit none
      include "../commons/activeflavours.h"
      integer nlmaxout,numaxout,ndmaxout
*
      nlmaxout = nlmax
      numaxout = numax
      ndmaxout = ndmax
*      
      return
      end
************************************************************************
      subroutine GetActiveFlavAem(nlmaxaemout,numaxaemout,ndmaxaemout)
      implicit none
      include "../commons/activeflavours.h"
      integer nlmaxaemout,numaxaemout,ndmaxaemout
*
      nlmaxaemout = nlmaxaem
      numaxaemout = numaxaem
      ndmaxaemout = ndmaxaem
*
      return
      end
************************************************************************      
      subroutine SetWaem(waemin)
      implicit none
      include "../commons/activeflavours.h"
      integer waemin
      waem = waemin
      return
      end
************************************************************************
      subroutine GetWaem(waemout)
      implicit none
      include "../commons/activeflavours.h"
      integer waemout
      waemout = waem
      return
      end
