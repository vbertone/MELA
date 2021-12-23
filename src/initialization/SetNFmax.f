************************************************************************
*
*     SetNFmax.f:
*
*     Sets the maximum number of active flavours to be used in the VFNS
*
************************************************************************
*     this set also nfmax for alpha
      subroutine SetNFmax(nlmax,numax,ndmax)
      implicit none
      include "../commons/nfmax.h"
      integer nlmax,numax,ndmax
*      
*     we always set nfmax to 8 i.e. we consider all thresholds
*
      nfmax = 8
      nfmaxalpha = nfmax
*      
      return
      end
************************************************************************
      subroutine SetNfmaxalpha(nlmaxaem,numaxaem,ndmaxaem)      
      implicit none
      include "../commons/nfmax.h"
      integer nlmaxaem,numaxaem,ndmaxaem
*
      nfmaxalpha = 8      
*
      return
      end      
************************************************************************      
      subroutine GetNFmax(nfmaxout)
      implicit none
      include "../commons/nfmax.h"      
      integer nfmaxout
      nfmaxout=nfmax
      return
      end
************************************************************************      
      subroutine GetNfmaxalpha(nfmaxout)
      implicit none
      include "../commons/nfmax.h"      
      integer nfmaxout
      nfmaxout=nfmaxalpha
      return
      end
