************************************************************************
*
*     SetNFmax.f:
*
*     Sets the maximum number of active flavours to be used in the VFNS
*
************************************************************************
      subroutine SetNFmax(NFmaxin)
*
      implicit none
*
      include "../commons/nfmax.h"
**
*     Input Variables
*
      integer NFmaxin
*
      NFmax = NFmaxin
      NFmaxalpha = NFmaxin
*
      return
      end
************************************************************************      
      subroutine SetNFmaxalpha(NFmaxalphain)
*
      implicit none
*
      include "../commons/nfmax.h"
**
*     Input Variables
*
      integer NFmaxalphain
*
      NFmaxalpha = NFmaxalphain
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
      subroutine GetNFmaxalpha(nfmaxout)
      implicit none
      include "../commons/nfmax.h"      
      integer nfmaxout
      nfmaxout=nfmaxalpha
      return
      end
