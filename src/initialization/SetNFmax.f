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
*
      return
      end
