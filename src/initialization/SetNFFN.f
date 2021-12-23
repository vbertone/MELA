************************************************************************
*
*     SetNFFN.f:
*
*     Sets the numer of flavours to be used when running in the FFNS
*
************************************************************************
*     this set also nffn for alpha
      subroutine SetNFFN(nlmax,numax,ndmax)      
      implicit none
      include "../commons/nffn.h"
      integer nlmax,numax,ndmax
*
*     we always set nfmax to 8 i.e. we consider all thresholds
*
      nffn = 8
      nffnalpha = nffn
*
      return
      end
************************************************************************      
      subroutine SetNffnalpha(nlmax,numax,ndmax)
      implicit none
      include "../commons/nffn.h"
      integer nlmax,numax,ndmax      
*
      nffnalpha = 8
*
      return
      end
************************************************************************
      subroutine GetNFFN(nffnout)
      implicit none
      include "../commons/nffn.h"      
      integer nffnout
      nffnout=nffn
      return
      end
************************************************************************
      subroutine GetNffnalpha(nffnalphaout)
      implicit none
      include "../commons/nffn.h"      
      integer nffnalphaout
      nffnalphaout=nffnalpha
      return
      end
      
