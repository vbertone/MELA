************************************************************************
*
*     SetNFFN.f:
*
*     Sets the numer of flavours to be used when running in the FFNS
*
************************************************************************
      subroutine SetNFFN(NFFNin)
*
      implicit none
*
      include "../commons/nffn.h"
**
*     Input Variables
*
      integer NFFNin
*
      NFFN = NFFNin
      NFFNalpha = NFFNin
*
      return
      end
************************************************************************      
      subroutine SetNFFNalpha(NFFNalphain)
*
      implicit none
*
      include "../commons/nffn.h"
**
*     Input Variables
*
      integer NFFNalphain
*
      NFFNalpha = NFFNalphain
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
      subroutine GetNFFNalpha(nffnalphaout)
      implicit none
      include "../commons/nffn.h"      
      integer nffnalphaout
      nffnalphaout=nffnalpha
      return
      end
      
