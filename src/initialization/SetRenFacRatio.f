************************************************************************
*
*     SetRenFacRatio.f:
*
*     This subroutine sets the ratio between renormalization and
*     factorization scales in the PDF evolution.
*
************************************************************************
      subroutine SetRenFacRatio(ratio)
*
      implicit none
*
      include "../commons/renfacscales.h"
*
*     Variables
*
      double precision ratio
*
      KRF = ratio
*
      return
      end
