************************************************************************
*
*     SetFlavourScheme.f:
*
*     Sets the evolution scheme (FFNS or VFNS)
*
************************************************************************
      subroutine SetFlavourScheme(nsin)
*
      implicit none
*
      include "../commons/ns.h"
**
*     Input Variables
*
      character*4 nsin
*
      ns = nsin
*
      return
      end
