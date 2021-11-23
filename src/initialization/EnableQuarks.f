************************************************************************
*
*     EnableQuarks.f:
*
*     Enables or disables quarks in the evolution
*
************************************************************************
      subroutine EnableQuarks(quarksin)
*
      implicit none
*
      include "../commons/activeflavours.h"
**
*     Input Variables
*
      logical quarksin
*
      quarks = quarksin
*
      return
      end
