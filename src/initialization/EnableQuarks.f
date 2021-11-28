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
      quarksalpha = quarksin
*
      return
      end
************************************************************************
      subroutine EnableQuarksalpha(quarksin)
*
      implicit none
*
      include "../commons/activeflavours.h"
**
*     Input Variables
*
      logical quarksin
*
      quarksalpha = quarksin
*
      return
      end
