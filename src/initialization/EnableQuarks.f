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
************************************************************************      
      subroutine GetEnableQuarks(quarksout)
      implicit none
      include "../commons/activeflavours.h"
      logical quarksout
      quarksout = quarks
      return
      end
************************************************************************
      subroutine GetEnableQuarksalpha(quarksout)
      include "../commons/activeflavours.h"
      logical quarksout
      quarksout = quarksalpha
      return
      end
