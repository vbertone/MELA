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
************************************************************************
*     Easy the C interface by using a int
*     0 = no quarks, 1 = yes quarks
************************************************************************      
      subroutine EnableQuarksInt(quarksinint)
*
      implicit none
*
      include "../commons/activeflavours.h"
**
*     Input Variables
*
      integer quarksinint
*     
      if (quarksinint.eq.0) then
         quarks = .false.
      elseif (quarksinint.eq.1) then
         quarks = .true.
      else
         write(6,*) "In EnableQuarks.f:"
         write(6,*) "Invalid value"
         call exit(-10)
      endif 
*
      return
      end
