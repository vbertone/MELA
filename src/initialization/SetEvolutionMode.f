************************************************************************
*
*     SetEvolutionMode.f:
*
*     This subroutine sets the evolution mode.
*
************************************************************************
      subroutine SetEvolutionMode(evmode)
*
      implicit none
*
      include "../commons/modev.h"
*
*     Variables
*
      character*3 evmode
*
      MODEV = evmode
*
      return
      end
