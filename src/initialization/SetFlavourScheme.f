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
************************************************************************
*     Easy the C interface by using a int
*     0 = FFNS, 1 = VFNS
************************************************************************      
      subroutine SetFlavourSchemeInt(nsinint)
*
      implicit none
*     
      include "../commons/ns.h"
**
*     Input Variables
*     
      integer nsinint
*     
      if (nsinint.eq.0) then
         ns = "FFNS"
      elseif (nsinint.eq.1) then
         ns = "VFNS"
      else
         write(6,*) "In SetFlavourScheme.f:"
         write(6,*) "Invalid value"
         call exit(-10)
      endif 
*
      return
      end

      
