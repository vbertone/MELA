************************************************************************
*
*     SetFactorisationScheme.f:
*
*     Sets the factorisation scheme (MSBAR or DELTA)
*
************************************************************************
      subroutine SetFactorisationScheme(fsin)
*
      implicit none
*
      include "../commons/facscheme.h"
**
*     Input Variables
*
      character*5 fsin
*
      facscheme = fsin
*
      return
      end
************************************************************************
*     Easy the C interface by using a int
*     0 = MSBAR, 1 = DELTA
************************************************************************      
      subroutine SetFactorisationSchemeInt(fsinint)
*
      implicit none
*
      include "../commons/facscheme.h"
**
*     Input Variables
*     
      integer fsinint
*     
      if (fsinint.eq.0) then
         facscheme = "MSBAR"
      elseif (fsinint.eq.1) then
         facscheme = "DELTA"
      else
         write(6,*) "In SetFactorisationScheme.f:"
         write(6,*) "Invalid value"
         call exit(-10)
      endif 
*
      return
      end

      
