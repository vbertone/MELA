************************************************************************
*
*     SetRenormalisationScheme.f:
*
*     Sets the renormalisation scheme (MSBAR, FIXED, ALPMZ)
*
************************************************************************
      subroutine SetRenormalisationScheme(rsin)
*
      implicit none
*
      include "../commons/renscheme.h"
**
*     Input Variables
*
      character*5 rsin
*
      renscheme = rsin
*
      return
      end
************************************************************************
*     Easy the C interface by using a int
*     0 = MSBAR, 1 = FIXED, 2 = ALPMZ
************************************************************************      
      subroutine SetRenormalisationSchemeInt(rsinint)
*
      implicit none
*
      include "../commons/renscheme.h"
**
*     Input Variables
*     
      integer rsinint
*     
      if (rsinint.eq.0) then
         renscheme = "MSBAR"
      elseif (rsinint.eq.1) then
         renscheme = "FIXED"
      elseif (rsinint.eq.2) then
         renscheme = "ALPMZ"
      else
         write(6,*) "In SetRenormalisationScheme.f:"
         write(6,*) "Invalid value"
         call exit(-10)
      endif 
*
      return
      end
************************************************************************
      subroutine GetRenormalisationSchemeInt(rsintout)
      implicit none
      include "../commons/renscheme.h"      
      integer rsintout
      if (renscheme.eq."MSBAR") then
         rsintout = 0
      elseif (renscheme.eq."FIXED") then
         rsintout = 1
      elseif (renscheme.eq."ALPMZ") then
         rsintout = 2
      endif
      return
      end

      
