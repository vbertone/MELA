************************************************************************
*
*     SetTecParam.f:
*
*     Sets the technical parameters of MELA
*
************************************************************************
      subroutine setnint(varin)
      implicit none
      include "../commons/tecparam.h"
      integer varin
      nint = varin
      return
      end
************************************************************************
      subroutine getnint(varout)
      implicit none
      include "../commons/tecparam.h"
      integer varout
      varout = nint
      return
      end
************************************************************************
      subroutine setnexp(varin)
      implicit none
      include "../commons/tecparam.h"
      integer varin
      nexp = varin
      return
      end
************************************************************************
      subroutine getnexp(varout)
      implicit none
      include "../commons/tecparam.h"
      integer varout
      varout = nexp
      return
      end
************************************************************************
      subroutine setnstep(varin)
      implicit none
      include "../commons/tecparam.h"
      integer varin
      nstep = varin
      return
      end
************************************************************************
      subroutine getnstep(varout)
      implicit none
      include "../commons/tecparam.h"
      integer varout
      varout = nstep
      return
      end
************************************************************************
      subroutine setminvmel(varin)
      implicit none
      include "../commons/tecparam.h"
      integer varin
      minvmel = varin
      return
      end
************************************************************************
      subroutine getminvmel(varout)
      implicit none
      include "../commons/tecparam.h"
      integer varout
      varout = minvmel
      return
      end
************************************************************************
      subroutine setrinvmel(varin)
      implicit none
      include "../commons/tecparam.h"
      integer varin
      rinvmel = varin
      return
      end
************************************************************************
      subroutine getrinvmel(varout)
      implicit none
      include "../commons/tecparam.h"
      integer varout
      varout = rinvmel
      return
      end
