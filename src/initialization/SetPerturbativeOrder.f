************************************************************************
*
*     SetPerturbativeOrder.f:
*
*     Sets the perturbative order of the evolution
*
************************************************************************
      subroutine SetPerturbativeOrder(iptin)
*
      implicit none
*
      include "../commons/ipt.h"
**
*     Input Variables
*
      integer iptin
*
      ipt = iptin
*
      return
      end
************************************************************************      
      subroutine GetPerturbativeOrder(iptout)
      implicit none
      include "../commons/ipt.h"
      integer iptout
      iptout = ipt
      return
      end
