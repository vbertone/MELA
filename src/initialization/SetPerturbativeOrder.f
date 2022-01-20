************************************************************************
*
*     SetPerturbativeOrder.f:
*
*     Sets the perturbative order of the evolution
*
************************************************************************
      subroutine SetPerturbativeOrder(iptin)
      implicit none
      include "../commons/ipt.h"
      integer iptin
      ipt = iptin
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
************************************************************************            
      subroutine SetPerturbativeOrderAlpha(iptalphain)
      implicit none
      include "../commons/ipt.h"
      integer iptalphain
      iptalpha = iptalphain
      return
      end
************************************************************************      
      subroutine GetPerturbativeOrderAlpha(iptalphaout)
      implicit none
      include "../commons/ipt.h"
      integer iptalphaout
      iptalphaout = iptalpha
      return
      end

      
