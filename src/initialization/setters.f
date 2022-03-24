************************************************************************
*
*     Setters.f
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
      subroutine SetPerturbativeOrderAlpha(iptalphain)
      implicit none
      include "../commons/ipt.h"
      integer iptalphain
      iptalpha = iptalphain
      return
      end
************************************************************************            
      subroutine SetFlavourScheme(nsin)
      implicit none
      include "../commons/ns.h"
      character*4 nsin
      ns = nsin
      return
      end
************************************************************************
      subroutine SetFlavourSchemeInt(nsinint)
      implicit none
      include "../commons/ns.h"
      integer nsinint
      if (nsinint.eq.0) then
         ns = "FFNS"
      elseif (nsinint.eq.1) then
         ns = "VFNS"
      else
         write(6,*) "In SetFlavourScheme.f:"
         write(6,*) "Invalid value"
         call exit(-10)
      endif 
      return
      end
************************************************************************
      subroutine SetFactorisationScheme(fsin)
      implicit none
      include "../commons/facscheme.h"
      character*5 fsin
      facscheme = fsin
      return
      end
************************************************************************
      subroutine SetFactorisationSchemeInt(fsinint)
      implicit none
      include "../commons/facscheme.h"
      integer fsinint
      if (fsinint.eq.0) then
         facscheme = "MSBAR"
      elseif (fsinint.eq.1) then
         facscheme = "DELTA"
      else
         write(6,*) "In SetFactorisationScheme.f:"
         write(6,*) "Invalid value"
         call exit(-10)
      endif 
      return
      end
************************************************************************      
      subroutine SetRenormalisationScheme(rsin)
      implicit none
      include "../commons/renscheme.h"
      character*5 rsin
      renscheme = rsin
      return
      end
************************************************************************
      subroutine SetRenormalisationSchemeInt(rsinint)
      implicit none
      include "../commons/renscheme.h"
      integer rsinint
      if (rsinint.eq.0) then
         renscheme = "MSBAR"
      elseif (rsinint.eq.1) then
         renscheme = "FIXED"
      elseif (rsinint.eq.2) then
         renscheme = "ALPMZ"
      elseif (rsinint.eq.3) then
         renscheme = "ALGMU"
      elseif (rsinint.eq.4) then
         renscheme = "AFAKE"
      else
         write(6,*) "In SetRenormalisationScheme.f:"
         write(6,*) "Invalid value"
         call exit(-10)
      endif 
      return
      end
************************************************************************      
      subroutine SetActiveFlav(nlmaxin,numaxin,ndmaxin)
      implicit none
      include "../commons/activeflavours.h"
      integer nlmaxin,numaxin,ndmaxin
      nlmax = nlmaxin
      numax = numaxin
      ndmax = ndmaxin
      return
      end
************************************************************************
      subroutine SetActiveFlavAem(nlmaxaemin,numaxaemin,ndmaxaemin)
      implicit none
      include "../commons/activeflavours.h"
      integer nlmaxaemin,numaxaemin,ndmaxaemin
      nlmaxaem = nlmaxaemin
      numaxaem = numaxaemin
      ndmaxaem = ndmaxaemin
      return
      end
************************************************************************      
      subroutine SetWaem(waemin)
      implicit none
      include "../commons/activeflavours.h"
      integer waemin
      waem = waemin
      return
      end
************************************************************************            
      subroutine SetAlpha(ain,Qin)
      implicit none
      include "../commons/alpha.h"
      double precision ain, Qin
      aref = ain
      Q2ref = Qin**2
      return
      end
************************************************************************
      subroutine SetThresholds(me,mu,md,ms,mm,mc,mt,mb,MW,MZ,mtp)
      implicit none
      include "../commons/massthrs.h"
      double precision me,mu,md,ms,mm,mc,mt,mb,MW,MZ,mtp
      Q2TH(1) = ME**2
      Q2TH(2) = MU**2
      Q2TH(3) = MD**2
      Q2TH(4) = MS**2
      Q2TH(5) = MM**2
      Q2TH(6) = MC**2
      Q2TH(7) = MT**2
      Q2TH(8) = MB**2
      Q2TH(9) = MW**2
      Q2TH(10)= MZ**2
      Q2TH(11) = MTP**2
      return
      end
************************************************************************
      subroutine setnint(varin)
      implicit none
      include "../commons/tecparam.h"
      integer varin
      nint = varin
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
      subroutine setnstepaem(varin)
      implicit none
      include "../commons/tecparam.h"
      integer varin
      nstepaem = varin
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
      subroutine setrinvmel(varin)
      implicit none
      include "../commons/tecparam.h"
      integer varin
      rinvmel = varin
      return
      end
************************************************************************      
      subroutine setalpxxsolint(intin)
      implicit none
      include "../commons/tecparam.h"
      integer intin
      if (intin.eq.0) then
         alpxxsol = "PATHOR"
      elseif (intin.eq.1) then
         alpxxsol = "MAGNUS"
      else
         write(6,*) "In setters.f:"
         write(6,*) "Invalid value"
         call exit(-10)
      endif 
      return
      end
