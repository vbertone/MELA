************************************************************************
*
*     Getters.f
*
************************************************************************
      subroutine GetPerturbativeOrder(iptout)
      implicit none
      include "../commons/ipt.h"
      integer iptout
      iptout = ipt
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
************************************************************************ 
      subroutine GetFlavourSchemeInt(nsintout)
      implicit none
      include "../commons/ns.h"
      integer nsintout
      if (ns.eq."FFNS") then
         nsintout = 0
      elseif (ns.eq."VFNS") then
         nsintout = 1
      endif
      return
      end
************************************************************************       
      subroutine GetFactorisationSchemeInt(fsintout)
      implicit none
      include "../commons/facscheme.h"      
      integer fsintout
      if (facscheme.eq."MSBAR") then
         fsintout = 0
      elseif (facscheme.eq."DELTA") then
         fsintout = 1
      endif
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
      elseif (renscheme.eq."ALGMU") then
         rsintout = 3
      elseif (renscheme.eq."AFAKE") then
         rsintout = 4         
      endif
      return
      end
************************************************************************      
      subroutine GetActiveFlav(nlmaxout,numaxout,ndmaxout)
      implicit none
      include "../commons/activeflavours.h"
      integer nlmaxout,numaxout,ndmaxout
      nlmaxout = nlmax
      numaxout = numax
      ndmaxout = ndmax
      return
      end
************************************************************************
      subroutine GetActiveFlavAem(nlmaxaemout,numaxaemout,ndmaxaemout)
      implicit none
      include "../commons/activeflavours.h"
      integer nlmaxaemout,numaxaemout,ndmaxaemout
      nlmaxaemout = nlmaxaem
      numaxaemout = numaxaem
      ndmaxaemout = ndmaxaem
      return
      end
************************************************************************
      subroutine GetWaem(waemout)
      implicit none
      include "../commons/activeflavours.h"
      integer waemout
      waemout = waem
      return
      end
************************************************************************ 
      subroutine GetAlphaRef(arefout)
      implicit none
      include "../commons/alpha.h"      
      double precision arefout
      arefout = aref
      return
      end
************************************************************************      
      subroutine GetAlphaQref(Qrefout)
      implicit none
      include "../commons/alpha.h"      
      double precision Qrefout
      Qrefout = dsqrt(Q2ref)
      return
      end
************************************************************************            
      subroutine GetThresholds2(q2thrs)
      implicit none
      include "../commons/massthrs.h"      
      double precision q2thrs(11)
      q2thrs=q2th
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
      subroutine getnminstep(varout)
      implicit none
      include "../commons/tecparam.h"
      integer varout
      varout = nminstep
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
      subroutine getnstepaem(varout)
      implicit none
      include "../commons/tecparam.h"
      integer varout
      varout = nstepaem
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
      subroutine getrinvmel(varout)
      implicit none
      include "../commons/tecparam.h"
      integer varout
      varout = rinvmel
      return
      end
************************************************************************      
      subroutine Getb0(b0)
      implicit none
      include "../commons/beta.h"
      include "../commons/consts.h"      
      double precision b0(1:11)
      b0(:) = -beta0(:)/4d0/pi
      return
      end      
************************************************************************
      subroutine Getb1(b1)
      implicit none
      include "../commons/beta.h"
      include "../commons/consts.h"      
      double precision b1(1:11)
      b1(:) = -beta1(:)/4d0/4d0/pi/pi
      return
      end      
************************************************************************
      subroutine GetC2(C2)
      implicit none
      include "../commons/nfsum.h"
      double precision C2(1:11)
      C2(:) = nfsum2(:)
      return
      end
************************************************************************
      subroutine GetC4(C4)
      implicit none
      include "../commons/nfsum.h"
      double precision C4(1:11)
      C4(:) = nfsum4(:)
      return
      end
************************************************************************
      subroutine GetMZ2(mz2out)      
      implicit none
      include "../commons/massthrs.h"      
      double precision mz2out
      mz2out=q2th(10)
      return
      end
************************************************************************
      subroutine GetMW2(mw2out)      
      implicit none
      include "../commons/massthrs.h"      
      double precision mw2out
      mw2out=q2th(9)
      return
      end
************************************************************************
      subroutine Geta0(a0)
      implicit none
      include "../commons/alpha.h"
      double precision a0
      a0 = ath(1)
      return
      end

      
      
