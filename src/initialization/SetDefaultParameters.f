************************************************************************
*
*     SetDefaultParameters.f:
*
*     It sets the default parameters of the evolution.
*
************************************************************************
      subroutine SetDefaultParameters
*
      implicit none
*
      include "../commons/ipt.h"
      include "../commons/ns.h"
      include "../commons/facscheme.h"
      include "../commons/renscheme.h"
      include "../commons/activeflavours.h"
      include "../commons/alpha.h"      
      include "../commons/massthrs.h"
      include "../commons/tecparam.h"
*
*     Default values
*
      IPT         = 1
      IPTALPHA    = 1
      NS          = "VFNS"
      FACSCHEME   = "MSBAR"
      RENSCHEME   = "MSBAR"
      NLMAX       = 3
      NUMAX       = 2
      NDMAX       = 3
      NLMAXAEM    = 3
      NUMAXAEM    = 2
      NDMAXAEM    = 3
      WAEM        = 0
      Q2REF       = 0.000510998928d0**2
      AREF        = 0.0072973525693d0     
      Q2TH(1)     = 0.000510998928d0**2
      Q2TH(2)     = 0.00216d0**2
      Q2TH(3)     = 0.00467D0**2
      Q2TH(4)     = 0.093d0**2
      Q2TH(5)     = 0.10566d0**2
      Q2TH(6)     = 1.27D0**2
      Q2TH(7)     = 1.77686d0**2
      Q2TH(8)     = 4.18d0**2
      Q2TH(9)     = 80.379d0**2
      Q2TH(10)    = 91.1876d0**2
      Q2TH(11)    = 172.76d0**2
      NINT        = 300
      NMINSTEP    = 50
      NEXP        = 30
      NSTEPAEM    = 10
      MINVMEL     = 51
      RINVMEL     = 20
C     ALPXXSOL    = "PATHOR"
      ALPXXSOL    = "MAGNUS"      
*
      return
      end
