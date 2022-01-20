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
      include "../commons/alpha.h"
      include "../commons/massthrs.h"
      include "../commons/nffn.h"
      include "../commons/nfmax.h"
      include "../commons/activeflavours.h"
      include "../commons/facscheme.h"
      include "../commons/renscheme.h"
      include "../commons/tecparam.h"
*
*     Default values
*
      IPT         = 1
      IPTALPHA    = 1
      NS          = "VFNS"
      FACSCHEME   = "MSBAR"
      RENSCHEME   = "MSBAR"
      NFMAX       = 8
      NFFN        = 8
      NLMAX       = 3
      NUMAX       = 2
      NDMAX       = 3      
      Q2REF       = 0.000510998928d0**2
      AREF        = 0.0072973525693d0
      NFMAXALPHA  = 8
      NFFNALPHA   = 8
      NLMAXAEM    = 3
      NUMAXAEM    = 2
      NDMAXAEM    = 3
      WAEM        = 0
      Q2TH(1)     = 0.000510998928d0**2
      Q2TH(2)     = 0.00216d0**2
      Q2TH(3)     = 0.00467D0**2
      Q2TH(4)     = 0.093d0**2
      Q2TH(5)     = 0.10566d0**2
      Q2TH(6)     = 1.27D0**2
      Q2TH(7)     = 1.77686d0**2
      Q2TH(8)     = 4.18d0**2
      Q2TH(9)     = 172.76d0**2
      NINT        = 300
C     NEXP        = 10
      NEXP        = 20
      NSTEP       = 10
C     NSTEP       = 20
      MINVMEL     = 51
C     MINVMEL     = 101
      RINVMEL     = 20
*
      return
      end
