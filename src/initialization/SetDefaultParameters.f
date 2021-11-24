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
*
*     Default values
*
      IPT     = 1
      NS      = "VFNS"
      NFMAX   = 9
      NFFN    = 8
      Q2REF   = 0.000510998928d0**2
      AREF    = 0.0072973525693d0
      QUARKS  = .TRUE.
      Q2TH(1) = 0.000510998928d0**2
      Q2TH(2) = 0.00216d0**2
      Q2TH(3) = 0.00467D0**2
      Q2TH(4) = 0.093d0**2
      Q2TH(5) = 0.10566d0**2
      Q2TH(6) = 1.27D0**2
      Q2TH(7) = 1.77686d0**2
      Q2TH(8) = 4.18d0**2
      Q2TH(9) = 172.76d0**2
*
      return
      end
