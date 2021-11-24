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
      include "../commons/modev.h"
      include "../commons/evol.h"
*
*     Default values
*
      IPT     = 1
      MODEV   = "PTH"
      NS      = "VFNS"      
      NFMAX   = 3
      NFFN    = 3
      EVOL    = "SPACE"
      Q2REF   = dcmplx(0.000510998928d0**2,0d0)
      AREF    = dcmplx(0.0072973525693d0,0d0)
      Q2TH(1) = 0.000510998928d0**2
      Q2TH(2) = 0.10566d0**2
      Q2TH(3) = 1.77686d0**2
*
      return
      end
