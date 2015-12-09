************************************************************************
*
*     ReadParameters.f:
*
*     It read the parameters of the evolution from the input card.
*
************************************************************************
      subroutine ReadParameters(card)
*
      implicit none
*
      include "../commons/ipt.h"
      include "../commons/modev.h"
      include "../commons/ns.h"
      include "../commons/alphas.h"
      include "../commons/hqmass.h"
      include "../commons/massthrs.h"
      include "../commons/nffn.h"
      include "../commons/nfmax.h"
      include "../commons/renfacscales.h"
      include "../commons/evol.h"
      include "../commons/process.h"
      include "../commons/distf.h"
      include "../commons/pol.h"
**
*     Input Variables
*
      character*100 card
**
*     Internal Variables
*
      integer lp,lu
      double precision mc,mb,mt
      double precision eps
      double precision ReQR,ImQR
      double precision ReASR,ImASR
      character*6  fhqmass
      character*50 str
      parameter(eps=1d-10)
*
*     Overwritten by the values read from the input card if found.
*
      IPT     = 2
      MODEV   = "ITE"
      NS      = "VFNS"
      NFMAX   = 6
      NFFN    = 3
      FHQMASS = "POLE"
      KRF     = 1d0
      EVOL    = "SPACE"
      POL     = "OFF"
      ReQR    = dsqrt(2d0) - eps
      ImQR    = 0d0
      ReASR   = 0.35d0
      ImASR   = 0d0
      MC      = dsqrt(2d0)
      MB      = 4.5d0
      MT      = 175d0
      PROC    = "EM"
      KRENQ   = 1d0
      KFACQ   = 1d0
      DISTF   = "internal"
*
c      open(unit=10,status="old",file="../run/"//card)
      open(unit=10,status="old",file=card)
      do
         read(10,"(a)",end=101) str
         if(str(1:1).ne."#")then
            lp = index(str," ") - 1
            lu = index(str,"=") + 1
            if(str(1:lp).eq."IPT")     read(str(lu:50),*) IPT
            if(str(1:lp).eq."MODEV")   read(str(lu:50),*) MODEV
            if(str(1:lp).eq."NS")      read(str(lu:50),*) NS
            if(str(1:lp).eq."NFMAX")   read(str(lu:50),*) NFMAX
            if(str(1:lp).eq."NFFN")    read(str(lu:50),*) NFFN
            if(str(1:lp).eq."HQMASS")  read(str(lu:50),*) FHQMASS
            if(str(1:lp).eq."KRF")     read(str(lu:50),*) KRF
            if(str(1:lp).eq."EVOL")    read(str(lu:50),*) EVOL
            if(str(1:lp).eq."POL")     read(str(lu:50),*) POL
            if(str(1:lp).eq."QREF")    read(str(lu:50),*) ReQR,ImQR
            if(str(1:lp).eq."ASREF")   read(str(lu:50),*) ReASR,ImASR
            if(str(1:lp).eq."MC")      read(str(lu:50),*) MC
            if(str(1:lp).eq."MB")      read(str(lu:50),*) MB
            if(str(1:lp).eq."MT")      read(str(lu:50),*) MT
            if(str(1:lp).eq."PROC")    read(str(lu:50),*) PROC
            if(str(1:lp).eq."KRENQ")   read(str(lu:50),*) KRENQ
            if(str(1:lp).eq."KFACQ")   read(str(lu:50),*) KFACQ
            if(str(1:lp).eq."DISTF")   read(str(lu:50),*) DISTF
         endif
      enddo
 101  close(10)
*
*     Some redefinitions
*
      Q2REF    = DCMPLX(ReQR,ImQR)**2d0
      ASREF    = DCMPLX(ReASR,ImASR)
      Q2TH(4)  = DCMPLX(MC**2d0,0d0)
      Q2TH(5)  = DCMPLX(MB**2d0,0d0)
      Q2TH(6)  = DCMPLX(MT**2d0,0d0)
*
      if(FHQMASS.EQ."POLE")  HQMASS = 0
      if(FHQMASS.EQ."MSBAR") HQMASS = 1
*
*     If KRENQ or KFACQ are different from 1, force KRF to be KRENQ/KFACQ
*
      if(KRENQ.ne.1d0.or.KFACQ.ne.1d0) KRF = KRENQ / KFACQ
*
      return
      end
