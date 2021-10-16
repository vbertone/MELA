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
      include "../commons/alpha.h"
      include "../commons/massthrs.h"
      include "../commons/nffn.h"
      include "../commons/nfmax.h"
      include "../commons/evol.h"
**
*     Input Variables
*
      character*100 card
**
*     Internal Variables
*
      integer lp,lu
      double precision me,mm,mt
      double precision ReQR,ImQR
      double precision ReAR,ImAR
      character*50 str
*
*     Overwritten by the values read from the input card if found.
*
      IPT   = 1
      MODEV = "ITE"
      NS    = "VFNS"
      NFMAX = 6
      NFFN  = 3
      EVOL  = "SPACE"
      ReQR  = 0.000510998928d0
      ImQR  = 0d0
      ReAR  = 0.0072973525693d0
      ImAR  = 0d0
      ME    = 0.000510998928d0
      MM    = 0.10566d0
      MT    = 1.777d0
*
c      open(unit=10,status="old",file="../run/"//card)
      open(unit = 10, status = "old", file = card)
      do
         read(10,"(a)",end=101) str
         if(str(1:1).ne."#")then
            lp = index(str," ") - 1
            lu = index(str,"=") + 1
            if(str(1:lp).eq."IPT")   read(str(lu:50),*) IPT
            if(str(1:lp).eq."MODEV") read(str(lu:50),*) MODEV
            if(str(1:lp).eq."NS")    read(str(lu:50),*) NS
            if(str(1:lp).eq."NFMAX") read(str(lu:50),*) NFMAX
            if(str(1:lp).eq."NFFN")  read(str(lu:50),*) NFFN
            if(str(1:lp).eq."EVOL")  read(str(lu:50),*) EVOL
            if(str(1:lp).eq."QREF")  read(str(lu:50),*) ReQR,ImQR
            if(str(1:lp).eq."AREF")  read(str(lu:50),*) ReAR,ImAR
            if(str(1:lp).eq."ME")    read(str(lu:50),*) ME
            if(str(1:lp).eq."MM")    read(str(lu:50),*) MM
            if(str(1:lp).eq."MT")    read(str(lu:50),*) MT
         endif
      enddo
 101  close(10)
*
*     Some redefinitions
*
      Q2REF   = DCMPLX(ReQR, ImQR)**2d0
      AREF    = DCMPLX(ReAR, ImAR)
      Q2TH(1) = DCMPLX(ME**2d0, 0d0)
      Q2TH(2) = DCMPLX(MM**2d0, 0d0)
      Q2TH(3) = DCMPLX(MT**2d0, 0d0)
*
      return
      end
