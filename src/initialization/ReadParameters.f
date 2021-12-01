************************************************************************
*
*     ReadParameters.f:
*
*     It reads the parameters of the evolution from the input card.
*
************************************************************************
      subroutine ReadParameters(card)
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
**
*     Input Variables
*
      character*100 card
**
*     Internal Variables
*
      integer lp,lu
      double precision me,mu,md,ms,mm,mc,mt,mb,mtp
      double precision QR
      character*50 str
*
*     Set default parameteres that are overwritten by the values read
*     from the input card if found.
*
      call SetDefaultParameters
*
      open(unit = 10, status = "old", file = card)
      do
         read(10,"(a)",end=101) str
         if(str(1:1).ne."#")then
            lp = index(str," ") - 1
            lu = index(str,"=") + 1
            if(str(1:lp).eq."IPT")       read(str(lu:50),*) IPT
            if(str(1:lp).eq."NS")        read(str(lu:50),*) NS
            if(str(1:lp).eq."FACSCHEME") read(str(lu:50),*) FACSCHEME
            if(str(1:lp).eq."NFMAX")     read(str(lu:50),*) NFMAX
            if(str(1:lp).eq."NFFN")      read(str(lu:50),*) NFFN
            if(str(1:lp).eq."QREF")      read(str(lu:50),*) QR
            if(str(1:lp).eq."AREF")      read(str(lu:50),*) AREF
            if(str(1:lp).eq."QUARKS")    read(str(lu:50),*) QUARKS
            if(str(1:lp).eq."ME")        read(str(lu:50),*) ME
            if(str(1:lp).eq."MU")        read(str(lu:50),*) MU
            if(str(1:lp).eq."MD")        read(str(lu:50),*) MD
            if(str(1:lp).eq."MS")        read(str(lu:50),*) MS
            if(str(1:lp).eq."MM")        read(str(lu:50),*) MM
            if(str(1:lp).eq."MC")        read(str(lu:50),*) MC
            if(str(1:lp).eq."MT")        read(str(lu:50),*) MT
            if(str(1:lp).eq."MB")        read(str(lu:50),*) MB
            if(str(1:lp).eq."MTP")       read(str(lu:50),*) MTP
         endif
      enddo
 101  close(10)
*
*     Some redefinitions
*
      Q2REF   = QR**2
      Q2TH(1) = ME**2
      Q2TH(2) = MU**2
      Q2TH(3) = MD**2
      Q2TH(4) = MS**2
      Q2TH(5) = MM**2
      Q2TH(6) = MC**2
      Q2TH(7) = MT**2
      Q2TH(8) = MB**2
      Q2TH(9) = MTP**2
*
      return
      end
