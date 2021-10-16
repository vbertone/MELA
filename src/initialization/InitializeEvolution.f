************************************************************************
*
*     InitializeEvolution.f:
*
*     Initialization for the production of the evolution tables.
*
************************************************************************
      subroutine InitializeEvolution
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
*
*     Initialize the coefficient of the beta functions in QED
*
      call initBeta
*
*     Alpha at Q20 and at the heavy lepton mass thresholds
*
      call initCouplings
*
*     Evolution mode
*
      write(6,*) "  "
      write(6,*) "======= Evolution parameters ======="
      write(6,*) "  "
      if(evol.eq."SPACE")then
         write(6,*) "Space-like evolution (PDFs)"
      elseif(evol(1:4).eq."TIME")then
         write(6,*) "Time-like evolution (FFs)"
      else
         write(6,*) "In InitializeEvolution.f:"
         write(6,*) "Invalid evolution mode, evol = ",evol
         call exit(-10)
      endif
*
*     Perturbative order
*
      if(ipt.eq.0)then
         write(6,*) "Perturbative order: LO"
      elseif(ipt.eq.1)then
         write(6,*) "Perturbative order: NLO"
      else
         write(6,*) "In InitializeEvolution.f:"
         write(6,*) "Invalid perturbative order, ipt = ",ipt
         call exit(-10)
      endif
*
*     Mode of solution of PT evolution equations
*
      if(MODEV.eq."TRN")then
         write(6,*) "Truncated solution chosen"
      elseif(MODEV.EQ."ITE")then
         write(6,*) "Iterated solution chosen"
      elseif(MODEV.EQ."PTH")then
         write(6,*) "Path-ordering solution chosen"
      elseif(MODEV.EQ."GFN")then
         write(6,*) "g-functions solution chosen"
      else
         write(6,*) "In InitializeEvolution.f:"
         write(6,*) "Unknown solution mode, MODEV = ",MODEV
         call exit(-10)
      endif
*
*     Evolution scheme
*
      if(NS.eq."FFNS")then
         if(NFFN.lt.0.or.NFFN.gt.3)then
            write(6,*)"In InitializeEvolution.f:"
            write(6,*)"NFFN out or range, NFFN =",NFFN
            call exit(-10)
         endif
         write(6,"(a,i1)") " Evolution scheme: FFNS with NF = ",NFFN
      elseif(NS.eq."VFNS")then
         write(6,*) "Evolution scheme: VFNS"
      else
         write(6,*) "In InitializeEvolution.f:"
         write(6,*) "Unknown mass scheme = ",NS
         call exit(-10)
      endif
*
*     Alpha reference values
*
      write(6,*) "Alpha reference value:"
      write(6,"(a,f14.9,a,f14.9,a)")"   Qref = (",real(sqrt(Q2REF)),",",
     1                             imag(sqrt(Q2REF))," ) GeV"
      write(6,"(a,f14.9,a,f14.9,a)")"   Alpha(Qref) = (",
     1                             real(AREF),",",imag(AREF)," )"
*
      write(6,*) "Lepton thresholds:"
      write(6,"(a,f14.9,a)") "   me =",dsqrt(real(q2th(1)))," GeV"
      write(6,"(a,f14.9,a)") "   mm =",dsqrt(real(q2th(2)))," GeV"
      write(6,"(a,f14.9,a)") "   mt =",dsqrt(real(q2th(3)))," GeV"
*
*     Maximum number of active flavours
*
      if(nfmax.ge.0.and.nfmax.le.3)then
         write(6,"(a,i1)") " Maximum number of active flavours = ",nfmax
      else
         write(6,*) "In InitializeEvolution.f:"
         write(6,*) "Invalid value for nffac =",nfmax
         call exit(-10)
      endif
*
      return
      end
