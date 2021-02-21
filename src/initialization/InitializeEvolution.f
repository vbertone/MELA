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
      include "../commons/alphas.h"
      include "../commons/hqmass.h"
      include "../commons/massthrs.h"
      include "../commons/nffn.h"
      include "../commons/nfmax.h"
      include "../commons/renfacscales.h"
      include "../commons/evol.h"
      include "../commons/process.h"
      include "../commons/coeffhqmellin.h"
      include "../commons/minimax.h"
      include "../commons/distf.h"
      include "../commons/pol.h"
*
*     Initialize the coefficient of the beta functions in QCD and QED 
*     and the coefficient of the QCD gamma function
*
      call initBeta
*
*     Alphas at Q20 and at the heavy quark mass thresholds
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
         if(distf(1:8).eq."internal")then
            write(6,*) "Initial conditions: toy Les Houches PDF",
     1                 " (Q0 = sqrt(2) GeV)"
         elseif(distf.eq."XFitter")then
            write(6,*) "Initial conditions: XFitter PDFs"
            write(6,*) "(These input distubutions can be used only ",
     1                 "with XFitter)"
         elseif(distf(1:9).eq."ZeroScale")then
            write(6,*) "Initial conditions: Zero scale PDFs"
         endif
      elseif(evol(1:4).eq."TIME")then
         write(6,*) "Time-like evolution (FFs)"
         if(distf(1:8).eq."internal")then
            write(6,*) "Initial conditions: HKNS 2007",
     1                 " for pi+ at NLO (Q0 = 1 GeV)"
         elseif(distf.eq."XFitter")then
            write(6,*) "Initial conditions: XFitter PDFs"
            write(6,*) "(These input distubutions can be used only ",
     1                 "with XFitter)"
         endif
      else
         write(6,*) "In InitializeEvolution.f:"
         write(6,*) "Invalid evolution mode, evol = ",evol
         call exit(-10)
      endif
*
*     Polarized or upolarized evolution
*
      if(pol.eq."OFF")then
         write(6,*) "Unpolarized evolution"
      elseif(pol(1:2).eq."ON")then
         write(6,*) "Polarized evolution"
         if(evol(1:4).eq."TIME")then
            write(6,*) "In InitializeEvolution.f:"
            write(6,*) "Time-like polarized evolution not",
     1                 " yet available"
            call exit(-10)
         endif
      else
         write(6,*) "In InitializeEvolution.f:"
         write(6,*) "Invalid polarization switch, pol = ",pol
         call exit(-10)
      endif
*
*     Perturbative order
*
      if(ipt.eq.0)then
         write(6,*) "Perturbative order: LO"
      elseif(ipt.eq.1)then
         write(6,*) "Perturbative order: NLO"
      elseif(ipt.eq.2)then
         write(6,*) "Perturbative order: NNLO"
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
         if(NFFN.lt.3.or.NFFN.gt.6)then
            write(6,*)"In InitializeEvolution.f:"
            write(6,*)"NFFN out or range, NFFN =",NFFN
            call exit(-10)
         endif
         write(6,"(a,i1)") " Evolution scheme: FFNS with NF = ",NFFN
      elseif(NS.eq."VFNS")then
         write(6,*) "Evolution scheme: VFNS"
         if(ipt.eq.2.and.evol(1:4).eq."TIME")then
            write(6,*) "In InitializeEvolution.f:"
            write(6,*) "Time-like evolution at NNLO in the VFNS",
     1                 " not available"
            call exit(-10)
         endif
      else
         write(6,*) "In InitializeEvolution.f:"
         write(6,*) "Unknown mass scheme = ",NS
         call exit(-10)
      endif
*
*     Alphas reference values
*
      write(6,*) "Alphas reference value:"
      write(6,"(a,f8.3,a,f8.3,a)") "   Qref = (",real(sqrt(Q2REF)),",",
     1                             imag(sqrt(Q2REF))," ) GeV"
      write(6,"(a,f7.3,a,f7.3,a)") "   Alpha_s(Qref) = (",
     1                             real(ASREF),",",imag(ASREF)," )"
*
*     MSbar or Pole masses reference values
*
      if(hqmass.eq.0)then
         write(6,*) "Pole masses chosen"
      elseif(hqmass.eq.1)then
         write(6,*) "MSbar masses chosen"
      else
         write(6,*) "In InitializeEvolution.f:"
         write(6,*) "Unknown value of HQMASS =",hqmass
         call exit(-10)
      endif
*
      write(6,*) "Heavy quark thresholds:"
      write(6,"(a,f8.3,a)") "   mc =",dsqrt(real(q2th(4)))," GeV"
      write(6,"(a,f8.3,a)") "   mb =",dsqrt(real(q2th(5)))," GeV"
      write(6,"(a,f8.3,a)") "   mt =",dsqrt(real(q2th(6)))," GeV"
*
*     Maximum number of active flavours
*
      if(nfmax.ge.3.and.nfmax.le.6)then
         write(6,"(a,i1)") " Maximum number of active flavours = ",nfmax
      else
         write(6,*) "In InitializeEvolution.f:"
         write(6,*) "Invalid value for nffac =",nfmax
         call exit(-10)
      endif
*
*     Factorization/renormalization scales (QCD and QED)
*
      write(6,"(a,f7.4)") " muR / muF = ",KRF
      write(6,*) "  "
*
      write(6,*) "========== DIS parameters =========="
      write(6,*) "  "
      if(proc.eq."EM")then
         write(6,*) "Electromagnetic process"
      elseif(proc.eq."NC")then
         write(6,*) "Neutral current process"
      else
         write(6,*) "In InitializeEvolution.f:"
         write(6,*) "Invalid DIS process, proc = ",proc
         call exit(-10)
      endif
*
      write(6,"(a,f7.4)") " muR / Q = ",KRENQ
      write(6,"(a,f7.4)") " muF / Q = ",KFACQ
      write(6,*) "  "
*
      return
      end
