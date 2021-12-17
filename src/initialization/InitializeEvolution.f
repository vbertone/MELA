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
      include "../commons/ns.h"
      include "../commons/alpha.h"
      include "../commons/massthrs.h"
      include "../commons/nffn.h"
      include "../commons/nfmax.h"
      include "../commons/activeflavours.h"
      include "../commons/facscheme.h"
      include "../commons/renscheme.h"
**
*     Internal variables
*
      integer nf
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
*     Factorisation scheme
*
      if(FACSCHEME.eq."MSBAR")then
         write(6,*) "Factorisation scheme: MSbar"
      elseif(FACSCHEME.eq."DELTA")then
         write(6,*) "Factorisation scheme: Delta"
      else
         write(6,*) "In InitializeEvolution.f:"
         write(6,*) "Unknown factorisation scheme = ",FACSCHEME
         call exit(-10)
      endif
*
*     Renormalisation scheme
*
      if(RENSCHEME.eq."MSBAR")then
         write(6,*) "Renormalisation scheme: MSbar"
      elseif(RENSCHEME.eq."FIXED")then
         write(6,*) "Renormalisation scheme: alpha fixed"
      elseif(RENSCHEME.eq."ALPMZ")then
         write(6,*) "Renormalisation scheme: alpha(MZ)"
      else
         write(6,*) "In InitializeEvolution.f:"
         write(6,*) "Unknown renormalisation scheme = ",RENSCHEME
         call exit(-10)
      endif
*
*     Evolution scheme
*
      if(NS.eq."FFNS")then
         if(NFFN.lt.0.or.NFFN.gt.9)then
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
*     Maximum number of active flavours
*
      if(nfmax.ge.0.and.nfmax.le.9)then
         write(6,"(a,i1)") " Maximum number of active flavours = ",nfmax
      else
         write(6,*) "In InitializeEvolution.f:"
         write(6,*) "Invalid value for nfmax =",nfmax
         call exit(-10)
      endif
**            
*     Thresholds
*
      write(6,*) "Fermion thresholds:"
      write(6,"(a,f14.9,a)") "   me  =",dsqrt(q2th(1))," GeV"
      write(6,"(a,f14.9,a)") "   mu  =",dsqrt(q2th(2))," GeV"
      write(6,"(a,f14.9,a)") "   md  =",dsqrt(q2th(3))," GeV"
      write(6,"(a,f14.9,a)") "   ms  =",dsqrt(q2th(4))," GeV"
      write(6,"(a,f14.9,a)") "   mm  =",dsqrt(q2th(5))," GeV"
      write(6,"(a,f14.9,a)") "   mc  =",dsqrt(q2th(6))," GeV"
      write(6,"(a,f14.9,a)") "   mt  =",dsqrt(q2th(7))," GeV"
      write(6,"(a,f14.9,a)") "   mb  =",dsqrt(q2th(8))," GeV"
      write(6,"(a,f14.9,a)") "   mtp =",dsqrt(q2th(9))," GeV"
*
*     Check that thresholds are ordered
*
      do nf = 2, 9
         if (dsqrt(q2th(nf - 1)).gt.dsqrt(q2th(nf))) then
            write(6,*) "In InitializeEvolution.f:"
            write(6,*) "Fermion masses are not ordered"
            call exit(-10)
         endif
      enddo
*
*     Quark are included in the evolution
*
      if (quarks) then
         write(6,*) "Quarks are included in the evolution"
      else
         write(6,*) "Quarks are NOT included in the evolution"
      endif
*      
*     Alpha reference values
*
      if(RENSCHEME.ne."MSBAR")then
         write(6,*) "Alpha FIXED with reference value:"
         write(6,"(a,f14.9)")   "   Alpha(Qref) = ",AREF         
      else
         write(6,*) "Alpha reference value:"
         write(6,"(a,f14.9,a)") "   Qref = ",sqrt(Q2REF)," GeV"
         write(6,"(a,f14.9)")   "   Alpha(Qref) = ",AREF
         if(NS.eq."FFNS")then
            if(NFFNalpha.lt.0.or.NFFNalpha.gt.9)then
               write(6,*)"In InitializeEvolution.f:"
               write(6,*)"NFFNalpha out or range, NFFNalpha =",
     .              NFFNalpha
               call exit(-10)
            endif
            write(6,"(a,i1)") " Alpha evolution: FFNS with NF = ",
     .           NFFNalpha
         elseif(NS.eq."VFNS")then
            write(6,*) "Alpha evolution: VFNS"
         else
            write(6,*) "In InitializeEvolution.f:"
            write(6,*) "Unknown mass scheme = ",NS
            call exit(-10)
         endif
         if(nfmaxalpha.ge.0.and.nfmaxalpha.le.9)then
            write(6,"(a,i1)") " Max.num. active flav in alpha = ",
     .           nfmaxalpha
         else
            write(6,*) "In InitializeEvolution.f:"
            write(6,*) "Invalid value for nfmaxalpha =",nfmaxalpha
            call exit(-10)
         endif
         if (quarksalpha) then
            write(6,*) "Quarks are included in the evolution of alpha"
         else
            write(6,*) "Quarks NOT included in the evolution of alpha"
         endif
      endif
*     
      return
      end
