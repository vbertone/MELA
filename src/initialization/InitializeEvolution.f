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
*     Alpha reference values
*
      write(6,*) "Alpha reference value:"
      write(6,"(a,f14.9,a)") "   Qref = ",sqrt(Q2REF)," GeV"
      write(6,"(a,f14.9)")   "   Alpha(Qref) = ",AREF
*
*     Quark are included in the evolution
*
      if (quarks) then
         write(6,*) "Quarks are included in the evolution"
      else
         write(6,*) "Quarks are NOT included in the evolution"
      endif
*
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
*     Maximum number of active flavours
*
      if(nfmax.ge.0.and.nfmax.le.9)then
         write(6,"(a,i1)") " Maximum number of active flavours = ",nfmax
      else
         write(6,*) "In InitializeEvolution.f:"
         write(6,*) "Invalid value for nffac =",nfmax
         call exit(-10)
      endif
*
      return
      end
