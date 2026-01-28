************************************************************************
*
*     NDistributions.f:
*
*     Routine that returns the evolved N-space distributions.
*
************************************************************************
      subroutine NDistributions(N,nQ,Q,xfevN)
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/evscale.h"
      include "../commons/nf.h"
      include "../commons/evol.h"
      include "../commons/distf.h"
      include "../commons/pol.h"
**
*     Input Variables
*
      integer nQ
      double complex Q(100)
      double complex N
**
*     Internal Variables
*
      integer iQ
      integer i,j
* KL mod
      double complex aQCD
      double complex zfuncns(3),zfuncsg(2,2)
      double complex zfuncns15(2),zfuncns24(2),zfuncns35(2)
      double complex zfuncnsv15,zfuncnsv24,zfuncnsv35
      double complex xfph0N(-6:6),xfev0N(13)
      external ExternalSetMELA
**
*     Output Variables
*
      double complex xfevN(13)
*
*     Check
*
      if(nQ.lt.2)then
         write(6,*) "In NDistributions.f:"
         write(6,*) "The number of scales must be bigger or equal to 2"
         write(6,*) "nQ =",nQ
         call exit(-10)
      endif
*
*     Call N-space distributions
*
      if(distf(1:8).eq."internal")then
         if(evol.eq."SPACE")then
            if(pol.eq."OFF")then
               call toyLHPDFsn(N,xfph0N)    ! ToyLH PDFs
            else
               call toyLHPDFsPoln(N,xfph0N) ! ToyLH PDFs
            endif
         elseif(evol(1:4).eq."TIME")then
            call HKNSFFsn(N,xfph0N) ! HKNS for pi+ at NLO
         endif
*      elseif(distf.eq."g3S11v4")then
*         call g3S11v4(N,Q(1),xfph0N)
*          write(6,*) "v4; alphas(Q)=",12.56637d0*asEvolIni(1)         
      elseif(distf.eq."deltagNLO")then
         call GluonDeltaNLO(N,Q(1),xfph0N)
*         write(6,*) "NLO function! Q=", Q(1)
          write(6,*) "NLO; alphas(Q)=",12.56637d0*asEvolIni(1)
      elseif(distf(1:6).eq."deltag")then
         call GluonDelta(N,Q,xfph0N)
*         write(6,*) "LO function! alphas(Q)=",12.56637d0*asEvolIni(1)
*        write(6,*) "LO function! alphas(Q)=",12.56637d0*asEvolFi(1)
      elseif(distf.eq."constant")then
         call Constant(N,Q(1),xfph0N)
*         write(6,*) "Constant function! Q=", Q(1)
      elseif(distf.eq."pwave0")then
         call P0wave(N,Q(1),xfph0N)
*         write(6,*) "J=0! Q=", Q(1)
      elseif(distf.eq."pwave1")then
         call P1wave(N,Q(1),xfph0N)
*         write(6,*) "J=1! Q=", Q(1)
      elseif(distf.eq."pwave2")then
         call P2wave(N,Q(1),xfph0N)
*         write(6,*) "J=2! Q=", Q(1)
      elseif(distf.eq."g3s18")then
         call g3s18NLO(N,Q(1),xfph0N)
         write(6,*) "alphas(Q)=",12.56637d0*asEvolIni(1)
         write(6,*) "Ma function! Q=", Q(1)
*         write(6,*) "LO function! alphas(Q)=",12.56637d0*asEvolIni(1)
      elseif(distf.eq."g3s18MA")then
         call Ma3S18(N,Q(1),xfph0N)
         write(6,*) "Ma function! Q=", Q(1)
         write(6,*) "alphas(Q)=",12.56637d0*asEvolIni(1)

*GLUON FF      
      elseif(distf.eq."d0g3S11")then
         call d0g3S11(N,Q(1),xfph0N) 
      elseif(distf.eq."d2g3S11")then
         call d2g3S11(N,Q(1),xfph0N) 
      elseif(distf.eq."d4g3S11")then
         call d4g3S11(N,Q(1),xfph0N) 
      elseif(distf.eq."g3S11c0")then
         call g3S11c0(N,Q(1),xfph0N)  
      elseif(distf.eq."g3S11c1")then
         call g3S11c1(N,Q(1),xfph0N)    
      elseif(distf.eq."FENGg3S11")then
         call FENGg3S11(N,Q(1),xfph0N)             
*CHARM FF      
      elseif(distf.eq."c3S11LO")then
         call c3S11LO(N,Q(1),xfph0N)       
      elseif(distf.eq."c3S11REL")then
         call c3S11REL(N,Q(1),xfph0N) 
      elseif(distf.eq."c3S11PQQ")then
         call c3S11PQQ(N,Q(1),xfph0N)       
      elseif(distf.eq."c3S11NLO")then
         call c3S11NLO(N,Q(1),xfph0N) 
*ALL together
      elseif(distf.eq."FF3S11")then
         call FF3S11(N,Q(1),xfph0N)  
      elseif(distf.eq."FF3S11v4")then
         call FF3S11v4(N,Q(1),xfph0N)   
*MISC      
      elseif(distf.eq."TEST")then
         call TEST(N,Q(1),xfph0N) 
      elseif(distf.eq."TESTF")then
         call TESTF(N,Q(1),xfph0N) 
      elseif(distf.eq."TESTS")then
         call TESTS(N,Q(1),xfph0N) 
      elseif(distf.eq."XFitter")then
         call XFitterParametrization(N-1d0,xfph0N) ! XFitter PDFs
      elseif(distf(1:9).eq."ZeroScale")then
         call ZeroScalePDFs(N-1d0,xfph0N)
      elseif(distf.eq."external")then
         call ExternalSetMELA(N,xfph0N)
      else
         write(6,*) "Unknown input distributions, distf = ",distf
      endif
*      write(6,*) "alphas(Q=", Q(1), "= ", 12.56637d0*asEvolIni(1)
*
*     Rotate distribution into the evolution basis
*
      call lha2evlnc(xfph0N,xfev0N)
*
*     Loop over the scales
*
      do iQ=1,nQ-1
*
*     Initial and final scales
*
         Q20 = Q(iQ)**2d0
         Q2  = Q(iQ+1)**2d0
*
*     Values of alphas at the initial and final scales
*
* from XDist KL      
         asEvolIni(1) = aQCD(Q(1)**2d0)
         asEvolFin(1) = aQCD(Q(2)**2d0)
  
         as0 = asEvolIni(iQ)
         asq = asEvolFin(iQ)
*         write(6,*) as0*12.56637d0
*         write(6,*) asq*12.56637d0
*
*     Call evolution kernels
*
         call zfunc(N,Q20,Q2,
     1              zfuncns,zfuncns15,zfuncns24,zfuncns35,
     2              zfuncsg,zfuncnsv15,zfuncnsv24,zfuncnsv35)
*
*     Convolute PDFs with evolution kernels
*
*     Singlet
         do i=1,2
            xfevN(i) = (0d0,0d0)
            do j=1,2
               xfevN(i) = xfevN(i) + zfuncsg(i,j) * xfev0N(j)
            enddo
         enddo
*     Valence
         xfevN(3)  = zfuncns(3) * xfev0N(3)
*     V3
         xfevN(4)  = zfuncns(2) * xfev0N(4)
*     V8
         xfevN(5)  = zfuncns(2) * xfev0N(5)
*     V15
         if(nfi.lt.4)then
            xfevN(6) = zfuncnsv15 * xfev0N(3)
         else
            xfevN(6) = zfuncnsv15 * xfev0N(6)
         endif
*     V24
         if(nfi.lt.5)then
            xfevN(7) = zfuncnsv24 * xfev0N(3)
         else
            xfevN(7) = zfuncnsv24 * xfev0N(7)
         endif
*     V35
         if(nfi.lt.6)then
            xfevN(8) = zfuncnsv35 * xfev0N(3)
         else
            xfevN(8) = zfuncnsv35 * xfev0N(8)
         endif
*     T3
         xfevN(9)  = zfuncns(1) * xfev0N(9)
*     T8
         xfevN(10) = zfuncns(1) * xfev0N(10)
*     T15
         if(nfi.lt.4)then
            xfevN(11) = (0d0,0d0)
            do j=1,2
               xfevN(11) = xfevN(11) + zfuncns15(j) * xfev0N(j)
            enddo
         else
            xfevN(11) = zfuncns15(1) * xfev0N(11)
         endif
*     T24
         if(nfi.lt.5)then
            xfevN(12) = (0d0,0d0)
            do j=1,2
               xfevN(12) = xfevN(12) + zfuncns24(j) * xfev0N(j)
            enddo
         else
            xfevN(12) = zfuncns24(1) * xfev0N(12)
         endif
*     T35
         if(nfi.lt.6)then
            xfevN(13) = (0d0,0d0)
            do j=1,2
               xfevN(13) = xfevN(13) + zfuncns35(j) * xfev0N(j)
            enddo
         else
            xfevN(13) = zfuncns35(1) * xfev0N(13)
         endif
*
*     Copy evolved distributions into the initial ones before
*     the next evolution.
*
         if(iQ.lt.nQ-1)then
            do i=1,13
               xfev0N(i) = xfevN(i)   
            enddo
         endif
      enddo

      return
      end
