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
      double complex zfuncns(3),zfuncsg(2,2)
      double complex zfuncns15(2),zfuncns24(2),zfuncns35(2)
      double complex zfuncnsv15,zfuncnsv24,zfuncnsv35
      double complex xfph0N(-6:6),xfev0N(13)
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
      elseif(distf.eq."XFitter")then
         call XFitterParametrization(N-1d0,xfph0N) ! XFitter PDFs
      elseif(distf(1:9).eq."ZeroScale")then
         call ZeroScalePDFs(N-1d0,xfph0N)
      else
         write(6,*) "Unknown input distributions, distf = ",distf
      endif
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
         as0 = asEvolIni(iQ)
         asq = asEvolFin(iQ)
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
*
      return
      end
