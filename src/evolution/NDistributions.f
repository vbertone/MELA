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
      include "../commons/alpha.h"
      include "../commons/evscale.h"
      include "../commons/nf.h"
      include "../commons/evol.h"
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
      double complex zfuncns3(2),zfuncns8(2)
      double complex zfuncnsv3,zfuncnsv8
      double complex xfph0N(-3:3),xfev0N(7)
      external ExternalSetMELA
**
*     Output Variables
*
      double complex xfevN(7)
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
      if(evol.eq."SPACE")then
         call electronPDFsn(N, xfph0N) ! Electron PDFs
      elseif(evol(1:4).eq."TIME")then
         call electronFFsn(N, xfph0N) ! Electron FFs
      endif
*
*     Rotate distribution into the evolution basis
*
      call lha2evlnc(xfph0N, xfev0N)
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
*     Values of alpha at the initial and final scales
*
         a0 = aEvolIni(iQ)
         aq = aEvolFin(iQ)
*
*     Call evolution kernels
*
         call zfunc(N,Q20,Q2,
     1              zfuncns,zfuncns3,zfuncns8,
     2              zfuncsg,zfuncnsv3,zfuncnsv8)
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
         if(nfi.lt.2)then
            xfevN(4) = zfuncnsv3 * xfev0N(3)
         else
            xfevN(4) = zfuncns(2) * xfev0N(4)
         endif
*     V8
         if(nfi.lt.3)then
            xfevN(5)  = zfuncnsv8 * xfev0N(3)
         else
            xfevN(5)  = zfuncns(2) * xfev0N(5)
         endif

*     T3
         if(nfi.lt.2)then
            xfevN(6) = (0d0,0d0)
            do j=1,2
               xfevN(6) = xfevN(6) + zfuncns3(j) * xfev0N(j)
            enddo
         else
            xfevN(6)  = zfuncns(1) * xfev0N(6)
         endif
*     T8
         if(nfi.lt.3)then
            xfevN(7) = (0d0,0d0)
            do j=1,2
               xfevN(7) = xfevN(7) + zfuncns8(j) * xfev0N(j)
            enddo
         else
            xfevN(7) = zfuncns(1) * xfev0N(7)
         endif
*
*     Copy evolved distributions into the initial ones before
*     the next evolution.
*
         if(iQ.lt.nQ-1)then
            do i=1,7
               xfev0N(i) = xfevN(i)
            enddo
         endif
      enddo
*
      return
      end
