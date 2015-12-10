************************************************************************
*
*     Timing.f:
*
*     code to measure the time taken by the evolution and the computation
*     of the structure functions.
*
************************************************************************
      program Timing
*
      implicit none
*
      integer nQ,nx
      integer ix,iQ
      double precision x,xmin,xmax,xstep
      double precision eps
      double precision t1,t2
      double complex Q(100),Qmin,Qmax,Qstep
      double complex SFx(3,0:6)
      double complex xf(-6:6)
      character*100 card
      parameter(eps=1d-10)
*
      write(6,*)
      write(6,*) "Type the name of the input card (e.g. Reference.ini)"
      read(5,*) card
*
*     Read parameters of the evolution from the card
*
      call ReadParameters(card)
*
*     Initialization of evolution parameters
*
      call InitializeEvolution
*
      nx    = 100
      xmin  = 1d-5
      xmax  = 1d0
      xstep = dexp( dlog( xmax / xmin ) / ( nx - 1 ) )
*
      Q(1) = dsqrt(2d0) - eps
*
      nQ    = 10
      Qmin  = dsqrt(3d0)
      Qmax  = dsqrt(10000d0)
      Qstep = exp( zlog( Qmax / Qmin ) / ( nQ - 1 ) )
*
      call cpu_time(t1)
      Q(2) = Qmin
      do iQ=1,nQ
         x = xmin
         do ix=1,nx
            call xDistributions(x,2,Q,xf)
            x = x * xstep
         enddo
         Q(2) = Q(2) * Qstep
      enddo
      write(*,*) "  "
      call cpu_time(t2)
      write(6,*) "Time take for the computation of ",nx*nQ," evolution",
     1     " points"
      write(6,*) nx," points in x for each of the ",nQ," points in Q"
      write(6,*) "equal to: ",t2-t1," s"
*
      call cpu_time(t1)
      Q(2) = Qmin
      do iQ=1,nQ
         x = xmin
         do ix=1,nx
            call xStructureFunctions(x,2,Q,SFx)
            x = x * xstep
         enddo
         Q(2) = Q(2) * Qstep
      enddo
      write(*,*) "  "
      call cpu_time(t2)
      write(6,*) "Time take for the computation of ",nx*nQ," structure",
     1     " functions"
      write(6,*) nx," points in x for each of the ",nQ," points in Q"
      write(6,*) "equal to: ",t2-t1," s"
*
      end
