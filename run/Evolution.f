************************************************************************
*
*     Evolution.f:
*
*     Driver for the generation of the evolution tables.
*
************************************************************************
      program Evolution
*
      implicit none
*
      integer ilha
      integer nQ
      double precision eps
      double precision xlha(11)
      double complex Q(100)
      double complex xf(-6:6)
      character*100 card

      parameter(eps=1d-10)
      data xlha / 1d-7, 1d-6, 1d-5, 1d-4, 1d-3, 1d-2,
     1            1d-1, 3d-1, 5d-1, 7d-1, 9d-1 /
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
      nQ = 2
      Q(1) = dsqrt(2d0) - eps
c      Q(2) = (10d0,10d0)
c      Q(3) = (20d0,-20d0)
c      Q(4) = (50d0,50d0)
c      Q(2) = Q(1)
      Q(2) = dsqrt(100d0)
*
      write(6,*) "Standard evolution:"
      write(6,*)
     1     "  x    ",
     2     " u-ubar     ",
     3     " d-dbar     ",
     4     " 2(ubr+dbr) ",
     5     " c+cbar     ",
     6     " gluon      "
      do ilha=3,11
*
*     Perform the N-space evolution
*
         call xDistributions(xlha(ilha),nQ,Q,xf)
*
         write(6,'(es7.1,5(es12.4))')
     1         xlha(ilha),
     2         dble(xf(2) - xf(-2)),
     3         dble(xf(1) - xf(-1)),
     4         dble(2d0 * ( xf(-1) + xf(-2) )),
     5         dble(xf(4) + xf(-4)),
     6         dble(xf(0))
      enddo
      write(*,*) "  "
*
      end
