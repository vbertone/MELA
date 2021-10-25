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
      double precision xlha(12)
      double precision pi
      double complex Q(100)
      double complex xf(-3:3)
      double complex aQED
      character*100 card

      parameter(eps = 1d-10)
      parameter(pi = 3.1415926535897932385d0)
      data xlha / 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0,
     1            0.8d0, 0.9d0, 0.95d0, 0.99d0, 0.999d0/;

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
      Q(1) = 0.000510998928d0
      Q(2) = 100d0
*
      write(6,*) "alpha(Q) = ", dble(aQED(Q(2)**2)) * 4d0 * pi
*
      write(6,*)
     1     "  x   ",
     2     "   e- + e+  ",
     3     "    photon  ",
     4     "   e- - e+  ",
     5     "  mu- + mu+ ",
     6     "  mu- - mu+ ",
     7     " tau- + tau+",
     8     " tau- - tau+"
      do ilha=1,12
*
*     Perform the N-space evolution
*
         call xDistributions(xlha(ilha),nQ,Q,xf)
*
         write(6,'(es8.2,7(es12.4))')
     1         xlha(ilha),
     2         dble(xf(1) + xf(-1)) / xlha(ilha),
     3         dble(xf(0)) / xlha(ilha),
     4         dble(xf(1) - xf(-1)) / xlha(ilha),
     5         dble(xf(2) + xf(-2)) / xlha(ilha),
     6         dble(xf(2) - xf(-2)) / xlha(ilha),
     7         dble(xf(3) + xf(-3)) / xlha(ilha),
     8         dble(xf(3) - xf(-3)) / xlha(ilha)
      enddo
      write(*,*) "  "
*
      end
