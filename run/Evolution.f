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
      double precision aQED
      double precision xlha(12)
      double precision Q
      double precision xf(-9:9)
      character*100 card
      data xlha / 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0,
     1            0.8d0, 0.9d0, 0.95d0, 0.99d0, 0.999d0/
*
c      write(6,*)
c      write(6,*) "Type the name of the input card (e.g. Reference.ini)"
c      read(5,*) card
      card = "Reference.ini"
*
*     Read parameters of the evolution from the card
*
      call ReadParameters(card)
*
*     Initialization of evolution parameters
*
      call InitializeEvolution
*
*     Final scale (initial scale assumed to be the electron mass)
*
      !Q = 100d0
      Q = 0.000510998928d0
*
      write(6,*) "alpha(Q) = ", aQED(Q**2) * 12.566370614359173d0
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
         call xDistributions(xlha(ilha),Q,xf)
*
         write(6,'(es8.2,7(es12.4))')
     1         xlha(ilha),
     2         ( xf(1) + xf(-1) ) / xlha(ilha),
     3         xf(0) / xlha(ilha),
     4         ( xf(1) - xf(-1) ) / xlha(ilha),
     5         ( xf(2) + xf(-2) ) / xlha(ilha),
     6         ( xf(2) - xf(-2) ) / xlha(ilha),
     7         ( xf(3) + xf(-3) ) / xlha(ilha),
     8         ( xf(3) - xf(-3) ) / xlha(ilha)
      enddo
      write(*,*) "  "
*
      end
