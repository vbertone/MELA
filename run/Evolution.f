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
      double precision pi
      double complex Q(100)
      double complex xf(-6:6)
      double complex aQED
      character*100 card

      parameter(eps = 1d-10)
      parameter(pi = 3.1415926535897932385d0)
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
      Q(2) = 10d0
*
      write(6,*) "alpha(Q) = ", dble(aQED(Q(2)**2)) * 4d0 * pi
*
      write(6,*)
     1     "  x   ",
     2     "      e-    ",
     3     "      e+    ",
     4     "     mu-    ",
     5     "     mu+    ",
     6     "     tau-   ",
     7     "     tau+   ",
     8     "    photon  "
      do ilha=3,11
*
*     Perform the N-space evolution
*
         call xDistributions(xlha(ilha),nQ,Q,xf)
*
         write(6,'(es7.1,7(es12.4))')
     1         xlha(ilha),
     2         dble(xf(1)),
     2         dble(xf(-1)),
     2         dble(xf(2)),
     2         dble(xf(-2)),
     2         dble(xf(3)),
     2         dble(xf(-3)),
     6         dble(xf(0))
      enddo
      write(*,*) "  "
*
      end
