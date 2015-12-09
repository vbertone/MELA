************************************************************************
*
*     StructureFunctions.f:
*
*     Driver for the generation of the evolution tables.
*
************************************************************************
      program StructureFunctions
*
      implicit none
*
      integer ilha
      integer nQ
      double precision eps
      double precision xlha(11)
      double complex Q(100)
      double complex SFx(3,0:6)
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
      Q(2) = dsqrt(10000d0)
*
      write(6,*) "MELA evolution:"
      write(6,*)
     1     "   x   ","    Q2    ",
     2     "   u-ubar   ",
     3     "   d-dbar   ",
     4     " 2(ubr+dbr) ",
     5     "   c+cbar   ",
     6     "   gluon    "
      do ilha=3,11
         call xDistibutions(xlha(ilha),nQ,Q,xf)
*
         write(6,'(es7.1,2x,es8.2,5(es12.4))')
     1         xlha(ilha),abs(Q(nQ)**2d0),
     2         dreal(xf(2) - xf(-2)),
     3         dreal(xf(1) - xf(-1)),
     4         dreal(2d0 * ( xf(-1) + xf(-2) )),
     5         dreal(xf(4) + xf(-4)),
     6         dreal(xf(0))
      enddo
      write(*,*) "  "
*
      write(6,*) "MELA structure functions:"
      write(6,*)
     1     "   x   ","    Q2    ",
     2     "  F2light   ",
     3     "  F2charm   ",
     4     "  FLlight   ",
     5     "  FLcharm   ",
     6     "  F3light   ",
     7     "  F3charm   "
      do ilha=3,11
         call xStructureFunctions(xlha(ilha),nQ,Q,SFx)
*
         write(6,'(es7.1,2x,es8.2,6(es12.4))')
     1        xlha(ilha),abs(Q(nQ)**2d0),
     2        dreal(SFx(1,1)+SFx(1,2)+SFx(1,3)),
     3        dreal(SFx(1,4)),
     4        dreal(SFx(2,1)+SFx(2,2)+SFx(2,3)),
     5        dreal(SFx(2,4)),
     6        dreal(SFx(3,1)+SFx(3,2)+SFx(3,3)),
     7        dreal(SFx(3,4))
      enddo
      write(*,*) "  "
*
      end
