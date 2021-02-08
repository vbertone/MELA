************************************************************************
*
*     SIDISStructureFunctions.f:
*
************************************************************************
      program SIDISStructureFunctions
*
      implicit none
*
      integer i, j, k
      double precision eps
      double precision xlha(7)
      double precision zlha(6)
      double complex Q0, Q, Qv(2)
      double complex SFx(2)
      double complex xf(-6:6), xd(-6:6), F2LO, f2lov
      double complex bq(6), dq(6), ch
      character*100 card

      parameter(eps=1d-10)
      data xlha / 1d-3, 1d-2, 1d-1, 3d-1, 5d-1, 7d-1, 9d-1 /
      data zlha / 1d-2, 1d-1, 3d-1, 5d-1, 7d-1, 9d-1 /
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
      Q0 = dsqrt(2d0) - eps
      Q  = 10
      Qv(1) = Q0
      Qv(2) = Q
*
      write(6,*) "MELA SIDIS structure functions:"
      write(6,*)
     1     "   x    ","   z   ","     Q    ",
     2     "   F2 (LO)   ","     F2     ","     FL     "
      do i = 1, 7
         do j = 1, 6
*     Compute LO F2 using PDFs and FFs
            f2lov = F2LO(xlha(i),zlha(j),Q0,Q)
*
            call SIDISxStructureFunctions(xlha(i), zlha(j), Q0, Q, SFx)
            write(6,'(2(es7.1,2x),es8.2,3(es12.4))')
     1           xlha(i),zlha(j),abs(Q),
     2           dreal(f2lov),
     3           dreal(SFx(1)),dreal(SFx(2))
         enddo
      enddo
      write(*,*) "  "
*
      end
