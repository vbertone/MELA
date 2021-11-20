************************************************************************
*
*     NDistributions.f:
*
*     Routine that returns the evolved N-space distributions.
*
************************************************************************
      subroutine NDistributions(N,Q,fn)
*
      implicit none
**
*     Input Variables
*
      double precision Q
      double complex N
**
*     Internal Variables
*
      double complex fn0(19)
      double complex efsg(4,4)
      double complex efV(3)
      double complex efT1(3,4),efT2(3,4)
**
*     Output Variables
*
      double complex fn(19)
*
*     Call N-space distributions
*
      call electronPDFsn(N, fn0) ! Electron PDFs
*
*     Call evolution kernels
*
      call zfunc(N,Q**2,efsg,efV,efT1,efT2)
*
*     Convolute PDFs with evolution kernels
*
*      1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
*     gm Sgl T1l T2l  Vl V1l V2l Sgu T1u T2u  Vu V1u V2u Sgd T1d T2d  Vd V1d V2d
*
*     Singlet
*
      fn(1)  = efsg(1,1) * fn0(1) + efsg(1,2) * fn0(2)
     1     + efsg(1,4) * fn0(8) + efsg(1,4) * fn0(14)
      fn(2)  = efsg(2,1) * fn0(1) + efsg(2,2) * fn0(2)
     1     + efsg(2,4) * fn0(8) + efsg(2,4) * fn0(14)
      fn(8)  = efsg(3,1) * fn0(1) + efsg(3,2) * fn0(2)
     1     + efsg(3,4) * fn0(8) + efsg(3,4) * fn0(14)
      fn(14) = efsg(4,1) * fn0(1) + efsg(4,2) * fn0(2)
     1     + efsg(4,4) * fn0(8) + efsg(4,4) * fn0(14)
*
*     Total valence
*
      fn(5)  = efV(1) * fn0(5)
      fn(11) = efV(2) * fn0(11)
      fn(17) = efV(3) * fn0(17)
*
*     T1
*
      fn(3)  = efT1(1,1) * fn0(1) + efT1(1,2) * fn0(2)
     1     + efT1(1,4) * fn0(8) + efT1(1,4) * fn0(14)
      fn(9)  = efT1(2,1) * fn0(1) + efT1(2,2) * fn0(2)
     1     + efT1(2,4) * fn0(8) + efT1(2,4) * fn0(14)
      fn(15) = efT1(3,1) * fn0(1) + efT1(3,2) * fn0(2)
     1     + efT1(3,4) * fn0(8) + efT1(3,4) * fn0(14)
*
*     T2
*
      fn(4)  = efT2(1,1) * fn0(1) + efT2(1,2) * fn0(2)
     1     + efT2(1,4) * fn0(8) + efT2(1,4) * fn0(14)
      fn(10) = efT2(2,1) * fn0(1) + efT2(2,2) * fn0(2)
     1     + efT2(2,4) * fn0(8) + efT2(2,4) * fn0(14)
      fn(16) = efT2(3,1) * fn0(1) + efT2(3,2) * fn0(2)
     1     + efT2(3,4) * fn0(8) + efT2(3,4) * fn0(14)
*
*     V1
*
      fn(6)  = efV(1) * fn0(5)
      fn(12) = efV(2) * fn0(11)
      fn(18) = efV(3) * fn0(17)
*
*     V2
*
      fn(7)  = efV(1) * fn0(5)
      fn(13) = efV(2) * fn0(11)
      fn(19) = efV(3) * fn0(17)
*
      return
      end
