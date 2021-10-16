************************************************************************
*     Convert from the LHA convention ordering:
*
*      -3   -2   -1   0   1   2   3
*      tau+ mu+  e+   gm  e-  mu- tau-
*
*     to the EVOLUTION convention:
*
*       1    2    3    4    5    6    7
*     Sigma  gm   V    V_3  V_8  T_3  T_8
*
************************************************************************
      subroutine lha2evln(fin,fout)
*
      implicit none
**
*     Input Variables
*
      double precision fin(-3:3)
**
*     Internal Variables
*
      integer jf
**
*     Output Variables
*
      double precision fout(7)
*
*     Singlet
      fout(1) = 0d0
      do jf=1,3
         fout(1) = fout(1) + ( fin(jf) + fin(-jf) )
      enddo
*
*     Gamma
      fout(2) = fin(0)
*
*     Total valence
      fout(3) = 0d0
      do jf=1,3
         fout(3) = fout(3) + ( fin(jf) - fin(-jf) )
      enddo
*
*     V3
      fout(4) = ( fin(1) - fin(-1) ) - ( fin(2) - fin(-2) )
*     V8
      fout(5) = ( fin(1) - fin(-1) ) + ( fin(2) - fin(-2) )
     1          - 2d0 * ( fin(3) - fin(-3) )
*     T3
      fout(6) = ( fin(1) + fin(-1) ) - ( fin(2) + fin(-2) )
*     T8
      fout(7) = ( fin(1) + fin(-1) ) + ( fin(2) + fin(-2) )
     1           - 2d0 * ( fin(3) + fin(-3) )
*
      return
      end
*
************************************************************************
*     Convert from the EVOLUTION convention ordering:
*
*       1    2    3    4    5    6    7
*     Sigma  gm   V    V_3  V_8  T_3  T_8
*
*     to the LHA convention:
*
*      -3   -2   -1   0   1   2   3
*      tau+ mu+  e+   gm  e-  mu- tau-
*
************************************************************************
      subroutine evln2lha(fin,fout)
*
      implicit none
**
*     Input Variables
*
      double precision fin(7)
**
*     Internal Variables
*
      integer i,j
      double precision rot(3,3)
      double precision ford(-3:3)
      double precision fp(3),fm(3)
**
*     Output Variables
*
      double precision fout(-3:3)
*
*     Gamma
*
      fout(0) = fin(2)
*
      rot(1,1) =    1d0 / 3d0
      rot(1,2) =    1d0 / 2d0
      rot(1,3) =    1d0 / 6d0
*
      rot(2,1) =    1d0 / 3d0
      rot(2,2) =  - 1d0 / 2d0
      rot(2,3) =    1d0 / 6d0
*
      rot(3,1) =    1d0 / 3d0
      rot(3,2) =    0d0
      rot(3,3) =  - 1d0 / 3d0
*
      ford(-1) = fin(3)
      ford(-2) = fin(4)
      ford(-3) = fin(5)
*
      ford(1)  = fin(1)
      ford(2)  = fin(6)
      ford(3)  = fin(7)
*
      do i=1,3
          fp(i) = 0d0
          fm(i) = 0d0
      enddo
*
      do i=1,3
         do j=1,3
            fp(i) = fp(i) + rot(i,j) * ford(j)
            fm(i) = fm(i) + rot(i,j) * ford(-j)
         enddo
         fout(i)  = ( fp(i) + fm(i) ) / 2d0
         fout(-i) = ( fp(i) - fm(i) ) / 2d0
      enddo
*
      return
      end
*
************************************************************************
*
*     Same rotation routines as above but for complex arguments
*
************************************************************************
      subroutine lha2evlnc(fin,fout)
*
      implicit none
**
*     Input Variables
*
      double complex fin(-3:3)
**
*     Internal Variables
*
      integer jf
**
*     Output Variables
*
      double complex fout(7)
*
*     Singlet
      fout(1) = (0d0, 0d0)
      do jf=1,3
         fout(1) = fout(1) + ( fin(jf) + fin(-jf) )
      enddo
*
*     Gamma
      fout(2) = fin(0)
*
*     Total valence
      fout(3) = (0d0, 0d0)
      do jf=1,3
         fout(3) = fout(3) + ( fin(jf) - fin(-jf) )
      enddo
*
*     V3
      fout(4) = ( fin(1) - fin(-1) ) - ( fin(2) - fin(-2) )
*     V8
      fout(5) = ( fin(1) - fin(-1) ) + ( fin(2) - fin(-2) )
     1          - 2d0 * ( fin(3) - fin(-3) )
*     T3
      fout(6) = ( fin(1) + fin(-1) ) - ( fin(2) + fin(-2) )
*     T8
      fout(7) = ( fin(1) + fin(-1) ) + ( fin(2) + fin(-2) )
     1           - 2d0 * ( fin(3) + fin(-3) )
*
      return
      end
*
************************************************************************
      subroutine evln2lhac(fin,fout)
*
      implicit none
**
*     Input Variables
*
      double complex fin(7)
**
*     Internal Variables
*
      integer i,j
      double precision rot(3,3)
      double complex ford(-3:3)
      double complex fp(3),fm(3)
**
*     Output Variables
*
      double complex fout(-3:3)
*
*     Gamma
*
      fout(0) = fin(2)
*
      rot(1,1) =    1d0 / 3d0
      rot(1,2) =    1d0 / 2d0
      rot(1,3) =    1d0 / 6d0
*
      rot(2,1) =    1d0 / 3d0
      rot(2,2) =  - 1d0 / 2d0
      rot(2,3) =    1d0 / 6d0
*
      rot(3,1) =    1d0 / 3d0
      rot(3,2) =    0d0
      rot(3,3) =  - 1d0 / 3d0
*
      ford(-1) = fin(3)
      ford(-2) = fin(4)
      ford(-3) = fin(5)
*
      ford(1)  = fin(1)
      ford(2)  = fin(6)
      ford(3)  = fin(7)
*
      do i=1,3
          fp(i) = (0d0, 0d0)
          fm(i) = (0d0, 0d0)
      enddo
*
      do i=1,3
         do j=1,3
            fp(i) = fp(i) + rot(i,j) * ford(j)
            fm(i) = fm(i) + rot(i,j) * ford(-j)
         enddo
         fout(i)  = ( fp(i) + fm(i) ) / 2d0
         fout(-i) = ( fp(i) - fm(i) ) / 2d0
      enddo
*
      return
      end
