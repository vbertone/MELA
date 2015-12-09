************************************************************************
*     Convert from the LHA convention ordering:
*
*      -6   -5   -4   -3   -2   -1   0   1   2   3   4   5   6 
*     tbar bbar cbar sbar ubar dbar  g   d   u   s   c   b   t
*
*     to the EVOLUTION convention:
*
*       1    2   3   4   5    6   7    8    9    10   11    12    13
*     Sigma  g   V  V_3 V_8 V_15 V_24 V_35  T_3  T_8  T_15  T_24  T_35
*
************************************************************************
      subroutine lha2evln(fin,fout)
*
      implicit none
**
*     Input Variables
*
      double precision fin(-6:6)
**
*     Internal Variables
*
      integer jf
**
*     Output Variables
*
      double precision fout(13)
*
      fout(1) = 0d0
      do jf=1,6
         fout(1) = fout(1) + ( fin(jf) + fin(-jf) )
      enddo
*
      fout(2) = fin(0)
*
      fout(3) = 0d0
      do jf=1,6
         fout(3) = fout(3) + ( fin(jf) - fin(-jf) )
      enddo
*
      fout(4)  = ( fin(2) - fin(-2) ) - ( fin(1) - fin(-1) )
      fout(5)  = ( fin(2) - fin(-2) ) + ( fin(1) - fin(-1) )
     1           - 2d0 * ( fin(3) - fin(-3) )
      fout(6)  = ( fin(2) - fin(-2) ) + ( fin(1) - fin(-1) )
     1           + ( fin(3) - fin(-3) ) 
     2           - 3d0 * ( fin(4) - fin(-4) )
      fout(7)  = ( fin(2) - fin(-2) ) + ( fin(1) - fin(-1) )
     1           + ( fin(3) - fin(-3) ) + ( fin(4) - fin(-4) )
     2           - 4d0 * ( fin(5) - fin(-5) )
      fout(8)  = ( fin(2) - fin(-2) ) + ( fin(1) - fin(-1) )
     1           + ( fin(3) - fin(-3) ) + ( fin(4) - fin(-4) )
     2           + ( fin(5) - fin(-5) ) 
     3           - 5d0 * ( fin(6) - fin(-6) )
      fout(9)  = ( fin(2) + fin(-2) ) - ( fin(1) + fin(-1) )
      fout(10) = ( fin(2) + fin(-2) ) + ( fin(1) + fin(-1) )
     1           - 2d0 * ( fin(3) + fin(-3) )
      fout(11) = ( fin(2) + fin(-2) ) + ( fin(1) + fin(-1) )
     1           + ( fin(3) + fin(-3) ) 
     2           - 3d0 * ( fin(4) + fin(-4) )
      fout(12) = ( fin(2) + fin(-2) ) + ( fin(1) + fin(-1) )
     1           + ( fin(3) + fin(-3) ) + ( fin(4) + fin(-4) )
     2           - 4d0 * ( fin(5) + fin(-5) )
      fout(13) = ( fin(2) + fin(-2) ) + ( fin(1) + fin(-1) )
     1           + ( fin(3) + fin(-3) ) + ( fin(4) + fin(-4) )
     2           + ( fin(5) + fin(-5) ) 
     3           - 5d0 * ( fin(6) + fin(-6) )
*
      return
      end
*
************************************************************************
*     Convert from the EVOLUTION convention ordering:
*
*       1    2   3   4   5    6   7    8    9    10   11    12    13
*     Sigma  g   V  V_3 V_8 V_15 V_24 V_35  T_3  T_8  T_15  T_24  T_35
*
*     to the LHA convention:
*
*      -6   -5   -4   -3   -2   -1   0   1   2   3   4   5   6 
*     tbar bbar cbar sbar ubar dbar  g   d   u   s   c   b   t
*
************************************************************************
      subroutine evln2lha(fin,fout)
*
      implicit none
**
*     Input Variables
*
      double precision fin(13)
**
*     Internal Variables
*
      integer i,j
      double precision rot(6,6)
      double precision ford(-6:6)
      double precision fp(6),fm(6)
**
*     Output Variables
*
      double precision fout(-6:6)
*
*     Gluone
*
      fout(0) = fin(2)
*
      rot(1,1) =    1d0 / 6d0
      rot(1,2) =    1d0 / 2d0
      rot(1,3) =    1d0 / 6d0
      rot(1,4) =    1d0 / 12d0
      rot(1,5) =    1d0 / 20d0
      rot(1,6) =    1d0 / 30d0
*
      rot(2,1) =    1d0 / 6d0
      rot(2,2) =  - 1d0 / 2d0
      rot(2,3) =    1d0 / 6d0
      rot(2,4) =    1d0 / 12d0
      rot(2,5) =    1d0 / 20d0
      rot(2,6) =    1d0 / 30d0
*
      rot(3,1) =    1d0 / 6d0
      rot(3,2) =    0d0
      rot(3,3) =  - 1d0 / 3d0
      rot(3,4) =    1d0 / 12d0
      rot(3,5) =    1d0 / 20d0
      rot(3,6) =    1d0 / 30d0
*
      rot(4,1) =    1d0 / 6d0
      rot(4,2) =    0d0
      rot(4,3) =    0d0
      rot(4,4) =  - 1d0 / 4d0
      rot(4,5) =    1d0 / 20d0
      rot(4,6) =    1d0 / 30d0
*
      rot(5,1) =    1d0 / 6d0
      rot(5,2) =    0d0
      rot(5,3) =    0d0
      rot(5,4) =    0d0
      rot(5,5) =  - 1d0 / 5d0
      rot(5,6) =    1d0 / 30d0
*
      rot(6,1) =    1d0 / 6d0
      rot(6,2) =    0d0
      rot(6,3) =    0d0
      rot(6,4) =    0d0
      rot(6,5) =    0d0
      rot(6,6) =  - 1d0 / 6d0
*
      ford(-1) = fin(3)
      ford(-2) = fin(4)
      ford(-3) = fin(5)
      ford(-4) = fin(6)
      ford(-5) = fin(7)
      ford(-6) = fin(8)
*
      ford(1)  = fin(1)
      ford(2)  = fin(9)
      ford(3)  = fin(10)
      ford(4)  = fin(11)
      ford(5)  = fin(12)
      ford(6)  = fin(13)
*
      do i=1,6
          fp(i) = 0d0
          fm(i) = 0d0
      enddo
*
      do i=1,6
         do j=1,6
            fp(i) = fp(i) + rot(i,j) * ford(j)
            fm(i) = fm(i) + rot(i,j) * ford(-j)
         enddo
      enddo
*
      fout(1)  = ( fp(2) + fm(2) ) / 2d0
      fout(-1) = ( fp(2) - fm(2) ) / 2d0
      fout(2)  = ( fp(1) + fm(1) ) / 2d0
      fout(-2) = ( fp(1) - fm(1) ) / 2d0
      fout(3)  = ( fp(3) + fm(3) ) / 2d0
      fout(-3) = ( fp(3) - fm(3) ) / 2d0
      fout(4)  = ( fp(4) + fm(4) ) / 2d0
      fout(-4) = ( fp(4) - fm(4) ) / 2d0
      fout(5)  = ( fp(5) + fm(5) ) / 2d0
      fout(-5) = ( fp(5) - fm(5) ) / 2d0
      fout(6)  = ( fp(6) + fm(6) ) / 2d0
      fout(-6) = ( fp(6) - fm(6) ) / 2d0
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
      double complex fin(-6:6)
**
*     Internal Variables
*
      integer jf
**
*     Output Variables
*
      double complex fout(13)
*
      fout(1) = (0d0,0d0)
      do jf=1,6
         fout(1) = fout(1) + ( fin(jf) + fin(-jf) )
      enddo
*
      fout(2) = fin(0)
*
      fout(3) = (0d0,0d0)
      do jf=1,6
         fout(3) = fout(3) + ( fin(jf) - fin(-jf) )
      enddo
*
      fout(4)  = ( fin(2) - fin(-2) ) - ( fin(1) - fin(-1) )
      fout(5)  = ( fin(2) - fin(-2) ) + ( fin(1) - fin(-1) )
     1           - 2d0 * ( fin(3) - fin(-3) )
      fout(6)  = ( fin(2) - fin(-2) ) + ( fin(1) - fin(-1) )
     1           + ( fin(3) - fin(-3) ) 
     2           - 3d0 * ( fin(4) - fin(-4) )
      fout(7)  = ( fin(2) - fin(-2) ) + ( fin(1) - fin(-1) )
     1           + ( fin(3) - fin(-3) ) + ( fin(4) - fin(-4) )
     2           - 4d0 * ( fin(5) - fin(-5) )
      fout(8)  = ( fin(2) - fin(-2) ) + ( fin(1) - fin(-1) )
     1           + ( fin(3) - fin(-3) ) + ( fin(4) - fin(-4) )
     2           + ( fin(5) - fin(-5) ) 
     3           - 5d0 * ( fin(6) - fin(-6) )
      fout(9)  = ( fin(2) + fin(-2) ) - ( fin(1) + fin(-1) )
      fout(10) = ( fin(2) + fin(-2) ) + ( fin(1) + fin(-1) )
     1           - 2d0 * ( fin(3) + fin(-3) )
      fout(11) = ( fin(2) + fin(-2) ) + ( fin(1) + fin(-1) )
     1           + ( fin(3) + fin(-3) ) 
     2           - 3d0 * ( fin(4) + fin(-4) )
      fout(12) = ( fin(2) + fin(-2) ) + ( fin(1) + fin(-1) )
     1           + ( fin(3) + fin(-3) ) + ( fin(4) + fin(-4) )
     2           - 4d0 * ( fin(5) + fin(-5) )
      fout(13) = ( fin(2) + fin(-2) ) + ( fin(1) + fin(-1) )
     1           + ( fin(3) + fin(-3) ) + ( fin(4) + fin(-4) )
     2           + ( fin(5) + fin(-5) ) 
     3           - 5d0 * ( fin(6) + fin(-6) )
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
      double complex fin(13)
**
*     Internal Variables
*
      integer i,j
      double precision rot(6,6)
      double complex ford(-6:6)
      double complex fp(6),fm(6)
**
*     Output Variables
*
      double complex fout(-6:6)
*
*     Gluone
*
      fout(0) = fin(2)
*
      rot(1,1) =    1d0 / 6d0
      rot(1,2) =    1d0 / 2d0
      rot(1,3) =    1d0 / 6d0
      rot(1,4) =    1d0 / 12d0
      rot(1,5) =    1d0 / 20d0
      rot(1,6) =    1d0 / 30d0
*
      rot(2,1) =    1d0 / 6d0
      rot(2,2) =  - 1d0 / 2d0
      rot(2,3) =    1d0 / 6d0
      rot(2,4) =    1d0 / 12d0
      rot(2,5) =    1d0 / 20d0
      rot(2,6) =    1d0 / 30d0
*
      rot(3,1) =    1d0 / 6d0
      rot(3,2) =    0d0
      rot(3,3) =  - 1d0 / 3d0
      rot(3,4) =    1d0 / 12d0
      rot(3,5) =    1d0 / 20d0
      rot(3,6) =    1d0 / 30d0
*
      rot(4,1) =    1d0 / 6d0
      rot(4,2) =    0d0
      rot(4,3) =    0d0
      rot(4,4) =  - 1d0 / 4d0
      rot(4,5) =    1d0 / 20d0
      rot(4,6) =    1d0 / 30d0
*
      rot(5,1) =    1d0 / 6d0
      rot(5,2) =    0d0
      rot(5,3) =    0d0
      rot(5,4) =    0d0
      rot(5,5) =  - 1d0 / 5d0
      rot(5,6) =    1d0 / 30d0
*
      rot(6,1) =    1d0 / 6d0
      rot(6,2) =    0d0
      rot(6,3) =    0d0
      rot(6,4) =    0d0
      rot(6,5) =    0d0
      rot(6,6) =  - 1d0 / 6d0
*
      ford(-1) = fin(3)
      ford(-2) = fin(4)
      ford(-3) = fin(5)
      ford(-4) = fin(6)
      ford(-5) = fin(7)
      ford(-6) = fin(8)
*
      ford(1)  = fin(1)
      ford(2)  = fin(9)
      ford(3)  = fin(10)
      ford(4)  = fin(11)
      ford(5)  = fin(12)
      ford(6)  = fin(13)
*
      do i=1,6
          fp(i) = (0d0,0d0)
          fm(i) = (0d0,0d0)
      enddo
*
      do i=1,6
         do j=1,6
            fp(i) = fp(i) + rot(i,j) * ford(j)
            fm(i) = fm(i) + rot(i,j) * ford(-j)
         enddo
      enddo
*
      fout(1)  = ( fp(2) + fm(2) ) / 2d0
      fout(-1) = ( fp(2) - fm(2) ) / 2d0
      fout(2)  = ( fp(1) + fm(1) ) / 2d0
      fout(-2) = ( fp(1) - fm(1) ) / 2d0
      fout(3)  = ( fp(3) + fm(3) ) / 2d0
      fout(-3) = ( fp(3) - fm(3) ) / 2d0
      fout(4)  = ( fp(4) + fm(4) ) / 2d0
      fout(-4) = ( fp(4) - fm(4) ) / 2d0
      fout(5)  = ( fp(5) + fm(5) ) / 2d0
      fout(-5) = ( fp(5) - fm(5) ) / 2d0
      fout(6)  = ( fp(6) + fm(6) ) / 2d0
      fout(-6) = ( fp(6) - fm(6) ) / 2d0
*
      return
      end

