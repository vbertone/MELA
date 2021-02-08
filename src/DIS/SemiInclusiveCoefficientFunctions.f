************************************************************************
*
*     The following routines return the Zero-Mass coefficient functions
*     for semi-nclusive DIS structure functions in Mellin
*     space. Expressions are taken from
*     https://arxiv.org/pdf/hep-ph/0107064.pdf.
*
************************************************************************
      DOUBLE COMPLEX FUNCTION CL1QQ_SIDIS(N, M)
*
      IMPLICIT NONE
*
      include "../commons/colfact.h"
*
*     Input variable
*
      DOUBLE COMPLEX N, M
*
      CL1QQ_SIDIS = 8D0 * CF / ( M + 1D0 ) / ( N + 1D0 ) 
*
      RETURN
      END
*
************************************************************************
      DOUBLE COMPLEX FUNCTION CL1GQ_SIDIS(N, M)
*
      IMPLICIT NONE
*
      include "../commons/colfact.h"
*
*     Input variable
*
      DOUBLE COMPLEX N, M
*
      CL1GQ_SIDIS = 8D0 * CF / M / ( M + 1D0 ) / ( N + 1D0 ) 
*
      RETURN
      END
*
************************************************************************
      DOUBLE COMPLEX FUNCTION CL1QG_SIDIS(N)
*
      IMPLICIT NONE
*
      include "../commons/colfact.h"
*
*     Input variable
*
      DOUBLE COMPLEX N
*
      CL1QG_SIDIS = 2D0 * TR / ( N + 1D0 ) / ( N + 2D0 ) 
*
      RETURN
      END
*
************************************************************************
c$$$
c$$$      subroutine c21p(cn,c2pls1n,c2g1n)
c$$$      IMPLICIT complex*16 (A-H,O-Z)
c$$$*
c$$$      real*8 dl,pi,ge,z2,z3,z4,z5,z6
c$$$
c$$$      pi = 3.1415926535897932385d0
c$$$      dl = dlog(2.0d0)
c$$$      ge = 5.7721566490153286061d-1
c$$$      z2 = 1.6449340668482264365d0
c$$$      z3 = 1.2020569031595942854d0
c$$$      z4 = 1.0823232337111381916d0
c$$$      z5 = 1.0369277551433699263d0
c$$$      z6 = 1.0173430619844491398d0
c$$$
c$$$      c2pls1n=+4.d0/3.d0
c$$$     &  *(-9.d0-2.d0*(cn+1)**(-2)+4.d0*(cn+1)**(-1)+2.d0*(cn)**(-2)+3.
c$$$     &  d0*(cn)**(-1)+2.d0*s1(1,1.d0+cn)*(cn+1)**(-1)+3.d0*s1(1,cn)-2.d0
c$$$     &    *s1(1,cn)*(cn)**(-1)+2.d0*s1(1,cn)**2-2.d0*s1(2,cn))
c$$$
c$$$      c2g1n=
c$$$     &    4.d0*(cn+2)**(-2)-16.d0*(cn+2)**(-1)-4.d0*(cn+1)**(-2)+16.d0
c$$$     &    *(cn+1)**(-1)+2.d0*(cn)**(-2)-2.d0*(cn)**(-1)+4.d0*s1(1,1.d0+
c$$$     &    cn)*(cn+1)**(-1)-4.d0*s1(1,2.d0+cn)*(cn+2)**(-1)-2.d0*s1(1,cn
c$$$     &    )*(cn)**(-1)
c$$$*
c$$$      return
c$$$      end
c$$$*
c$$$************************************************************************
c$$$      subroutine c22p(cn,nf,c2pls2n,c2min2n,c2ps2n,c2g2n)
c$$$      IMPLICIT complex*16 (A-H,O-Z)
c$$$*
c$$$      integer nf
c$$$      real*8 dl,pi,ge,z2,z3,z4,z5,z6
c$$$
c$$$      pi = 3.1415926535897932385d0
c$$$      dl = dlog(2.0d0)
c$$$      ge = 5.7721566490153286061d-1
c$$$      z2 = 1.6449340668482264365d0
c$$$      z3 = 1.2020569031595942854d0
c$$$      z4 = 1.0823232337111381916d0
c$$$      z5 = 1.0369277551433699263d0
c$$$      z6 = 1.0173430619844491398d0
c$$$
c$$$      c2pls2n = 
c$$$     &     + nf * ( + 216883.d0/4629.d0 - 8113.d0/1000.d0*z2*cn**(-1)
c$$$     &     + 40.d0/9.d0*cn**(-3) - 16.d0/3.d0*cn**(-2)
c$$$     &     - 73101.d0/50000.d0*cn**(-1) + 111.d0/100.d0*(cn+1.d0)**(-4)
c$$$     &     - 891.d0/50.d0*(cn+1.d0)**(-1)
c$$$     &     - 1297.d0/100.d0*(cn+2.d0)**(-1) - 79361.d0/12500.d0*s1(1,cn)
c$$$     &     + 8113.d0/1000.d0*s1(1,cn)*cn**(-2)
c$$$     &     - 43949.d0/2700.d0*s1(1,cn)*cn**(-1)
c$$$     &     - 116.d0/27.d0*s1(1,cn)**2 + 5.d0/18.d0*s1(1,cn)**2*cn**(-1)
c$$$     &     - 16.d0/27.d0*s1(1,cn)**3 - 16.d0/9.d0*s1(1,cn)*s1(2,cn)
c$$$     &     - 116.d0/27.d0*s1(2,cn) + 75517.d0/9000.d0*s1(2,cn)*cn**(-1)
c$$$     &     - 32.d0/27.d0*s1(3,cn) )
c$$$      c2pls2n = c2pls2n
c$$$     &     - 338513.d0/1000.d0 + 151.d0/2.d0*z2*cn**(-2)
c$$$     &     - 2187.d0/10.d0*z3*cn**(-1) + 160.d0/9.d0*cn**(-4)
c$$$     &     - 207.d0/5.d0*cn**(-3) + 3548.d0/125.d0*cn**(-2)
c$$$     &     + 7641.d0/1000.d0*cn**(-1) + 2157.d0/125.d0*(cn+1.d0)**(-5)
c$$$     &     - 8067.d0/10.d0*(cn+1.d0)**(-1) - 188641.d0/1000.d0*s1(1,cn)
c$$$     &     - 1471.d0/5.d0*s1(1,cn)*z2*cn**(-1)
c$$$     &     - 151.d0/2.d0*s1(1,cn)*cn**(-3)
c$$$     &     + 1649763.d0/2500.d0*s1(1,cn)*cn**(-1)
c$$$     &     - 77763.d0/5000.d0*s1(1,cn)**2
c$$$     &     + 1471.d0/10.d0*s1(1,cn)**2*cn**(-2)
c$$$      c2pls2n = c2pls2n
c$$$     &     + 818.d0/75.d0*s1(1,cn)**2*cn**(-1) + 184.d0/9.d0*s1(1,cn)**3
c$$$     &     + 1583.d0/450.d0*s1(1,cn)**3*cn**(-1)
c$$$     &     + 32.d0/9.d0*s1(1,cn)**4 + 64.d0/3.d0*s1(1,cn)**2*s1(2,cn)
c$$$     &     + 184.d0/3.d0*s1(1,cn)*s1(2,cn)
c$$$     &     + 45713.d0/150.d0*s1(1,cn)*s1(2,cn)*cn**(-1)
c$$$     &     + 256.d0/9.d0*s1(1,cn)*s1(3,cn) - 77763.d0/5000.d0*s1(2,cn)
c$$$     &     + 358.d0/5.d0*s1(2,cn)*cn**(-2)
c$$$     &     + 818.d0/75.d0*s1(2,cn)*cn**(-1) + 32.d0/3.d0*s1(2,cn)**2
c$$$     &     + 368.d0/9.d0*s1(3,cn) + 101581.d0/450.d0*s1(3,cn)*cn**(-1)
c$$$     &     + 64.d0/3.d0*s1(4,cn)
c$$$
c$$$      c2ps2n =
c$$$     &     + ( - 1013.d0/25.d0*z2*cn**(-2)
c$$$     &     + 3393.d0/50.d0*z2*(cn+1.d0)**(-2)
c$$$     &     - 1013.d0/25.d0*z3*cn**(-1)
c$$$     &     + 3393.d0/50.d0*z3*(cn+1.d0)**(-1) - 80.d0/3.d0*cn**(-4)
c$$$     &     - 1053.d0/500.d0*cn**(-3) - 1084.d0/25.d0*cn**(-2)
c$$$     &     + 88831.d0/10000.d0*cn**(-1)
c$$$     &     + 52903.d0/10000.d0*(cn-1.d0)**(-1)
c$$$     &     - 17001.d0/500.d0*(cn+1.d0)**(-4)
c$$$     &     + 1084.d0/25.d0*(cn+1.d0)**(-2)
c$$$     &     - 185767.d0/5000.d0*(cn+1.d0)**(-1)
c$$$     &     + 3639.d0/100.d0*(cn+2.d0)**(-1)
c$$$     &     - 1341.d0/100.d0*(cn+3.d0)**(-1) )
c$$$      c2ps2n = c2ps2n
c$$$     &     + ( - 3393.d0/50.d0*s1(1,cn+1.d0)*(cn+1.d0)**(-3)
c$$$     &     - 32.d0/3.d0*s1(1,cn+1.d0)*(cn+1.d0)**(-1)
c$$$     &     - 8.d0/3.d0*s1(1,cn+1.d0)**2*(cn+1.d0)**(-1)
c$$$     &     + 4.d0/25.d0*s1(1,cn+1.d0)**3*(cn+1.d0)**(-1)
c$$$     &     + 12.d0/25.d0*s1(1,cn+1.d0)*s1(2,cn+1.d0)*(cn+1.d0)**(-1)
c$$$     &     - 2.d0/25.d0*s1(1,cn+2.d0)**3*(cn+2.d0)**(-1)
c$$$     &     - 6.d0/25.d0*s1(1,cn+2.d0)*s1(2,cn+2.d0)*(cn+2.d0)**(-1)
c$$$     &     + 1013.d0/25.d0*s1(1,cn)*cn**(-3)
c$$$     &     + 32.d0/3.d0*s1(1,cn)*cn**(-1)
c$$$     &     + 8.d0/3.d0*s1(1,cn)**2*cn**(-1)
c$$$     &     - 2.d0/25.d0*s1(1,cn)**3*cn**(-1)
c$$$     &     - 6.d0/25.d0*s1(1,cn)*s1(2,cn)*cn**(-1)
c$$$     &     - 3393.d0/50.d0*s1(2,cn+1.d0)*(cn+1.d0)**(-2)
c$$$     &     - 8.d0/3.d0*s1(2,cn+1.d0)*(cn+1.d0)**(-1) )
c$$$      c2ps2n = c2ps2n
c$$$     &     + ( + 1013.d0/25.d0*s1(2,cn)*cn**(-2)
c$$$     &     + 8.d0/3.d0*s1(2,cn)*cn**(-1)
c$$$     &     - 3377.d0/50.d0*s1(3,cn+1.d0)*(cn+1.d0)**(-1)
c$$$     &     - 4.d0/25.d0*s1(3,cn+2.d0)*(cn+2.d0)**(-1)
c$$$     &     + 1009.d0/25.d0*s1(3,cn)*cn**(-1) )
c$$$
c$$$      c2g2n = 
c$$$     &     + ( + 2086.d0*z2*cn**(-2) - 1469.d0/5.d0*z2*cn**(-1)
c$$$     &     - 3556.d0/5.d0*z2*(cn+1.d0)**(-1) + 2086.d0*z3*cn**(-1)
c$$$     &     - 140.d0/3.d0*cn**(-4) - 7109.d0/500.d0*cn**(-3)
c$$$     &     - 578.d0/5.d0*cn**(-2) - 63973.d0/10000.d0*cn**(-1)
c$$$     &     + 119033.d0/10000.d0*(cn-1.d0)**(-1)
c$$$     &     + 2408.d0*(cn+1.d0)**(-3) - 36761.d0/50.d0*(cn+1.d0)**(-1)
c$$$     &     + 7603.d0/10.d0*(cn+2.d0)**(-1)
c$$$     &     + 3556.d0/5.d0*s1(1,cn+1.d0)*(cn+1.d0)**(-2)
c$$$     &     - 593.d0/20.d0*s1(1,cn+1.d0)**3*(cn+1.d0)**(-1) )
c$$$      c2g2n = c2g2n
c$$$     &     + ( - 1779.d0/20.d0*s1(1,cn+1.d0)*s1(2,cn+1.d0)*
c$$$     &    (cn+1.d0)**(-1) - 2086.d0*s1(1,cn)*cn**(-3)
c$$$     &     + 1469.d0/5.d0*s1(1,cn)*cn**(-2)
c$$$     &     + 872.d0/25.d0*s1(1,cn)*cn**(-1) - 24.d0*s1(1,cn)**2*cn**(-1)
c$$$     &     + 4177.d0/180.d0*s1(1,cn)**3*cn**(-1)
c$$$     &     + 4177.d0/60.d0*s1(1,cn)*s1(2,cn)*cn**(-1)
c$$$     &     + 3556.d0/5.d0*s1(2,cn+1.d0)*(cn+1.d0)**(-1)
c$$$     &     - 2086.d0*s1(2,cn)*cn**(-2) + 1349.d0/5.d0*s1(2,cn)*cn**(-1)
c$$$     &     - 593.d0/10.d0*s1(3,cn+1.d0)*(cn+1.d0)**(-1)
c$$$     &     - 183563.d0/90.d0*s1(3,cn)*cn**(-1) )
c$$$
c$$$      call C22PSGNS(cn,nf,c2min2n)
c$$$*
c$$$      return
c$$$      end
c$$$*
c$$$************************************************************************
c$$$      subroutine cl1p(cn,clpls1n,clg1n)
c$$$      IMPLICIT complex*16 (A-H,O-Z)
c$$$*
c$$$      clpls1n=
c$$$     &  + 16.d0/3.d0*(cn+1)**(-1)
c$$$
c$$$      clg1n=( 8.d0*(cn+1)**(-1) - 8.d0*(cn+2)**(-1) )
c$$$*
c$$$      return
c$$$      end
c$$$*
c$$$************************************************************************
c$$$      subroutine cl2p(cn,nf,clpls2n,clmin2n,clps2n,clg2n)
c$$$      IMPLICIT complex*16 (A-H,O-Z)
c$$$*
c$$$      integer nf
c$$$      real*8 dl,pi,ge,z2,z3,z4,z5,z6
c$$$
c$$$      pi = 3.1415926535897932385d0
c$$$      dl = dlog(2.0d0)
c$$$      ge = 5.7721566490153286061d-1
c$$$      z2 = 1.6449340668482264365d0
c$$$      z3 = 1.2020569031595942854d0
c$$$      z4 = 1.0823232337111381916d0
c$$$      z5 = 1.0369277551433699263d0
c$$$      z6 = 1.0173430619844491398d0
c$$$
c$$$      clpls2n=
c$$$     &  + nf * ( 32.d0/9.d0*cn**(-1) + 64.d0/9.d0*(cn+1)**(-2) - 400.d0/
c$$$     &    27.d0*(cn+1)**(-1) - 32.d0/9.d0*s1(1,1.d0 + cn)*(cn+1)**(-1)
c$$$     &     )
c$$$      clpls2n=clpls2n
c$$$     &     - 3.d0/250.d0 + 42047.d0/500.d0*z2*cn**(-1) + 128.d0
c$$$     &    /9.d0*cn**(-2) - 18669.d0/500.d0*cn**(-1) + 1841.d0/50.d0*
c$$$     &    (cn+1)**(-3) - 329.d0/10.d0*(cn+1)**(-2) + 8953.d0/100.d0*
c$$$     &    (cn+1)**(-1) + 1691.d0/50.d0*(cn+2)**(-1) + 93.d0/2.d0*s1(1,1.
c$$$     &    d0 + cn)*(cn+1)**(-1) + 128.d0/9.d0*s1(1,1.d0 + cn)**2*
c$$$     &    (cn+1)**(-1) - 42047.d0/500.d0*s1(1,cn)*cn**(-2) + 128.d0/9.d0
c$$$     &    *s1(2,1.d0 + cn)*(cn+1)**(-1) - 42047.d0/500.d0*s1(2,cn)*
c$$$     &    cn**(-1)
c$$$
c$$$      clps2n=
c$$$     &  + ( 421.d0/500.d0*cn**(-3) - 2809.d0/100.d0*cn**(-2) + 
c$$$     &    1319.d0/50.d0*cn**(-1) - 237.d0/100.d0*(cn-1)**(-1) + 76.d0/
c$$$     &    25.d0*(cn+1)**(-3) + 2809.d0/100.d0*(cn+1)**(-2) - 1623.d0/25.
c$$$     &    d0*(cn+1)**(-1) + 3009.d0/50.d0*(cn+2)**(-1) - 1927.d0/100.d0
c$$$     &    *(cn+3)**(-1) + 9273.d0/250.d0*s1(1,1.d0 + cn)*(cn+1)**(-1)
c$$$     &     - 6591.d0/250.d0*s1(1,2.d0 + cn)*(cn+2)**(-1) + 1303.d0/250.d
c$$$     &    0*s1(1,3.d0 + cn)*(cn+3)**(-1) - 797.d0/50.d0*s1(1,cn)*
c$$$     &    cn**(-1) )
c$$$
c$$$      clg2n=
c$$$     &  + (  - 1161.d0*z2*(cn+1)**(-1) - 1983.d0/50.d0*cn**(-2) + 
c$$$     &    5333.d0/1000.d0*cn**(-1) - 5333.d0/1000.d0*(cn-1)**(-1) + 
c$$$     &    3003.d0/25.d0*(cn+1)**(-3) + 1983.d0/50.d0*(cn+1)**(-2) + 
c$$$     &    1161.d0*s1(1,1.d0 + cn)*(cn+1)**(-2) + 4324.d0/5.d0*s1(1,1.d0
c$$$     &     + cn)*(cn+1)**(-1) - 7197.d0/50.d0*s1(1,1.d0 + cn)**2*
c$$$     &    (cn+1)**(-1) + 246.d0/5.d0*s1(1,2.d0 + cn)**2*(cn+2)**(-1) - 
c$$$     &    4324.d0/5.d0*s1(1,cn)*cn**(-1) + 4737.d0/50.d0*s1(1,cn)**2*
c$$$     &    cn**(-1) + 50853.d0/50.d0*s1(2,1.d0 + cn)*(cn+1)**(-1) + 246.d
c$$$     &    0/5.d0*s1(2,2.d0 + cn)*(cn+2)**(-1) + 4737.d0/50.d0*s1(2,cn)*
c$$$     &    cn**(-1) )
c$$$
c$$$      call CL2PSGNS(cn,nf,clmin2n)
c$$$*
c$$$      return
c$$$      end
c$$$*
c$$$************************************************************************
c$$$*
c$$$* ..harmonic sums and the like
c$$$*
c$$$************************************************************************
c$$$      complex*16 function s1(k,cn)
c$$$      IMPLICIT complex*16 (A-H,O-Z)
c$$$      integer k
c$$$
c$$$      real*8 dl,pi,ge,z2,z3,z4,z5,z6
c$$$
c$$$      pi = 3.1415926535897932385d0
c$$$      dl = dlog(2.0d0)
c$$$      ge = 5.7721566490153286061d-1
c$$$      z2 = 1.6449340668482264365d0
c$$$      z3 = 1.2020569031595942854d0
c$$$      z4 = 1.0823232337111381916d0
c$$$      z5 = 1.0369277551433699263d0
c$$$      z6 = 1.0173430619844491398d0
c$$$
c$$$      if (k.eq.1) then
c$$$        cn1=cn+1
c$$$        call psi0(cn1,res)
c$$$        s1=ge + res
c$$$      end if
c$$$      if (k.eq.2) then
c$$$        cn1=cn+1
c$$$        call psi1(cn1,res)
c$$$        s1=z2 - res
c$$$      end if
c$$$      if (k.eq.3) then
c$$$        cn1=cn+1
c$$$        call psi2(cn1,res)
c$$$        s1=z3 + res/2.d0
c$$$      end if
c$$$      if (k.eq.4) then
c$$$        cn1=cn+1
c$$$        call psi3(cn1,res)
c$$$        s1=z4 - res/6.d0
c$$$      end if
c$$$      if (k.eq.5) then
c$$$        cn1=cn+1
c$$$        call psi4(cn1,res)
c$$$        s1=z5 + res/24.d0
c$$$      end if
c$$$      if (k.eq.6) then
c$$$        cn1=cn+1
c$$$        call psi5(cn1,res)
c$$$        s1=z6 - res/120.d0
c$$$      end if
c$$$      return
c$$$      end
c$$$*
c$$$************************************************************************
c$$$*
c$$$* ..psi(cn) for complex argument
c$$$*
c$$$************************************************************************
c$$$      subroutine psi0(cn,res)
c$$$      real*8 r
c$$$      complex*16 cn,t,t0,y,y2,res
c$$$
c$$$      t=dcmplx(0.d0,0.d0)
c$$$2     r=dsqrt(dreal(cn)**2+dimag(cn)**2)
c$$$      if(r.gt.16.d0) goto 1
c$$$      t=t-1.d0/cn
c$$$      cn=cn+1.d0
c$$$      goto 2
c$$$
c$$$1     y=1.d0/cn
c$$$      y2=y*y
c$$$
c$$$      t0 = 
c$$$     & (-1.d0/2.d0+(-1.d0/12.d0+(1.d0/120.d0+(-1.d0/252.d0+(1.d0/240.d0
c$$$     & +(-1.d0/132.d0+(691.d0/32760.d0+(-1.d0/12.d0+(3617.d0/8160.d0
c$$$     & +(-43867.d0/14364.d0+(174611.d0/6600.d0+(-77683.d0/276.d0
c$$$     & +(236364091.d0/65520.d0+(-657931.d0/12.d0+3392780147.d0/3480.d0
c$$$     & *y2)*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y)*y-log(y)
c$$$
c$$$      res=t+t0
c$$$
c$$$      return
c$$$      end
c$$$*
c$$$************************************************************************
c$$$*
c$$$* ..psi'(cn) (1st derivative) for complex argument
c$$$*
c$$$************************************************************************
c$$$      subroutine psi1(cn,res)
c$$$      real*8 r
c$$$      complex*16 cn,t,t0,y,y2,res
c$$$
c$$$      t=dcmplx(0.d0,0.d0)
c$$$2     r=dsqrt(dreal(cn)**2+dimag(cn)**2)
c$$$      if(r.gt.16.d0) goto 1
c$$$      t=t+1.d0/cn**2
c$$$      cn=cn+1.d0
c$$$      goto 2
c$$$1     y=1.d0/cn
c$$$      y2=y*y
c$$$
c$$$      t0 = 
c$$$     & (1.d0+(1.d0/2.d0+(1.d0/6.d0+(-1.d0/30.d0+(1.d0/42.d0
c$$$     & +(-1.d0/30.d0+(5.d0/66.d0+(-691.d0/2730.d0+(7.d0/6.d0
c$$$     & +(-3617.d0/510.d0+(43867.d0/798.d0+(-174611.d0/330.d0
c$$$     & +(854513.d0/138.d0+(-236364091.d0/2730.d0+(8553103.d0/6
c$$$     & -23749461029.d0/870*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y2)
c$$$     & *y2)*y2)*y)*y)*y
c$$$
c$$$      res=t+t0
c$$$
c$$$      return
c$$$      end
c$$$*
c$$$************************************************************************
c$$$*
c$$$* ..psi''(cn) (2nd derivative) for complex argument
c$$$*
c$$$************************************************************************
c$$$      subroutine psi2(cn,res)
c$$$      real*8 r
c$$$      complex*16 cn,y,y2,t,t0,res
c$$$
c$$$      t=dcmplx(0.d0,0.d0)
c$$$2     r=dsqrt(dreal(cn)**2+dimag(cn)**2)
c$$$      if(r.gt.16.d0) goto 1
c$$$      t=t-2.d0/cn**3
c$$$      cn=cn+1.d0
c$$$      goto 2
c$$$1     y=1.d0/cn
c$$$      y2=y*y
c$$$
c$$$      t0 = 
c$$$     & (-1.d0+(-1.d0+(-1.d0/2.d0+(1.d0/6.d0+(-1.d0/6.d0+(3.d0/10.d0
c$$$     & +(-5.d0/6.d0+(691.d0/210.d0+(-35.d0/2.d0+(3617.d0/30.d0
c$$$     & +(-43867.d0/42.d0+(1222277.d0/110.d0+(-854513.d0/6.d0
c$$$     & +(1181820455.d0/546-76977927.d0/2*y2)*y2)*y2)*y2)*y2)*y2)*y2)
c$$$     & *y2)*y2)*y2)*y2)*y2)*y)*y)*y2;
c$$$
c$$$      res=t+t0
c$$$
c$$$      return
c$$$      end
c$$$*
c$$$************************************************************************
c$$$*
c$$$* ..psi'''(cn) (3rd derivative) for complex argument
c$$$*
c$$$************************************************************************
c$$$      subroutine psi3(cn,res)
c$$$      real*8 r
c$$$      complex*16 cn,y,y2,t,t0,res
c$$$
c$$$      t=dcmplx(0.d0,0.d0)
c$$$2     r=dsqrt(dreal(cn)**2+dimag(cn)**2)
c$$$      if(r.gt.16.d0) goto 1
c$$$      t=t+6.d0/cn**4
c$$$      cn=cn+1.d0
c$$$      goto 2
c$$$1     y=1.d0/cn
c$$$      y2=y*y
c$$$
c$$$      t0 = 
c$$$     & (2.d0+(3.d0+(2.d0+(-1.d0+(4.d0/3.d0+(-3.d0+(10.d0+(-691.d0/15.d0
c$$$     & +(280.d0+(-10851.d0/5.d0+(438670.d0/21.d0+(-1222277.d0/5.d0
c$$$     & +(3418052.d0+(-1181820455.d0/21.d0+1077690978*y2)*y2)*y2)*y2)
c$$$     & *y2)*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y)*y)*y2*y;
c$$$
c$$$      res=t+t0
c$$$
c$$$      return
c$$$      end
c$$$*
c$$$************************************************************************
c$$$*
c$$$* ..psi''''(cn) (4th derivative) for complex argument
c$$$*
c$$$************************************************************************
c$$$      subroutine psi4(cn,res)
c$$$      real*8 r
c$$$      complex*16 cn,y,y2,t,t0,res
c$$$
c$$$      t=dcmplx(0.d0,0.d0)
c$$$2     r=dsqrt(dreal(cn)**2+dimag(cn)**2)
c$$$      if(r.gt.16.d0) goto 1
c$$$      t=t-24.d0/cn**5
c$$$      cn=cn+1.d0
c$$$      goto 2
c$$$1     y=1.d0/cn
c$$$      y2=y*y
c$$$
c$$$      t0 = 
c$$$     & (-6.d0+(-12.d0+(-10.d0+(7.d0+(-12.d0+(33.d0+(-130.d0+(691.d0
c$$$     & +(-4760.d0+(206169.d0/5.d0+(-438670.d0+(28112371.d0/5.d0
c$$$     & +(-85451300.d0+10636384095.d0/7.d0*y2)*y2)*y2)*y2)*y2)*y2)
c$$$     & *y2)*y2)*y2)*y2)*y2)*y)*y)*y2*y2;
c$$$
c$$$      res=t+t0
c$$$
c$$$      return
c$$$      end
c$$$*
c$$$************************************************************************
c$$$*
c$$$* ..psi'''''(cn) (5th derivative) for complex argument
c$$$*
c$$$************************************************************************
c$$$      subroutine psi5(cn,res)
c$$$      real*8 r
c$$$      complex*16 cn,y,y2,t,t0,res
c$$$
c$$$      t=dcmplx(0.d0,0.d0)
c$$$2     r=dsqrt(dreal(cn)**2+dimag(cn)**2)
c$$$      if(r.gt.16.d0) goto 1
c$$$      t=t+120.d0/cn**6
c$$$      cn=cn+1.d0
c$$$      goto 2
c$$$1     y=1.d0/cn
c$$$      y2=y*y
c$$$
c$$$      t0 = 
c$$$     & (24.d0+(60.d0+(60.d0+(-56.d0+(120.d0+(-396.d0+(1820.d0
c$$$     & +(-11056.d0+(85680.d0+(-824676.d0+(9650740.d0
c$$$     & +(-674696904.d0/5.d0+(2221733800.d0-42545536380.d0*y2)*y2)*y2)
c$$$     & *y2)*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y)*y)*y2*y2*y;
c$$$
c$$$      res=t+t0
c$$$
c$$$      return
c$$$      end
c$$$*
c$$$************************************************************************
c$$$*
c$$$*     C31P and C32P return the N-space coefficient functions for F3
c$$$*     at NLO and NNLO respectively
c$$$*
c$$$*     C22P returnS the N-space non-singlet coefficient functions for F2
c$$$*     at NNLO
c$$$*
c$$$*     CL2P returnS the N-space non-singlet coefficient functions for FL
c$$$*     at NNLO
c$$$*
c$$$*     Reference: hep-ph/9907472
c$$$*
c$$$************************************************************************
c$$$      SUBROUTINE C31P(N,C31PNS,C31G)
c$$$*
c$$$      IMPLICIT NONE
c$$$*
c$$$      include "../commons/consts.h"
c$$$      include "../commons/colfact.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      DOUBLE COMPLEX N
c$$$**
c$$$*     Internal Variables
c$$$*
c$$$      DOUBLE COMPLEX N1,NS
c$$$      DOUBLE COMPLEX S1,S2
c$$$      DOUBLE COMPLEX PSI,DPSI
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      DOUBLE COMPLEX C31PNS,C31G
c$$$*
c$$$      N1 = N + 1D0
c$$$      NS = N * N
c$$$*
c$$$      S1 = EMC + PSI(N1)
c$$$      S2 = ZETA2 - DPSI(N1,1)
c$$$*
c$$$      C31PNS = CF * ( 2.D0 * S1**2.D0 - 2.D0 * S2 + 3.D0 * S1
c$$$     1              - 2.D0 * S1 / ( N * N1 ) + 3.D0 / N + 4D0 / N1 
c$$$     2              + 2.D0 / NS - 9.D0 - ( 4.D0 * N + 2.D0 ) 
c$$$     3              / ( N * N1 ) )
c$$$*
c$$$      C31G = (0D0,0D0)
c$$$*
c$$$      RETURN
c$$$      END
c$$$*
c$$$************************************************************************
c$$$c      SUBROUTINE C32P(N,NF,C32PNSP,C32PNSM)
c$$$      SUBROUTINE C32P(N,NF,C32PNSP,C32PNSM,C32PS,C32G)
c$$$*
c$$$      IMPLICIT NONE
c$$$*
c$$$      include "../commons/consts.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      INTEGER NF
c$$$      DOUBLE COMPLEX N
c$$$**
c$$$*     Internal Variables
c$$$*
c$$$      DOUBLE COMPLEX N1,NS,NC,NQ
c$$$      DOUBLE COMPLEX S1,S2,S3,S4
c$$$      DOUBLE COMPLEX PSI,DPSI
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      DOUBLE COMPLEX C32PNSP,C32PNSM,C32PS,C32G
c$$$*
c$$$      N1 = N + 1D0
c$$$      NS = N * N
c$$$      NC = NS * N
c$$$      NQ = NC * N
c$$$*
c$$$      S1 = EMC + PSI(N1)
c$$$      S2 = ZETA2 - DPSI(N1,1)
c$$$      S3 = ZETA3 + DPSI(N1,2) / 2D0
c$$$      S4 = ZETA4 - DPSI(N1,3) / 6D0
c$$$*
c$$$      C32PNSP = - 0.59259D0 * ( 2D0 * S3 + 3D0 * S2 * S1 
c$$$     1        + S1**3D0 ) - 4.2963D0 * ( S2 + S1**2 ) 
c$$$     2        - 6.3489D0 * S1 - 0.042D0 * ( 2D0 * S3 
c$$$     3        + 3D0 * S2 * S1 + S1**3D0 ) / N 
c$$$     4        + 0.96978D0 * S1**2D0 / N  + ( 9.684D0 / NS 
c$$$     5        - 16.4074D0 / N ) * S1 + 10.6538D0 * S2 / N 
c$$$     6        + 4.414D0 / NC - 8.683D0 / NS - 15.9177D0 / N 
c$$$     7        - 14.97D0 / N1 + 46.856D0
c$$$      C32PNSP = NF * C32PNSP
c$$$     1        + 3.55555D0 * ( 6D0 * S4 + 8D0 * S3 * S1
c$$$     2        + 3D0 * S2**2D0 + 6D0 * S2 * S1**2D0 + S1**4D0 )
c$$$     3        + 20.4444D0 * ( 2D0 * S3 + 3D0 * S2 *S1 + S1**3D0 )
c$$$     4        - 15.5525D0 * ( S2 + S1**2D0 ) - 188.64D0 * S1
c$$$     5        + 186.816D0 * S3 / N + ( 92.43D0 / NS 
c$$$     6        + 33.2767D0 / N ) * S2 + 187.793D0 * S1 * S2 / N
c$$$     7        + 0.9778D0 * S1**3D0 / N + ( 92.42D0 / NS 
c$$$     8        + 33.2767D0 / N ) * S1**2D0 + 123.121D0 * S1 / N
c$$$     9        + 18.294D0 / NQ - 60.28D0 / NC + 79.14D0 / NS
c$$$     1        - 276.473D0 / N - 467.2D0 / N1 - 338.681D0
c$$$*
c$$$      C32PNSM = - 0.59259D0 * ( 2D0 * S3 + 3D0 * S2 * S1 
c$$$     2        + S1**3D0 ) - 4.2963D0 * ( S2 + S1**2 ) 
c$$$     3        - 6.3489D0 * S1 - 0.042D0 * ( 2D0 * S3 
c$$$     4        + 3D0 * S2 * S1 + S1**3D0 ) / N 
c$$$     5        + 0.96978D0 * S1**2D0 / N  + ( 9.684D0 / NS 
c$$$     6        - 16.4074D0 / N ) * S1 + 10.6538D0 * S2 / N 
c$$$     7        + 4.414D0 / NC - 8.683D0 / NS - 15.9177D0 / N 
c$$$     8       - 14.97D0 / N1 + 46.856D0
c$$$      C32PNSM = NF * C32PNSM
c$$$     1        + 3.55555D0 * ( 6D0 * S4 + 8D0 * S3 * S1
c$$$     1        + 3D0 * S2**2D0 + 6D0 * S2 * S1**2D0 + S1**4D0 )
c$$$     2        + 20.4444D0 * ( 2D0 * S3 + 3D0 * S2 *S1 + S1**3D0 )
c$$$     3        - 15.5525D0 * ( S2 + S1**2D0 ) - 188.64D0 * S1
c$$$     4        + 297.756D0 * S3 / N + ( 147.9D0 / NS 
c$$$     5        + 33.2767D0 / N ) * S2 + 298.733D0 * S1 * S2 / N
c$$$     6        + 0.9778D0 * S1**3D0 / N + ( 147.9D0 / NS 
c$$$     7        + 33.2767D0 / N ) * S1**2D0 - 45.8683D0 * S1 / N
c$$$     8        + 23.532D0 / NQ - 66.62D0 / NC + 67.6D0 / NS
c$$$     9        - 373.029D0 / N - 576.8D0 / N1 - 338.625D0
c$$$*
c$$$      C32PS = (0D0,0D0)
c$$$      C32G  = (0D0,0D0)
c$$$*
c$$$      RETURN
c$$$      END
c$$$*
c$$$************************************************************************
c$$$      SUBROUTINE C22PSGNS(N,NF,C22PNSM)
c$$$*
c$$$      IMPLICIT NONE
c$$$*
c$$$      include "../commons/consts.h"
c$$$      include "../commons/colfact.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      INTEGER NF
c$$$      DOUBLE COMPLEX N
c$$$**
c$$$*     Internal Variables
c$$$*
c$$$      DOUBLE COMPLEX N1,NS,NC,NQ
c$$$      DOUBLE COMPLEX S1,S2,S3,S4
c$$$      DOUBLE COMPLEX PSI,DPSI
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      DOUBLE COMPLEX C22PNSM
c$$$*
c$$$      N1  = N + 1D0
c$$$      NS  = N * N
c$$$      NC  = NS * N
c$$$      NQ  = NC * N
c$$$*
c$$$      S1   = EMC + PSI(N1)
c$$$      S2   = ZETA2 - DPSI(N1,1)
c$$$      S3   = ZETA3 + DPSI(N1,2) / 2D0
c$$$      S4   = ZETA4 - DPSI(N1,3) / 6D0
c$$$*
c$$$*     C22PNSM is the the non-singlet minus coefficient functions, respectively.
c$$$*     It has been taken from Appendix A of hep-ph/9907472.
c$$$*
c$$$      C22PNSM = 3.55555D0 * ( 6D0 * S4 + 8D0 * S3 * S1
c$$$     1        + 3D0 * S2**2D0 + 6D0 * S2 * S1**2D0 + S1**4D0 )
c$$$     2        + 20.4444D0 * ( 2D0 * S3 + 3D0 * S2 *S1 + S1**3D0 )
c$$$     3        - 15.5525D0 * ( S2 + S1**2D0 ) - 188.64D0 * S1
c$$$     4        + 229.916D0 * S3 / N + ( 31.58D0 / NS
c$$$     5        + 9.7467D0 / N ) * S2 + 393.703D0 * S1 * S2 / N
c$$$     6        + 2.9678D0 * S1**3D0 / N + ( 192.4D0 / NS 
c$$$     7        + 9.7467D0 / N ) * S1**2D0 - ( 160.82D0 / NC 
c$$$     8        - 61.1321D0 / N ) * S1
c$$$     9        + 22.488D0 / NQ - 39.12D0 / NC + 265.774D0 / NS
c$$$     1        - 164.777D0 / N - 1010D0 / N1 - 337.992D0
c$$$     2 + NF * ( - 0.59259D0 * ( 2D0 * S3 + 3D0 * S2 * S1 
c$$$     3        + S1**3D0 ) - 4.2963D0 * ( S2 + S1**2 ) 
c$$$     4        - 6.3489D0 * S1 - 6.072 * S3 / N - ( 6.072D0 / NS 
c$$$     5        - 18.0408D0 / N ) * S2 - ( 6.072D0 / NC - 17.97D0 / NS 
c$$$     6        + 14.3574D0 / N ) * S1 + 0.07078D0 * S1**2D0 / N
c$$$     7        + 4.488D0 / NC + 4.21808D0 / NS - 21.6028D0 / N
c$$$     8        - 37.91D0 / N1 + 46.8406D0 )
c$$$*
c$$$      RETURN
c$$$      END
c$$$*
c$$$************************************************************************
c$$$      SUBROUTINE CL2PSGNS(N,NF,CL2PNSM)
c$$$*
c$$$      IMPLICIT NONE
c$$$*
c$$$      include "../commons/consts.h"
c$$$      include "../commons/colfact.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      INTEGER NF
c$$$      DOUBLE COMPLEX N
c$$$**
c$$$*     Internal Variables
c$$$*
c$$$      DOUBLE COMPLEX N1,NS,N1S,NC,N1C
c$$$      DOUBLE COMPLEX S1,S2
c$$$      DOUBLE COMPLEX PSI,DPSI
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      DOUBLE COMPLEX CL2PNSM
c$$$*
c$$$      N1  = N + 1D0
c$$$      NS  = N * N
c$$$      N1S = N1 * N1
c$$$      NC  = NS * N
c$$$      N1C = N1S * N1
c$$$*
c$$$      S1 = EMC + PSI(N1)
c$$$      S2 = ZETA2 - DPSI(N1,1)
c$$$*
c$$$*     C22PNSM is the the non-singlet minus coefficient functions, respectively.
c$$$*     It has been taken from Appendix A of hep-ph/9907472.
c$$$*
c$$$      CL2PNSM = - 128.4D0 * S2 / N + 13.3D0 * S1**2D0 / N
c$$$     1        + ( 59.12D0 / N - 141.7D0 / NS ) * S1 - 0.086D0 / NC
c$$$     2        + 22.21D0 / NS + 180.818D0 / N + 46.58D0 / N1C
c$$$     3        + 100.8D0 / N1 - 0.150D0
c$$$     4 + NF * 16D0 * ( - 6D0 * S1 / N1 + 6D0 / N + 6D0 / N1S 
c$$$     5                 - 25D0 / N1 ) / 27D0
c$$$*
c$$$      RETURN
c$$$      END
