************************************************************************
*
*     The following routines return the netral current coefficient
*     functions in Mellin space in the Zero-Mass Variable-Flavor-Number
*     scheme.
*
************************************************************************
*
* ..File: c2lmom.f
*
*
* ..The subroutines c21p, c22p, c23p, cl1p, cl2p and cl3p
*    for the (complex) moments of the parametrized non-singlet 
*    and singlet MS(bar) coefficient functions for the unpolarized 
*    deep-inelastic structure functions F_2 and F_L, mu_r = mu_f.
*
* ..The QCD colour factors have been hard-wired in the parametrizations.
*    The coupling constant is normalized as  a_s = alpha_s/(4*pi). 
*    
* ..The required routines for simple harmonic sums (expressed through 
*    Euler's psi function and its derivatives) are appended at the end 
*    of the file.
* 
* ..References: hep-ph/0411112 = Phys. Lett. B606 (2005) 123 
*               hep-ph/0504242 = Nucl. Phys. B724 (2005) 3 
*
************************************************************************
      subroutine c21p(cn,c2pls1n,c2g1n)
      IMPLICIT complex*16 (A-H,O-Z)
*
      real*8 dl,pi,ge,z2,z3,z4,z5,z6

      pi = 3.1415926535897932385d0
      dl = dlog(2.0d0)
      ge = 5.7721566490153286061d-1
      z2 = 1.6449340668482264365d0
      z3 = 1.2020569031595942854d0
      z4 = 1.0823232337111381916d0
      z5 = 1.0369277551433699263d0
      z6 = 1.0173430619844491398d0

      c2pls1n=+4.d0/3.d0
     &  *(-9.d0-2.d0*(cn+1)**(-2)+4.d0*(cn+1)**(-1)+2.d0*(cn)**(-2)+3.
     &  d0*(cn)**(-1)+2.d0*s1(1,1.d0+cn)*(cn+1)**(-1)+3.d0*s1(1,cn)-2.d0
     &    *s1(1,cn)*(cn)**(-1)+2.d0*s1(1,cn)**2-2.d0*s1(2,cn))

      c2g1n=
     &    4.d0*(cn+2)**(-2)-16.d0*(cn+2)**(-1)-4.d0*(cn+1)**(-2)+16.d0
     &    *(cn+1)**(-1)+2.d0*(cn)**(-2)-2.d0*(cn)**(-1)+4.d0*s1(1,1.d0+
     &    cn)*(cn+1)**(-1)-4.d0*s1(1,2.d0+cn)*(cn+2)**(-1)-2.d0*s1(1,cn
     &    )*(cn)**(-1)
*
      return
      end
*
************************************************************************
      subroutine c22p(cn,nf,c2pls2n,c2min2n,c2ps2n,c2g2n)
      IMPLICIT complex*16 (A-H,O-Z)
*
      integer nf
      real*8 dl,pi,ge,z2,z3,z4,z5,z6

      pi = 3.1415926535897932385d0
      dl = dlog(2.0d0)
      ge = 5.7721566490153286061d-1
      z2 = 1.6449340668482264365d0
      z3 = 1.2020569031595942854d0
      z4 = 1.0823232337111381916d0
      z5 = 1.0369277551433699263d0
      z6 = 1.0173430619844491398d0

      c2pls2n = 
     &     + nf * ( + 216883.d0/4629.d0 - 8113.d0/1000.d0*z2*cn**(-1)
     &     + 40.d0/9.d0*cn**(-3) - 16.d0/3.d0*cn**(-2)
     &     - 73101.d0/50000.d0*cn**(-1) + 111.d0/100.d0*(cn+1.d0)**(-4)
     &     - 891.d0/50.d0*(cn+1.d0)**(-1)
     &     - 1297.d0/100.d0*(cn+2.d0)**(-1) - 79361.d0/12500.d0*s1(1,cn)
     &     + 8113.d0/1000.d0*s1(1,cn)*cn**(-2)
     &     - 43949.d0/2700.d0*s1(1,cn)*cn**(-1)
     &     - 116.d0/27.d0*s1(1,cn)**2 + 5.d0/18.d0*s1(1,cn)**2*cn**(-1)
     &     - 16.d0/27.d0*s1(1,cn)**3 - 16.d0/9.d0*s1(1,cn)*s1(2,cn)
     &     - 116.d0/27.d0*s1(2,cn) + 75517.d0/9000.d0*s1(2,cn)*cn**(-1)
     &     - 32.d0/27.d0*s1(3,cn) )
      c2pls2n = c2pls2n
     &     - 338513.d0/1000.d0 + 151.d0/2.d0*z2*cn**(-2)
     &     - 2187.d0/10.d0*z3*cn**(-1) + 160.d0/9.d0*cn**(-4)
     &     - 207.d0/5.d0*cn**(-3) + 3548.d0/125.d0*cn**(-2)
     &     + 7641.d0/1000.d0*cn**(-1) + 2157.d0/125.d0*(cn+1.d0)**(-5)
     &     - 8067.d0/10.d0*(cn+1.d0)**(-1) - 188641.d0/1000.d0*s1(1,cn)
     &     - 1471.d0/5.d0*s1(1,cn)*z2*cn**(-1)
     &     - 151.d0/2.d0*s1(1,cn)*cn**(-3)
     &     + 1649763.d0/2500.d0*s1(1,cn)*cn**(-1)
     &     - 77763.d0/5000.d0*s1(1,cn)**2
     &     + 1471.d0/10.d0*s1(1,cn)**2*cn**(-2)
      c2pls2n = c2pls2n
     &     + 818.d0/75.d0*s1(1,cn)**2*cn**(-1) + 184.d0/9.d0*s1(1,cn)**3
     &     + 1583.d0/450.d0*s1(1,cn)**3*cn**(-1)
     &     + 32.d0/9.d0*s1(1,cn)**4 + 64.d0/3.d0*s1(1,cn)**2*s1(2,cn)
     &     + 184.d0/3.d0*s1(1,cn)*s1(2,cn)
     &     + 45713.d0/150.d0*s1(1,cn)*s1(2,cn)*cn**(-1)
     &     + 256.d0/9.d0*s1(1,cn)*s1(3,cn) - 77763.d0/5000.d0*s1(2,cn)
     &     + 358.d0/5.d0*s1(2,cn)*cn**(-2)
     &     + 818.d0/75.d0*s1(2,cn)*cn**(-1) + 32.d0/3.d0*s1(2,cn)**2
     &     + 368.d0/9.d0*s1(3,cn) + 101581.d0/450.d0*s1(3,cn)*cn**(-1)
     &     + 64.d0/3.d0*s1(4,cn)

      c2ps2n =
     &     + ( - 1013.d0/25.d0*z2*cn**(-2)
     &     + 3393.d0/50.d0*z2*(cn+1.d0)**(-2)
     &     - 1013.d0/25.d0*z3*cn**(-1)
     &     + 3393.d0/50.d0*z3*(cn+1.d0)**(-1) - 80.d0/3.d0*cn**(-4)
     &     - 1053.d0/500.d0*cn**(-3) - 1084.d0/25.d0*cn**(-2)
     &     + 88831.d0/10000.d0*cn**(-1)
     &     + 52903.d0/10000.d0*(cn-1.d0)**(-1)
     &     - 17001.d0/500.d0*(cn+1.d0)**(-4)
     &     + 1084.d0/25.d0*(cn+1.d0)**(-2)
     &     - 185767.d0/5000.d0*(cn+1.d0)**(-1)
     &     + 3639.d0/100.d0*(cn+2.d0)**(-1)
     &     - 1341.d0/100.d0*(cn+3.d0)**(-1) )
      c2ps2n = c2ps2n
     &     + ( - 3393.d0/50.d0*s1(1,cn+1.d0)*(cn+1.d0)**(-3)
     &     - 32.d0/3.d0*s1(1,cn+1.d0)*(cn+1.d0)**(-1)
     &     - 8.d0/3.d0*s1(1,cn+1.d0)**2*(cn+1.d0)**(-1)
     &     + 4.d0/25.d0*s1(1,cn+1.d0)**3*(cn+1.d0)**(-1)
     &     + 12.d0/25.d0*s1(1,cn+1.d0)*s1(2,cn+1.d0)*(cn+1.d0)**(-1)
     &     - 2.d0/25.d0*s1(1,cn+2.d0)**3*(cn+2.d0)**(-1)
     &     - 6.d0/25.d0*s1(1,cn+2.d0)*s1(2,cn+2.d0)*(cn+2.d0)**(-1)
     &     + 1013.d0/25.d0*s1(1,cn)*cn**(-3)
     &     + 32.d0/3.d0*s1(1,cn)*cn**(-1)
     &     + 8.d0/3.d0*s1(1,cn)**2*cn**(-1)
     &     - 2.d0/25.d0*s1(1,cn)**3*cn**(-1)
     &     - 6.d0/25.d0*s1(1,cn)*s1(2,cn)*cn**(-1)
     &     - 3393.d0/50.d0*s1(2,cn+1.d0)*(cn+1.d0)**(-2)
     &     - 8.d0/3.d0*s1(2,cn+1.d0)*(cn+1.d0)**(-1) )
      c2ps2n = c2ps2n
     &     + ( + 1013.d0/25.d0*s1(2,cn)*cn**(-2)
     &     + 8.d0/3.d0*s1(2,cn)*cn**(-1)
     &     - 3377.d0/50.d0*s1(3,cn+1.d0)*(cn+1.d0)**(-1)
     &     - 4.d0/25.d0*s1(3,cn+2.d0)*(cn+2.d0)**(-1)
     &     + 1009.d0/25.d0*s1(3,cn)*cn**(-1) )

      c2g2n = 
     &     + ( + 2086.d0*z2*cn**(-2) - 1469.d0/5.d0*z2*cn**(-1)
     &     - 3556.d0/5.d0*z2*(cn+1.d0)**(-1) + 2086.d0*z3*cn**(-1)
     &     - 140.d0/3.d0*cn**(-4) - 7109.d0/500.d0*cn**(-3)
     &     - 578.d0/5.d0*cn**(-2) - 63973.d0/10000.d0*cn**(-1)
     &     + 119033.d0/10000.d0*(cn-1.d0)**(-1)
     &     + 2408.d0*(cn+1.d0)**(-3) - 36761.d0/50.d0*(cn+1.d0)**(-1)
     &     + 7603.d0/10.d0*(cn+2.d0)**(-1)
     &     + 3556.d0/5.d0*s1(1,cn+1.d0)*(cn+1.d0)**(-2)
     &     - 593.d0/20.d0*s1(1,cn+1.d0)**3*(cn+1.d0)**(-1) )
      c2g2n = c2g2n
     &     + ( - 1779.d0/20.d0*s1(1,cn+1.d0)*s1(2,cn+1.d0)*
     &    (cn+1.d0)**(-1) - 2086.d0*s1(1,cn)*cn**(-3)
     &     + 1469.d0/5.d0*s1(1,cn)*cn**(-2)
     &     + 872.d0/25.d0*s1(1,cn)*cn**(-1) - 24.d0*s1(1,cn)**2*cn**(-1)
     &     + 4177.d0/180.d0*s1(1,cn)**3*cn**(-1)
     &     + 4177.d0/60.d0*s1(1,cn)*s1(2,cn)*cn**(-1)
     &     + 3556.d0/5.d0*s1(2,cn+1.d0)*(cn+1.d0)**(-1)
     &     - 2086.d0*s1(2,cn)*cn**(-2) + 1349.d0/5.d0*s1(2,cn)*cn**(-1)
     &     - 593.d0/10.d0*s1(3,cn+1.d0)*(cn+1.d0)**(-1)
     &     - 183563.d0/90.d0*s1(3,cn)*cn**(-1) )

      call C22PSGNS(cn,nf,c2min2n)
*
      return
      end
*
************************************************************************
      subroutine cl1p(cn,clpls1n,clg1n)
      IMPLICIT complex*16 (A-H,O-Z)
*
      clpls1n=
     &  + 16.d0/3.d0*(cn+1)**(-1)

      clg1n=( 8.d0*(cn+1)**(-1) - 8.d0*(cn+2)**(-1) )
*
      return
      end
*
************************************************************************
      subroutine cl2p(cn,nf,clpls2n,clmin2n,clps2n,clg2n)
      IMPLICIT complex*16 (A-H,O-Z)
*
      integer nf
      real*8 dl,pi,ge,z2,z3,z4,z5,z6

      pi = 3.1415926535897932385d0
      dl = dlog(2.0d0)
      ge = 5.7721566490153286061d-1
      z2 = 1.6449340668482264365d0
      z3 = 1.2020569031595942854d0
      z4 = 1.0823232337111381916d0
      z5 = 1.0369277551433699263d0
      z6 = 1.0173430619844491398d0

      clpls2n=
     &  + nf * ( 32.d0/9.d0*cn**(-1) + 64.d0/9.d0*(cn+1)**(-2) - 400.d0/
     &    27.d0*(cn+1)**(-1) - 32.d0/9.d0*s1(1,1.d0 + cn)*(cn+1)**(-1)
     &     )
      clpls2n=clpls2n
     &     - 3.d0/250.d0 + 42047.d0/500.d0*z2*cn**(-1) + 128.d0
     &    /9.d0*cn**(-2) - 18669.d0/500.d0*cn**(-1) + 1841.d0/50.d0*
     &    (cn+1)**(-3) - 329.d0/10.d0*(cn+1)**(-2) + 8953.d0/100.d0*
     &    (cn+1)**(-1) + 1691.d0/50.d0*(cn+2)**(-1) + 93.d0/2.d0*s1(1,1.
     &    d0 + cn)*(cn+1)**(-1) + 128.d0/9.d0*s1(1,1.d0 + cn)**2*
     &    (cn+1)**(-1) - 42047.d0/500.d0*s1(1,cn)*cn**(-2) + 128.d0/9.d0
     &    *s1(2,1.d0 + cn)*(cn+1)**(-1) - 42047.d0/500.d0*s1(2,cn)*
     &    cn**(-1)

      clps2n=
     &  + ( 421.d0/500.d0*cn**(-3) - 2809.d0/100.d0*cn**(-2) + 
     &    1319.d0/50.d0*cn**(-1) - 237.d0/100.d0*(cn-1)**(-1) + 76.d0/
     &    25.d0*(cn+1)**(-3) + 2809.d0/100.d0*(cn+1)**(-2) - 1623.d0/25.
     &    d0*(cn+1)**(-1) + 3009.d0/50.d0*(cn+2)**(-1) - 1927.d0/100.d0
     &    *(cn+3)**(-1) + 9273.d0/250.d0*s1(1,1.d0 + cn)*(cn+1)**(-1)
     &     - 6591.d0/250.d0*s1(1,2.d0 + cn)*(cn+2)**(-1) + 1303.d0/250.d
     &    0*s1(1,3.d0 + cn)*(cn+3)**(-1) - 797.d0/50.d0*s1(1,cn)*
     &    cn**(-1) )

      clg2n=
     &  + (  - 1161.d0*z2*(cn+1)**(-1) - 1983.d0/50.d0*cn**(-2) + 
     &    5333.d0/1000.d0*cn**(-1) - 5333.d0/1000.d0*(cn-1)**(-1) + 
     &    3003.d0/25.d0*(cn+1)**(-3) + 1983.d0/50.d0*(cn+1)**(-2) + 
     &    1161.d0*s1(1,1.d0 + cn)*(cn+1)**(-2) + 4324.d0/5.d0*s1(1,1.d0
     &     + cn)*(cn+1)**(-1) - 7197.d0/50.d0*s1(1,1.d0 + cn)**2*
     &    (cn+1)**(-1) + 246.d0/5.d0*s1(1,2.d0 + cn)**2*(cn+2)**(-1) - 
     &    4324.d0/5.d0*s1(1,cn)*cn**(-1) + 4737.d0/50.d0*s1(1,cn)**2*
     &    cn**(-1) + 50853.d0/50.d0*s1(2,1.d0 + cn)*(cn+1)**(-1) + 246.d
     &    0/5.d0*s1(2,2.d0 + cn)*(cn+2)**(-1) + 4737.d0/50.d0*s1(2,cn)*
     &    cn**(-1) )

      call CL2PSGNS(cn,nf,clmin2n)
*
      return
      end
*
************************************************************************
*
* ..harmonic sums and the like
*
************************************************************************
      complex*16 function s1(k,cn)
      IMPLICIT complex*16 (A-H,O-Z)
      integer k

      real*8 dl,pi,ge,z2,z3,z4,z5,z6

      pi = 3.1415926535897932385d0
      dl = dlog(2.0d0)
      ge = 5.7721566490153286061d-1
      z2 = 1.6449340668482264365d0
      z3 = 1.2020569031595942854d0
      z4 = 1.0823232337111381916d0
      z5 = 1.0369277551433699263d0
      z6 = 1.0173430619844491398d0

      if (k.eq.1) then
        cn1=cn+1
        call psi0(cn1,res)
        s1=ge + res
      end if
      if (k.eq.2) then
        cn1=cn+1
        call psi1(cn1,res)
        s1=z2 - res
      end if
      if (k.eq.3) then
        cn1=cn+1
        call psi2(cn1,res)
        s1=z3 + res/2.d0
      end if
      if (k.eq.4) then
        cn1=cn+1
        call psi3(cn1,res)
        s1=z4 - res/6.d0
      end if
      if (k.eq.5) then
        cn1=cn+1
        call psi4(cn1,res)
        s1=z5 + res/24.d0
      end if
      if (k.eq.6) then
        cn1=cn+1
        call psi5(cn1,res)
        s1=z6 - res/120.d0
      end if
      return
      end
*
************************************************************************
*
* ..psi(cn) for complex argument
*
************************************************************************
      subroutine psi0(cn,res)
      real*8 r
      complex*16 cn,t,t0,y,y2,res

      t=dcmplx(0.d0,0.d0)
2     r=dsqrt(dreal(cn)**2+dimag(cn)**2)
      if(r.gt.16.d0) goto 1
      t=t-1.d0/cn
      cn=cn+1.d0
      goto 2

1     y=1.d0/cn
      y2=y*y

      t0 = 
     & (-1.d0/2.d0+(-1.d0/12.d0+(1.d0/120.d0+(-1.d0/252.d0+(1.d0/240.d0
     & +(-1.d0/132.d0+(691.d0/32760.d0+(-1.d0/12.d0+(3617.d0/8160.d0
     & +(-43867.d0/14364.d0+(174611.d0/6600.d0+(-77683.d0/276.d0
     & +(236364091.d0/65520.d0+(-657931.d0/12.d0+3392780147.d0/3480.d0
     & *y2)*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y)*y-log(y)

      res=t+t0

      return
      end
*
************************************************************************
*
* ..psi'(cn) (1st derivative) for complex argument
*
************************************************************************
      subroutine psi1(cn,res)
      real*8 r
      complex*16 cn,t,t0,y,y2,res

      t=dcmplx(0.d0,0.d0)
2     r=dsqrt(dreal(cn)**2+dimag(cn)**2)
      if(r.gt.16.d0) goto 1
      t=t+1.d0/cn**2
      cn=cn+1.d0
      goto 2
1     y=1.d0/cn
      y2=y*y

      t0 = 
     & (1.d0+(1.d0/2.d0+(1.d0/6.d0+(-1.d0/30.d0+(1.d0/42.d0
     & +(-1.d0/30.d0+(5.d0/66.d0+(-691.d0/2730.d0+(7.d0/6.d0
     & +(-3617.d0/510.d0+(43867.d0/798.d0+(-174611.d0/330.d0
     & +(854513.d0/138.d0+(-236364091.d0/2730.d0+(8553103.d0/6
     & -23749461029.d0/870*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y2)
     & *y2)*y2)*y)*y)*y

      res=t+t0

      return
      end
*
************************************************************************
*
* ..psi''(cn) (2nd derivative) for complex argument
*
************************************************************************
      subroutine psi2(cn,res)
      real*8 r
      complex*16 cn,y,y2,t,t0,res

      t=dcmplx(0.d0,0.d0)
2     r=dsqrt(dreal(cn)**2+dimag(cn)**2)
      if(r.gt.16.d0) goto 1
      t=t-2.d0/cn**3
      cn=cn+1.d0
      goto 2
1     y=1.d0/cn
      y2=y*y

      t0 = 
     & (-1.d0+(-1.d0+(-1.d0/2.d0+(1.d0/6.d0+(-1.d0/6.d0+(3.d0/10.d0
     & +(-5.d0/6.d0+(691.d0/210.d0+(-35.d0/2.d0+(3617.d0/30.d0
     & +(-43867.d0/42.d0+(1222277.d0/110.d0+(-854513.d0/6.d0
     & +(1181820455.d0/546-76977927.d0/2*y2)*y2)*y2)*y2)*y2)*y2)*y2)
     & *y2)*y2)*y2)*y2)*y2)*y)*y)*y2;

      res=t+t0

      return
      end
*
************************************************************************
*
* ..psi'''(cn) (3rd derivative) for complex argument
*
************************************************************************
      subroutine psi3(cn,res)
      real*8 r
      complex*16 cn,y,y2,t,t0,res

      t=dcmplx(0.d0,0.d0)
2     r=dsqrt(dreal(cn)**2+dimag(cn)**2)
      if(r.gt.16.d0) goto 1
      t=t+6.d0/cn**4
      cn=cn+1.d0
      goto 2
1     y=1.d0/cn
      y2=y*y

      t0 = 
     & (2.d0+(3.d0+(2.d0+(-1.d0+(4.d0/3.d0+(-3.d0+(10.d0+(-691.d0/15.d0
     & +(280.d0+(-10851.d0/5.d0+(438670.d0/21.d0+(-1222277.d0/5.d0
     & +(3418052.d0+(-1181820455.d0/21.d0+1077690978*y2)*y2)*y2)*y2)
     & *y2)*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y)*y)*y2*y;

      res=t+t0

      return
      end
*
************************************************************************
*
* ..psi''''(cn) (4th derivative) for complex argument
*
************************************************************************
      subroutine psi4(cn,res)
      real*8 r
      complex*16 cn,y,y2,t,t0,res

      t=dcmplx(0.d0,0.d0)
2     r=dsqrt(dreal(cn)**2+dimag(cn)**2)
      if(r.gt.16.d0) goto 1
      t=t-24.d0/cn**5
      cn=cn+1.d0
      goto 2
1     y=1.d0/cn
      y2=y*y

      t0 = 
     & (-6.d0+(-12.d0+(-10.d0+(7.d0+(-12.d0+(33.d0+(-130.d0+(691.d0
     & +(-4760.d0+(206169.d0/5.d0+(-438670.d0+(28112371.d0/5.d0
     & +(-85451300.d0+10636384095.d0/7.d0*y2)*y2)*y2)*y2)*y2)*y2)
     & *y2)*y2)*y2)*y2)*y2)*y)*y)*y2*y2;

      res=t+t0

      return
      end
*
************************************************************************
*
* ..psi'''''(cn) (5th derivative) for complex argument
*
************************************************************************
      subroutine psi5(cn,res)
      real*8 r
      complex*16 cn,y,y2,t,t0,res

      t=dcmplx(0.d0,0.d0)
2     r=dsqrt(dreal(cn)**2+dimag(cn)**2)
      if(r.gt.16.d0) goto 1
      t=t+120.d0/cn**6
      cn=cn+1.d0
      goto 2
1     y=1.d0/cn
      y2=y*y

      t0 = 
     & (24.d0+(60.d0+(60.d0+(-56.d0+(120.d0+(-396.d0+(1820.d0
     & +(-11056.d0+(85680.d0+(-824676.d0+(9650740.d0
     & +(-674696904.d0/5.d0+(2221733800.d0-42545536380.d0*y2)*y2)*y2)
     & *y2)*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y)*y)*y2*y2*y;

      res=t+t0

      return
      end
*
************************************************************************
*
*     C31P and C32P return the N-space coefficient functions for F3
*     at NLO and NNLO respectively
*
*     C22P returnS the N-space non-singlet coefficient functions for F2
*     at NNLO
*
*     CL2P returnS the N-space non-singlet coefficient functions for FL
*     at NNLO
*
*     Reference: hep-ph/9907472
*
************************************************************************
      SUBROUTINE C31P(N,C31PNS,C31G)
*
      IMPLICIT NONE
*
      include "../commons/consts.h"
      include "../commons/colfact.h"
**
*     Input Variables
*
      DOUBLE COMPLEX N
**
*     Internal Variables
*
      DOUBLE COMPLEX N1,NS
      DOUBLE COMPLEX S1,S2
      DOUBLE COMPLEX PSI,DPSI
**
*     Output Variables
*
      DOUBLE COMPLEX C31PNS,C31G
*
      N1 = N + 1D0
      NS = N * N
*
      S1 = EMC + PSI(N1)
      S2 = ZETA2 - DPSI(N1,1)
*
      C31PNS = CF * ( 2.D0 * S1**2.D0 - 2.D0 * S2 + 3.D0 * S1
     1              - 2.D0 * S1 / ( N * N1 ) + 3.D0 / N + 4D0 / N1 
     2              + 2.D0 / NS - 9.D0 - ( 4.D0 * N + 2.D0 ) 
     3              / ( N * N1 ) )
*
      C31G = (0D0,0D0)
*
      RETURN
      END
*
************************************************************************
c      SUBROUTINE C32P(N,NF,C32PNSP,C32PNSM)
      SUBROUTINE C32P(N,NF,C32PNSP,C32PNSM,C32PS,C32G)
*
      IMPLICIT NONE
*
      include "../commons/consts.h"
**
*     Input Variables
*
      INTEGER NF
      DOUBLE COMPLEX N
**
*     Internal Variables
*
      DOUBLE COMPLEX N1,NS,NC,NQ
      DOUBLE COMPLEX S1,S2,S3,S4
      DOUBLE COMPLEX PSI,DPSI
**
*     Output Variables
*
      DOUBLE COMPLEX C32PNSP,C32PNSM,C32PS,C32G
*
      N1 = N + 1D0
      NS = N * N
      NC = NS * N
      NQ = NC * N
*
      S1 = EMC + PSI(N1)
      S2 = ZETA2 - DPSI(N1,1)
      S3 = ZETA3 + DPSI(N1,2) / 2D0
      S4 = ZETA4 - DPSI(N1,3) / 6D0
*
      C32PNSP = - 0.59259D0 * ( 2D0 * S3 + 3D0 * S2 * S1 
     1        + S1**3D0 ) - 4.2963D0 * ( S2 + S1**2 ) 
     2        - 6.3489D0 * S1 - 0.042D0 * ( 2D0 * S3 
     3        + 3D0 * S2 * S1 + S1**3D0 ) / N 
     4        + 0.96978D0 * S1**2D0 / N  + ( 9.684D0 / NS 
     5        - 16.4074D0 / N ) * S1 + 10.6538D0 * S2 / N 
     6        + 4.414D0 / NC - 8.683D0 / NS - 15.9177D0 / N 
     7        - 14.97D0 / N1 + 46.856D0
      C32PNSP = NF * C32PNSP
     1        + 3.55555D0 * ( 6D0 * S4 + 8D0 * S3 * S1
     2        + 3D0 * S2**2D0 + 6D0 * S2 * S1**2D0 + S1**4D0 )
     3        + 20.4444D0 * ( 2D0 * S3 + 3D0 * S2 *S1 + S1**3D0 )
     4        - 15.5525D0 * ( S2 + S1**2D0 ) - 188.64D0 * S1
     5        + 186.816D0 * S3 / N + ( 92.43D0 / NS 
     6        + 33.2767D0 / N ) * S2 + 187.793D0 * S1 * S2 / N
     7        + 0.9778D0 * S1**3D0 / N + ( 92.42D0 / NS 
     8        + 33.2767D0 / N ) * S1**2D0 + 123.121D0 * S1 / N
     9        + 18.294D0 / NQ - 60.28D0 / NC + 79.14D0 / NS
     1        - 276.473D0 / N - 467.2D0 / N1 - 338.681D0
*
      C32PNSM = - 0.59259D0 * ( 2D0 * S3 + 3D0 * S2 * S1 
     2        + S1**3D0 ) - 4.2963D0 * ( S2 + S1**2 ) 
     3        - 6.3489D0 * S1 - 0.042D0 * ( 2D0 * S3 
     4        + 3D0 * S2 * S1 + S1**3D0 ) / N 
     5        + 0.96978D0 * S1**2D0 / N  + ( 9.684D0 / NS 
     6        - 16.4074D0 / N ) * S1 + 10.6538D0 * S2 / N 
     7        + 4.414D0 / NC - 8.683D0 / NS - 15.9177D0 / N 
     8       - 14.97D0 / N1 + 46.856D0
      C32PNSM = NF * C32PNSM
     1        + 3.55555D0 * ( 6D0 * S4 + 8D0 * S3 * S1
     1        + 3D0 * S2**2D0 + 6D0 * S2 * S1**2D0 + S1**4D0 )
     2        + 20.4444D0 * ( 2D0 * S3 + 3D0 * S2 *S1 + S1**3D0 )
     3        - 15.5525D0 * ( S2 + S1**2D0 ) - 188.64D0 * S1
     4        + 297.756D0 * S3 / N + ( 147.9D0 / NS 
     5        + 33.2767D0 / N ) * S2 + 298.733D0 * S1 * S2 / N
     6        + 0.9778D0 * S1**3D0 / N + ( 147.9D0 / NS 
     7        + 33.2767D0 / N ) * S1**2D0 - 45.8683D0 * S1 / N
     8        + 23.532D0 / NQ - 66.62D0 / NC + 67.6D0 / NS
     9        - 373.029D0 / N - 576.8D0 / N1 - 338.625D0
*
      C32PS = (0D0,0D0)
      C32G  = (0D0,0D0)
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE C22PSGNS(N,NF,C22PNSM)
*
      IMPLICIT NONE
*
      include "../commons/consts.h"
      include "../commons/colfact.h"
**
*     Input Variables
*
      INTEGER NF
      DOUBLE COMPLEX N
**
*     Internal Variables
*
      DOUBLE COMPLEX N1,NS,NC,NQ
      DOUBLE COMPLEX S1,S2,S3,S4
      DOUBLE COMPLEX PSI,DPSI
**
*     Output Variables
*
      DOUBLE COMPLEX C22PNSM
*
      N1  = N + 1D0
      NS  = N * N
      NC  = NS * N
      NQ  = NC * N
*
      S1   = EMC + PSI(N1)
      S2   = ZETA2 - DPSI(N1,1)
      S3   = ZETA3 + DPSI(N1,2) / 2D0
      S4   = ZETA4 - DPSI(N1,3) / 6D0
*
*     C22PNSM is the the non-singlet minus coefficient functions, respectively.
*     It has been taken from Appendix A of hep-ph/9907472.
*
      C22PNSM = 3.55555D0 * ( 6D0 * S4 + 8D0 * S3 * S1
     1        + 3D0 * S2**2D0 + 6D0 * S2 * S1**2D0 + S1**4D0 )
     2        + 20.4444D0 * ( 2D0 * S3 + 3D0 * S2 *S1 + S1**3D0 )
     3        - 15.5525D0 * ( S2 + S1**2D0 ) - 188.64D0 * S1
     4        + 229.916D0 * S3 / N + ( 31.58D0 / NS
     5        + 9.7467D0 / N ) * S2 + 393.703D0 * S1 * S2 / N
     6        + 2.9678D0 * S1**3D0 / N + ( 192.4D0 / NS 
     7        + 9.7467D0 / N ) * S1**2D0 - ( 160.82D0 / NC 
     8        - 61.1321D0 / N ) * S1
     9        + 22.488D0 / NQ - 39.12D0 / NC + 265.774D0 / NS
     1        - 164.777D0 / N - 1010D0 / N1 - 337.992D0
     2 + NF * ( - 0.59259D0 * ( 2D0 * S3 + 3D0 * S2 * S1 
     3        + S1**3D0 ) - 4.2963D0 * ( S2 + S1**2 ) 
     4        - 6.3489D0 * S1 - 6.072 * S3 / N - ( 6.072D0 / NS 
     5        - 18.0408D0 / N ) * S2 - ( 6.072D0 / NC - 17.97D0 / NS 
     6        + 14.3574D0 / N ) * S1 + 0.07078D0 * S1**2D0 / N
     7        + 4.488D0 / NC + 4.21808D0 / NS - 21.6028D0 / N
     8        - 37.91D0 / N1 + 46.8406D0 )
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE CL2PSGNS(N,NF,CL2PNSM)
*
      IMPLICIT NONE
*
      include "../commons/consts.h"
      include "../commons/colfact.h"
**
*     Input Variables
*
      INTEGER NF
      DOUBLE COMPLEX N
**
*     Internal Variables
*
      DOUBLE COMPLEX N1,NS,N1S,NC,N1C
      DOUBLE COMPLEX S1,S2
      DOUBLE COMPLEX PSI,DPSI
**
*     Output Variables
*
      DOUBLE COMPLEX CL2PNSM
*
      N1  = N + 1D0
      NS  = N * N
      N1S = N1 * N1
      NC  = NS * N
      N1C = N1S * N1
*
      S1 = EMC + PSI(N1)
      S2 = ZETA2 - DPSI(N1,1)
*
*     C22PNSM is the the non-singlet minus coefficient functions, respectively.
*     It has been taken from Appendix A of hep-ph/9907472.
*
      CL2PNSM = - 128.4D0 * S2 / N + 13.3D0 * S1**2D0 / N
     1        + ( 59.12D0 / N - 141.7D0 / NS ) * S1 - 0.086D0 / NC
     2        + 22.21D0 / NS + 180.818D0 / N + 46.58D0 / N1C
     3        + 100.8D0 / N1 - 0.150D0
     4 + NF * 16D0 * ( - 6D0 * S1 / N1 + 6D0 / N + 6D0 / N1S 
     5                 - 25D0 / N1 ) / 27D0
*
      RETURN
      END
