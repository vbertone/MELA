************************************************************************
*
*     Evaluation of the integral in eq. (97) of ArXiv.1001.2312 needed 
*     for the Adler sum rule to be true
*
************************************************************************
      SUBROUTINE ADLERSR_COEFF(Q)
*
      IMPLICIT NONE
*
      include "../commons/adlersr.h"
      include "../commons/massthrs.h"
**
*     Input Variables
*
      DOUBLE COMPLEX Q
**
*     Internal Variables
*
      INTEGER I
      DOUBLE PRECISION DGAUSSMELA,A,B,EPS
      DOUBLE PRECISION QQ2,MM2
      COMMON/SCALES/QQ2,MM2

      DOUBLE PRECISION L22Q
      EXTERNAL L22Q
*
      QQ2 = ABS(Q*Q)
*
      A = 1D-6
      B = 1D0
      EPS = 1D-5
*
      DO I=4,6
         MM2 = ABS(Q2TH(I))
         ADLERSR(I) = DGAUSSMELA(L22Q,A,B,EPS)
      ENDDO
*
      RETURN
      END
*
************************************************************************
*
*     Reference: Appendix A of hep-ph/9601302
*
************************************************************************
      function l22q(z)
*
      implicit none
*
      include "../commons/colfact.h"
**
*     Input Variables
*
      double precision z
**
*     Internal Variables
*
      double precision xi,zmax
      double precision sq1,sq2,l1,l2,l3,dil1,dil2,dil3,dil4
      double precision ddilogmela

      double precision qq2,mm2
      common/scales/qq2,mm2
**
*     Output Variables
*
      double precision l22q
*
      xi = qq2 / mm2
      zmax = 1d0/(1d0+4d0/xi)
      if(z.ge.zmax)then
         l22q = 0d0
      else
         sq1 = dsqrt( 1d0 - 4d0 * z / xi / ( 1d0 - z ) )
         sq2 = dsqrt( 1d0 - 4d0 * z / xi )
*
         l1 = dlog( ( 1d0 + sq1 ) / ( 1d0 - sq1 ) )
         l2 = dlog( ( 1d0 + sq2 ) / ( 1d0 - sq2 ) )
         l3 = dlog( ( sq2 + sq1 ) / ( sq2 - sq1 ) )
*
         dil1 = ddilogmela( ( 1d0 - z ) * (1d0 + sq1 ) / ( 1d0 + sq2 ) )
         dil2 = ddilogmela( ( 1d0 - sq2 ) / ( 1d0 + sq1 ) )
         dil3 = ddilogmela( ( 1d0 - sq1 ) / ( 1d0 + sq2 ) )
         dil4 = ddilogmela( ( 1d0 + sq1 ) / ( 1d0 + sq2 ) )
*
         l22q =CF*TR*((4d0/3d0*(1d0+z**2d0)/(1d0-z)
     1     -16d0/(1d0-z)*(z/xi)**2d0*(1d0-9d0*z+9d0*z**2d0))
     2     *(dlog((1d0-z)/z**2d0)*l1+l1*l2+2d0*(-dil1+dil2+dil3-dil4))
     3     +(-8d0/3d0+4d0/(1d0-z)+(z/(1d0-z)/xi)**2d0*(128d0-432d0*z
     4     +288*z**2d0-8d0/(1d0-z)))*l1+(88d0/9d0+136d0/9d0*z
     5     -152d0/9d0/(1d0-z)+(z/(1d0-z)/xi)*(464d0/9d0-512d0/3d0*z
     6     +2048d0/9d0*z**2d0)+(z/(1d0-z)/xi)**2d0*(-832d0/9d0
     7     +6208d0/9d0*z-11392d0/9d0*z**2d0+6016d0/9d0*z**3d0))*l3/sq2
     8     +(-272d0/27d0-1244d0/27d0*z+718d0/27d0/(1d0-z)+(z/(1d0-z)/xi)
     9     *(-3424d0/27d0+15608d0/27d0*z-4304d0/9d0*z**2d0
     1     +20d0/27d0/(1d0-z)))*sq1)
      endif
*
      return
      end
