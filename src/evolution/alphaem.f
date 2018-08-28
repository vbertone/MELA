****************************************************************
*
*     alphaem.f:
*
*     Routine that returns alpha_QED / 4 / pi .
*
*     Note that the only argument of this routine is the 
*     factorization scale mu2F which is convereted into the 
*     renormalization scale inside the routine itself.
*
****************************************************************
      function aQED(mu2f)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/alphaem.h"
      include "../commons/renfacscales.h"
**
*     Input Variables
*
      double precision mu2F
**
*     Internal Variables
*
      double precision mur2,mur20
      double precision aem0
      double precision alphaqedevMELA
**
*     Output Variables
*
      double precision aQED
*
      aem0 = aemref / 4d0 / pi
*
      mur20 = q2emref
      mur2  = krf * mu2F
*
      aQED = alphaqedevMELA(0,1,mur2,mur20,aem0)
*
      return
      end
*
****************************************************************
      FUNCTION ALPHAQEDEVMELA(NF,NL,Q2F,Q2I,ALPHAREF)
*
      IMPLICIT NONE
**
*     Input Variables
*
      INTEGER NF,NL
      DOUBLE PRECISION Q2F,Q2I
      DOUBLE PRECISION ALPHAREF
**
*     Internal Variables
*
      DOUBLE PRECISION BETA0QEDMELA,BETA0
      DOUBLE PRECISION L
**
*     Output Variables
*
      DOUBLE PRECISION ALPHAQEDEVMELA
*
      BETA0 = BETA0QEDMELA(NF,NL)
      L = DLOG( Q2F / Q2I )
*
      ALPHAQEDEVMELA = ALPHAREF / ( 1D0 + ALPHAREF * BETA0 * L )
*
      RETURN
      END
*
****************************************************************
      function beta0qedMELA(nf,nl)
*
      implicit none
**
*     Input Variables
*
      integer nf,nl
**
*     Internal Variables
*
      integer nc
      double precision sumch2(0:6)
**
*     Output Variables
*
      double precision beta0qedMELA
*
*     Number of colours
*
      nc = 3
*
*     Sum of the first 3, 4, 5 and 6 squared electric charges
*
      sumch2(0) = 0d0
      sumch2(1) = 1d0 / 9d0
      sumch2(2) = 5d0 / 9d0
      sumch2(3) = 2d0 / 3d0
      sumch2(4) = 10d0 / 9d0
      sumch2(5) = 11d0 / 9d0
      sumch2(6) = 5d0 / 3d0
*
      beta0qedMELA = - 4d0 / 3d0 * ( nc * sumch2(nf) + nl )
*
      return
      end
*
****************************************************************************
*
*     QED beta function.
*
****************************************************************************
      function fbetaQEDMELA(a,nf,nl,ipt)
*
      implicit none
**
*     Input Variables
*
      double precision a
      integer nf,nl,ipt
**
*     Internal Variables
*
      double precision beta0qedMELA
**
*     Output Variables
*
      double precision fbetaQEDMELA
*
      fbetaQEDMELA = - a**2 * beta0qedMELA(nf,nl)
*
      return
      end
