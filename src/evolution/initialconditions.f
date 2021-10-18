************************************************************************
*
*     initialconditions.f:
*
*     These routines return PDFs/FFs at the intitial scale in N space.
*
************************************************************************
      subroutine electronPDFsn(N, npdf)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/massthrs.h"
      include "../commons/evscale.h"
      include "../commons/alpha.h"
      include "../commons/ipt.h"
**
*     Input Variables
*
      double complex N
**
*     Internal Variables
*
      integer ipdf
      double complex Nm, N1
      double complex s1, s2
      double complex psi, dpsi
      double complex deN, dgmN
      double complex logfact
      double complex sg, gm
**
*     Output Variables
*
      double complex npdf(-3:3)
*
*     Useful definitions
*
      Nm = N - 1d0
      N1 = N + 1d0
      s1 = emc + psi(N1)
      s2 = zeta2 - dpsi(N1, 1)
*
      deN  = 2d0 * ( - ( 1d0 / N - 1d0 / N1 - 2d0 * s1 + 3d0 / 2d0 )
     1     - 2d0 * ( s1**2 + s2 - s1 / N + s1 / N1 + 1d0 / N1**2
     2     - 7d0 / 4d0 ) )
      dgmN = 2d0 * ( - ( 2d0 / Nm - 2d0 / N + 1d0 / N1 )
     1     - 2d0 * ( - 2d0 / Nm**2 + 2d0 / N**2 - 1d0 / N1**2 ) )
*
      logfact = zlog(Q20 / q2th(1));
*
      sg = 1d0
      gm = 0d0
      if (ipt.gt.0) then
         sg = sg + a0 * ( logfact * 2d0
     1        * ( 1d0 / N - 1d0 / N1 - 2d0 * s1 + 3. / 2. ) + deN )
         gm = gm + a0 * ( logfact * 2d0
     1        * ( 2d0 / Nm - 2d0 / N + 1d0 / N1 ) + dgmN )
      endif
*
*     Initialize PDFs to zero
*      
      do ipdf = -3,3
         npdf(ipdf) = (0d0,0d0)
      enddo
*
*     Set electron and photon PDFs
*
      npdf(0) = gm
      npdf(1) = sg
*
      return
      end
*
************************************************************************
      subroutine electronFFsn(N, nff)
*
      implicit none
**
*     Input Variables
*
      double complex N
**
*     Internal Variables
*
      integer ipdf
**
*     Output Variables
*
      double complex nff(-3:3)
*
*     Initialize FFs to zero
*
      do ipdf = -3,3
         nff(ipdf) = (0d0,0d0)
      enddo
*
      return
      end
