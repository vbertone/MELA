************************************************************************
*
*     xDistributions.f:
*
*     Routine that returns the evolved x-space distributions.
*
************************************************************************
      subroutine xDistributions(x,Q,xfph)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/alpha.h"
      include "../commons/massthrs.h"
      include "../commons/ns.h"
      include "../commons/nffn.h"
      include "../commons/nfmax.h"
      include "../commons/tecparam.h"
**
*     Input Variables
*
      double precision x
      double precision Q
**
*     Internal Variables
*
      integer nfi,nff
      integer ip
      double precision aQED
      double precision xfev(19)
      double precision integrand
      external integrand
      integer j
      double precision r
      double precision theta,sigma,t
      double complex s,tmp(19)
      double complex xfN(19)
      double precision sq2thnfi,sq2thnffp1      
**
*     Output Variables
*
      double precision xfph(-9:9)
*
c$$$      double precision dgauss,dtrap      
c$$$      double precision xval
c$$$      integer nfival,nffval,ipval,mval
c$$$      common/integvar/xval,nfival,nffval,ipval,mval
*
*     Determine number of active flavours at the final scale
*
      if(ns.eq."FFNS")then
         nfi = nffn
         nff = nffn
      elseif(ns.eq."VFNS")then
         nfi = 1
         do nff = 9, 1, -1
            if (Q**2.ge.q2th(nff)) exit
         enddo
         if(nff.gt.nfmax) nff = nfmax
      endif
*
      if(aemfix)then
         sq2thnfi = q2th(nfi)
         sq2thnffp1 = q2th(nff+1)      
         q2th(nfi)  = q2th(1)
         q2th(nff+1) = Q**2d0
      else
         ath(nfi)   = ath(1)
         ath(nff+1) = aQED(Q**2d0)
      endif
*     
*     trapezioidal integration with Talbot path
      do ip=1,19
         tmp(ip) = (0d0,0d0)
      enddo      
      t = - dlog(x)
      r = 2d0 * rinvmel / 5d0 / t
      do j=1,minvmel-1
         theta = - pi + dble(j) * ( 2d0 * pi / dble(minvmel) )
         sigma = theta + (theta/tan(theta)-1d0)/tan(theta)
         s     = r * theta * dcmplx(1d0/tan(theta),1d0) + 1d0
         call NDistributions(s,nfi,nff,xfN)
         do ip=1,19
            tmp(ip) = tmp(ip)
     1              + exp( t * s ) * dcmplx(1d0,sigma) * xfN(ip)
         enddo
      enddo
      do ip=1,19
         xfev(ip) = x * r * dble(tmp(ip)) / minvmel
      enddo
*
*     Try a different integrator. Slower because xfN is evaluated
*     (m-2)*19 times, whereas above only (m-2) times, since the x points
*     where the integrand is evaluated are the same for all ip.     
c$$$*     Fill the common
c$$$      xval = x
c$$$      nfival = nfi
c$$$      nffval = nff
c$$$      mval = 83
c$$$*     Integrate
c$$$      do ip=1,19
c$$$         ipval = ip
c$$$c        xfev(ip) = x * dtrap(integrand,-pi,pi,mval)/(2d0*pi)
c$$$         xfev(ip) = x * dgauss(integrand,0d0,pi,1d-7)/pi
c$$$      enddo
c$$$*            
*
*     Rotate back to the physical basis
*
      call evln2lha(xfev,xfph)
*
*     Restore thresholds or ath(nff+1) to proper values
      if(aemfix)then
         q2th(nfi) = sq2thnfi
         q2th(nff+1) = sq2thnffp1
      else
         call initCouplings
      endif         
*     
      return
      end
************************************************************************
      double precision function integrand(theta)
      implicit none
      include "../commons/consts.h"      
*
      double precision theta
      double precision t,r
      double precision sigma
      double complex s
      double complex xfN(19)      
*
      double precision xval
      integer nfival,nffval,ipval,mval
      common/integvar/xval,nfival,nffval,ipval,mval
*      
      t = - dlog(xval)
      r = 2d0 * 43 / 5d0 / t
      sigma = theta + (theta/tan(theta)-1d0)/tan(theta)
      s     = r * theta * dcmplx(1d0/tan(theta),1d0) + 1d0
      call NDistributions(s,nfival,nffval,xfN)
      integrand = r * dble(exp(t*s) * dcmplx(1d0,sigma) * xfN(ipval))
c     write(*,*) "ipval",ipval,"theta/pi",theta/pi,"integrand",integrand
*      
      return 
      end
************************************************************************
      double precision function dtrap(f,a,b,n)
      implicit none
*     Input
      double precision f,a,b
      integer n      
*     Internal
      integer j
      double precision deltax
*
      dtrap = 0d0
      deltax = (b-a)/n
      do j=1,n-1
         dtrap = dtrap + f(a+j*deltax)
*     N.B. We are ignoring the additional piece (f(a)+f(b))/2
*     which for our function is null, when a=-pi and b=+pi,
*     but problematic to compute         
c        dtrap = dtrap + (f(a)+f(b))/2d0
      enddo
      dtrap = deltax*dtrap
*      
      return 
      end
************************************************************************
c$$$      DOUBLE PRECISION FUNCTION DGAUSS(F,A,B,EPS)
c$$$      DOUBLE PRECISION W(12),X(12),A,B,EPS,DELTA,CONST,AA,BB,Y,C1,C2,S8,
c$$$     1                 S16,U,F
c$$$      DATA CONST /1.0D-25/
c$$$      DATA W
c$$$     1       / 0.10122 85362 90376 25915 25313 543D0,
c$$$     2         0.22238 10344 53374 47054 43559 944D0,
c$$$     3         0.31370 66458 77887 28733 79622 020D0,
c$$$     4         0.36268 37833 78361 98296 51504 493D0,
c$$$     5         0.02715 24594 11754 09485 17805 725D0,
c$$$     6         0.06225 35239 38647 89286 28438 370D0,
c$$$     7         0.09515 85116 82492 78480 99251 076D0,
c$$$     8         0.12462 89712 55533 87205 24762 822D0,
c$$$     9         0.14959 59888 16576 73208 15017 305D0,
c$$$     A         0.16915 65193 95002 53818 93120 790D0,
c$$$     B         0.18260 34150 44923 58886 67636 680D0,
c$$$     C         0.18945 06104 55068 49628 53967 232D0 /
c$$$      DATA X
c$$$     1       / 0.96028 98564 97536 23168 35608 686D0,
c$$$     2         0.79666 64774 13626 73959 15539 365D0,
c$$$     3         0.52553 24099 16328 98581 77390 492D0,
c$$$     4         0.18343 46424 95649 80493 94761 424D0,
c$$$     5         0.98940 09349 91649 93259 61541 735D0,
c$$$     6         0.94457 50230 73232 57607 79884 155D0,
c$$$     7         0.86563 12023 87831 74388 04678 977D0,
c$$$     8         0.75540 44083 55003 03389 51011 948D0,
c$$$     9         0.61787 62444 02643 74844 66717 640D0,
c$$$     A         0.45801 67776 57227 38634 24194 430D0,
c$$$     B         0.28160 35507 79258 91323 04605 015D0,
c$$$     C         0.09501 25098 37637 44018 53193 354D0 /
c$$$      DELTA=CONST*DABS(A-B)
c$$$      DGAUSS=0.
c$$$      AA=A
c$$$    5 Y=B-AA
c$$$      IF(DABS(Y) .LE. DELTA) RETURN
c$$$    2 BB=AA+Y
c$$$      C1=0.5D0*(AA+BB)
c$$$      C2=C1-AA
c$$$      S8=0.
c$$$      S16=0.
c$$$      DO 1 I = 1,4
c$$$      U=X(I)*C2
c$$$    1 S8=S8+W(I)*(F(C1+U)+F(C1-U))
c$$$      DO 3 I = 5,12
c$$$      U=X(I)*C2
c$$$    3 S16=S16+W(I)*(F(C1+U)+F(C1-U))
c$$$      S8=S8*C2
c$$$      S16=S16*C2
c$$$      IF(DABS(S16-S8) .GT. EPS*(1.0D0+DABS(S16))) GO TO 4
c$$$      DGAUSS=DGAUSS+S16
c$$$      AA=BB
c$$$      GO TO 5
c$$$    4 Y=0.5D0*Y
c$$$      IF(DABS(Y) .GT. DELTA) GO TO 2
c$$$      PRINT 7
c$$$      DGAUSS=0.
c$$$      RETURN
c$$$    7 FORMAT(1X,'DGAUSS ... TOO HIGH ACCURACY REQUIRED')
c$$$      END
