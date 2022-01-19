************************************************************************
*
*     NCoefficientFunctions.f:
*
*     Routine that returns the evolved N-space distributions.
*
************************************************************************
      subroutine NCoefficientFunctions(N,Q,CF)
*
      implicit none
*
      include "../commons/ipt.h"
      include "../commons/massthrs.h"
      include "../commons/nffn.h"
      include "../commons/nfmax.h"
      include "../commons/ns.h"
      include "../commons/renfacscales.h"
      include "../commons/beta.h"
      include "../commons/alphas.h"
**
*     Input Variables
*
      double complex Q
      double complex N
**
*     Internal Variables
*
      integer isf,id
      integer nf
      integer j,k
      double complex Q2
      double complex a(2)
      double complex Cg(3,6),Cps(3,6),Cnsp(3,6),Cnsm(3,6)
      double complex c0nsp(3,6)
      double complex c1nsp(3,6),c1g(3,6)
      double complex c2nsp(3,6),c2nsm(3,6),c2ps(3,6),c2g(3,6)
      double complex tR,tF,tF2h
      double complex P0ns,P0sg(2,2)
      double complex P1ns(3),P1sg(2,2)
      double complex c2gh_ffns_nc,clgh_ffns_nc
      double complex c2qh_ffns_nc,clqh_ffns_nc
      double complex c2qgrh_ffns_nc,clqgrh_ffns_nc
**
*     Output Variables
*
      double complex CF(3,6,13)
*
      Q2   = Q * Q
*
*     Powers of alphas
*
      a(1) = asDIS
      a(2) = a(1) * a(1)
*
*     number of active flavours
*
      if(ns.eq."FFNS")then
         nf = nffn
         nf = nffn
      elseif(ns.eq."VFNS")then
         if(abs(Q2).ge.abs(q2th(6)))then
            nf = 6
         elseif(abs(Q2).ge.abs(q2th(5)))then
            nf = 5
         elseif(abs(Q2).ge.abs(q2th(4)))then
            nf = 4
         else
            nf = 3
         endif
         if(nf.gt.nfmax) nf = nfmax
      endif
*
*     Initialize coefficient functions
*
      do isf=1,3
         do k=1,6
            Cg(isf,k)   = (0d0,0d0)
            Cps(isf,k)  = (0d0,0d0)
            Cnsp(isf,k) = (0d0,0d0)
            Cnsm(isf,k) = (0d0,0d0)
         enddo
      enddo
*
*     LO
*
      do isf=1,3
         do k=1,6
            c0nsp(isf,k) = (0d0,0d0)
         enddo
      enddo
*     F3 (always ZM)
      do k=1,6
         c0nsp(3,k) = (1d0,0d0)
      enddo
*     F2
      if(ns.eq."VFNS")then
         do k=1,6
            c0nsp(1,k) = (1d0,0d0)
         enddo
      elseif(ns.eq."FFNS")then
         do k=1,3
            c0nsp(1,k) = (1d0,0d0)
         enddo
      endif
*
      do isf=1,3
         do k=1,6
            Cnsp(isf,k) = Cnsp(isf,k) + c0nsp(isf,k)
            Cnsm(isf,k) = Cnsm(isf,k) + c0nsp(isf,k)
         enddo
      enddo
*     
*     NNLO
*     
      if(ipt.ge.2)then
         if(ns.eq."VFNS")then
            call c22p(N,nf,c2nsp(1,1),c2nsm(1,1),c2ps(1,1),c2g(1,1))
            call cl2p(N,nf,c2nsp(2,1),c2nsm(2,1),c2ps(2,1),c2g(2,1))
            call c32p(N,nf,c2nsp(3,1),c2nsm(3,1),c2ps(3,1),c2g(3,1))
            do isf=1,3
               do k=2,6
                  c2nsp(isf,k) = c2nsp(isf,1)
                  c2nsm(isf,k) = c2nsm(isf,1)
                  c2ps(isf,k)  = c2ps(isf,1)
                  c2g(isf,k)   = c2g(isf,1)
               enddo
            enddo
         elseif(ns.eq."FFNS")then
            call c22p(N,nf,c2nsp(1,1),c2nsm(1,1),c2ps(1,1),c2g(1,1))
            call cl2p(N,nf,c2nsp(2,1),c2nsm(2,1),c2ps(2,1),c2g(2,1))
            call c32p(N,nf,c2nsp(3,1),c2nsm(3,1),c2ps(3,1),c2g(3,1))
            do isf=1,3
               do k=2,3
                  c2nsp(isf,k) = c2nsp(isf,1)
                  c2nsm(isf,k) = c2nsm(isf,1)
                  c2ps(isf,k)  = c2ps(isf,1)
                  c2g(isf,k)   = c2g(isf,1)
               enddo
            enddo
            do k=1,3
               c2nsp(1,k) = c2nsp(1,k)
     1              + c2qgrh_ffns_nc(2,N,Q2,q2th(4))
     2              + c2qgrh_ffns_nc(2,N,Q2,q2th(5))
     3              + c2qgrh_ffns_nc(2,N,Q2,q2th(6))
               c2nsp(2,k) = c2nsp(2,k)
     1              + clqgrh_ffns_nc(2,N,Q2,q2th(4))
     2              + clqgrh_ffns_nc(2,N,Q2,q2th(5))
     3              + clqgrh_ffns_nc(2,N,Q2,q2th(6))
            enddo
            do k=4,6
               c2nsp(1,k) = (0d0,0d0)
               c2nsp(2,k) = (0d0,0d0)
               c2nsp(3,k) = c2nsp(3,1)

               c2nsm(1,k) = (0d0,0d0)
               c2nsm(2,k) = (0d0,0d0)
               c2nsm(3,k) = c2nsm(3,1)

               c2ps(1,k)  = c2qh_ffns_nc(2,N,Q2,q2th(k))
               c2ps(2,k)  = clqh_ffns_nc(2,N,Q2,q2th(k))
               c2ps(3,k)  = c2ps(3,1)

               c2g(1,k)   = c2gh_ffns_nc(2,N,Q2,q2th(k),0)
               c2g(2,k)   = clgh_ffns_nc(2,N,Q2,q2th(k),0)
               c2g(3,k)   = c2g(3,1)
            enddo
         endif
*
*     Rescale coefficient functions appropriately in the presence of
*     scale variations.
*
         if(krenQ.ne.1d0.or.kfacQ.ne.1d0)then
            call andim_nlo(N,nf,P1ns,P1sg)
*
            do isf=1,3
               do k=1,6
*
*     The O(as^2) massive coefficient functions already contain the
*     factorization = renormalization scale variation (set tF and tF2h to zero).
*
                  if(ns.eq."FFNS".and.k.ge.4.and.isf.ne.3)then
                     tF   = 0d0
                     tF2h = 0d0
                     tR   = 2d0 * dlog( krenQ / kfacQ )
                  else
                     tF   = 2d0 * dlog(kfacQ)
                     tF2h = tF * tF / 2d0
                     tR   = 2d0 * dlog(krenQ)
                  endif
*     C+
                  c2nsp(isf,k) = c2nsp(isf,k)
     1                 + tR * beta0(nf)
     2                 * ( c1nsp(isf,k) + tF * P0ns * c0nsp(isf,k) )
     3                 - tF * ( c1nsp(isf,k)
     4                 + tF * P0ns * c0nsp(isf,k) ) * P0ns
     5                 + c0nsp(isf,k)
     6                  * ( tF2h * ( P0ns * P0ns - beta0(nf) * P0ns )
     7                 - tF * ( P1ns(1)
     8                 - ( tF - tR ) * beta0(nf) * P0ns ) )
*     C-
                  c2nsm(isf,k) = c2nsm(isf,k)
     1                 + tR * beta0(nf)
     2                 * ( c1nsp(isf,k) + tF * P0ns * c0nsp(isf,k) )
     3                 - tF * ( c1nsp(isf,k)
     4                 + tF * P0ns * c0nsp(isf,k) ) * P0ns
     5                 + c0nsp(isf,k)
     6                 * ( tF2h * ( P0ns * P0ns - beta0(nf) * P0ns )
     7                 - tF * ( P1ns(2)
     8                 - ( tF - tR ) * beta0(nf) * P0ns ) )
*     Cps
                  c2ps(isf,k) = c2ps(isf,k)
     1                 - tF * ( c1g(isf,k)
     2                 + tF * P0sg(1,2) * c0nsp(isf,k) / nf )
     3                 * P0sg(2,1)
     4                 + c0nsp(isf,k)
     5                 * ( tF2h * P0sg(1,2) * P0sg(2,1) / nf
     6                 - tF * ( P1sg(1,1) - P1ns(1) ) / nf )
*     Cg
                  c2g(isf,k) = c2g(isf,k)
     1                 + tR * beta0(nf) * ( c1g(isf,k)
     2                 + tF * P0sg(1,2) * c0nsp(isf,k) / nf )
     3                 - tF * ( ( c1nsp(isf,k)
     4                 + tF * P0ns * c0nsp(isf,k) ) * P0sg(1,2) / nf
     5                 + ( c1g(isf,k)
     6                 + tF * P0sg(1,2) * c0nsp(isf,k) / nf )
     7                 * P0sg(2,2) ) + c0nsp(isf,k) 
     8                 * ( tF2h / nf * ( P0sg(1,1) * P0sg(1,2)
     9                 + P0sg(1,2) * P0sg(2,2) - beta0(nf) * P0sg(1,2) )
     1                 - tF * ( P1sg(1,2)
     2                 - ( tF - tR ) * beta0(nf) * P0sg(1,2) ) / nf )
               enddo
            enddo
         endif
*
         do isf=1,3
            do k=1,6
               Cg(isf,k)   = Cg(isf,k)   + a(2) * c2g(isf,k)
               Cps(isf,k)  = Cps(isf,k)  + a(2) * c2ps(isf,k)
               Cnsp(isf,k) = Cnsp(isf,k) + a(2) * c2nsp(isf,k)
               Cnsm(isf,k) = Cnsm(isf,k) + a(2) * c2nsm(isf,k)
            enddo
         enddo
      endif
*
*     NLO
*
      if(ipt.ge.1)then
         if(ns.eq."VFNS")then
            call c21p(N,c1nsp(1,1),c1g(1,1))
            call cl1p(N,c1nsp(2,1),c1g(2,1))
            call c31p(N,c1nsp(3,1),c1g(3,1))
            do isf=1,3
               do k=2,6
                  c1nsp(isf,k) = c1nsp(isf,1)
                  c1g(isf,k)   = c1g(isf,1)
               enddo
            enddo
         elseif(ns.eq."FFNS")then
            call c21p(N,c1nsp(1,1),c1g(1,1))
            call cl1p(N,c1nsp(2,1),c1g(2,1))
            call c31p(N,c1nsp(3,1),c1g(3,1))
            do isf=1,3
               do k=2,3
                  c1nsp(isf,k) = c1nsp(isf,1)
                  c1g(isf,k)   = c1g(isf,1)
               enddo
            enddo
            do k=4,6
               c1nsp(1,k) = (0d0,0d0)
               c1nsp(2,k) = (0d0,0d0)
               c1nsp(3,k) = c1nsp(3,1)
               c1g(1,k)   = c2gh_ffns_nc(1,N,Q2,q2th(k),0)
               c1g(2,k)   = clgh_ffns_nc(1,N,Q2,q2th(k),0)
               c1g(3,k)   = c1g(3,1)
            enddo
         endif
*
*     Rescale coefficient functions appropriately in the presence of
*     scale variations.
*
         if(kfacQ.ne.1d0)then
            call andim_lo(N,nf,P0ns,P0sg)
*
            tF = 2d0 * dlog(kfacQ)
            do isf=1,3
               do k=1,6
                  c1nsp(isf,k) = c1nsp(isf,k) - tF * P0ns * c0nsp(isf,k)
                  c1g(isf,k)   = c1g(isf,k)
     1                         - tF * P0sg(1,2) * c0nsp(isf,k) / nf
               enddo
            enddo
         endif
*
         do isf=1,3
            do k=1,6
               Cg(isf,k)   = Cg(isf,k)   + a(1) * c1g(isf,k)
               Cnsp(isf,k) = Cnsp(isf,k) + a(1) * c1nsp(isf,k)
               Cnsm(isf,k) = Cnsm(isf,k) + a(1) * c1nsp(isf,k)
            enddo
         enddo
      endif
*
*     F2 and FL
*
      do isf=1,2
         do k=1,6
            do id=1,13
               CF(isf,k,id) = (0d0,0d0)
            enddo
*     Singlet
            CF(isf,k,1) = bq(k) * ( Cps(isf,k) + Cnsp(isf,k) / nf )
*     Gluon
            CF(isf,k,2) = bq(k) * Cg(isf,k)
*     T3, T8, T15, T24, T35
            if(k.ge.2)then
               CF(isf,k,7+k) = - bq(k) * Cnsp(isf,k) / k
            endif
            do j=k+1,nf
               CF(isf,k,7+j) = bq(k) * Cnsp(isf,k) / j / ( j - 1 )
            enddo
         enddo
      enddo
*
*     F3
*
      do k=1,6
         do id=1,13
            CF(3,k,id) = (0d0,0d0)
         enddo
         if(k.le.nf)then
*     Valence
            CF(3,k,3) = dq(k) * Cnsm(3,k) / nf
*     V3, V8, V15, V24, V35
            if(k.ge.2)then
               CF(3,k,2+k) = - dq(k) * Cnsm(3,k) / k
            endif
            do j=k+1,nf
               CF(3,k,2+j) = dq(k) * Cnsm(3,k) / j / ( j - 1 )
            enddo
         endif
      enddo
*
      return
      end
