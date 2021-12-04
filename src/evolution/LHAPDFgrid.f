************************************************************************
*
*     LHAPDFgrid.f:
*
************************************************************************
      subroutine LHAPDFgrid(Qin,fname)
*
      implicit none
*
      include "../commons/LHAgrid.h"
      include "../commons/ipt.h"
      include "../commons/alphas.h"
      include "../commons/massthrs.h"
      include "../commons/renfacscales.h"
      include "../commons/ns.h"
      include "../commons/nfmax.h"
      include "../commons/nffn.h"
**
*     Input Variables
*
      double precision Qin
      character*50 fname
**
*     Intenal Variables
*
      integer ln
      integer ix,iq2,iq2c,ipdf,ilep,krep,iq2in,iq2fi
      integer nfin,nffi
      integer isg,nQ(3:7)
      double precision xbLHA(nxmax),q2LHA(nq2max),as(nxmax)
      double precision lnQmin,lnQmax
      double precision eps,eps2
      double complex aQCD
      double complex Q(2)
      double complex xf(-6:6)
      character*30  str
      character*100 pdfsetbkp
      character*3   ids(-6:6)
      parameter(eps=1d-4)
      parameter(eps2=1d-8)
*
*     Initialize APFEL
*
      call initializeEvolution
*
*     Report grid parameters
*
      write(6,*) "Report of the LHAPDF grid parameters:"
      write(6,"(a,i3,a,f11.9,a,f6.4,a)") " - x-space grid: ",nxLHA,
     1     " points in [",xminLHA," :",xmaxLHA,"]"
      write(6,"(a,f6.4)") "    transition from log to lin in x = ",xmLHA
      write(6,"(a,i3,a,f5.2,a,f13.1,a)") " - Q2-space grid: ",nq2LHA,
     1     " points in [",q2minLHA,":",q2maxLHA,"] GeV^2"
*
*     Compute x-space grid
*
      do ix=1,nxLHA
         if(ix.le.nxmLHA)then
            xbLHA(ix) = xminLHA * ( xmLHA / xminLHA )
     1                **( 2d0 * dble( ix-1 ) / dble( nxLHA - 1 ) )
         else
            xbLHA(ix) = xmLHA + ( xmaxLHA - xmLHA )
     1                * ( dble( ix - nxmLHA - 1 )
     2                / dble( nxLHA - nxmLHA - 1 ) )
         endif
      enddo
*
*     Compute Q2 grid.
*     Use a distribution of the Q2 nodes uniform in ln(ln(Q2/Lambda2))
*     and that has nodes on the heavy quark thresholds.
*
      lnQmin = dlog( q2minLHA / Lambda2 )
      lnQmax = dlog( q2maxLHA / Lambda2 )
*
*     Initialize number of points per subgrid
*
      do isg=3,7
         nQ(isg) = 0
      enddo
*
      if(Ns.eq."VFNS")then
         if(q2minLHA.gt.dble(q2th(6)))then
            nfin = 6
         elseif(q2minLHA.gt.dble(q2th(5)))then
            nfin = 5
         elseif(q2minLHA.gt.dble(q2th(4)))then
            nfin = 4
         else
            nfin = 3
         endif
         if(nfin.gt.nfmax) nfin = nfmax
*
         if(q2maxLHA.gt.dble(q2th(6)))then
            nffi = 6
         elseif(q2maxLHA.gt.dble(q2th(5)))then
            nffi = 5
         elseif(q2maxLHA.gt.dble(q2th(4)))then
            nffi = 4
         else
            nffi = 3
         endif
         if(nffi.gt.nfmax) nffi = nfmax
*
         isg = nfin
         do iq2=1,nq2LHA
            q2LHA(iq2) = Lambda2 * dexp( lnQmin
     1                 * dexp( dble( iq2 - 1 ) / dble( nq2LHA - 1 )
     2                 * dlog( lnQmax / lnQmin ) ) )
            if(q2LHA(iq2).lt.dble(q2th(isg+1))+eps)then
               nQ(isg) = nQ(isg) + 1
            else
               isg = isg + 1
               if(isg.gt.nfmax) isg = nfmax
               nQ(isg) = nQ(isg) + 1
            endif
         enddo
*
*     Make sure that all subgrids have at least two points.
*
         do isg=nfin,nffi-1
            if(nQ(isg).lt.2)then
               nQ(isg) = nQ(isg) + 1
               nQ(isg+1) = nQ(isg+1) - 1
            endif
         enddo
*
*     Redefine the grid with the subgrids
*
         lnQmin = dlog( ( q2minLHA ) / Lambda2 )
         if(nfin.eq.nffi)then
            lnQmax = dlog( q2maxLHA / Lambda2 )
         else
            lnQmax = dlog( ( dble(q2th(nfin+1)) ) / Lambda2 )
         endif
*
         iq2c = 0
         do isg=nfin,nffi
            do iq2=1,nQ(isg)
               iq2c = iq2c + 1
               q2LHA(iq2c) = Lambda2 * dexp( lnQmin
     1              * dexp( dble( iq2 - 1 )
     2              / dble( nQ(isg) - 1 )
     3              * dlog( lnQmax / lnQmin ) ) )
            enddo
            lnQmin = dlog( ( q2LHA(iq2c) ) / Lambda2 )
            if(isg.eq.nffi-1)then
               lnQmax = dlog( q2maxLHA / Lambda2 )
            else
               lnQmax = dlog( ( dble(q2th(isg+2)) ) / Lambda2 )
            endif
         enddo
*
         if(iq2c.ne.nq2LHA)then
            write(6,*) "In LHAPDFgrid.f:"
            write(6,*) "Mismatch in the Number of Q2 nodes"
            write(6,*) "- Expected = ",nq2LHA
            write(6,*) "- Found = ",iq2c
            call exit(-10)
         endif
      else
         nfin = nffn
         nffi = nffn
*
         nQ(nffn) = nq2LHA
*
         do iq2=1,nq2LHA
            q2LHA(iq2) = Lambda2 * dexp( lnQmin
     1                 * dexp( dble( iq2 - 1 ) / dble( nq2LHA - 1 )
     2                 * dlog( lnQmax / lnQmin ) ) )
         enddo
      endif
*
*     Compute alphas on the grid being careful with the grids
*
      iq2in = 1
      iq2fi = nQ(nfin)
      do isg=nfin,nffi
         do iq2=iq2in,iq2fi
            if(iq2.eq.iq2in.and.isg.ne.nfin)then
               as(iq2) = 12.566370614359173d0
     1              * dble(aQCD(q2LHA(iq2)+eps2))
            elseif(iq2.eq.iq2fi.and.isg.ne.nffi)then
               as(iq2) = 12.566370614359173d0
     1              * dble(aQCD(q2LHA(iq2)-eps2))
            else
               as(iq2) = 12.566370614359173d0 * dble(aQCD(q2LHA(iq2)))
            endif
         enddo
         iq2in = iq2in + nQ(isg)
         iq2fi = iq2fi + nQ(isg+1)
      enddo
*
*     Define quark IDs
*
      ids(-6) = " -6"
      ids(-5) = " -5"
      ids(-4) = " -4"
      ids(-3) = " -3"
      ids(-2) = " -2"
      ids(-1) = " -1"
      ids(0)  = " 21"
      ids(1)  = "  1"
      ids(2)  = "  2"
      ids(3)  = "  3"
      ids(4)  = "  4"
      ids(5)  = "  5"
      ids(6)  = "  6"
*
*     LHAPDF6 output
*
      ln = index(fname,char(0)) - 1
      if(ln.eq.-1) ln = index(fname,char(32)) - 1
*     creating main folder
      call mkdir(fname(1:ln))
*     creating info file
      open(unit=13,status="unknown",file=fname(1:ln)//"/"
     1                    //fname(1:ln)//".info")
*
*     Write header
*
      write(13,*) 'SetDesc: "set generated with MELA - ',
     1     fname(1:ln),'"'
      write(13,*) "Authors: V. Bertone"
      write(13,*) "Reference: ArXiv:xxxx.xxxxx"
      write(13,*) "Format: lhagrid1"
      write(13,*) "DataVersion: 1"
      write(13,*) "NumMembers: 1"
      write(13,*) "Particle: 0000"
      write(13,*) "Flavors: [",
     1     (ids(ipdf),",",ipdf=-nfmax,nfmax-1),
     2     ids(nfmax),"]"
      write(13,*) "OrderQCD:",ipt
      write(13,*) "FlavorScheme: variable"
      write(13,*) "NumFlavors: ",nfmax
      write(13,*) "ErrorType: replicas"
      write(13,*) "XMin:",xminLHA
      write(13,*) "XMax:",xmaxLHA
      write(13,*) "QMin:",dsqrt(q2minLHA)
      write(13,*) "QMax:",dsqrt(q2maxLHA)
      write(13,*) "MZ:",dsqrt(dble(q2ref))
      write(13,*) "MUp: 0"
      write(13,*) "MDown: 0"
      write(13,*) "MStrange: 0"
      write(13,*) "MCharm:",dsqrt(dble(q2th(4)))
      write(13,*) "MBottom:",dsqrt(dble(q2th(5)))
      write(13,*) "MTop:",dsqrt(dble(q2th(6)))
      write(13,*) "AlphaS_MZ:",asref
      write(13,*) "AlphaS_OrderQCD:",ipt
      write(13,*) "AlphaS_Type: ipol"
      write(13,*) "AlphaS_Qs: [",(dsqrt(q2LHA(iq2)),",",
     1     iq2=1,nq2LHA-1),dsqrt(q2LHA(nq2LHA)),"]"
      write(13,*) "AlphaS_Vals: [",(as(iq2),",",
     1     iq2=1,nq2LHA-1),as(nq2LHA),"]"
      close(13)
*
*     Now Loop over all replicas and print to file
*
      open(unit=13,status="unknown",file=fname(1:ln)//"/"
     1     //fname(1:ln)//"_0000.dat")

      write(13,"(a)") "PdfType: central"
      write(13,"(a)") "Format: lhagrid1"
      write(13,"(a)") "---"
*
      iq2in = 1
      iq2fi = nQ(nfin)
      Q(1)  = dcmplx(Qin, 0d0)
      do isg=nfin,nffi
         write(13,*) (xbLHA(ix),ix=1,nxLHA)
         write(13,*) (dsqrt(q2LHA(iq2)),iq2=iq2in,iq2fi)
         write(13,*) (ids(ipdf),ipdf=-nfmax,nfmax)
         do ix=1,nxLHA
            do iq2=iq2in,iq2fi
               Q(2) = dcmplx(dsqrt(q2LHA(iq2)), 0d0)
               if(iq2.eq.iq2in.and.isg.ne.nfin)then
                  Q(2) = dcmplx(dsqrt(q2LHA(iq2)) + eps2, 0d0)
               elseif(iq2.eq.iq2fi.and.isg.ne.nffi)then
                  Q(2) = dcmplx(dsqrt(q2LHA(iq2)) - eps2, 0d0)
               endif
               call xDistributions(xbLHA(ix),2,Q,xf)
               write(13,40) (dble(xf(ipdf)), ipdf=-nfmax,nfmax)
            enddo
         enddo
         write(13,"(a)") "---"
         iq2in = iq2in + nQ(isg)
         iq2fi = iq2fi + nQ(isg+1)
      enddo
      close(13)
*
      write(6,*) achar(27)//"[1;32m"
      write(6,*) "File ",fname(1:ln)," grid produced"
      write(6,*) achar(27)//"[0m"
*
 40   format(13(es14.7,1x))
*
      return
      end
