************************************************************************
*
*     toyLHPDFs.f:
*
*     These routines return some PDFs/FFs at the intitial scale in N
*     space.
*
************************************************************************
      subroutine toyLHPDFsn(N,npdf)
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
      double complex N_uv,auv,buv,N_dv,adv,bdv,N_g,ag
      double complex bg,N_db,adb,bdb,fs
      double complex nuv,ndv,ng,ndbar,nubar,ns,nsbar
      double complex betac
**
*     Output Variables
*
      double complex npdf(-6:6)
*
*     Parameters of the User defined PDFs
*
      N_uv = dcmplx(5.107200d0,0d0)
      auv  = dcmplx(0.8d0,0d0) - dcmplx(1d0,0d0)
      buv  = dcmplx(3d0,0d0)
      N_dv = dcmplx(3.064320d0,0d0)
      adv  = dcmplx(0.8d0,0d0) - dcmplx(1d0,0d0)
      bdv  = dcmplx(4d0,0d0)
      N_g  = dcmplx(1.7d0,0d0)
      ag   = dcmplx(-0.1d0,0d0) - dcmplx(1d0,0d0)
      bg   = dcmplx(5d0,0d0)
      N_db = dcmplx(0.1939875d0,0d0)
      adb  = dcmplx(-0.1d0,0d0) - dcmplx(1d0,0d0)
      bdb  = dcmplx(6d0,0d0)
      fs   = dcmplx(0.2d0,0d0)
*
*     User defined PDFs
*
      nuv   = N_uv * betac(N+auv,buv+dcmplx(1d0,0d0))
      ndv   = N_dv * betac(N+adv,bdv+dcmplx(1d0,0d0))
      ng    = N_g  * betac(N+ag ,bg +dcmplx(1d0,0d0))
      ndbar = N_db * betac(N+adb,bdb+dcmplx(1d0,0d0))
      nubar = N_db * betac(N+adb,bdb+dcmplx(2d0,0d0))
      ns    = fs * ( ndbar + nubar )
      nsbar = ns
*
*     Initialize PDFs to zero
*
      do ipdf=-6,6
         npdf(ipdf) = (0d0,0d0)
      enddo
*
      npdf(3)  = ns
      npdf(2)  = nuv + nubar
      npdf(1)  = ndv + ndbar
      npdf(0)  = ng
      npdf(-1) = ndbar
      npdf(-2) = nubar
      npdf(-3) = nsbar
*
      return
      end
*
************************************************************************
      subroutine toyLHPDFsPoln(N,npdf)
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
      double complex N_uv,auv,buv,N_dv,adv,bdv,N_g,ag
      double complex bg,N_db,adb,bdb,fs
      double complex nuv,ndv,ng,ndbar,nubar,ns,nsbar
      double complex betac
**
*     Output Variables
*
      double complex npdf(-6:6)
*
*     Parameters of the User defined PDFs
*
      N_uv = dcmplx(1.3d0,0d0)
      auv  = dcmplx(0.7d0,0d0) - dcmplx(1d0,0d0)
      buv  = dcmplx(3d0,0d0)
      N_dv = dcmplx(-0.5d0,0d0)
      adv  = dcmplx(0.7d0,0d0) - dcmplx(1d0,0d0)
      bdv  = dcmplx(4d0,0d0)
      N_g  = dcmplx(1.5d0,0d0)
      ag   = dcmplx(0.5d0,0d0) - dcmplx(1d0,0d0)
      bg   = dcmplx(5d0,0d0)
      N_db = dcmplx(-0.05d0,0d0)
      adb  = dcmplx(0.3d0,0d0) - dcmplx(1d0,0d0)
      bdb  = dcmplx(7d0,0d0)
      fs   = dcmplx(0.5d0,0d0)
*
*     User defined PDFs
*
      nuv   = N_uv * ( betac(N+auv,buv+dcmplx(1d0,0d0))
     1      + 3d0 * betac(N+auv+1d0,buv+dcmplx(1d0,0d0)) )
      ndv   = N_dv * ( betac(N+adv,bdv+dcmplx(1d0,0d0))
     1      + 4d0 * betac(N+adv+1d0,bdv+dcmplx(1d0,0d0)) )
      ng    = N_g  * betac(N+ag ,bg +dcmplx(1d0,0d0))
      ndbar = N_db * betac(N+adb,bdb+dcmplx(1d0,0d0))
      nubar = ndbar
      ns    = fs * ndbar
      nsbar = ns
*
*     Initialize PDFs to zero
*
      do ipdf=-6,6
         npdf(ipdf) = (0d0,0d0)
      enddo
*
      npdf(3)  = ns
      npdf(2)  = nuv + nubar
      npdf(1)  = ndv + ndbar
      npdf(0)  = ng
      npdf(-1) = ndbar
      npdf(-2) = nubar
      npdf(-3) = nsbar
*
      return
      end
*
************************************************************************
*
*     Kretzer's parametrization at Q2 = 0.4 GeV^2 of the light partons
*     for pi+ taken at NLO from hep-ph/0003177.
*
************************************************************************
      subroutine KretzerFFsn(N,nff)
*
      implicit none
**
*     Input Variables
*
      double complex N
**
*     Internal Variables
*
      integer iff
      double complex ag,bg,N_g
      double complex as,bs,N_s
      double complex al,bl,N_l
      double complex betac
**
*     Output Variables
*
      double complex nff(-6:6)
*
*     Parameters of the User defined PDFs
*
      al  = dcmplx(-0.829d0,0d0)
      bl  = dcmplx(0.949d0,0d0)
      N_l = dcmplx(0.264d0,0d0) 
     1    / betac(al+dcmplx(2d0,0d0),bl+dcmplx(1d0,0d0))
      as  = al
      bs  = bl + dcmplx(1d0,0d0)
      N_s = dcmplx(0.165d0,0d0) 
     1    / betac(as+dcmplx(2d0,0d0),bs+dcmplx(1d0,0d0))
      ag  = dcmplx(4.374d0,0d0)
      bg  = dcmplx(9.778d0,0d0)
      N_g = dcmplx(0.215d0,0d0) 
     1    / betac(ag+dcmplx(2d0,0d0),bg+dcmplx(1d0,0d0))
*
*     Initialize PDFs to zero
*
      do iff=-6,6
         nff(iff) = (0d0,0d0)
      enddo
*
      nff(3)  = N_s * betac(N+as,bs+dcmplx(1d0,0d0))
      nff(2)  = N_l * betac(N+al,bl+dcmplx(1d0,0d0))
      nff(1)  = nff(3)
      nff(0)  = N_g * betac(N+ag,bg+dcmplx(1d0,0d0))
      nff(-1) = nff(2)
      nff(-2) = nff(1)
      nff(-3) = nff(3)
*
      return
      end
*
************************************************************************
*
*     HKNS parametrization at Q2 = 1 GeV^2 of the light partons
*     for pi+ taken at NLO from hep-ph/0702250.
*
************************************************************************
      subroutine HKNSFFsn(N,nff)
*
      implicit none
**
*     Input Variables
*
      double complex N
**
*     Internal Variables
*
      integer iff
      double complex ag,bg,N_g
      double complex as,bs,N_s
      double complex al,bl,N_l
      double complex betac
**
*     Output Variables
*
      double complex nff(-6:6)
*
*     Parameters of the User defined PDFs
*
      al  = dcmplx(-0.963d0,0d0)
      bl  = dcmplx(1.370d0,0d0)
      N_l = dcmplx(0.401d0,0d0) 
     1    / betac(al+dcmplx(2d0,0d0),bl+dcmplx(1d0,0d0))
      as  = dcmplx(0.718d0,0d0)
      bs  = dcmplx(6.266d0,0d0)
      N_s = dcmplx(0.094d0,0d0) 
     1    / betac(as+dcmplx(2d0,0d0),bs+dcmplx(1d0,0d0))
      ag  = dcmplx(1.943d0,0d0)
      bg  = dcmplx(8.000d0,0d0)
      N_g = dcmplx(0.238d0,0d0) 
     1    / betac(ag+dcmplx(2d0,0d0),bg+dcmplx(1d0,0d0))
*
*     Initialize PDFs to zero
*
      do iff=-6,6
         nff(iff) = (0d0,0d0)
      enddo
*
      nff(3)  = N_s * betac(N+as,bs+dcmplx(1d0,0d0))
      nff(2)  = N_l * betac(N+al,bl+dcmplx(1d0,0d0))
      nff(1)  = nff(3)
      nff(0)  = N_g * betac(N+ag,bg+dcmplx(1d0,0d0))
      nff(-1) = nff(2)
      nff(-2) = nff(1)
      nff(-3) = nff(3)
*
      return
      end
*
************************************************************************
      subroutine ZeroScalePDFs(N,npdf)
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
      double complex npdf(-6:6)
*
*     Initialize PDFs to zero
*
      do ipdf=-6,6
         npdf(ipdf) = (0d0,0d0)
      enddo
*
      npdf(1)  = zexp( - ( N - 1d0 ) * dlog(3d0) )
      npdf(2)  = 2d0 * npdf(1)
*
      return
      end
*
************************************************************************
*
*     HERAFitter paratrization.
*     There parameters are passed by means of a common block.
*
************************************************************************
      subroutine HERAFitterParametrization(N,npdf)
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
      double complex a2pN,a3p1,a8pN,a9p1
      double complex betac
      double complex ubar,dbar,uval,dval,glue,strg
*
      double precision ubarMELA(10)
      double precision dbarMELA(10)
      double precision uvalMELA(10)
      double precision dvalMELA(10)
      double precision glueMELA(10)
      double precision fsMELA,fcMELA
      common / HERAFitterParametersMELA / ubarMELA,dbarMELA,
     1     uvalMELA,dvalMELA,glueMELA,fsMELA,fcMELA
**
*     Output Variables
*
      double complex npdf(-6:6)
*
*     Initialize PDFs to zero
*
      a2pN = ubarMELA(2) + N
      a3p1 = ubarMELA(3) + 1d0
      a8pN = ubarMELA(8) + N
      a9p1 = ubarMELA(9) + 1d0
      ubar = ubarMELA(1) * betac(a2pN,a3p1)
     1     + ubarMELA(1) * ubarMELA(4)  * betac(a2pN+1d0,a3p1)
     2     + ubarMELA(1) * ubarMELA(5)  * betac(a2pN+2d0,a3p1)
     3     + ubarMELA(1) * ubarMELA(6)  * betac(a2pN+3d0,a3p1)
     4     + ubarMELA(1) * ubarMELA(10) * betac(a2pN+0.5d0,a3p1)
     5     + ubarMELA(7) * betac(a8pN,a9p1)
      ubar = ubar / ( 1d0 - fcMELA )
*
      a2pN = dbarMELA(2) + N
      a3p1 = dbarMELA(3) + 1d0
      a8pN = dbarMELA(8) + N
      a9p1 = dbarMELA(9) + 1d0
      dbar = dbarMELA(1) * betac(a2pN,a3p1)
     1     + dbarMELA(1) * dbarMELA(4)  * betac(a2pN+1d0,a3p1)
     2     + dbarMELA(1) * dbarMELA(5)  * betac(a2pN+2d0,a3p1)
     3     + dbarMELA(1) * dbarMELA(6)  * betac(a2pN+3d0,a3p1)
     4     + dbarMELA(1) * dbarMELA(10) * betac(a2pN+0.5d0,a3p1)
     5     + dbarMELA(7) * betac(a8pN,a9p1)
*
      a2pN = uvalMELA(2) + N
      a3p1 = uvalMELA(3) + 1d0
      a8pN = uvalMELA(8) + N
      a9p1 = uvalMELA(9) + 1d0
      uval = uvalMELA(1) * betac(a2pN,a3p1)
     1     + uvalMELA(1) * uvalMELA(4)  * betac(a2pN+1d0,a3p1)
     2     + uvalMELA(1) * uvalMELA(5)  * betac(a2pN+2d0,a3p1)
     3     + uvalMELA(1) * uvalMELA(6)  * betac(a2pN+3d0,a3p1)
     4     + uvalMELA(1) * uvalMELA(10) * betac(a2pN+0.5d0,a3p1)
     5     + uvalMELA(7) * betac(a8pN,a9p1)
*
      a2pN = dvalMELA(2) + N
      a3p1 = dvalMELA(3) + 1d0
      a8pN = dvalMELA(8) + N
      a9p1 = dvalMELA(9) + 1d0
      dval = dvalMELA(1) * betac(a2pN,a3p1)
     1     + dvalMELA(1) * dvalMELA(4)  * betac(a2pN+1d0,a3p1)
     2     + dvalMELA(1) * dvalMELA(5)  * betac(a2pN+2d0,a3p1)
     3     + dvalMELA(1) * dvalMELA(6)  * betac(a2pN+3d0,a3p1)
     4     + dvalMELA(1) * dvalMELA(10) * betac(a2pN+0.5d0,a3p1)
     5     + dvalMELA(7) * betac(a8pN,a9p1)
*
      a2pN = glueMELA(2) + N
      a3p1 = glueMELA(3) + 1d0
      a8pN = glueMELA(8) + N
      a9p1 = glueMELA(9) + 1d0
      glue = glueMELA(1) * betac(a2pN,a3p1)
     1     + glueMELA(1) * glueMELA(4)  * betac(a2pN+1d0,a3p1)
     2     + glueMELA(1) * glueMELA(5)  * betac(a2pN+2d0,a3p1)
     3     + glueMELA(1) * glueMELA(6)  * betac(a2pN+3d0,a3p1)
     4     + glueMELA(1) * glueMELA(10) * betac(a2pN+0.5d0,a3p1)
     5     + glueMELA(7) * betac(a8pN,a9p1)
*
      strg = fsMELA * dbar
      dbar = ( 1d0 - fsMELA ) * dbar
*
      npdf(-6) = (0d0,0d0)
      npdf(-5) = (0d0,0d0)
      npdf(-4) = (0d0,0d0)
      npdf(-3) = strg
      npdf(-2) = ubar
      npdf(-1) = dbar
      npdf(0)  = glue
      npdf(1)  = dval + dbar
      npdf(2)  = uval + ubar
      npdf(3)  = strg
      npdf(4)  = (0d0,0d0)
      npdf(5)  = (0d0,0d0)
      npdf(6)  = (0d0,0d0)
*
      return
      end
*
************************************************************************
*
*     Set the parameters of the HERAFitter parametrization
*
************************************************************************
      subroutine SetHERAFitterParametersMELA(ubar,dbar,
     1                                       uval,dval,
     2                                       glue,
     3                                       fs,fc)
*
      implicit none
**
*     Input parameters
*
      double precision ubar(10)
      double precision dbar(10)
      double precision uval(10)
      double precision dval(10)
      double precision glue(10)
      double precision fs,fc
**
*     Internal parameters
*
      double precision ubarMELA(10)
      double precision dbarMELA(10)
      double precision uvalMELA(10)
      double precision dvalMELA(10)
      double precision glueMELA(10)
      double precision fsMELA,fcMELA
      common / HERAFitterParametersMELA / ubarMELA,dbarMELA,
     1     uvalMELA,dvalMELA,glueMELA,fsMELA,fcMELA
*
      ubarMELA = ubar
      dbarMELA = dbar
      uvalMELA = uval
      dvalMELA = dval
      glueMELA = glue
      fsMELA   = fs
      fcMELA   = fc
*
      return
      end
