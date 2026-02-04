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
*     XFitter paratrization.
*     There parameters are passed by means of a common block.
*
************************************************************************
      subroutine XFitterParametrization(N,npdf)
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
      common / XFitterParametersMELA / ubarMELA,dbarMELA,
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
*     Set the parameters of the XFitter parametrization
*
************************************************************************
      subroutine SetXFitterParametersMELA(ubar,dbar,
     1                                    uval,dval,
     2                                    glue,
     3                                    fs,fc)
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
      common / XFitterParametersMELA / ubarMELA,dbarMELA,
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
*
************************************************************************
*
*     Gluon delta function multiplied by alpha_s
*     KL 5/11/24
*
************************************************************************
      subroutine GluonDelta(N,Q,nff)
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
**
*     Internal Variables
*
      double complex N,Q
      integer iff
      double complex alphas
      double precision mc,LDME
**
*     Output Variables
*
      double complex nff(-6:6)
*
*     Initialize PDFs to zero
*
      LDME=1d0
      mc=1.5d0
      alphas=PI*4d0*asEvolIni(1) 
      do iff=-6,6
         nff(iff) = (0d0,0d0)
      enddo
      nff(0) = alphas*PI*LDME/(mc**3*24d0)
*     write(*,*) "KL inside function !nff(0)", nff(0)
*
      return
      end
************************************************************************
*
*     constant
*     KL 5/11/24
*
************************************************************************
      subroutine Constant(N,Q,nff)
*
      implicit none
*
      include "../commons/alphas.h"
**
*     Internal Variables
*
      double complex N,Q
      integer iff
**
*     Output Variables
*
      double complex nff(-6:6)
*
*     Initialize PDFs to zero
*
      do iff=-6,6
         nff(iff) = (0d0,0d0)
      enddo
      nff(0) = 1/N
*      write(*,*) "CONSTANT KL inside function !nff(0)", nff(0)
*
      return
      end
      
************************************************************************
*
*     Gluon delta function at NLO KL 08/10/24
*
************************************************************************
      subroutine GluonDeltaNLO(N,Q,nff)
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
      
      double complex DPSI,PSI,S1,S2,dL,dLD,dT,dTLO,dTNLO,FdT
      double complex SHarm,S1N,S1N1,S1N2,S1N22,S2N2,PGG,FRAC
**
*     Internal Variables
*
      integer iff
      double precision mc,A,beta
      double complex N,Q,logscale,alphas
**
*     Output Variables
*
      double complex nff(-6:6)
      alphas=PI*4d0*asEvolIni(1) 
      mc = 1.5d0
      beta = (11d0*CA - 4d0*TR*NUMFL)/3d0
*
* 	In the file "andim_nlo.f" the harmonic functions are defined
*	S[L=1,N]	S1 = EMC + PSI(N+1d0)
*      	S[L=2,N]	S2 = ZETA2 - DPSI(N+1d0,1)
*      
      S1 = EMC + PSI(N+1d0)
      S1N = (EMC + PSI(N+1d0))/N
      S1N1 = (EMC + PSI(N+2d0))/(N+1d0)
      S1N2 = (EMC + PSI(N+3d0))/(N+2d0)
      S1N22 = (EMC + PSI(N-1d0))**2
      S2N2 = ZETA2 - DPSI(N-1d0,1)
      
*	Longitudinal Component
      dLD = 8d0 * mc**3 * N * (N-1d0)
      dL = (alphas)**2 * dLD**(-1)
*     
*	Transverse component
      FdT =  alphas * PI / (48d0 * mc**3)
      dTLO = FdT * (1d0, 0d0)
      logscale = log(Q/(2*mc))
      A = beta * (logscale + 13d0/6d0) + 4.15975994110758d0 
      SHarm = -6d0*(2*S1N - S1N1 + S1N2) - 3d0*(S1N22  + S2N2 )
      FRAC = 2d0 * (N**2 + N + 1d0)/(N * (N**2 - 1d0) * (N+2d0) )
      PGG = (logscale -1d0/2d0)*6d0*(-S1 +beta/6d0 + FRAC)
      dTNLO = SHarm + A + PGG
      dT = dTLO + dTNLO*FdT*asEvolIni(1) * 12.566370614359173d0/PI
*     Initialize PDFs to zero
*
      do iff=-6,6
         nff(iff) = (0d0,0d0)
      enddo
*    
      nff(0) = 2*dT + dL

      return
      end
      
************************************************************************
*
*     P-wave KL 14/03/25
*     J=0
*
************************************************************************
      subroutine P0wave(N,Q,nff)
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
      
      double complex DPSI,PSI,S1
      double complex T1,T2,T3
**
*     Internal Variables
*
      integer iff
      double precision Nc,QJ
      double precision mc,mulambda
      double complex N,Q,logscale,alphas
**
*     Output Variables
*
      double complex nff(-6:6)
      alphas=PI*4d0*asEvolIni(1)
      Nc=3d0 
      QJ=1/4d0
      mc = 1.5d0
      S1 = EMC + PSI(N+1d0)
      mulambda=mc
      T1=QJ-0.5d0*log(mulambda**2/(4*mc**2))
      T2=-S1
      T3=85/8*1/(N+1)-26/8*1/(N+2)+45/4*1/(N)-27/4*1/(N+1)
*
      do iff=-6,6
         nff(iff) = (0d0,0d0)
      enddo
*    
      nff(0) = (4d0,0d0)/(9d0*Nc)*(T1+T2+T3)*alphas**2/mc**5

      return
      end
************************************************************************
*
*     P-wave KL 14/03/25
*     J=1
*
************************************************************************
      subroutine P1wave(N,Q,nff)
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
      
      double complex DPSI,PSI,S1
      double complex T1,T2,T3
**
*     Internal Variables
*
      integer iff
      double precision Nc,QJ
      double precision mc,mulambda
      double complex N,Q,logscale,alphas
**
*     Output Variables
*
      double complex nff(-6:6)
      alphas=PI*4d0*asEvolIni(1)
      Nc=3d0 
      QJ=3/8d0/(3d0)
      mc = 1.5d0
      S1 = EMC + PSI(N+1d0)
      mulambda=mc
      T1=QJ-0.5d0*log(mulambda**2/(4*mc**2))
      T2=-S1
      T3=-3/4*(1/(N+1)+4/(N+2))/3d0
*
      do iff=-6,6
         nff(iff) = (0d0,0d0)
      enddo
*    
      nff(0) = (4d0,0d0)/(9d0*Nc)*(T1+T2+T3)*alphas**2/mc**5

      return
      end
************************************************************************
*
*     P-wave KL 14/03/25
*     J=2
*
************************************************************************
      subroutine P2wave(N,Q,nff)
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
      
      double complex DPSI,PSI,S1,S2
      double complex T1,T2,T3
**
*     Internal Variables
*
      integer iff
      double precision Nc,QJ
      double precision mc,mulambda
      double complex N,Q,logscale,alphas
**
*     Output Variables
*
      double complex nff(-6:6)
      alphas=PI*4d0*asEvolIni(1)
      Nc=3d0 
      QJ=7/8d0/(5d0)
      mc = 1.5d0
      S1 = EMC + PSI(N+1d0)
      S1 = EMC + PSI(N+2d0)
      mulambda=mc
      T1=QJ-0.5d0*log(mulambda**2/(4*mc**2))
      T2=-S1
      T3=5/4*(11/(N+1)-4/(N+2)) +9*(-2*S1/N + S2/(N+1)) 
*
      do iff=-6,6
         nff(iff) = (0d0,0d0)
      enddo
*    
      nff(0) = (4d0,0d0)/(9d0*Nc)*(T1+T2+T3/5d0)*alphas**2/mc**5

      return
      end
************************************************************************
*
*     Gluon delta function at NLO KL 05/06/25
*
************************************************************************
      subroutine g3s18NLO(N,Q,nff)
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
      
      double complex DPSI,PSI,S1,S2,S3,S4,S5
      double complex PGG,S,LOGR
      double complex T1,T2,T3,T4,T5
**
*     Internal Variables
*
      integer iff
      double precision mc,beta,NC,alphas
      double complex N,Q
**
*     Output Variables
*
      double complex nff(-6:6)
*
* 	In the file "andim_nlo.f" the harmonic functions are defined
*	S[L=1,N]	S1 = EMC + PSI(N+1d0)
*      	S[L=2,N]	S2 = ZETA2 - DPSI(N+1d0,1)
*     
      NC=3d0
      alphas=PI*4d0*asEvolIni(1) 
      mc = 1.5d0
      beta = (11d0*CA - 4d0*TR*3d0)/3d0
      LOGR=log(Q**2/(4d0*mc**2))
      
      T1=PI*alphas/24d0
      
      T2=beta/NC*(LOGR+13d0/3d0)+4d0/NC**2-PI**2/3d0+16d0/3d0*log(2d0)
      
      S=EMC + PSI(N+1d0)
      PGG=2*NC*(-S+1/(N-1d0)-1/(N)+1/(N+1d0)-1/(N+2d0)+beta/(2*NC))
      T3=1/NC*PGG*(LOGR-1d0)
      
      T4=2d0/(N-1d0)-2d0/(N)
 
      S1 = (EMC + PSI(N-1d0))**2 + ZETA2 - DPSI(N-1d0,1)
      S2 = (EMC + PSI(N))**2 + ZETA2 - DPSI(N,1)
      S3 = (EMC + PSI(N+1d0))**2 + ZETA2 - DPSI(N+1d0,1)
      S4 = (EMC + PSI(N+2d0))**2 + ZETA2 - DPSI(N+2d0,1)
      S5 = (EMC + PSI(N+3d0))**2 + ZETA2 - DPSI(N+3d0,1)
      T5= -2d0*(S1-2d0*S2+3d0*S3-2d0*S4+S5)
*     Initialize PDFs to zero
*
      do iff=-6,6
         nff(iff) = (0d0,0d0)
      enddo
*    
      nff(0)=(T1+alphas**2/12d0/CF*(T2+T3+T4+T5))/mc**3

      return
      end

************************************************************************
*
*     TEST FUNC
*    
*
************************************************************************
      subroutine TEST(N,Q,nff)
*
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
      
**
*     Internal Variables
*
      integer iff
      double complex N,Q
      double complex PSI,DPSI
      double complex S1N,S1Np1,S1Np2,S1Np3
      double complex Talpha
      double complex Ta0,Ta1,Ta2,Ta3
      double precision alpha1,alpha2,alpha3
**
*     Output Variables
*
      double complex nff(-6:6)
      do iff=-6,6
         nff(iff) = (0d0,0d0)
      enddo
*
*    define params from table 3
*    
      alpha1=-4.9866d-3
      alpha2=9.8448d-3
      alpha3=1.9512d-2
      
* In the file "andim_nlo.f" the harmonic functions are defined
*	S[L=1,N]	S1 = EMC + PSI(N+1d0)
*      	S[L=2,N]	S2 = ZETA2 - DPSI(N+1d0,1)
*   
*--- Harmonic sums
*
      S1N   = EMC + PSI(N+1d0)
      S1Np1 = EMC + PSI(N+2d0)
      S1Np2 = EMC + PSI(N+3d0)
      S1Np3 = EMC + PSI(N+4d0)
*
*  co-ef of alpha part
      Ta0=-S1N/N
      Ta1=-S1Np1/(N+1d0)
      Ta2=-S1Np2/(N+2d0)
      Ta3=-S1Np3/(N+3d0)

      
*  making alpha, beta, etc.     
      Talpha=(alpha1+alpha2+alpha3)*Ta0
     & +(-alpha1-2d0*alpha2-3d0*alpha3)*Ta1
     & +(alpha2+3d0*alpha3)*Ta2
     & -alpha3*Ta3
      
      
      nff(0)=Talpha

      return
      end
************************************************************************
*
*     TEST FUNC 2
*    
*
************************************************************************
      subroutine TESTS(N,Q,nff)
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
**
*     Internal Variables
*
      integer iff
      double complex SMKL8
      double complex N,Q
      double precision m
**
**
*     Output Variables
*
      double complex nff(-6:6)
      do iff=-6,6
         nff(iff) = (0d0,0d0)
      enddo
      m=0d0
      nff(0)=SMKL8(N,m) 
      
      return
      end   
************************************************************************
*
*     g->3S18 NLO from MA: Phys. Rev. D, 89(9):094029, 2014
*    
*
************************************************************************
      subroutine Ma3S18(N,Q,nff)
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
**
*     Internal Variables
*
      integer iff
      double complex N,Q,LL
      double complex mc,alphas,beta0
      double complex PSI,S,PGG
      double complex D1,D2,D21,D22,D23,D24
**
**
*     Output Variables
*
      double complex nff(-6:6)
      double complex GH
      external GH
      
      do iff=-6,6
         nff(iff) = (0d0,0d0)
      enddo
      mc=1.5d0
      LL=Q**2/(4d0*mc**2)
*     beta def with nf=3      
      beta0=(11d0*3d0-2d0*3d0)/6d0
      alphas=PI*4d0*asEvolIni(1) 
*	S[L=1,N]	S1 = EMC + PSI(N+1d0)      
      S=EMC + PSI(N+1d0)
      PGG=6d0*(-S+1d0/(N-1d0)-1d0/N+1d0/(N+1d0)
     & -1d0/(N+2d0)+beta0/6d0)
      D1=(1d0,0d0)/(24d0)
      D21=beta0/3d0*(log(LL)+13/3)
     & +4d0/3d0**2-PI**2/3d0+16d0/3d0*log(2d0)
      D22=1d0/3d0*PGG*(log(LL)-1d0)
      D23=2d0/(N-1d0)-2d0/N
      D24=-2*(GH(N+3d0)-2*GH(N+2d0)+3*GH(N+1d0)-2*GH(N)+GH(N-1d0))
      D2=(1d0,0d0)*(D21+D22+D23+D24)/(12d0*CF)
*      D2=(1d0,0d0)/(12d0*CF)*(D22)	
	
      nff(0)=((PI*alphas*D1) + (alphas**2*D2)) / (mc**3)
*      nff(0)=(PI*alphas*D1)/mc**3
      	
*      nff(0)=(alphas**2*D2)/mc**3
      return
      end   
C=====================================================================
C   Function GHARM(N): computes S[L=1,N]**2 + S[L=2,N]
C=====================================================================
*    *
* In the file "andim_nlo.f" the harmonic functions are defined
*	S[L=1,N]	S1 = EMC + PSI(N+1d0)
*      	S[L=2,N]	S2 = ZETA2 - DPSI(N+1d0,1)
*	S[L=1,N-1]	S1 = EMC + PSI(N)
*      	S[L=2,N-1]	S2 = ZETA2 - DPSI(N,1)
*   
      double complex function GH(X)
      implicit none
      double complex X
      double complex PSI, DPSI
      include "../commons/consts.h"
      include "../commons/alphas.h"
      include "../commons/colfact.h"

      GH = (EMC + PSI(X))**2 + (ZETA2 - DPSI(X,1))

      return
      end
************************************************************************
*
*     v4 gluon divergent part 
*     JHEP11 2012 020, eq. 7.7
*     Does not include finite piece
*
************************************************************************
      double complex function G4G3S11DIST(N,muL)
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
*
      double complex DPSI,PSI
      double complex S1N,S1Nm1,S1Np1,S1Np2,S1Nm2
      double complex S2N,S2Nm1,S2Np1
      double complex PolyN,PolyNp1
      double complex bracket
      double complex N
      double complex T1,T2,T3,T4,T5,T6,T7,T8,T9
      double complex T91,T92,T93,T93i,T93ii
      double precision mc,NC,dimd
      double complex muL,Lmu
      double precision prefacN,prefacD,prefac
*
*--- parameters
*
      NC   = 3d0
      mc   = 1.5d0
      dimd = 4d0
*
*
      Lmu = log(muL/(2d0*mc))


* In the file "andim_nlo.f" the harmonic functions are defined
*	S[L=1,N]	S1 = EMC + PSI(N+1d0)
*      	S[L=2,N]	S2 = ZETA2 - DPSI(N+1d0,1)
*
      S1N   = EMC + PSI(N+1d0)
      S1Nm1 = EMC + PSI(N    )
      S1Np1 = EMC + PSI(N+2d0)
      S1Np2 = EMC + PSI(N+3d0)
      S1Nm2 = EMC + PSI(N-1d0)
*
      S2N   = ZETA2 - DPSI(N+1d0,1)
      S2Np1 = ZETA2 - DPSI(N+2d0,1)
      S2Nm1 = ZETA2 - DPSI(N    ,1)
*
*--- Polygamma combinations
*

      PolyN=(S1N**2-S2N+2d0*ZETA2)/N-2d0*S1N/N**2
      PolyNp1=(S1Np1**2-S2Np1+2d0*ZETA2)/(N+1d0)-2d0*S1Np1/(N+1d0)**2

*
*--- build bracket
*
      T1=1d0/24d0 - ZETA2 - Lmu/3d0 + Lmu*Lmu
 
      T2= - S1Nm1*(1d0/3d0 - 2d0*Lmu)
      T3= S1Nm1*S1Nm1 + S2Nm1

      T4=(-1d0/24d0)*(104d0/N - 29d0/(N+1d0) - 10d0/(N+2d0))

      T5=(7d0/2d0)*(1d0/(N-1d0) - S1Nm2/(N-2d0) - S1Nm1/(N-1d0))

      T6=(8d0/N - 21d0/(N+1d0) + 10d0/(N+2d0))*Lmu/4d0
     
      T7=1/4d0*(5d0*S1N/(N)+36d0*S1Nm1/(N-1d0)-25d0*S1Np1/(N+1d0)
     * + 6d0*S1Np2/(N+2d0))

      T8=39d0/(4d0*(N+1d0)**2) - 3d0/(2d0*(N+2d0)**2)  
     &     -  2d0*(S2N - ZETA2)

      T91= (13d0/2d0)*PolyN - (7d0/2d0)*PolyNp1
      T92= Lmu*((13d0/2d0)*S1N/N - (7d0/2d0)*S1Np1/(N+1d0))
      T93i=ZETA2/N -S1N/(N**2)
      T93ii=ZETA2/(N+1d0) -S1Np1/(N+1d0)**2
      T93=-13d0/2d0*T93i+7d0/2d0*T93ii
*
*--- prefactor
*
      prefacN = 2d0*(NC*NC - 4d0) 
      prefacD = (3d0*PI*(dimd-1d0)**3*NC**3)
      prefac = prefacN/prefacD
      T9 = T91+T92+T93

      bracket=T1+T2+T3+T4+T5+T6+T7+T8+T9

      G4G3S11DIST=prefac*bracket
*
      return
      end
************************************************************************
*
*     function v0 gluon paramterised FF kernel
*     JHEP11 2012 020, eq. B2
*
************************************************************************
      double complex function G0G3S11(N)
*
      implicit none
*
      double complex N
      double complex PSI,DPSI
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
*
*     Internal variables
*
      double complex S1N,S1Np1,S1Np2,S1Np3
      double complex S2N,S2Np1,S2Np2
      double complex Talpha,Tbeta,Tmu,Tnu,Tomega
      double complex Ta0,Ta1,Ta2,Ta3
      double complex Tb0,Tb1,Tb2
      double complex Tnu1,Tnu2,Tnu3,Tnu4
      double complex To11,To12
      double precision alpha1,alpha2,alpha3,beta1,beta2
      double precision mu1,nu1,nu2,nu3,nu4,omega11,omega12
*
*     Parameters (Table 3)
*
      alpha1=-4.9866d-3
      alpha2= 9.8448d-3
      alpha3= 1.9512d-2
      beta1 =-5.1697d-6
      beta2 = 1.0462d-2
      mu1   =-1.8921d-3
      nu1   = 1.2154d-3
      nu2   = 1.3039d-3
      nu3   =-2.7246d-3
      nu4   =-1.4814d-3
      omega11=-1.6910d-2
      omega12= 3.8110d-2
*
*--- Harmonic sums
*
      S1N   = EMC + PSI(N+1d0)
      S1Np1 = EMC + PSI(N+2d0)
      S1Np2 = EMC + PSI(N+3d0)
      S1Np3 = EMC + PSI(N+4d0)

      S2N   = ZETA2 - DPSI(N+1d0,1)
      S2Np1 = ZETA2 - DPSI(N+2d0,1)
      S2Np2 = ZETA2 - DPSI(N+3d0,1)
*
*--- Coefficients
*
      Ta0=-S1N/N
      Ta1=-S1Np1/(N+1d0)
      Ta2=-S1Np2/(N+2d0)
      Ta3=-S1Np3/(N+3d0)

      Tb0=(S1N*S1N+S2N)/N
      Tb1=(S1Np1*S1Np1+S2Np1)/(N+1d0)
      Tb2=(S1Np2*S1Np2+S2Np2)/(N+2d0)

      Tnu1=2d0/(N+1d0)**3
      Tnu2=2d0/(N+2d0)**3
      Tnu3=2d0/(N+3d0)**3
      Tnu4=2d0/(N+4d0)**3

      To11=1d0/(N+1d0)-1d0/(N+2d0)
      To12=1d0/(N+1d0)-2d0/(N+2d0)+1d0/(N+3d0)
*
*--- Assemble pieces
*
      Talpha=(alpha1+alpha2+alpha3)*Ta0
     & +(-alpha1-2d0*alpha2-3d0*alpha3)*Ta1
     & +(alpha2+3d0*alpha3)*Ta2
     & -alpha3*Ta3

      Tbeta=(beta1+beta2)*Tb0
     & +(-beta1-2d0*beta2)*Tb1
     & + beta2*Tb2

      Tmu=-mu1/(N+1d0)**2
      Tnu=nu1*Tnu1+nu2*Tnu2+nu3*Tnu3+nu4*Tnu4
      Tomega=omega11*To11+omega12*To12
*
*--- Final result
*
      G0G3S11 = Talpha + Tbeta + Tmu + Tnu + Tomega

      return
      end
************************************************************************
*
*     Subroutine calling gluon FF kernel
*
************************************************************************
      subroutine d0g3S11(N,Q,nff)
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
*
*     Function declaration
*
      double complex G0G3S11
      external G0G3S11
*
*     Internal variables
*
      integer iff
      double complex N,Q
*
*     Output
*
      double complex nff(-6:6)

      do iff=-6,6
         nff(iff) = (0d0,0d0)
      enddo
*
*     Call function
*
      nff(0) = G0G3S11(N)

      return
      end
************************************************************************
*
*     function v2 gluon paramterised FF kernel
*     JHEP11 2012 020, eq. B2
*
************************************************************************
      double complex function G2G3S11(N)
*
      implicit none
*
      double complex N
      double complex PSI,DPSI
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
*
*     Internal variables
*
      double complex S1N,S1Np1,S1Np2,S1Np3
      double complex S2N,S2Np1,S2Np2
      double complex Tb,Talpha,Tbeta,Tmu,Tnu,Tomega
      double complex Ta0,Ta1
      double complex Tb0,Tb1
      double complex Tnu1
      double complex To11,To12,To13,To14
      double precision b,alpha1,beta1
      double precision mu1,mu2,nu1,omega11,omega12,omega13,omega14
*
*     Parameters (Table 3)
*
      b=1.4710d-2
      alpha1=-1.8127d-2
      beta1 =-1.2282d-2
      mu1   =4.1069d-3
      mu2   =1.8341d-2
      nu1   =-1.3653d-3
      omega11=-2.1832d-2
      omega12=4.1531d-3
      omega13=8.4949d-4
      omega14=-3.8207d-3
*
*--- Harmonic sums
*
      S1N   = EMC + PSI(N+1d0)
      S1Np1 = EMC + PSI(N+2d0)
      S1Np2 = EMC + PSI(N+3d0)
      S1Np3 = EMC + PSI(N+4d0)

      S2N   = ZETA2 - DPSI(N+1d0,1)
      S2Np1 = ZETA2 - DPSI(N+2d0,1)
      S2Np2 = ZETA2 - DPSI(N+3d0,1)
*
*--- Coefficients
*
      Ta0=-S1N/N
      Ta1=-S1Np1/(N+1d0)

      Tb0=(S1N*S1N+S2N)/N
      Tb1=(S1Np1*S1Np1+S2Np1)/(N+1d0)

      Tnu1=2d0/(N+1d0)**3

      To11=1d0/(N+1d0)-1d0/(N+2d0)
      To12=1d0/(N+1d0)-2d0/(N+2d0)+1d0/(N+3d0)
      To13 = 1d0/(N+1d0) - 3d0/(N+2d0) + 3d0/(N+3d0) - 1d0/(N+4d0)
      To14=1d0/(N+1d0)-4d0/(N+2d0)+6d0/(N+3d0)-4d0/(N+4d0)+1d0/(N+5d0)

*
*--- Assemble pieces
*
      Tb=b/(N+1d0)
      
      Talpha=alpha1*Ta0-alpha1*Ta1

      Tbeta=beta1*Tb0-beta1*Tb1

      Tmu=-mu1/(N+1d0)**2-mu2/(N+2d0)**2
      
      Tnu=nu1*Tnu1
      Tomega=omega11*To11+omega12*To12+omega13*To13+omega14*To14
*
*--- Final result
*
      G2G3S11 = Tb + Talpha + Tbeta + Tmu + Tnu + Tomega

      return
      end
************************************************************************
*
*     Subroutine calling gluon FF kernel
*
************************************************************************
      subroutine d2g3S11(N,Q,nff)
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
*
*     Function declaration
*
      double complex G2G3S11
      external G2G3S11
*
*     Internal variables
*
      integer iff
      double complex N,Q
*
*     Output
*
      double complex nff(-6:6)

      do iff=-6,6
         nff(iff) = (0d0,0d0)
      enddo
*
*     Call function
*
      nff(0) = G2G3S11(N)

      return
      end
************************************************************************
*
*     function v4 gluon paramterised FF kernel [Finite part]
*     JHEP11 2012 020, eq. B2
*
************************************************************************
      double complex function G4G3S11(N)
*
      implicit none
*
      double complex N
      double complex PSI,DPSI
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
*
*     Internal variables
*
      double complex S1N,S1Np1,S1Np2,S1Np3
      double complex S2N,S2Np1,S2Np2
      double complex Tb,Talpha,Tbeta,Tmu,Tnu,Tomega
      double complex Ta0,Ta1,Ta2
      double complex Tb0,Tb1
      double complex Tnu1
      double complex To11,To12,To13,To14
      double precision b,alpha1,alpha2,beta1
      double precision mu1,mu2,mu3,mu4,nu1,omega11
*
*     Parameters (Table 3)
*
      b=-3.6531d-2
      alpha1=5.2157d-2
      alpha2=2.4588d-3
      beta1 =3.6020d-2
      mu1   =-7.3565d-3
      mu2   =-2.5387d-2
      mu3   =-9.3477d-3
      mu4   =-3.6926d-3
      nu1   =1.3839d-3
      omega11=7.7264d-2
*
*--- Harmonic sums
*
      S1N   = EMC + PSI(N+1d0)
      S1Np1 = EMC + PSI(N+2d0)
      S1Np2 = EMC + PSI(N+3d0)
      S1Np3 = EMC + PSI(N+4d0)

      S2N   = ZETA2 - DPSI(N+1d0,1)
      S2Np1 = ZETA2 - DPSI(N+2d0,1)
      S2Np2 = ZETA2 - DPSI(N+3d0,1)
*
*--- Coefficients
*
      Ta0=-S1N/N
      Ta1=-S1Np1/(N+1d0)
      Ta2=-S1Np2/(N+2d0)

      Tb0=(S1N*S1N+S2N)/N
      Tb1=(S1Np1*S1Np1+S2Np1)/(N+1d0)

      Tnu1=2d0/(N+1d0)**3

      To11=1d0/(N+1d0)-1d0/(N+2d0)

*
*--- Assemble pieces
*
      Tb=b/(N+1d0)
      Talpha=(alpha1+alpha2)*Ta0
     & +(-alpha1-2d0*alpha2)*Ta1
     & +alpha2*Ta2
     

      Tbeta=beta1*Tb0-beta1*Tb1

      Tmu=-mu1/(N+1d0)**2-mu2/(N+2d0)**2
     &  -mu3/(N+3d0)**2
     &  -mu4/(N+4d0)**2
      
      Tnu=nu1*Tnu1
      Tomega=omega11*To11
*
*--- Final result
*
      G4G3S11 = Tb + Talpha + Tbeta + Tmu + Tnu + Tomega

      return
      end
************************************************************************
*
*     c->3S11 at v^4 KL (Function) 04/02/26
*
************************************************************************
      double complex function Fv4c3S11(N)
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
      
      double complex SMKL10
      double complex N
      double complex T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11
      double precision m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11
      double precision a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11

*     powers (m)
      m1=1d0
      m2=2d0
      m3=3d0
      m4=4d0
      m5=5d0
      m6=6d0
      m7=7d0
      m8=8d0
      m9=9d0
      m10=10d0
      m11=11d0
      
*     coefs of powers (a)
      a1=1639680d0
      a2=-12648960d0
      a3=53962240d0
      a4=-141141120d0
      a5=225693536d0
      a6=-222599040d0
      a7=136706160d0
      a8=-52504200d0
      a9=12461885d0 
      a10=-1668066d0 
      a11=97885d0      
*    
      T1=a1*SMKL10(N,m1)
      T2=a2*SMKL10(N,m2)
      T3=a3*SMKL10(N,m3)
      T4=a4*SMKL10(N,m4)
      T5=a5*SMKL10(N,m5)
      T6=a6*SMKL10(N,m6) 
      T7=a7*SMKL10(N,m7)       
      T8=a8*SMKL10(N,m8)       
      T9=a9*SMKL10(N,m9) 
      T10=a10*SMKL10(N,m10)       
      T11=a11*SMKL10(N,m11) 
      Fv4c3S11=2d0/164025d0*(T1+T2+T3+T4+T5+T6+T7+T8+T9+T10+T11)
*
      return
      end      
************************************************************************
*
*     Subroutine calling gluon FF kernel
*
************************************************************************
      subroutine d4g3S11(N,Q,nff)
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
*
*     Function declaration
*
      double complex G4G3S11,G4G3S11DIST
      external G4G3S11,G4G3S11DIST
*
*     Internal variables
*
      integer iff
      double complex N,Q,muL
      double precision mc
*
*     Output
*
      double complex nff(-6:6)
      
      mc=1.5d0
      muL=(1.5d0,0d0)
      
      do iff=-6,6
         nff(iff) = (0d0,0d0)
      enddo
*
*     Call function
*
      nff(0) = G4G3S11(N)+G4G3S11DIST(N,muL) 

      return
      end
      
      
************************************************************************
*
*     Subroutine calling v4 corr to charm AND gluon
*
************************************************************************
      subroutine d4gc3S11(N,Q,nff)
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
*
*     Function declaration
*
      double complex G4G3S11,G4G3S11DIST,Fv4c3S11
      external G4G3S11,G4G3S11DIST,Fv4c3S11
*
*     Internal variables
*
      integer iff
      double complex N,Q,muL
      double precision mc,alphas      
*
*     Output
*
      double complex nff(-6:6)
      alphas=0.31494871292d0
      mc=1.5d0
      muL=(1.5d0,0d0)
      
      do iff=-6,6
         nff(iff) = (0d0,0d0)
      enddo
*
*     Call function
*
      nff(0) = alphas**3*(G4G3S11(N)+G4G3S11DIST(N,muL)) 
      nff(4) = alphas**2*(Fv4c3S11(N))
      nff(-4) = nff(4) 
      return
      end      
************************************************************************
*
*     testing function which works
*
************************************************************************
      double complex function TESTFUNC(X)
      implicit none
      double complex X
      double complex PSI, DPSI
      include "../commons/consts.h"
      include "../commons/alphas.h"
      include "../commons/colfact.h"

      TESTFUNC = 1/X

      return
      end
*
*   subroutine calling above function
*
      subroutine TESTF(N,Q,nff)
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
      
**
*     Internal Variables
*
      integer iff
      double complex N,Q
      double complex TESTFUNC
      external TESTFUNC
**
**
*     Output Variables
*
      double complex nff(-6:6)
      do iff=-6,6
         nff(iff) = (0d0,0d0)
      enddo

      nff(0)=TESTFUNC(N)
      
      return
      end
************************************************************************
*
*     c->3S11 at LO KL (Function F) 12/05/25
*
************************************************************************
      double complex function Fc3S11(N)
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
      
      double complex SMKL6
      double complex N,Q
      double complex T1,T2,T3,T4,T5,T6,T7
      double precision m1,m2,m3,m4,m5,m6,m7
      double precision a1,a2,a3,a4,a5,a6,a7

*     powers (m)
      m1=1d0
      m2=2d0
      m3=3d0
      m4=4d0
      m5=5d0
      m6=6d0
      m7=7d0

*     coefs of powers (a)
      a1=16d0
      a2=-64d0
      a3=152d0
      a4=-208d0
      a5=141d0
      a6=-42d0
      a7=5d0
*    
      T1=a1*SMKL6(N,m1)
      T2=a2*SMKL6(N,m2)
      T3=a3*SMKL6(N,m3)
      T4=a4*SMKL6(N,m4)
      T5=a5*SMKL6(N,m5)
      T6=a6*SMKL6(N,m6) 
      T7=a7*SMKL6(N,m7) 
      Fc3S11=16d0/243d0*(T1+T2+T3+T4+T5+T6+T7)
*
      return
      end
*
************************************************************************
*
*     subroutine calling above function
*
************************************************************************
      subroutine c3S11LO(N,Q,nff)
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
*
*     Internal variables
*
      integer iff
      double complex N, Q
      double complex Fc3S11
      external Fc3S11
*
*     Output
*
      double complex nff(-6:6)
*
      do iff = -6, 6
         nff(iff) = (0d0,0d0)
      end do
*
      nff(4)  = Fc3S11(N)
      nff(-4) = nff(4)
*
      return
      end
************************************************************************
*
*     c->3S11 at v^2 KL (Function G) 12/05/25
*
************************************************************************
      double complex function Gc3S11(N)
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
      
      double complex SMKL8
      double complex N
      double complex T1,T2,T3,T4,T5,T6,T7,T8,T9
      double precision m1,m2,m3,m4,m5,m6,m7,m8,m9
      double precision a1,a2,a3,a4,a5,a6,a7,a8,a9

*     powers (m)
      m1=1d0
      m2=2d0
      m3=3d0
      m4=4d0
      m5=5d0
      m6=6d0
      m7=7d0
      m8=8d0
      m9=9d0

*     coefs of powers (a)
      a1=-10752d0
      a2=62976d0
      a3=-209024d0
      a4=417536d0
      a5=-477536d0
      a6=304928d0
      a7=-106328d0
      a8=19664d0
      a9=-1464d0      
*    
      T1=a1*SMKL8(N,m1)
      T2=a2*SMKL8(N,m2)
      T3=a3*SMKL8(N,m3)
      T4=a4*SMKL8(N,m4)
      T5=a5*SMKL8(N,m5)
      T6=a6*SMKL8(N,m6) 
      T7=a7*SMKL8(N,m7)       
      T8=a8*SMKL8(N,m8)       
      T9=a9*SMKL8(N,m9) 
      Gc3S11=1d0/2187d0*(T1+T2+T3+T4+T5+T6+T7+T8+T9)
*
      return
      end
************************************************************************
*
*     subroutine calling above function
*
************************************************************************
      subroutine c3S11REL(N,Q,nff)
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
*
*     Internal variables
*
      integer iff
      double complex N, Q
      double complex Gc3S11
      external Gc3S11
*
*     Output
*
      double complex nff(-6:6)
*
      do iff = -6, 6
         nff(iff) = (0d0,0d0)
      end do
*
      nff(4)  = Gc3S11(N)
      nff(-4) = nff(4)
*
      return
      end
************************************************************************
*
*     c->3S11 NLO (PQQ) KL 12/05/25
*
************************************************************************
      double complex function PQQc3S11(N)
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
      
      double complex PSI
      double complex N
      double complex S1Np1,S1Nm1,DELTA
      
*	S[L=1,N]	S1 = EMC + PSI(N+1d0)
      
      S1Nm1 = EMC + PSI(N)
      S1Np1 = EMC + PSI(N+2d0)
      DELTA = (1.5d0,0d0)
      PQQc3S11=2d0*CF*(-S1Nm1-S1Np1+DELTA)
*
      return
      end
************************************************************************
*
*     subroutine calling above function
*
************************************************************************
      subroutine c3S11PQQ(N,Q,nff)
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
*
*     Internal variables
*
      integer iff
      double complex N, Q
      double complex PQQc3S11,Fc3S11,P,F
      external PQQc3S11,Fc3S11
*
*     Output
*
      double complex nff(-6:6)
*
      do iff = -6, 6
         nff(iff) = (0d0,0d0)
      end do
*
      P=PQQc3S11(N)
      F=Fc3S11(N)
      nff(4)  = P*F
      nff(-4) = nff(4)
*
      return
      end
************************************************************************
*
*     c->3S11 NLO (POLY) KL 12/05/25
*
************************************************************************
      double complex function POLYc3S11(N)
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
      
      double complex N
      double complex FACT,F
      
      FACT=0.6981317008d0
      F=-9.01726/(N+10d0)+18.22777/(N+9d0)
      F=F+16.11858/(N+8d0)-82.54936/(N+7d0) 
      F=F+106.57565/(N+6d0)-72.30107/(N+5d0)
      F=F+28.85798/(N+4d0)-6.70607/(N+3d0)
      F=F+0.84950/(N+2d0)-0.05376/(N+1d0)-0.00205/N
      POLYc3S11=F*FACT
*
      return
      end
************************************************************************
*
*     subroutine calling above function
*
************************************************************************
      subroutine c3S11NLO(N,Q,nff)
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
*
*     Internal variables
      integer iff
      double complex N, Q
      double complex muR0C,mc,beta0
      double complex F,cAS3
      double complex POLYc3S11,Fc3S11,PQQc3S11
      external POLYc3S11,Fc3S11,PQQc3S11
*
*     Output
*
      double complex nff(-6:6)
*
      do iff = -6, 6
         nff(iff) = (0d0,0d0)
      end do
*
 
      beta0=(11d0*CA - 4d0*TR*3d0)/3d0
      mc = 1.5d0
      muR0C = 2*mc
      F = Fc3S11(N)
      cAS3 = ( 
     *           +F/(2d0*PI)*beta0*log(muR0C**2/(2d0*mc)**2)
     *           +PQQc3S11(N)*F/(4d0*PI)*log(Q**2/(3d0*mc)**2) 
     *           +POLYc3S11(N)  
     *           )
      nff(4)  = cAS3
      nff(-4) = nff(4)
*
      return
      end      
************************************************************************
*
*     g->3S11 FENG 12/05/25
*
************************************************************************
      double complex function c0g3S11(N)
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
      
      double complex N
      double precision FACT
      double precision a0,a1,a2,a3,a4,a5,a6,a7
      double complex T20,T21,T22,T23,T24,T25,T26,T27
      double complex T1,T2,T3
      double complex SMKL5,LOG2KL
      external SMKL5,LOG2KL
      
      
      FACT = 1/(PI*4d0*3645d0)
*  a's for (z-2)^-5 term (T2)
      a0=-6112d0
      a1=22632d0
      a2=-35176d0
      a3=30160d0
      a4=-15685d0
      a5=5032d0
      a6=-927d0
      a7=76d0

* Terms for T2
      T20=a0*SMKL5(N,0d0)
      T21=a1*SMKL5(N,1d0)
      T22=a2*SMKL5(N,2d0)
      T23=a3*SMKL5(N,3d0)
      T24=a4*SMKL5(N,4d0)
      T25=a5*SMKL5(N,5d0)
      T26=a6*SMKL5(N,6d0)
      T27=a7*SMKL5(N,7d0)

*  Terms            
      T1=-480d0*(LOG2KL(N,2d0)-19d0*LOG2KL(N,1d0)+36d0*LOG2KL(N,0d0))
      T2=-64d0*(T20+T21+T22+T23+T24+T25+T26+T27)      
      T3=480d0*(-1d0/(N+2d0)**2-1d0/(N+1d0)**2)
      
      c0g3S11=FACT*(T1+T2+T3)
*
      return
      end
      
      double complex function c1g3S11(N)
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
      
      double complex N
      double complex F
      double precision a0,a1,a2,a3,a4,a5,a6,a7,a8
      
      a0=0.068789941142313182
      a1=-0.28294758545234799
      a2=0.43755027838051319
      a3=-0.24403801281005144
      a4=-0.12154015467967838
      a5=0.27499430533498526
      a6=-0.20288020418956876
      a7=0.092142697889357805
      a8=-0.022071291547035798

      F=a8/(N+8d0)+a7/(N+7d0)+a6/(N+6d0)+a5/(N+5d0)
      F=F+a4/(N+4d0)+a3/(N+3d0)+a2/(N+2d0)+a1/(N+1d0)+a0/N
      c1g3S11=F
*
      return
      end
*************************************************************
*   c1,c0 subroutines
************************************************************      
      subroutine g3S11c0(N,Q,nff)
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
      
**
*     Internal Variables
*
      integer iff
      double complex N,Q
      double complex c0g3S11
      external c0g3S11
**
**
*     Output Variables
*
      double complex nff(-6:6)
      do iff=-6,6
         nff(iff) = (0d0,0d0)
      enddo

      nff(0)=c0g3S11(N)
      
      return
      end

      subroutine g3S11c1(N,Q,nff)
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
      
**
*     Internal Variables
*
      integer iff
      double complex N,Q
      double complex c1g3S11
      external c1g3S11
**
**
*     Output Variables
*
      double complex nff(-6:6)
      do iff=-6,6
         nff(iff) = (0d0,0d0)
      enddo

      nff(0)=c1g3S11(N)
      
      return
      end

      subroutine FENGg3S11(N,Q,nff)
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
      
**
*     Internal Variables
*
      integer iff
      double complex N,Q
      double complex c1g3S11,c0g3S11
      external c1g3S11, c0g3S11
**
**
*     Output Variables
*
      double complex nff(-6:6)
      do iff=-6,6
         nff(iff) = (0d0,0d0)
      enddo

      nff(0)=c0g3S11(N)*log(Q/1.5d0) + c1g3S11(N)
      
      return
      end
      
*************************************************************
*   FULL sigma up to v^2
************************************************************      
      subroutine FF3S11(N,Q,nff)
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
      
**
*     Internal Variables
*
      integer iff
      double complex N,Q
      double complex gFENG,gLO,g
      double complex CAS2,CAS3,CAS2v2,c
      
      double precision mc,v2,beta0
      double precision alphasG,alphasCLO,alphasCNLO,muR0G,muR0C
      
      double complex c0g3S11,c1g3S11,G0G3S11,G2G3S11
      double complex F,Fc3S11,POLYc3S11,Gc3S11,PQQc3S11
      external c0g3S11,c1g3S11,G0G3S11,G2G3S11
      external Fc3S11,POLYc3S11,Gc3S11,PQQc3S11
**
**
*     Output Variables
*
      double complex nff(-6:6)
      
      
* set scales and alphas
 
      beta0=(11d0*CA - 4d0*TR*3d0)/3d0
      
      mc = 1.5d0
      muR0G=2d0*mc
      muR0C=2d0*mc
      
      v2=0.5d0
      
      alphasG=0.31494871292
      alphasCLO=0.31494871292
      alphasCNLO=0.251822
      
* gluon induced part

      gFENG = alphasG**3*(c0g3S11(N)*log(Q/mc) + c1g3S11(N))
      gLO = alphasG**3*(G0G3S11(N) + v2*G2G3S11(N))
      g = (gLO + gFENG)/mc**3

* charm induced part

      F = Fc3S11(N)
      cAS2 = alphasCNLO**2*F
      cAS3 = alphasCNLO**3*( 
     *           +F/(2d0*PI)*beta0*log(muR0C**2/(2d0*mc)**2)
     *           +PQQc3S11(N)*F/(4d0*PI)*log(Q**2/(3d0*mc)**2) 
     *           +POLYc3S11(N)  
     *           )
      CAS2v2 =Gc3S11(N)*v2*alphasCLO**2
      c = (CAS2+CAS3+CAS2v2)/mc**3 

      do iff=-6,6
         nff(iff) = (0d0,0d0)
      enddo

      nff(0)=g
      nff(4)=c
      nff(-4)=c
      
      return
      end
*************************************************************
*   FULL sigma up to v4
************************************************************      
      subroutine FF3S11v4(N,Q,nff)
*
      implicit none
*
      include "../commons/alphas.h"
      include "../commons/consts.h"
      include "../commons/colfact.h"
      
**
*     Internal Variables
*
      integer iff
      double complex N,Q,muL
      double complex gFENG,gLO,g
      double complex CAS2,CAS3,CAS2v2,c
      
      double precision mc,v2,beta0
      double precision alphasG,alphasCLO,alphasCNLO,muR0G,muR0C
      
      double complex c0g3S11,c1g3S11,G0G3S11,G2G3S11
      double complex F,Fc3S11,POLYc3S11,Gc3S11,PQQc3S11
      double complex G4G3S11,G4G3S11DIST,gv4
      external c0g3S11,c1g3S11,G0G3S11,G2G3S11
      external Fc3S11,POLYc3S11,Gc3S11,PQQc3S11
      external G4G3S11,G4G3S11DIST
**
**
*     Output Variables
*
      double complex nff(-6:6)
      
      
* set scales and alphas
 
      beta0=(11d0*CA - 4d0*TR*3d0)/3d0
      
      mc = 1.5d0
      muL = (1.5d0,0d0)
      muR0G=2d0*mc
      muR0C=2d0*mc
      
      v2=0.5d0
      
      alphasG=0.31494871292
      alphasCLO=0.31494871292
      alphasCNLO=0.251822
      
* gluon induced part

      gFENG = alphasG**3*(c0g3S11(N)*log(Q/mc) + c1g3S11(N))
      gLO = alphasG**3*(G0G3S11(N) + v2*G2G3S11(N))
      gv4 =  alphasG**3*v2**2*(G4G3S11(N)+G4G3S11DIST(N,muL))
      g = (gLO + gFENG + gv4)/mc**3

* charm induced part

      F = Fc3S11(N)
      cAS2 = alphasCNLO**2*F
      cAS3 = alphasCNLO**3*( 
     *           +F/(2d0*PI)*beta0*log(muR0C**2/(2d0*mc)**2)
     *           +PQQc3S11(N)*F/(4d0*PI)*log(Q**2/(3d0*mc)**2) 
     *           +POLYc3S11(N)  
     *           )
      CAS2v2 =Gc3S11(N)*v2*alphasCLO**2
      c = (CAS2+CAS3+CAS2v2)/mc**3 

      do iff=-6,6
         nff(iff) = (0d0,0d0)
      enddo

      nff(0)=g
      nff(4)=c
      nff(-4)=c
      
      return
      end
