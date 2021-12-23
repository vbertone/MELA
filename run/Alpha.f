************************************************************************
*
*     Evolution.f:
*
*     Driver for the generation of the evolution tables.
*
************************************************************************
      program Evolution
*
      implicit none
*
      integer i
      character*4 ns
      double precision aQED
      double precision Q
      double precision thrs(9)
      double precision a0,Qref,L
      double precision aana
      double precision beta0(9), beta1(9)
      data thrs / 0.000510998928d0,0.00216d0,0.00467d0,
     1     0.093d0,0.10566d0,1.27d0,1.77686d0,4.18d0,172.76d0 /
*
      ns = "VFNS"
      a0 = 0.0078152650036d0/12.566370614359173d0
      Qref = 91.187600000000d0
*
*     Start by setting default values. This will always need to be
*     called at first to set all parameters and avoid
*     misbehaviours. Parameters can be adjusted after wards using the
*     setting functions.
*
      call SetDefaultParameters
*
*     Set custom parameters
*
      call SetPerturbativeOrder(1)
      call SetFlavourScheme(ns)
      call SetNFmax(8)
      call SetNFFN(8)
      call SetAlpha(a0*12.566370614359173d0, Qref)
      call SetThresholds(thrs(1), thrs(2), thrs(3), thrs(4), thrs(5),
     1     thrs(6), thrs(7), thrs(8), thrs(9))
*
*     Initialization of evolution parameters
*
      call InitializeEvolution
*
*     Final scale (initial scale assumed to be the electron mass)
*
      Q = thrs(1)
      thrs(9) = Qref
*
      write(6,*) "Q = ", Q, " GeV"
      write(6,*) "Numerical 1/alpha(Q) = ",1d0
     1     / ( aQED(Q**2) * 12.566370614359173d0 )
*
      aana = 0
      call BetaCoeffs(beta0,beta1)
      if (ns.eq."FFNS") then
         L = 2 * dlog(Qref / Q)
         aana = a0 + a0**2 * beta0(8) * L
     1        + a0**3 * ( beta0(8)**2 * L**2 + beta1(8) * L )
      elseif (ns.eq."VFNS") then
         do i = 8, 1, -1
            L = 2 * dlog(thrs(i+1) / thrs(i))
            aana = a0 + a0**2 * beta0(i) * L
     1           + a0**3 * ( beta0(i)**2 * L**2 + beta1(i) * L )
            a0 = aana
         enddo
      endif
*
      write(6,*) "Analytic 1/alpha(Q) = ",1d0
     1     / ( aana * 12.566370614359173d0 )
*
      end
*
************************************************************************
      subroutine BetaCoeffs(beta0,beta1)
*
      implicit none
**
*     Internal Variables
*
      integer i,nc,nl(0:9),nu(0:9),nd(0:9)
      double precision el2,el4,eu2,eu4,ed2,ed4
      double precision beta0(9), beta1(9)
*
*     Initialise numver of active lepton, up-type, and down-type
*     flavours.
*
      nl(0) = 0
      nl(1) = 1
      nl(2) = 1
      nl(3) = 1
      nl(4) = 1
      nl(5) = 2
      nl(6) = 2
      nl(7) = 3
      nl(8) = 3
      nl(9) = 3
*
      nu(0) = 0
      nu(1) = 0
      nu(2) = 1
      nu(3) = 1
      nu(4) = 1
      nu(5) = 1
      nu(6) = 2
      nu(7) = 2
      nu(8) = 2
      nu(9) = 3
*
      nd(0) = 0
      nd(1) = 0
      nd(2) = 0
      nd(3) = 1
      nd(4) = 2
      nd(5) = 2
      nd(6) = 2
      nd(7) = 2
      nd(8) = 3
      nd(9) = 3
*
*     Number of colours
*
      nc = 3
*
*     electric charges
*
      el2 = 1d0
      eu2 = 4d0 / 9d0
      ed2 = 1d0 / 9d0
*
      el4 = el2**2
      eu4 = eu2**2
      ed4 = ed2**2
*
*     Beta function and mass anomalous dimension coefficients
*
      do i = 1, 9
         beta0(i) = - 4d0 * ( nl(i)
     1        + nc * ( ed2 * nd(i) + eu2 * nu(i) ) ) / 3d0
         beta1(i) = - 4d0 * ( nl(i)
     1        + nc * ( ed4 * nd(i) + eu4 * nu(i) ) )
      enddo
*
      return
      end
