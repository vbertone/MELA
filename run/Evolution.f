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
      integer ilha,i
      double precision aQED
      double precision xlha(12)
      double precision Q
      double precision xf(-9:9)
      double complex fn(19),fnph(-9:9),sr1
      data xlha / 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0,
     1            0.8d0, 0.9d0, 0.95d0, 0.99d0, 0.999d0/
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
c     call SetPerturbativeOrder(0)
      call SetFlavourScheme("VFNS")
      !call SetFlavourScheme("FFNS")
      call SetAlpha(0.007815265003645828d0, 91.1876d0)
      call SetFactorisationScheme("MSBAR")
c     call SetRenormalisationScheme("FIXED")
      call SetRenormalisationScheme("ALPMZ")
      !call SetAlpha(0.007755543505006d0, 91.1876d0)
      !call SetRenormalisationScheme("ALPMZ")
*
*     Initialization of evolution parameters
*
      call InitializeEvolution
*
*     Final scale (initial scale assumed to be the electron mass)
*
      Q = 50d0
*
      write(6,*) "alpha(Q) = ", aQED(Q**2) * 12.566370614359173d0
      write(6,*)
     1     "  x    ",
     2     "   e- + e+  ",
     3     "    photon  ",
     4     "   e- - e+  ",
     5     "  mu- + mu+ ",
     6     "  mu- - mu+ ",
     7     " tau- + tau+",
     8     " tau- - tau+"
      do ilha=1,12
*
*     Perform the N-space evolution
*
         call xDistributions(xlha(ilha),Q,xf)
*
         write(6,'(es8.2,7(es12.4))')
     1         xlha(ilha),
     2         ( xf(1) + xf(-1) ) / xlha(ilha),
     3         xf(0) / xlha(ilha),
     4         ( xf(1) - xf(-1) ) / xlha(ilha),
     5         ( xf(2) + xf(-2) ) / xlha(ilha),
     6         ( xf(2) - xf(-2) ) / xlha(ilha),
     7         ( xf(3) + xf(-3) ) / xlha(ilha),
     8         ( xf(3) - xf(-3) ) / xlha(ilha)
      enddo
      write(*,*) "  "
*
      call NDistributionsEv(dcmplx(2d0, 0d0),Q,fn)
      call evln2lhac(fn,fnph)
      sr1 = 0d0
      do i = -9, 9
         sr1 = sr1 + fnph(i)
      enddo
      write(6, *) "Momentum sum rule = ", sr1
*
      call NDistributionsEv(dcmplx(1d0, 0.00001d0),Q,fn)
      call evln2lhac(fn,fnph)
      do i = 1, 9
         write(6, *) "Valence sum rule (", i ,") = ",fnph(i) - fnph(-i)
      enddo
*
      end
