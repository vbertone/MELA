*     -*-fortran-*-
*
*     Parameters of the LHAPDF evolution grid     
*
      integer nxmax,nq2max
      integer nxLHA,nxmLHA,nq2LHA
      double precision q2minLHA,q2maxLHA
      double precision xminLHA,xmLHA,xmaxLHA
      double precision Lambda2
*
      parameter(Lambda2 = 0.0625d0)
      parameter(nxmax   = 300)
      parameter(nq2max  = 200)
*
      common / LHgridParamMELA / xminLHA,xmLHA,xmaxLHA,q2minLHA,
     1                           q2maxLHA,nxLHA,nxmLHA,nq2LHA
