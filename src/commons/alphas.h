*     -*-fortran-*-
*
*     Values of alphas at the reference scale, final scale and at the heavy quark mass thresholds
*
      double complex asref,Q2ref
      double complex as0
      double complex asq
      double complex asc,ascm,asb,asbm,ast,astm
      double complex asEvolIni(100),asEvolFin(100)
      double complex asDIS
      double complex bq(6),dq(6)
*
      common / alphasrefMELA / asref,Q2ref
      common / alphasiniMELA / as0
      common / alphasfinMELA / asq
      common / alphasthrMELA / asc,ascm,asb,asbm,ast,astm
      common / alphasEvoMELA / asEvolIni,asEvolFin
      common / alphasDISMELA / asDIS
      common / alphasCouMELA / bq,dq
