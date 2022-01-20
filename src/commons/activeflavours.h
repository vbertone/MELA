*     -*-fortran-*-
*
*     Active flavours of leptons, up-type and down-type quarks
*
      integer nlmax, numax, ndmax
      integer nlmaxaem, numaxaem, ndmaxaem
      integer nl(0:9), nu(0:9), nd(0:9)
      integer waem
*
      common / activeFlavoursMELA / nl,nu,nd
      common / nmaxMELA / nlmax,numax,ndmax,nlmaxaem,numaxaem,ndmaxaem
      common / waemMELA / waem
