*     -*-fortran-*-
*
*     Active flavours of leptons, up-type and down-type quarks
*
      integer nlmax, numax, ndmax
      integer nlmaxaem, numaxaem, ndmaxaem
      integer nl(1:11), nu(1:11), nd(1:11)
      integer waem
*
      common / activeFlavoursMELA / nl,nu,nd
      common / nmaxMELA / nlmax,numax,ndmax,nlmaxaem,numaxaem,ndmaxaem
      common / waemMELA / waem
