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
      character*100 card
*
      card = "Reference.ini"
*
*     Read parameters of the evolution from the card
*
      call ReadParameters(card)
*
      call LHAPDFgrid(2.8d0, "DeltaGluonFF")
*
      end
