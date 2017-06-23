noinst_LTLIBRARIES = libMELA.la

INIT = initialization/initBeta.f initialization/initCouplings.f \
	initialization/InitializeEvolution.f initialization/ReadParameters.f

EVOL = evolution/alphas.f evolution/andim_lo_tl.f evolution/andim_nlo_tl.f evolution/andim_nnlo_tl.f \
	evolution/crossing.f evolution/matchcoef2ns.f evolution/NDistributions.f evolution/toyLHPDFs.f \
	evolution/xDistributions.f evolution/andim_lo.f evolution/andim_nlo.f evolution/andim_nnlo.f \
	evolution/BasisTransformations.f evolution/evolfactn.f evolution/matchcoef2sg.f \
	evolution/path_ordering.f evolution/umatrix.f evolution/zfunc.f evolution/andim_lo_pol.f \
        evolution/andim_nlo_pol.f evolution/externalsetmela.cc

DIS = DIS/adlersr_coeff.f DIS/MassiveCoefficientFunctions.f DIS/NStructureFunctions.f \
	DIS/ZeroMassCoefficientFunctions.f DIS/ComputeChargesDIS.f DIS/NCoefficientFunctions.f \
	DIS/xStructureFunctions.f

CORE = core/betac.f core/dgauss.f core/gammal.f core/hypgeo.f core/mmult.f core/psi.f core/ddilog.f

EXTRA_DIST = commons

libMELA_la_SOURCES = $(INIT) $(EVOL) $(DIS) $(CORE)
