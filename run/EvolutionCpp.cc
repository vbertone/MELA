#include "MELA/MELA.h"

#include <iostream>
#include <vector>
#include <math.h>

int main()
{
  // Start by setting default values. This will always need to be
  // called at first to set all parameters and avoid
  // misbehaviours. Parameters can be adjusted after wards using the
  // setting functions.
  MELA::SetDefaultParameters();

  // Set custom parameters
  MELA::SetPerturbativeOrder(1);
  MELA::SetFlavourScheme("VFNS");
  //MELA::SetFactorisationScheme("DELTA")
  MELA::SetNFmax(8);
  MELA::SetNFFN(8);
  MELA::SetAlpha(0.007815265003645828, 91.1876);
  MELA::EnableQuarks(true);
  MELA::SetThresholds(0.000510998928,0.00216,0.00467,0.093,0.10566,1.27,1.77686,4.18,172.76);

  // Initialization of evolution parameters
  MELA::InitializeEvolution();

  // Final scale (initial scale assumed to be the electron mass)
  const double Q = 100;

  // Test values of x
  const std::vector<double> xlha{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999};

  // Print results
  std::cout << std::scientific;
  std::cout << "alpha(Q) = " << MELA::aQED(pow(Q, 2)) * 12.566370614359173 << std::endl;
  std::cout << "  x     " << "\t" <<
          "   e- + e+  " << "\t" <<
          "    photon  " << "\t" <<
          "   e- - e+  " << "\t" <<
          "  mu- + mu+ " << "\t" <<
          "  mu- - mu+ " << "\t" <<
          " tau- + tau+" << "\t" <<
          " tau- - tau+" << std::endl;
  for (double x : xlha)
    {
      const std::map<int, double> xf = MELA::xDistributions(x, Q);
      std::cout.precision(2);
      std::cout << x << "\t";
      std::cout.precision(4);
      std::cout << ( xf.at(1) + xf.at(-1) ) / x << "\t"
		<< xf.at(0) / x << "\t"
		<< ( xf.at(1) - xf.at(-1) ) / x << "\t"
		<< ( xf.at(2) + xf.at(-2) ) / x << "\t"
		<< ( xf.at(2) - xf.at(-2) ) / x << "\t"
		<< ( xf.at(3) + xf.at(-3) ) / x << "\t"
		<< ( xf.at(3) - xf.at(-3) ) / x << "\t"
		<< std::endl;
    }

  return 0;
}
