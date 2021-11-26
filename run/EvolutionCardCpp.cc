#include "MELA/MELA.h"

#include <iostream>
#include <vector>
#include <math.h>
#include <cstring>

int main(int argc, char* argv[])
{
  // Check that the input is correct otherwise stop the code
  if (argc < 2 || strcmp(argv[1], "--help") == 0)
    {
      std::cout << "\nInvalid Parameters:" << std::endl;
      std::cout << "Syntax: ./EvolutionCardCpp <input card>\n" << std::endl;
      exit(-10);
    }
  // Start by setting default values. This will always need to be
  // called at first to set all parameters and avoid
  // misbehaviours. Parameters can be adjusted after wards using the
  // setting functions.
  MELA::SetDefaultParameters();

  // Read parameters of the evolution from the card
  MELA::ReadParameters(std::string(argv[1]));

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