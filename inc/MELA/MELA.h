#pragma once

#include <string>
#include <map>
#include <vector>

/// Namespace containing all the MELA wrapper functions.
namespace MELA {
  /*
   * Main Methods
   */ 
  /// Initialize the library
  void InitializeEvolution(void);

  /// Read parameters from card
  void ReadParameters(std::string const& fcard);

  /// Set default parameters 
  void SetDefaultParameters(void);

  /// Enable quarks in the evolution
  void EnableQuarks(int);
  void EnableQuarksalpha(int);

  /// Set reference values for alpha
  void SetAlpha(double ain, double Qin);

  /// Set the flavour scheme to be used in the evolution as a string
  void SetFlavourScheme(std::string const& fnsin);

  /// Set the flavour scheme to be used in the evolution as a int
  void SetFlavourSchemeInt(int fnsin);

  /// Set the number of fermions in the FFNS
  void SetNFFN(int);

  /// Set maximum number of fermions in the VFNS
  void SetNFmax(int);
  int  GetNFmax();
  
  void SetNFFNalpha(int);
  void SetNFmaxalpha(int);
  
  /// Set the perturbative order
  void SetPerturbativeOrder(int iptin);

  /// Set fermion thresholds
  void SetThresholds(double me, double mu, double md, double ms, double mm, double mc, double mt, double mb, double mtp);

  /// Coupling
  double aQED4pi(double q2);  
  double aQED(double q2);
  double aQEDinit();
  
  /// PDFs
  std::map<int, double> xDistributions(double x, double Q);

  std::vector<double> GetThresholds();
  int GetRegionMU2(double mu2);  
  double GetC2(int region);
  double GetC4(int region);
  double Getb0(int region);
  double Getb1(int region);
}
