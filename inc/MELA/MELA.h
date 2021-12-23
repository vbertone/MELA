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
  /// Alpha/4*Pi at Q2
  double aQED4pi(double q2);
  /// Alpha at Q2  
  double aQED(double q2);
  /// Alpha at the electron mass  
  double aQEDinit();
  /// PDFs multiplied by x
  std::map<int, double> xDistributions(double x, double Q);
  /// PDFs multiplied by x (evolution basis)
  std::map<int, double> xDistributionsEv(double x, double Q);

  
  /*
   * SETTERS AND GETTERS
   */ 
  /// Perturbative order
  void SetPerturbativeOrder(int iptin);
  int GetPerturbativeOrder();

  /// Flavour scheme
  void SetFlavourScheme(std::string const& fnsin);
  void SetFlavourSchemeInt(int fnsin);
  int GetFlavourSchemeInt();

  void SetRenormalisationScheme(std::string const& rdin);
  void SetRenormalisationSchemeInt(int rnsinint);
  int GetRenormalisationSchemeInt();

  /// Factorisation scheme
  void SetFactorisationScheme(std::string const& fsin);
  void SetFactorisationSchemeInt(int fsin);
  int GetFactorisationSchemeInt();

  /// Maximum number of fermions in the VFNS
  void SetNFmax(int);
  int GetNFmax();

  /// Number of fermions in the FFNS
  void SetNFFN(int);
  int GetNFFN();

  /* /// Enable quarks in the evolution */
  /* void SetEnableQuarks(int); */
  /* bool GetEnableQuarks(); */

  /* void SetEvBasis(int); */
  /* bool GetEvBasis(); */
  
  /// Reference values for alpha
  void SetAlpha(double ain, double Qin);
  double GetAlphaRef();
  double GetAlphaQref();
  
  /* /// Analogous quantities but only for alpha evolution */
  /* void SetEnableQuarksalpha(int); */
  /* bool GetEnableQuarksalpha(); */
  
  void SetNFFNalpha(int);
  int GetNFFNalpha();
  void SetNFmaxalpha(int);
  int GetNFmaxalpha();
  
  /// Fermion thresholds
  void SetThresholds(double me, double mu, double md, double ms, double mm,
		     double mc, double mt, double mb, double mtp);
  void SetThresholds(std::vector<double> thresholds);
  std::vector<double> GetThresholds();
  std::vector<double> GetThresholds2();  

  double GetMZ2();
  
  /// Used in the analytic expressions
  int GetRegionMU2(double mu2);  
  double GetC2(int region);  
  double GetC4(int region);
  double Getb0(int region);
  double Getb1(int region);

  double GetDk(int region);

  /// Technical parameters
  void SetNint(int);
  int GetNint();
  void SetNexp(int);
  int GetNexp();
  void SetNstep(int);
  int GetNstep();
  void SetMinvmel(int);
  int GetMinvmel();
  void SetRinvmel(int);
  int GetRinvmel();

  void SetActiveFlavours(int nlmax, int numax, int ndmax);
  void SetActiveFlavoursAlpha(int nlmaxaem, int numaxaem, int ndmaxaem);  
}
