#pragma once

#include <string>
#include <map>
#include <vector>

/// Namespace containing all the MELA wrapper functions.
namespace MELA {

  //////////////////////////////////////////////////////////////////////
  /// Main methods
  //////////////////////////////////////////////////////////////////////  
  /// Set default parameters 
  void SetDefaultParameters(void);  
  /// Initialize the library
  void InitializeEvolution(void);

  /// Alpha at Q2  
  double aQED(double q2);
  /// Alpha/(4*Pi) at Q2
  double aQED4pi(double q2);
  /// Alpha at the electron mass  
  double aQEDinit();
  
  /// PDFs multiplied by x
  std::map<int, double> xDistributions(double x, double Q);
  /// PDFs multiplied by x (evolution basis)
  std::map<int, double> xDistributionsEv(double x, double Q);

  //////////////////////////////////////////////////////////////////////
  /// Setters
  //////////////////////////////////////////////////////////////////////  
  void SetPerturbativeOrder(int);
  void SetPerturbativeOrderAlpha(int);
  void SetFlavourScheme(std::string const&);
  void SetFlavourSchemeInt(int);
  void SetFactorisationScheme(std::string const&);
  void SetFactorisationSchemeInt(int);  
  void SetRenormalisationScheme(std::string const&);
  void SetRenormalisationSchemeInt(int);
  void SetActiveFlavours(int nlmax, int numax, int ndmax);
  void SetActiveFlavoursAlpha(int nlmaxaem, int numaxaem, int ndmaxaem);
  void SetActiveFlavours(std::vector<int> active);
  void SetActiveFlavoursAlpha(std::vector<int> active);
  void SetWalpha(int);
  void SetAlpha(double ain, double Qin);
  void SetThresholds(double me, double mu, double md, double ms, double mm,
		     double mc, double mt, double mb,
		     double MW, double MZ, double mtp);
  void SetThresholds(std::vector<double> thresholds);
  void SetNint(int);
  void SetNminstep(int);  
  void SetNexp(int);
  void SetNstepAem(int);
  void SetMinvmel(int);
  void SetRinvmel(int);
  void SetAlphaXXSolInt(int);

  //////////////////////////////////////////////////////////////////////
  /// Getters
  //////////////////////////////////////////////////////////////////////  
  int GetPerturbativeOrder();
  int GetPerturbativeOrderAlpha();
  int GetFlavourSchemeInt();
  int GetFactorisationSchemeInt();
  int GetRenormalisationSchemeInt();
  std::vector<int> GetActiveFlavours();
  std::vector<int> GetActiveFlavoursAlpha();
  int GetWalpha();      
  double GetAlphaRef();
  double GetAlphaQref();  
  std::vector<double> GetThresholds();
  std::vector<double> GetThresholds2();    
  int GetNint();
  int GetNminstep();
  int GetNexp();
  int GetNstepAem();  
  int GetMinvmel();
  int GetRinvmel();

  //////////////////////////////////////////////////////////////////////
  /// Used in the analytic expressions
  //////////////////////////////////////////////////////////////////////  
  // Return the number of the region of mu2
  // i.e. 1 if me2 <= mu2 < mu2
  //      2 if mu2 <= mu2 < md2
  //      ...
  //      8 if mb2 <= mu2 < mw2
  //      9 if mw2 <= mu2 < mz2
  //     10 if mz2 <= mu2 < mtp2
  //     11 if mtp2 <= mu2
  int GetRegionMU2(double mu2);

  // Get each of the expressions in the relevant region
  double GetC2(int region);  
  double GetC4(int region);
  double Getb0(int region);
  double Getb1(int region);
  double GetDk(int region);

  double aQEDinit();
  double GetMZ2();
  double GetMW2();  
}
