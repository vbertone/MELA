#include "MELA/MELA.h"

#include <vector>

/// Namespace containing all the MELA wrapper functions.
namespace MELA {
  extern"C"
  {
    void initializeevolution_(void);
    void readparameters_(char* fcard);
    void setdefaultparameters_(void);
    void enablequarks_(int*);
    void setalpha_(double* ain, double* qin);
    void setflavourscheme_(char* fnsin);
    void setflavourschemeint_(int* fnsin);
    void setnffn_(int* nffnin);
    void setnfmax_(int* nfmaxin);
    void setperturbativeorder_(int* iptin);
    void setthresholds_(double* me, double* mu, double* md, double* ms, double* mm, double* mc, double* mt, double* mb, double* mtp);
    double aqed_(double* q2);
    void xdistributions_(double* x, double* q, double* xfph);
  }

  /// Initialize the library
  void InitializeEvolution(void)
  {
    initializeevolution_();
  };

  /// Read parameters from card
  void ReadParameters(std::string const& fcard)
  {
    std::vector<char> cstr(fcard.c_str(), fcard.c_str() + fcard.size() + 1);
    readparameters_(cstr.data());
  };

  /// Set default parameters 
  void SetDefaultParameters(void) { setdefaultparameters_(); };

  /// Enable quarks in the evolution
  void EnableQuarks(int eq)
  {
    enablequarks_(&eq);
  };

  /// Set reference values for alpha
  void SetAlpha(double ain, double Qin)
  {
    setalpha_(&ain, &Qin);
  };

  /// Set the flavour scheme to be used in the evolution as a string
  void SetFlavourScheme(std::string const& fnsin)
  {
    std::vector<char> cstr(fnsin.c_str(), fnsin.c_str() + fnsin.size() + 1);
    setflavourscheme_(cstr.data());
  };

  /// Set the flavour scheme to be used in the evolution as a int
  void SetFlavourSchemeInt(int fnsin)
  {
    setflavourschemeint_(&fnsin);
  };

  /// Set the number of fermions in the FFNS
  void SetNFFN(int NFFNin)
  {
    setnffn_(&NFFNin);
  };

  /// Set maximum number of fermions in the VFNS
  void SetNFmax(int NFmaxin)
  {
    setnfmax_(&NFmaxin);
  };

  /// Set the perturbative order
  void SetPerturbativeOrder(int iptin)
  {
    setperturbativeorder_(&iptin);
  };

  /// Set fermion thresholds
  void SetThresholds(double me, double mu, double md, double ms, double mm, double mc, double mt, double mb, double mtp)
  {
    setthresholds_(&me, &mu, &md, &ms, &mm, &mc, &mt, &mb, &mtp);
  };

  /// Coupling
  double aQED(double q2)
  {
    return aqed_(&q2);
  };

  /// PDFs
  std::map<int, double> xDistributions(double x, double Q)
  {
    double* xfph = new double[19];
    xdistributions_(&x, &Q, xfph);
    std::map<int, double> xfout;
    for (int i = -9; i < 10; i++)
      xfout.insert({i, xfph[i+9]});
    return xfout;
  };
}
