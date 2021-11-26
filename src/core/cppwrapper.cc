#include "MELA/MELA.h"

#include <vector>
#include <math.h>

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
    void getc2_(double* C2);
    void getc4_(double* C4);
    void getb0_(double* b0);
    void getb1_(double* b1);
    void getthresholds_(double* q2thrs);
    void geta0_(double* a0);
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
  double aQED4pi(double q2)
  {
    return aqed_(&q2);
  };
  double aQED(double q2)
  {
    return aqed_(&q2)*4.0*M_PI;
  };
  
  /// Coupling at the initial scale
  double aQEDinit()
  {
    double a0;
    geta0_(&a0);
    return a0*4.0*M_PI;
  };
  
  /// PDFs
  std::map<int, double> xDistributions(double x, double Q)
  {
    double* xfph = new double[19];
    xdistributions_(&x, &Q, xfph);
    std::map<int, double> xfout;
    for (int i = -9; i < 10; i++)
      xfout.insert({i, xfph[i+9]});
    delete[] xfph;
    return xfout;
  };

  std::vector<double> GetThresholds()
  {
    double* q2thrsf = new double[9];
    getthresholds_(q2thrsf);
    std::vector<double> q2thrs;
    for (int i = 0; i < 8; i++)
      q2thrs.push_back(q2thrsf[i]);
    delete[] q2thrsf;
    return q2thrs;
  }

  int GetRegionMU2(double mu2)
  {
    std::vector<double> q2thrs = GetThresholds();
    int region = 8;
    while (region >= 0) {
      if (mu2 > q2thrs[region]) {
	break;
      } else {
	region -= 1;
      }
    }
    return region;
  }
  
  double GetC2(int region)
  {
    double* C2f = new double[8];
    getc2_(C2f);
    double res =  C2f[region];
    delete[] C2f;
    return res;
  }

  double GetC4(int region)
  {
    double* C4f = new double[8];
    getc4_(C4f);
    double res = C4f[region];
    delete[] C4f;
    return res;
  }

  double Getb0(int region)
  {
    double* b0f = new double[8];
    getb0_(b0f);
    double res = b0f[region];
    delete[] b0f;
    return res;
  }

  double Getb1(int region)
  {
    double* b1f = new double[8];
    getb1_(b1f);
    double res = b1f[region];
    delete[] b1f;
    return res;
  }
}
