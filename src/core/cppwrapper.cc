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
    double aqed_(double* q2);
    void geta0_(double* a0);    
    void xdistributions_(double* x, double* q, double* xfph);
    void xdistributionsev_(double* x, double* q, double* xfev);    

    void setperturbativeorder_(int* iptin);
    void getperturbativeorder_(int* iptout);

    void setflavourscheme_(char* fnsin);
    void setflavourschemeint_(int* fnsin);
    void getflavourschemeint_(int* fscheme);

    void setfactorisationscheme_(char* fsin);
    void setfactorisationschemeint_(int* fsinint);
    void getfactorisationschemeint_(int* fsintout);

    void setrenormalisationscheme_(char* rsin);
    void setrenormalisationschemeint_(int* rsinint);
    void getrenormalisationschemeint_(int* fsintout);

    void setnfmax_(int* nfmaxin);
    void getnfmax_(int* nfmax);
    void setnffn_(int* nffnin);
    void getnffn_(int* nffn);
    
    //void enablequarks_(int*);
    //void getenablequarks_(int*);
    
    void setalpha_(double* ain, double* qin);
    void getalpharef_(double*);
    void getalphaqref_(double*);

    void setactiveflav_(int* nlmax, int* numax, int* ndmax);
    void getactiveflav_(int* nlmax, int* numax, int* ndmax);
    void setactiveflavaem_(int* nlmax, int* numax, int* ndmax);
    void getactiveflavaem_(int* nlmax, int* numax, int* ndmax);
    
    //void enablequarksalpha_(int*);
    //void getenablequarksalpha_(int*);
    
    void setnffnalpha_(int*);
    void getnffnalpha_(int*);
    void setnfmaxalpha_(int*);
    void getnfmaxalpha_(int*);

    void setthresholds_(double* me, double* mu, double* md, double* ms, double* mm,
			double* mc, double* mt, double* mb, double* mtp);
    void getthresholds2_(double* q2thrs);
    void getmz2_(double*);    

    void getc2_(double* C2);
    void getc4_(double* C4);
    void getb0_(double* b0);
    void getb1_(double* b1);

    void getdk_(int* region, double* dk);
    
    void setnint_(int*);
    void getnint_(int*);
    void setnexp_(int*);
    void getnexp_(int*);
    void setnstep_(int*);
    void getnstep_(int*);
    void setminvmel_(int*);
    void getminvmel_(int*);
    void setrinvmel_(int*);
    void getrinvmel_(int*);
  }

  void InitializeEvolution(void)
  {
    initializeevolution_();
  };

  void ReadParameters(std::string const& fcard)
  {
    std::vector<char> cstr(fcard.c_str(), fcard.c_str() + fcard.size() + 1);
    readparameters_(cstr.data());
  };

  void SetDefaultParameters(void) { setdefaultparameters_(); };

  double aQED4pi(double q2)
  {
    return aqed_(&q2);
  };
  double aQED(double q2)
  {
    return aqed_(&q2)*4.0*M_PI;
  };
  
  double aQEDinit()
  {
    double a0;
    geta0_(&a0);
    return a0*4.0*M_PI;
  };

  /// PDFs (from -9 to 9)
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

  /// PDFs (from 1 to 19)
  std::map<int, double> xDistributionsEv(double x, double Q)
  {
    double* xfev = new double[19];
    xdistributionsev_(&x, &Q, xfev);
    std::map<int, double> xfout;
    for (int i = 1; i < 20; i++)
      xfout.insert({i, xfev[i-1]});
    delete[] xfev;
    return xfout;
  };
  
 
  void SetPerturbativeOrder(int iptin)
  {
    setperturbativeorder_(&iptin);
  };

  int GetPerturbativeOrder()
  {
    int pto;
    getperturbativeorder_(&pto);
    return pto;
  }

  void SetFlavourScheme(std::string const& fnsin)
  {
    std::vector<char> cstr(fnsin.c_str(), fnsin.c_str() + fnsin.size() + 1);
    setflavourscheme_(cstr.data());
  };

  void SetFlavourSchemeInt(int fnsin)
  {
    setflavourschemeint_(&fnsin);
  };

  int GetFlavourSchemeInt()
  {
    int fscheme;
    getflavourschemeint_(&fscheme);
    return fscheme;
  }

  void SetFactorisationScheme(std::string const& fsin)
  {
    std::vector<char> cstr(fsin.c_str(), fsin.c_str() + fsin.size() + 1);
    setfactorisationscheme_(cstr.data());
  };

  void SetFactorisationSchemeInt(int fnsinint)
  {
    setfactorisationschemeint_(&fnsinint);
  };

  int GetFactorisationSchemeInt()
  {
    int fscheme;
    getfactorisationschemeint_(&fscheme);
    return fscheme;
  }

  void SetRenormalisationScheme(std::string const& rdin)
  {
    std::vector<char> cstr(rdin.c_str(), rdin.c_str() + rdin.size() + 1);
    setrenormalisationscheme_(cstr.data());
  };

  void SetRenormalisationSchemeInt(int rnsinint)
  {
    setrenormalisationschemeint_(&rnsinint);
  };

  int GetRenormalisationSchemeInt()
  {
    int rdcheme;
    getrenormalisationschemeint_(&rdcheme);
    return rdcheme;
  }

  void SetNFmax(int NFmaxin)
  {
    setnfmax_(&NFmaxin);
  };
    
  int GetNFmax()
  {
    int nfmax;
    getnfmax_(&nfmax);
    return nfmax;
  }

  void SetNFFN(int NFFNin)
  {
    setnffn_(&NFFNin);
  };
  
  int GetNFFN()
  {
    int nffn;
    getnffn_(&nffn);
    return nffn;
  }
  
  // void SetEnableQuarks(int eq)
  // {
  //   enablequarks_(&eq);
  // };

  // bool GetEnableQuarks()
  // {
  //   int quarks;
  //   getenablequarks_(&quarks);
  //   return quarks;
  // }

  
  /// Set reference values for alpha
  void SetAlpha(double ain, double Qin)
  {
    setalpha_(&ain, &Qin);
  };

  double GetAlphaRef()
  {
    double res;
    getalpharef_(&res);
    return res;
  }

  double GetAlphaQref()
  {
    double res;
    getalphaqref_(&res);
    return res;
  }

  // void SetEnableQuarksalpha(int eq)
  // {
  //   enablequarksalpha_(&eq);
  // };
  
  // bool GetEnableQuarksalpha()
  // {
  //   int quarks;
  //   getenablequarksalpha_(&quarks);
  //   return quarks;
  // }

  void SetNFFNalpha(int NFFNalphain)
  {
    setnffnalpha_(&NFFNalphain);
  };
  int GetNFFNalpha()
  {
    int nffnalpha;
    getnffnalpha_(&nffnalpha);
    return nffnalpha;
  }
  
  void SetNFmaxalpha(int NFmaxalphain)
  {
    setnfmaxalpha_(&NFmaxalphain);
  };
  int GetNFmaxalpha()
  {
    int nffnalpha;
    getnfmaxalpha_(&nffnalpha);
    return nffnalpha;
  }
  
    
  /// Set fermion thresholds
  void SetThresholds(double me, double mu, double md, double ms, double mm,
		     double mc, double mt, double mb, double mtp)
  {
    setthresholds_(&me, &mu, &md, &ms, &mm, &mc, &mt, &mb, &mtp);
  };

  void SetThresholds(std::vector<double> thrs)
  {
    double & me = thrs[0];
    double & mu = thrs[1];
    double & md = thrs[2];
    double & ms = thrs[3];
    double & mm = thrs[4];
    double & mc = thrs[5];
    double & mt = thrs[6];
    double & mb = thrs[7];
    double & mtp = thrs[8];
    
    setthresholds_(&me, &mu, &md, &ms, &mm, &mc, &mt, &mb, &mtp);
  };
  

  std::vector<double> GetThresholds()
  {
    double* q2thrsf = new double[9];
    getthresholds2_(q2thrsf);
    std::vector<double> thrs;
    for (int i = 0; i < 9; i++)
      thrs.push_back(sqrt(q2thrsf[i]));
    delete[] q2thrsf;
    return thrs;
  }  

  std::vector<double> GetThresholds2()
  {
    double* q2thrsf = new double[9];
    getthresholds2_(q2thrsf);
    std::vector<double> q2thrs;
    for (int i = 0; i < 9; i++)
      q2thrs.push_back(q2thrsf[i]);
    delete[] q2thrsf;
    return q2thrs;
  }  

  
  int GetRegionMU2(double mu2)
  {
    std::vector<double> q2thrs = GetThresholds2();
    int region;
    // Initialize region to nfmax-1
    getnfmax_(&region);
    region = region - 1;
    while (region >= 0) {
      if (mu2 > q2thrs[region]) {
	break;
      } else {
	region = region - 1;
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

  double GetDk(int region)
  {
    // C++ and Fortran indices are shifted by one
    region += 1;
    double dk;
    getdk_(&region, &dk);
    return dk;
  }
  
  void SetNint(int res)
  {
    setnint_(&res);
  }

  int GetNint()
  {
    int res;
    getnint_(&res);
    return res;
  }

  void SetNexp(int res)
  {
    setnexp_(&res);
  }

  int GetNexp()
  {
    int res;
    getnexp_(&res);
    return res;
  }

  void SetNstep(int res)
  {
    setnstep_(&res);
  }

  int GetNstep()
  {
    int res;
    getnstep_(&res);
    return res;
  }

  void SetMinvmel(int res)
  {
    setminvmel_(&res);
  }

  int GetMinvmel()
  {
    int res;
    getminvmel_(&res);
    return res;
  }

  void SetRinvmel(int res)
  {
    setrinvmel_(&res);
  }

  int GetRinvmel()
  {
    int res;
    getrinvmel_(&res);
    return res;
  }

  void SetActiveFlavours(int nlmax, int numax, int ndmax)
  {
    setactiveflav_(&nlmax, &numax, &ndmax);    
  }

  void SetActiveFlavoursAlpha(int nlmaxaem, int numaxaem, int ndmaxaem)
  {
    setactiveflavaem_(&nlmaxaem, &numaxaem, &ndmaxaem);    
  }

  double GetMZ2()
  {
    double mz2;
    getmz2_(&mz2);
    return mz2;
  }
  
}
