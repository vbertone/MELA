// ePDF_private
#include <ePDF/epdf.h>

#include <iostream>
#include <vector>
#include <iomanip>

extern "C" void readparameters_(char* card);
extern "C" void initializeevolution_();
extern "C" void xdistributionsreal_(double*, double*, double*, double*);

int main() {
  double Q0 = 0.000510998928;
  double Q  = 100;

  // MELA nf = 1
  readparameters_("ReferenceNf3.ini");
  initializeevolution_();
  double* xfm = new double[7];

  // ePDF
  ePDF::epdf pdfNLL;
  pdfNLL.SetBase("Evolution");

  std::vector<double> xv(1000);

  const int nx = 500;
  double xmin = 0.1;
  double xmax = 0.9;
  double xstp = exp( log( xmax / xmin ) / ( nx - 1 ) );
  double x = xmin;
  for (int ix = 0; ix < nx; ix++)
    {
      xv[ix] = x;
      x *= xstp;
    }
  xmin = 0.9;
  xmax = 1 - 1e-5;
  xstp = ( xmax - xmin ) / ( nx - 1 );
  x = xmin;
  for (int ix = 500; ix < nx + 500; ix++)
    {
      xv[ix] = x;
      x += xstp;
    }

  //std::vector<double> xv{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999, 0.9999, 0.99999};

  std::cout << std::scientific;
  for (auto& x : xv)
    {
      xdistributionsreal_(&x, &Q0, &Q, xfm);
      const std::vector<double> xfe = pdfNLL.pdfxQ(x, Q);
      std::cout << std::setprecision(8) << x << "\t";
      std::cout << std::setprecision(8) << xfm[3] / x / xfe[1] << "\t"
		<< std::setprecision(8) << ( xfm[4] + xfm[2] ) / x / xfe[0] << "\t"
		<< std::setprecision(8) << ( xfm[4] - xfm[2] ) / x / xfe[2] << "\t"
		//<< std::setprecision(4) << ( xfm[5] + xfm[1] ) / x << "\t"
		//<< std::setprecision(4) << ( xfm[5] - xfm[1] ) / x << "\t"
		//<< std::setprecision(4) << ( xfm[6] + xfm[0] ) / x << "\t"
		//<< std::setprecision(4) << ( xfm[6] - xfm[0] ) / x << "\t"
		//<< std::setprecision(4) << xfe[1] << "\t"
		//<< std::setprecision(4) << xfe[0] << "\t"
		//<< std::setprecision(4) << xfe[2] << "\t"
		<< std::endl;
    }
  std::cout << "\n";
/*
  // Construct chi2 object with the irep-th replica
  MontBlanc::AnalyticChiSquare chi2{DSVect, new MontBlanc::LHAPDFparameterisation{FFsp, pg}};
  std::cout << "ID        Q [GeV]            x               z          MontBlanc/MELA" << std::endl;
  for (int iexp = 0; iexp < (int) DSVect.size(); iexp++)
    {
      // Get binning and kinematics
      const std::vector<NangaParbat::DataHandler::Binning> bins = DSVect[iexp].first->GetBinning();
      const NangaParbat::DataHandler::Kinematics kins = DSVect[iexp].first->GetKinematics();

      // Get predictions 
      const std::vector<double> preds = DSVect[iexp].second->GetPredictions([](double const &, double const &, double const &) -> double { return 0; });

      // Print results
      for (int i = 0; i < (int) bins.size(); i++)
	{
	  const NangaParbat::DataHandler::Binning b = bins[i];

	  // Inelasticity (Q^2 = sxy ==> y = Q^2 / (sx))
	  const double y  = pow(b.Qav / kins.Vs, 2) / b.xav;
	  const double y2 = y * y;
	  const double yp = 1 + pow(1 - y, 2);
	  double x = b.xav;
	  double z = b.zav;
	  double Q = b.Qav;

	  xstructurefunctionsreal_(&x, &Q0, &Q, isf);
	  sidisxstructurefunctions_(&x, &z, &Q0, &Q, sisf);
	  const double mela = ( yp * sisf[0] - y2 * sisf[1] ) / ( yp * isf[0] - y2 * isf[1] ) / z;

	  std::cout << std::scientific << i << "\t" << b.Qav << "\t" << b.xav << "\t" << b.zav << "\t" << preds[i] / mela << std::endl;
	}
    }
*/
  return 0;
}
