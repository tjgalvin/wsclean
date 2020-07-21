#include "../nlplfitter.h"

#include "../model/powerlawsed.h"

#include <aocommon/uvector.h>

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(nlplfitter)

BOOST_AUTO_TEST_CASE(first_order) {
  const double xFact = 1e1;

  NonLinearPowerLawFitter fitter;
  aocommon::UVector<double> terms{1.0, -0.7};

  for (size_t x = 1; x != 10; ++x) {
    double y = NonLinearPowerLawFitter::Evaluate(x * xFact, terms);
    fitter.AddDataPoint(x * xFact, y);
  }

  double e = 0.0, fact = 0.0;
  fitter.Fit(e, fact);
  BOOST_CHECK_SMALL(std::fabs(1.0 - fact), 1e-6);
  BOOST_CHECK_SMALL(std::fabs(-0.7 - e), 1e-6);
}

BOOST_AUTO_TEST_CASE(second_order_zero) {
  const double xFact = 1e1;

  NonLinearPowerLawFitter fitter;
  aocommon::UVector<double> terms{1.0, -0.7};

  for (size_t x = 1; x != 10; ++x) {
    double y = NonLinearPowerLawFitter::Evaluate(x * xFact, terms);
    fitter.AddDataPoint(x * xFact, y);
  }
  aocommon::UVector<double> fitted;
  fitter.Fit(fitted, 3);
  BOOST_CHECK_SMALL(std::fabs(1.0 - fitted[0]), 1e-6);
  BOOST_CHECK_SMALL(std::fabs(-0.7 - fitted[1]), 1e-6);
  BOOST_CHECK_SMALL(std::fabs(0.0 - fitted[2]), 1e-6);
}

BOOST_AUTO_TEST_CASE(first_order_stability) {
  const double xFact = 1e1;

  NonLinearPowerLawFitter fitter;
  aocommon::UVector<double> terms{1.0, -0.7, -0.01};
  for (size_t x = 1; x != 10; ++x) {
    double y = NonLinearPowerLawFitter::Evaluate(x * xFact, terms);
    fitter.AddDataPoint(x * xFact, y);
  }
  double e = 0.0, fact = 0.0;
  fitter.Fit(e, fact);
  BOOST_CHECK_SMALL(std::fabs(terms[0] - fact), 0.1);
  BOOST_CHECK_SMALL(std::fabs(terms[1] - e), 0.1);
}

BOOST_AUTO_TEST_CASE(second_order_nonzero) {
  const double xFact = 1e1;

  NonLinearPowerLawFitter fitter;
  aocommon::UVector<double> terms{1.0, -0.7, -0.01};
  for (size_t x = 1; x != 10; ++x) {
    double y = NonLinearPowerLawFitter::Evaluate(x * xFact, terms);
    fitter.AddDataPoint(x * xFact, y);
  }
  aocommon::UVector<double> fitted;
  fitter.Fit(fitted, 3);
  BOOST_CHECK_SMALL(std::fabs(terms[0] - fitted[0]), 1e-5);
  BOOST_CHECK_SMALL(std::fabs(terms[1] - fitted[1]), 1e-5);
  BOOST_CHECK_SMALL(std::fabs(terms[2] - fitted[2]), 1e-5);
}

BOOST_AUTO_TEST_CASE(third_order) {
  const double xFact = 1e1;

  NonLinearPowerLawFitter fitter;
  aocommon::UVector<double> terms{1.0, -0.7, -0.01, 0.05};
  for (size_t x = 1; x != 10; ++x) {
    double y = NonLinearPowerLawFitter::Evaluate(x * xFact, terms);
    fitter.AddDataPoint(x * xFact, y);
  }
  aocommon::UVector<double> fitted;
  fitter.Fit(fitted, 4);
  BOOST_CHECK_SMALL(std::fabs(terms[0] - fitted[0]), 1e-2);
  BOOST_CHECK_SMALL(std::fabs(terms[1] - fitted[1]), 1e-2);
  BOOST_CHECK_SMALL(std::fabs(terms[2] - fitted[2]), 1e-2);
  BOOST_CHECK_SMALL(std::fabs(terms[3] - fitted[3]), 1e-2);
}

/**
 * TODO
 */
/*
BOOST_AUTO_TEST_CASE(sed)
{
        PowerLawSED sed;
        double brightness3c196c1[4] = {22.6653, 0.0, 0.0, 0.0};
        aocommon::UVector<double> si3c196c1;
        si3c196c1.push_back(-0.840);
        si3c196c1.push_back(-0.288);
        sed.SetData(150, brightness3c196c1, si3c196c1);
        NonLinearPowerLawFitter fitter4;
        fitter4.AddDataPoint(140.0/150.0, sed.FluxAtFrequencyFromIndex(140, 0));
        fitter4.AddDataPoint(150.0/150.0, sed.FluxAtFrequencyFromIndex(150, 0));
        fitter4.AddDataPoint(160.0/150.0, sed.FluxAtFrequencyFromIndex(160, 0));
        fitter4.AddDataPoint(170.0/150.0, sed.FluxAtFrequencyFromIndex(170, 0));
        aocommon::UVector<double> terms3c196;
        fitter4.Fit(terms3c196, 3);
        std::cout << "3c196 fit: " << terms3c196[0] << " " << terms3c196[1] << "
" << terms3c196[2] << '\n';

        NonLinearPowerLawFitter fitter5;
        double refFreq5 = 5008.0;
        aocommon::UVector<double> data5, freqs5;
        data5.push_back(-0.00489664); freqs5.push_back(4880.0);
        data5.push_back(-0.00493774); freqs5.push_back(5008.0);
        data5.push_back(-0.00476208); freqs5.push_back(5136.0);
        data5.push_back(-0.00366475); freqs5.push_back(5264.0);
        fitter5.AddDataPoint(freqs5[0]/refFreq5, data5[0]);
        fitter5.AddDataPoint(freqs5[1]/refFreq5, data5[1]);
        fitter5.AddDataPoint(freqs5[2]/refFreq5, data5[2]);
        fitter5.AddDataPoint(freqs5[3]/refFreq5, data5[3]);
        aocommon::UVector<double> negTerms;
        fitter5.Fit(negTerms, 2);
        std::cout << "Negative fit: " << negTerms[0] << ' ' << negTerms[1] <<
'\n'; for(size_t i=0; i!=4; ++i)
        {
                std::cout << data5[i] << "->" <<
NonLinearPowerLawFitter::Evaluate(freqs5[i], negTerms, refFreq5) << ' ';
        }
        std::cout << '\n';

        NonLinearPowerLawFitter fitter6;
        fitter6.AddDataPoint(1.0, -1.0);
        fitter6.AddDataPoint(2.0, 1.0);
        fitter6.AddDataPoint(3.0, 0.5);
        fitter6.AddDataPoint(4.0, -0.1);
        aocommon::UVector<double> terms6;
        fitter6.Fit(terms6, 4);
        std::cout << "Pos/neg fit: " << terms6[0] << ' ' << terms6[1] << '\n';
        for(size_t i=0; i!=4; ++i)
                std::cout << '[' << i << "]=" <<
NonLinearPowerLawFitter::Evaluate(i+1, terms6) << ' '; std::cout << '\n';

        NonLinearPowerLawFitter fitter7;
        double refFreq7 = 5008.0;
        aocommon::UVector<double> data7, freqs7;
        data7.push_back(  0.00340717); freqs7.push_back(4880.0);
        data7.push_back(-0.000555695); freqs7.push_back(5008.0);
        data7.push_back(  0.00131385); freqs7.push_back(5136.0);
        data7.push_back( 0.000235965); freqs7.push_back(5264.0);
        fitter7.AddDataPoint(freqs7[0]/refFreq7, data7[0]);
        fitter7.AddDataPoint(freqs7[1]/refFreq7, data7[1]);
        fitter7.AddDataPoint(freqs7[2]/refFreq7, data7[2]);
        fitter7.AddDataPoint(freqs7[3]/refFreq7, data7[3]);
        aocommon::UVector<double> terms7;
        fitter7.Fit(terms7, 2);
        std::cout << "Small pos/neg fit: " << terms7[0] << ' ' << terms7[1] <<
'\n'; for(size_t i=0; i!=4; ++i)
        {
                std::cout << data7[i] << "->" <<
NonLinearPowerLawFitter::Evaluate(freqs7[i], terms7, refFreq7) << ' ';
        }
        std::cout << '\n';

}
*/

BOOST_AUTO_TEST_SUITE_END()
