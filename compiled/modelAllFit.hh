// new version Sept 10 2025
// time is in nanoseconds
#include "TString.h"
#include "TF1.h"
#include <TNtuple.h>

// fit parameters
enum
{
  NORM = 0,
  SFRAC,
  PPM,
  TAU3,
  TAUM,
  BKGCONST,
  BKGTAU,
  NPARS
};

enum
{
  MAXSAMPLE = 7500
};

// return fitted function
static double fitWave[MAXSAMPLE];

TNtuple *ntScan = new TNtuple("ntScan", "ntScan", "ppm:fx:f");
// units are nanoseconds
static double tResolution = 7.0;
static double tTriplet0 = 1600.0; // 2100.0;
static double tSinglet0 = 7.0;
static double tMix0 = 4700.;
static double tXe0 = 20.0;
// from paper
static double kUnit0 = 1.0E-4;       // unit to convert to inverse nanoseconds
static double kqZero = 1.3 * kUnit0; // kq in the paper collisiona de-excitation quenching rate
static double kxZero = 2.9 * kUnit0; // kx in the paper diffusion limited reaction rate /[PPM]
//
static double LY = 25.6;           //  photons/kev LEGEND number , ref see Doke
static double nPhotons = 60. * LY; // 60 keV gamma
//
static double distanceLevel[4];
static double xTrigger = 1411.; // peak fit

// NOT USING THIS
/* Ion-beam excitation of liquid argon M. Hofmann et al.  Eur. Phys. J. C (2013) 73:2618 */
/*
Electron transport and electron–ion recombination in liquid argon simulation based on the Cohen–Lekner theory
  Mariusz Wojcik a, Tachiya b doi.org/10.1016/S0009-2614(02)01177-6
*/
// static double trecomb = 2.976081E-03; // ns
// NOT USING THIS

static double buff[13][7500];    // buffer to store light curve data
static double lpar[NPARS];       // pass parameters to light model
static TString lparNames[NPARS]; // parameter names

// effiecienies
static double SiPMQE128Ham = 0.15;
static double SiPMQE150 = 0.238;
static double SiPMQE175 = 0.238;
// static double PMTQE175 = 0.38;
static double PMTQE150 = 0.01;
static double PMTQE175 = 0.38;
static double PMTQE400 = 0.35;

static void setParNames() // tousif
{
  lparNames[NORM] = TString("norm");
  lparNames[SFRAC] = TString("sfrac");
  lparNames[PPM] = TString("ppm");
  lparNames[TAU3] = TString("tau3");
  lparNames[TAUM] = TString("taumix");
  lparNames[BKGCONST] = TString("bkgconst");
  lparNames[BKGTAU] = TString("bkgtau");
}

static double Absorbtion(double ppm, double dist)
{
  // Calculate absorption as a function of distance and xenon concentration.%
  // Taken from fits to Neumeier data at 0.1 PPM and scaled;
  double A = 0.615;
  ppm = max(1.0E-9, ppm);
  double lambda1 = 12.7 * 0.1 / ppm;
  double lambda2 = 740 * 0.1 / ppm;
  double Tr128 = A * exp(-dist / lambda1) + (1 - A) * exp(-dist / lambda2);
  return 1. - Tr128;
}

static double expGaus(double x, double tau)
{
  x += 7.2; // compensate for shift in mean due to smearing of 10 percent
  double arg1 = (tResolution * tResolution / tau - 2. * x) / 2. / tau;
  double arg2 = (tResolution * tResolution / tau - x) / sqrt(2) / tResolution;
  double f = 0.5 * TMath::Exp(arg1) * TMath::Erfc(arg2);
  return f;
}

/*
fcn is required by Minuit to have exactly these argements
returns likelihood value for some set of parameters
*/
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // pack parameters into static array onto lightModel
  // for (int k = 0; 8 < npar; ++k)
  //  lpar[k] = par[k];

  f = 0;          // return value
  double bw = 2.; // ns
  double ppm = par[PPM];
  double norm = par[NORM];
  double tTriplet = par[TAU3];
  double sfrac = par[SFRAC];
  double tMix = par[TAUM];
  double bkg = par[BKGCONST];

  distanceLevel[0] = 11.6;
  distanceLevel[1] = 23.2;
  distanceLevel[2] = 34.8;
  distanceLevel[3] = 36.0;
  // loop over channels
  for (int ic = 0; ic < 13; ++ic)
  {
    if (ic == 5 || ic == 6 || ic == 8 || ic == 3 || ic == 9 || ic == 10 || ic == 11)
      continue;

    // Fit only the PMT
    if (ic != 12)
      continue;

    // level
    int ilevel = -1;
    if (ic == 6 || ic == 7 || ic == 8)
      ilevel = 0;
    else if (ic == 3 || ic == 4 || ic == 5)
      ilevel = 1;
    else if (ic == 0 || ic == 1 || ic == 2)
      ilevel = 2;
    else if (ic == 12)
      ilevel = 3;
    /** dist */
    double dist = distanceLevel[ilevel];
    // SiPMQ128
    double SiPMQ128 = SiPMQE128Ham;

    // absorption
    double lambda1 = 12.7 * 0.1 / ppm;
    double lambda2 = 740 * 0.1 / ppm;
    double Tr128 = 0.615 * exp(-dist / lambda1) + (1 - 0.615) * exp(-dist / lambda2);
    double ab = 1. - Tr128;

    /* geometric efficiencies */
    double fourPi = 12.566371;
    double effGeo = pow(0.6, 2.) / fourPi / pow(dist, 2.);
    double aPmt = TMath::Pi() / 4.0 * pow(6.40, 2); // R11410-20  Effective area : 64 mm dia units here are cm
    if (ic == 12)
      effGeo = aPmt / fourPi / pow(dist, 2.);

    int ilow = 0;
    int ihigh = 7500;
    // loop over bins to fit
    for (int j = ilow; j < ihigh; ++j) // 7500 is total
    {
      double xbin = bw * (double(j) + 0.5); // bin center convert to ns mutiplying by bin width
      double x = xbin - xTrigger;
      double alpha1 = sfrac * bw * norm * effGeo;        // singlet norm N1 in paper
      double alpha3 = (1. - sfrac) * bw * norm * effGeo; // triplet norm N3 in paper
      double kx = kxZero * ppm;                          // rate of tansfer to mixed state
      double kxPrime = kqZero + kx + 1. / tMix;          // k_x^\prime in paper
      double tkxPrime = 1. / kxPrime;                    // corresponding time

      // convenient rates lambda in paper
      double l1 = 1. / tSinglet0 + kqZero + kx;
      double l3 = 1. / tTriplet + kqZero + kx;
      double lX = 1. / tXe0;

      // corresponding times
      double t1 = 1. / l1;
      double t3 = 1. / l3;

      // C constants in paper
      double c1 = kx + ab / tSinglet0;
      double c3 = kx + ab / tTriplet;

      // model emission components terms in equation 6
      // convoute with resolution using expGaus
      double fs = (1. - ab) * alpha1 / tSinglet0 * expGaus(x, t1); // singlet
      double ft = (1. - ab) * alpha3 / tTriplet * expGaus(x, t3);  // triplet

      // xenenon emission x_i terms in paper
      double xterm1 = c1 * kx * alpha1 / (l1 - kxPrime) * ((expGaus(x, tkxPrime) - expGaus(x, tXe0)) / (lX - kxPrime) - (expGaus(x, t1) - expGaus(x, tXe0)) / (lX - l1));

      double xterm3 = c1 * kx * alpha1 / (l3 - kxPrime) * ((expGaus(x, tkxPrime) - expGaus(x, tXe0)) / (lX - kxPrime) - (expGaus(x, t3) - expGaus(x, tXe0)) / (lX - l3));
      double fx = (xterm1 + xterm3) / tXe0; // xenon

      // mixed component
      double mterm1 = alpha1 * c1 / (l1 - kxPrime) * (expGaus(x, tkxPrime) - expGaus(x, t1));
      double mterm3 = alpha3 * c3 / (l3 - kxPrime) * (expGaus(x, tkxPrime) - expGaus(x, t3));
      double fm = (mterm1 + mterm3) / tMix; // mixed

      // multiply by efficiencies
      fs = fs * SiPMQ128;
      ft = ft * SiPMQ128;
      fm = fm * SiPMQE150;

      // additional factors depend on SIPM channel
      if (ic == 1 || ic == 3) // glass covered sees only  175
      {
        fs = 0;
        ft = 0;
        fm = 0;
      }
      else if (ic == 12) // PMT
      {
        fx = fx * PMTQE175;
        fm = fm * PMTQE150;
        fs = 0;
        ft = 0;
      }
      else
      { // all other sipms
        fx = fx * SiPMQE175;
      }

      // total light for channel
      double mval = fs + ft + fx + fm + bkg;

      // fitted function
      fitWave[j] = mval;

      /*******/
      if (mval < 0)
      {
        printf("line226 NEGATIVE model value ibin %i x %E set to 1 \n", j, x);
        mval = 1.;
      }

      // build the likelihood f to minimize
      double y = buff[ic][j]; // observed
      // we do an NLL
      double yterm = 0.0; // in this case Prob=1 so log=0
      if (y > 0)
        yterm = y - y * log(y);

      f += mval - y * log(mval) - yterm;

      ntScan->Fill(par[PPM], fx, f);
      // leave warnning printout
      if (isnan(f))
      {
        printf("line243  ibin f is NAN %i x = %E y = %E \n", j, x, y);
      }
    }
  } // loop over channels
}
