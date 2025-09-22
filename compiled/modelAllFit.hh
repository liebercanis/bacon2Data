// new version Sept 10 2025
// arXiv:2009.10755v4 [physics.ins-det] 18 Jul 2022
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
  NCHAN = 13,
  MAXSAMPLE = 7500
};

// return fitted function
static double buff[NCHAN][MAXSAMPLE]; // buffer to store light curve data
static double fitWave[NCHAN][MAXSAMPLE];

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
static double distanceLevel[5];
static int iTrigger = 695;

// NOT USING THIS
/* Ion-beam excitation of liquid argon M. Hofmann et al.  Eur. Phys. J. C (2013) 73:2618 */
/*
Electron transport and electron–ion recombination in liquid argon simulation based on the Cohen–Lekner theory
  Mariusz Wojcik a, Tachiya b doi.org/10.1016/S0009-2614(02)01177-6
*/
// static double trecomb = 2.976081E-03; // ns
// NOT USING THIS

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
  // protect against very small value
  // f = max(f, 1.0E-20);
  if (isnan(f))
    printf(" !!!! expGaus NAN!! x %f tau %f args %f %f f%E m\n", x, tau, arg1, arg2, f);
  return f;
}

/*
fcn is required by Minuit to have exactly these argements
returns likelihood value for some set of parameters
*/
static void printModel(int ibin, Double_t *par, double *fsChan, double *ftChan)
{
  double x = double(ibin - iTrigger); // subract trigger sample
  double bw = 2.;                     // ns
  double ppm = par[PPM];
  double norm = par[NORM];
  double tTriplet = par[TAU3];
  double sfrac = par[SFRAC];
  double tMix = par[TAUM];
  double bkg = par[BKGCONST];

  distanceLevel[0] = 1.062; // 11.6;
  distanceLevel[1] = 11.48; // 11.6;
  distanceLevel[2] = 21.41; // 23.2;
  distanceLevel[3] = 31.35; // 34.8;
  distanceLevel[4] = 42.70; // 36.0;

  double kx = kxZero * ppm;                 // rate of tansfer to mixed state
  double kxPrime = kqZero + kx + 1. / tMix; // k_x^\prime in paper
  double tkxPrime = 1. / kxPrime;           // corresponding time

  // convenient rates lambda in paper
  double l1 = 1. / tSinglet0 + kqZero + kx;
  double l3 = 1. / tTriplet + kqZero + kx;
  double lX = 1. / tXe0;

  // corresponding times
  double t1 = 1. / l1;
  double t3 = 1. / l3;

  double effChan[NCHAN];
  double abChan[NCHAN];

  double alpha1Chan[NCHAN];
  double alpha3Chan[NCHAN];

  /* geometric efficiencies */
  double fourPi = 2. * TMath::TwoPi();
  for (int ic = 0; ic < NCHAN; ++ic)
  {
    // level
    int ilevel = 0; // triggger 9,10,11
    if (ic == 6 || ic == 7 || ic == 8)
      ilevel = 1;
    if (ic == 3 || ic == 4 || ic == 5)
      ilevel = 2;
    if (ic == 0 || ic == 1 || ic == 2)
      ilevel = 3;
    if (ic == 12)
      ilevel = 4;
    // ilevel

    /** dist */
    double dist = distanceLevel[ilevel];
    effChan[ic] = pow(0.6, 2.) / fourPi / pow(dist, 2.);
    double aPmt = TMath::Pi() / 4.0 * pow(6.40, 2); // R11410-20  Effective area : 64 mm dia units here are cm
    if (ic == 12)
      effChan[ic] = aPmt / fourPi / pow(dist, 2.);

    // absorption
    double lambda1 = 12.7 * 0.1 / ppm;
    double lambda2 = 740 * 0.1 / ppm;
    double Tr128 = 0.615 * exp(-dist / lambda1) + (1 - 0.615) * exp(-dist / lambda2);
    double ab = 1. - Tr128;
    abChan[ic] = ab;

    double alpha1 = sfrac * bw * norm * effChan[ic];        // singlet norm N1 in paper
    double alpha3 = (1. - sfrac) * bw * norm * effChan[ic]; // triplet norm N3 in paper
    alpha1Chan[ic] = alpha1;
    alpha3Chan[ic] = alpha3;
    // printf("ic %i alpha1 %E alpha3 %E \n", ic, alpha1, alpha3);

    // C constants in paper
    double c1 = kx + ab / tSinglet0;
    double c3 = kx + ab / tTriplet;

    // model emission components terms in equation 6
    // convoute with resolution using expGaus
    fsChan[ic] = (1. - ab) * alpha1 / tSinglet0 * expGaus(x, t1); // singlet
    printf("xxxx ic %i x %f ab %f alpha1 %f t1 %f  exp %E \n", ic, x, ab, alpha1, t1, expGaus(x, t1));
    ftChan[ic] = (1. - ab) * alpha3 / tTriplet * expGaus(x, t3); // triplet
  }

  printf(" \n\n >>> modelFit start parameters\n");
  for (int ii = 0; ii < NPARS; ++ii)
  {
    printf("\t  param %i %s %.4E  \n", ii, lparNames[ii].Data(), par[ii]);
  }

  printf("SiPMQE128Ham %.3f tSinglet0 %E kx %E kxPrime %E l1 %E l3 %E lX %E \n", SiPMQE128Ham, tSinglet0, kx, kxPrime, l1, l3, lX);

  printf("printModel ppm %.2f sample %.0f  \n", ppm, x);
  for (int ic = 0; ic < NCHAN; ++ic)
    printf("chan ic %i effGeo %E abs %f fs %E ft %E alpha1 %E alpha3 %E \n", ic, effChan[ic], abChan[ic], fsChan[ic], ftChan[ic], alpha1Chan[ic], alpha3Chan[ic]);
}
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

  distanceLevel[0] = 1.062; // 11.6;
  distanceLevel[1] = 11.48; // 11.6;
  distanceLevel[2] = 21.41; // 23.2;
  distanceLevel[3] = 31.35; // 34.8;
  distanceLevel[4] = 42.70; // 36.0;

  // loop over channels
  double chanList[5] = {9, 8, 5, 0, 12};
  for (int ichan = 0; ichan < 5; ++ichan)
  {
    int ic = chanList[ichan];
    // if (ic == 5 || ic == 6 || ic == 8 || ic == 3 || ic == 9 || ic == 10 || ic == 11)
    //   continue;

    // **** level
    int ilevel = 0; // triggger 9,10,11
    if (ic == 6 || ic == 7 || ic == 8)
      ilevel = 1;
    else if (ic == 3 || ic == 4 || ic == 5)
      ilevel = 2;
    else if (ic == 0 || ic == 1 || ic == 2)
      ilevel = 3;
    else if (ic == 12)
      ilevel = 4;
    //  **** ilevel
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
    double fourPi = 2. * TMath::TwoPi();
    double effGeo = pow(0.6, 2.) / fourPi / pow(dist, 2.);
    double aPmt = TMath::Pi() / 4.0 * pow(6.40, 2); // R11410-20  Effective area : 64 mm dia units here are cm
    if (ic == 12)
      effGeo = aPmt / fourPi / pow(dist, 2.);

    // values in samples
    int ilow = 650;        //
    int ihigh = MAXSAMPLE; // singlet MAXSAMPLE;
    // singlet region
    // ihigh = 1500;
    // loop over bins to fit
    for (int j = ilow; j < ihigh; ++j) // 7500 is total
    {
      /* skip dip region */
      bool dip = j > 725 && j < 800;
      if (dip)
        continue;
      double x = bw * (double(j - iTrigger) + 0.5);      // bin center convert to ns mutiplying by bin width
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
      fs = fs * SiPMQE128Ham;
      ft = ft * SiPMQE128Ham;
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

      // output fitted function
      fitWave[ic][j] = mval;
      /*if (ic == 9 && (j > 1000 && j < 2000))
        printf("xxxx %i %i %E ", ic, j, fitWave[ic][j]);
        */

      /*******/
      if (mval <= 0)
      {
        // printf("line226 NEGATIVE model value chan %i sample %i x %E fs %E ft %E fx %E fm %E bkg %E set to 0 \n", ic, j, x, fs, ft, fx, fm, bkg);
        continue;
      }

      // build the likelihood f to minimize
      double y = buff[ic][j]; // observe)d
      // we do an NLL
      double yterm = 0.0; // in this case Prob=1 so log=0
      if (y <= 0)         // skip zero bins
        continue;

      yterm = y - y * log(y);
      if (isnan(yterm))
      {
        printf("line255  ibin YTERM is NAN chan %i sample  %i f=%E  y = %E fs %E ft %E fx %E fm %E bkg %E \n", ic, j, f, y, fs, ft, fx, fm, bkg);
      }
      else
      {
        f += mval - y * log(mval) - yterm;
      }

      // ntScan->Fill(par[PPM], fx, f);
      //  leave warnning printout
      if (isnan(f))
      {
        printf("line265  ibin F is NAN chan %i sample %i f=%E mval = %E x = %E y = %E yterm %E effGeo  %E t1 %E t3 %E fs %E ft %E fx %E fm %E bkg %E \n", ic, j, f, mval, x, y, yterm, effGeo, t1, t3, fs, ft, fx, fm, bkg);
      }
    }
  } // loop over channels
}
