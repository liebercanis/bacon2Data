/******************************************************
 file with fit fcn
 new version Sept 10 2025
 the model:
    arXiv:2009.10755v4 [physics.ins-det] 18 Jul 2022
 time is in nanoseconds
 add in radiative component Aug 28 2026
 Agnes et al arXiv:2410.22863v2
*******************************************************/
#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include <TNtuple.h>
#include "distanceLevels.hh"

// light components
enum
{
  SINGLETCOMP,
  TRIPLETCOMP,
  XENONCOMP,
  MIXEDCOMP,
  BKGCOMP,
  RADIATIVECOMP,
  VISIBLECOMP,
  NUMCOMP
};
// fit parameters
enum
{
  NORM = 0,
  TRIGSTART,
  SFRAC,
  PPM,
  TAU3,
  TAUM,
  BKGCONST,
  KXCONST,
  R2CONST,
  R3CONST,
  C1CONST,
  ABSORB1CONST,
  THECHANNEL,
  NPARS
};

enum
{
  NCHAN = 13,
  NCHANPMT = 12,
  MAXSAMPLE = 7500
};

// return fitted function
static double buff[NCHAN][MAXSAMPLE]; // buffer to store light curve data
static double fitWave[NCHAN][MAXSAMPLE];
static double fitComp[NCHAN][NUMCOMP][MAXSAMPLE];
static double lateBkg[NCHAN]; // no longer used

/***** units are nanoseconds ****/
static double shift = 9.;
double tResolution = 3.7 / 2.0; // 6.0;
// static double tResolution = 20.0;
static double tTriplet0 = 944.0; // from fit in range document in fitting 1600.0; // 2100.0;
// static double tTriplet0 = 811.0;
//  static double tTriplet0 = 1600.0; // 2100.0;
static double tSinglet0 = 7.0;
static double tMix0 = 4700.;
static double tXe0 = 20.0;
// from paper
static double kUnit0 = 1.0E-4; // unit to convert to inverse nanoseconds
// static double kqZero = 1.3 * kUnit0; // kq in the paper collisiona de-excitation quenching rate
static double kqZero = 1.3 * kUnit0; // kq in the paper collisiona de-excitation quenching rate
static double kxZero = 2.9 * kUnit0; // kx in the paper diffusion limited reaction rate /[PPM]
//
static double LY = 41.; // Doke, April 2009 https://arxiv.org/abs/0910.4956v1
// LEGEND value 25.6;                              //  photons/kev LEGEND number , ref see Doke
static double nPhotons = 60. * LY;                  // 60 keV gamma
static double visibleYield = 0.2;                   // photons/keV
static double nPhotonsVisible = 60. * visibleYield; //
/* from Bondar arXiv:2202.09*/
static double r2ConstDefault = 3.76E-5;
static double r3ConstDefault = 1.95E-4;
static double C1ConstDefault = 0.09;
static double absorb1ConstDefault = 2568.0;
//
static int iTrigger = 729;

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
static TString compNames[NPARS]; // parameter names
static bool badChannel[NCHAN];
static bool isCovered[NCHAN];

// effiecienies
static double SiPMQE128Ham = 0.15; // 0.15;
static double SiPMQE150 = 0.238;
static double SiPMQE175 = 0.238;
// static double PMTQE175 = 0.38;
static double PMTQE150 = 0.01;
static double PMTQE175 = 0.38;
static double PMTQE400 = 0.35;

static bool isBadChannel(int ichan)
{
  return badChannel[ichan];
}

static void setBadChannels(std::vector<unsigned> list)
{
  for (unsigned ic = 0; ic < NCHAN; ++ic)
  {
    badChannel[ic] = false;
  }
  for (unsigned ib = 0; ib < list.size(); ++ib)
  {
    badChannel[list[ib]] = true;
  }
  printf("modelAllFit::setBadChannels \n");
  for (unsigned ic = 0; ic < NCHAN; ++ic)
  {
    if (badChannel[ic])
      printf("\t bad channel %u  \n", ic);
  }
}

static void setParNames() // tousif
{
  lparNames[NORM] = TString("norm");
  lparNames[TRIGSTART] = TString("trigStart");
  lparNames[SFRAC] = TString("sfrac");
  lparNames[PPM] = TString("ppm");
  lparNames[TAU3] = TString("tau3");
  lparNames[TAUM] = TString("taumix");
  lparNames[BKGCONST] = TString("bkgconst");
  lparNames[KXCONST] = TString("kxconst");
  lparNames[R2CONST] = TString("r2const");
  lparNames[R3CONST] = TString("r3const");
  lparNames[C1CONST] = TString("c1onst");
  lparNames[ABSORB1CONST] = TString("absorb1const");
  lparNames[THECHANNEL] = TString("theChannel");
}

static void setCompNames() // tousif
{
  compNames[SINGLETCOMP] = TString("singletComp");
  compNames[TRIPLETCOMP] = TString("tripletCpmp");
  compNames[XENONCOMP] = TString("xenonComp");
  compNames[MIXEDCOMP] = TString("mixedComp");
  compNames[BKGCOMP] = TString("bkgComp");
  compNames[RADIATIVECOMP] = TString("radiativeComp");
  compNames[VISIBLECOMP] = TString("visibleComp");
}

static void setupModelAllFit()
{
  /* set geomegtry version for modelAllFit,hh */
  // ntScan = new TNtuple("ntScan", "ntScan", "ppm:x:f");

  geoVersionOld = false;
  setDistanceLevels(geoVersionOld);
  printf("setParNames and setCompNames\n");
  setParNames();
  setCompNames();

  /* set bad channels used in fit modelAllFit.hh */
  std::vector<unsigned> badList;
  /*
  badList.push_back(1);
  badList.push_back(2);
  badList.push_back(3);
  badList.push_back(4);
  badList.push_back(0);
  badList.push_back(8);
  */
  setBadChannels(badList);

  // covered channels
  for (unsigned ic = 0; ic < NCHAN; ++ic)
    isCovered[ic] = false;

  // covered set to true
  isCovered[0] = true;
  isCovered[8] = true;
  // isCovered[3] = true;
}

// level
static int getLevel(int ichan)
{
  int ilevel = 0; // triggger 9,10,11
  if (ichan == 6 || ichan == 7 || ichan == 8)
    ilevel = 1;
  if (ichan == 3 || ichan == 4 || ichan == 5)
    ilevel = 2;
  if (ichan == 0 || ichan == 1 || ichan == 2)
    ilevel = 3;
  if (ichan == 12)
    ilevel = 4;
  return ilevel;
}

static double effGeoFunc(int ichan)
{
  /*
  Area of SiPMs is 6.0mm x 6.0mm

      Channels 6, 7, and 8 are at 11.6 cm
      from the source Channels 3, 4, and 5 are at 23.2 cm
      from the source Channels 0, 1, and 2 are at 34.8 cm from the source Channel 12 is at 36 cm from the source.
  */
  int ilevel = getLevel(ichan);
  double aPmt = TMath::Pi() / 4.0 * pow(6.4, 2); // R11410-20  Effective area : 64 mm dia
  double a = pow(0.6, 2.);
  if (ichan == 12)
    a = aPmt;

  double b = 4.0 * TMath::Pi();
  double e = a / b / pow(distanceLevel[ilevel], 2.);
  return e;
}

static double Absorption(double ppm, double dist)
{
  // Calculate absorption as a function of distance and xenon concentration.%
  // Taken from fits to Neumeier data at 0.1 PPM and scaled;
  /* corrected for ppm by mass
  double A = 0.615;
  ppm = max(1.0E-9, ppm);
  double lambda1 = 12.7 * .3 * 0.1 / ppm;
  double lambda2 = 740 * .3 * 0.1 / ppm;
  */
  /* new fit from Doug August 19*/
  ppm = max(1.0E-9, ppm);
  double A = 0.0228;
  double C = 0.643;
  double lambda1 = absorb1ConstDefault * 0.1 / ppm;
  double lambda2 = 3.38 * 0.1 / ppm;
  double lambda3 = 50.5 * 0.1 / ppm;
  double Tr128 = A * exp(-dist / lambda1) + C * exp(-dist / lambda2) + (1 - A - C) * exp(-dist / lambda3);
  // return 1. - Tr128;
  return 0;
}

/** exponential convolution with  gaussian time resolution ***/
static double expGaus(double x, double tau)
{
  x -= shift; // compensate for shift in mean due to smearing of 10 percent
  double arg1 = (tResolution * tResolution / tau - 2. * x) / 2. / tau;
  double arg2 = (tResolution * tResolution / tau - x) / sqrt(2) / tResolution;
  double f = 0.5 * TMath::Exp(arg1) * TMath::Erfc(arg2);
  // protect against very small value
  // f = max(f, 1.0E-20);
  if (isnan(f))
  {
    printf(" !!!! expGaus NAN!! x %f tau %f args %f %f f%E set to 0 \n", x, tau, arg1, arg2, f);
    return 0;
  }
  return f;
}

/*
fcn is required by Minuit to have exactly these argements
returns likelihood value for some set of parameters
*/
static void printModel(int ibin, Double_t *par)
{

  printf(" \n\n >>> printModel modelFit parameters\n");
  for (int ii = 0; ii < NPARS; ++ii)
  {
    printf("\t  param %i %s %.4E  \n", ii, lparNames[ii].Data(), par[ii]);
  }

  printf(" printModel ppm %.2f sample %i \n", par[PPM], ibin);

  double xTrigger = par[TRIGSTART];
  double x = double(ibin) - xTrigger; // subract trigger sample
  double bw = 2.;                     // ns
  double ppm = max(1.0E-9, par[PPM]);
  double norm = par[NORM]; // start is 60*LY;
  double tTriplet = par[TAU3];
  double sfrac = par[SFRAC];
  double tMix = par[TAUM];
  double rad3 = par[R3CONST];
  double rad2 = par[R2CONST];
  double crad = par[C1CONST];
  // double bkg = par[BKGCONST];

  double kx = par[KXCONST] * kxZero * ppm;  // rate of tansfer to mixed state
  double kxPrime = kqZero + kx + 1. / tMix; // k_x^\prime in paper
  double tkxPrime = 1. / kxPrime;           // corresponding time

  // convenient rates lambda in paper
  double l1 = 1. / tSinglet0 + kqZero + kx;
  double l3 = 1. / tTriplet + kqZero + kx;
  double lX = 1. / tXe0;

  // corresponding times
  double t1 = 1. / l1;
  double t3 = 1. / l3;

  double fsChan[NCHAN];
  double ftChan[NCHAN];
  // double fvisChan[NCHAN];
  double mVal[NCHAN];
  double effChan[NCHAN];
  double abChan[NCHAN];
  double alpha1Chan[NCHAN];
  double alpha3Chan[NCHAN];

  /* geometric efficiencies */
  double fourPi = 2. * TMath::TwoPi();
  /*
    distanceLevel[0] = 1.251; // trigger
    distanceLevel[1] = 11.789;
    distanceLevel[2] = 21.723;
    distanceLevel[3] = 31.668;
    distanceLevel[4] = 42.017; // PMT
  */
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
    effChan[ic] = pow(0.6, 2.) / fourPi / pow(distanceLevel[ilevel], 2.);
    double aPmt = TMath::Pi() / 4.0 * pow(6.40, 2); // R11410-20  Effective area : 64 mm dia units here are cm
    if (ic == 12)
      effChan[ic] = aPmt / fourPi / pow(distanceLevel[ilevel], 2.);

    // absorption
    double ab = 1.0;
    if (ppm > 1.0E-3)
    {
      double lambda1 = 12.7 * 0.1 / ppm;
      double lambda2 = 740 * 0.1 / ppm;
      double Tr128 = 0.615 * exp(-distanceLevel[ilevel] / lambda1) + (1 - 0.615) * exp(-distanceLevel[ilevel] / lambda2);
      ab = 1. - Tr128;
    }
    abChan[ic] = ab;

    double alpha1 = sfrac * bw * norm * effChan[ic];        // singlet norm N1 in paper
    double alpha3 = (1. - sfrac) * bw * norm * effChan[ic]; // triplet norm N3 in paper
    double visRatio = visibleYield / LY;
    double alpha1Vis = visRatio * alpha1; // visible norm Nvis in paper
    double alpha2Vis = visRatio * alpha1; // visible norm Nvis in paper
    alpha1Chan[ic] = alpha1;
    alpha3Chan[ic] = alpha3;
    // printf("ic %i alpha1 %E alpha3 %E \n", ic, alpha1, alpha3);

    // C constants in paper
    double c1 = kx + ab / tSinglet0;
    double c3 = kx + ab / tTriplet;

    // model emission components terms in equation 6
    // convoute with resolution using expGaus
    fsChan[ic] = (1. - ab) * alpha1 / tSinglet0 * expGaus(x, t1); // singlet
    // printf("xxxx ic %i x %f ab %f alpha1 %f t1 %f  exp %E \n", ic, x, ab, alpha1, t1, expGaus(x, t1));
    ftChan[ic] = (1. - ab) * alpha3 / tTriplet * expGaus(x, t3); // triplet

    // xenenon emission x_i terms in paper
    double xterm1 = c1 * kx * alpha1 / (l1 - kxPrime) * ((expGaus(x, tkxPrime) - expGaus(x, tXe0)) / (lX - kxPrime) - (expGaus(x, t1) - expGaus(x, tXe0)) / (lX - l1));
    /* fix was c1 not c3 Nov 14 2025*/
    double xterm3 = c3 * kx * alpha3 / (l3 - kxPrime) * ((expGaus(x, tkxPrime) - expGaus(x, tXe0)) / (lX - kxPrime) - (expGaus(x, t3) - expGaus(x, tXe0)) / (lX - l3));
    double fx = (xterm1 + xterm3) / tXe0; // xenon
    // double fx = (xterm1 + xterm3) / tXe0; // xenon
    //  mixed component
    double mterm1 = alpha1 * c1 / (l1 - kxPrime) * (expGaus(x, tkxPrime) - expGaus(x, t1));
    double mterm3 = alpha3 * c3 / (l3 - kxPrime) * (expGaus(x, tkxPrime) - expGaus(x, t3));
    double fm = (mterm1 + mterm3) / tMix; // mixed

    // total light for channel
    mVal[ic] = (fsChan[ic] + ftChan[ic] + fx + fm) * SiPMQE128Ham;
    printf("plotMoedl chan %i ibin %i time %f eff %.2E alpha1 %.2E alpha3 %.2E c1 %.2E c3 %.2E t1  %.2E fs %.2E ft %.2E ft.2E fx %.2E fm %.2E mval%.2E \n", ic, ibin, x, effChan[ic], alpha1, alpha3, c1, c3, t1, fsChan[ic], ftChan[ic], fx, fm, mVal[ic]);
    printf("printxxxx alpha1 %.2E sfrac %.2E bw %.2E norm %.2E eff %.2E \n", alpha1, sfrac, bw, norm, effChan[ic]);
  }

  for (int ic = 9; ic < NCHAN; ++ic)
    printf("chan ic %i mval %.2E effGeo %E abs %f fs %E ft %E alpha1 %E alpha3 %E \n", ic, mVal[ic], effChan[ic], abChan[ic], fsChan[ic], ftChan[ic], alpha1Chan[ic], alpha3Chan[ic]);

  printf("time %i = %.3f  siPMQE128Ham %.3f tSinglet0 %E kx %E kxPrime %E l1 %E l3 %E lX %E \n", ibin, x, SiPMQE128Ham, tSinglet0, kx, kxPrime, l1, l3, lX);

  printf("DENOMINATORS kxPrime %E lx %E l1  %E  l3  %E\n", kxPrime, lX, l1, l3);
  printf("DENOMINATORS lx - kxPrime %E lx - l1 %E lx - l3 %E\n", lX - kxPrime, lX - l1, lX - l3);
}
/*
 **************************************************************************
 **********************   fcn function fit by minuit **************************
 **********************   Minuit specifies function arguments
 **********************   returned f is NLL value ************************
 **************************************************************************
 */
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // pack parameters into static array onto lightModel
  // for (int k = 0; 8 < npar; ++k)
  //  lpar[k] = par[k];

  f = 0;          // return value
  double bw = 2.; // ns
  double ppm = max(1.0E-9, par[PPM]);

  double norm = par[NORM];
  double tTriplet = par[TAU3];
  double sfrac = par[SFRAC];
  double tMix = par[TAUM];
  double rad3 = par[R3CONST];
  double rad2 = par[R2CONST];
  double c1rad = par[C1CONST];
  double abs1Const = par[ABSORB1CONST];

  // double bkg = par[BKGCONST];

  /****************
   **************** loop over channels *************
   * **************/
  // double chanList[3] = {8, 5, 0};
  for (int ic = 0; ic < NCHAN; ++ic)
  {
    /* for fitting single channel
     ************  par[THECHANNEL] =-1 for all
     */
    if (par[THECHANNEL] > 0 && ic != par[THECHANNEL])
    {
      continue;
    }

    if (par[THECHANNEL] == -2)
    { // dont use trigger sipms
      if (ic > 8)
        continue;
    }

    /* skip bad channels */
    if (isBadChannel(ic))
    {
      printf("fcn skipping bad channel %i \n", ic);
      continue;
    }

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
    double dist = distanceLevel[ilevel]; // from header distanceLevels.hh
    // SiPMQ128
    double SiPMQ128 = SiPMQE128Ham;

    /****
     * absorption as a function of distance and xenon concentration.%
     ****/
    double ab = 0.0;
    double A = 0.0228;
    double C = 0.643;
    double lambda1 = abs1Const * 0.1 / ppm;
    double lambda2 = 3.38 * 0.1 / ppm;
    double lambda3 = 50.5 * 0.1 / ppm;
    double Tr128 = A * exp(-dist / lambda1) + C * exp(-dist / lambda2) + (1 - A - C) * exp(-dist / lambda3);
    ab = 1. - Tr128;
    ab = 0;

    /* geometric efficiencies */
    double fourPi = 2. * TMath::TwoPi();
    double effGeo = pow(0.6, 2.) / fourPi / pow(dist, 2.);
    double aPmt = TMath::Pi() / 4.0 * pow(6.40, 2); // R11410-20  Effective area : 64 mm dia units here are cm
    if (ic == 12)
      effGeo = aPmt / fourPi / pow(dist, 2.);
    /*
        set to 1 using geo normalized data
    */
    effGeo = 1.;

    // values in samples
    int ilow = 650; //
    int ihigh = 7500;
    // MAXSAMPLE; // singlet MAXSAMPLE;
    //  singlet region
    //  ihigh = 1500;
    /****
     ****   loop over bins to fit
     */
    for (int j = ilow; j < ihigh; ++j) // 7500 is total samples
    {
      /* skip dip region for trigger sipms */
      bool dip = j > 1400 / 2 && j < 1700 / 2;
      // if (dip && ilevel == 0) //
      //   continue;
      /* set sample starting based on TRIGSTART*/
      double xTrigger = par[TRIGSTART];
      double x = (double(j - int(xTrigger / 2.0)) + 0.5); // bin center convert to ns mutiplying by bin width
      double alpha1 = sfrac * bw * norm * effGeo;         // singlet norm N1 in paper
      double alpha3 = (1. - sfrac) * bw * norm * effGeo;  // triplet norm N3 in paper
      double kx = par[KXCONST] * kxZero * ppm;            // rate of tansfer to mixed state
      double kxPrime = kqZero + kx + 1. / tMix;           // k_x^\prime in paper
      double tkxPrime = 1. / kxPrime;                     // corresponding time

      /* visible yield */
      double visRatio = visibleYield / LY;
      double alpha1Vis = visRatio * alpha1; // visible norm Nvis in paper
      double alpha2Vis = visRatio * alpha1; // visible norm Nvis in paper

      // convenient rates lambda in paper
      double l1 = 1. / tSinglet0 + kqZero + kx; // lamda 1 in paper
      double l3 = 1. / tTriplet + kqZero + kx;  // lamda 3 in paper
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
      // visible component
      double fvis = 0.; // visibleYield * norm * effGeo;

      // radiative conpontent
      double krad3 = rad3 + 1. / tMix;
      double krad23 = rad2 - krad3;
      // here time x is nanoseconds
      double tnano = 2 * x;
      double pRadiative = rad2 * krad3 / krad23 * TMath::Exp(-krad3 * tnano) * (1. - TMath::Exp(-krad23 * tnano));
      pRadiative = max(pRadiative, 1.E-9);
      double frad = c1rad * bw * norm * effGeo * abs(pRadiative);

      // xenenon emission x_i terms in paper
      double xterm1 = c1 * kx * alpha1 / (l1 - kxPrime) * ((expGaus(x, tkxPrime) - expGaus(x, tXe0)) / (lX - kxPrime) - (expGaus(x, t1) - expGaus(x, tXe0)) / (lX - l1));
      /* fix was c1 not c3 Nov 14 2025*/
      double xterm3 = c3 * kx * alpha3 / (l3 - kxPrime) * ((expGaus(x, tkxPrime) - expGaus(x, tXe0)) / (lX - kxPrime) - (expGaus(x, t3) - expGaus(x, tXe0)) / (lX - l3));
      // divide by Xenon lifetime from equation 6
      double fx = (xterm1 + xterm3) / tXe0;
      //  mixed component
      double mterm1 = alpha1 * c1 / (l1 - kxPrime) * (expGaus(x, tkxPrime) - expGaus(x, t1));
      double mterm3 = alpha3 * c3 / (l3 - kxPrime) * (expGaus(x, tkxPrime) - expGaus(x, t3));
      double fm = (mterm1 + mterm3) / tMix; // mixed

      // additional factors depend on SIPM channel
      if (ic == 0) // poly carbon covered sees only above 4000 nm
      {
        fs = 0;
        ft = 0;
        fm = 0;
        fx = 0;
        frad = 0;
      }
      else if (ic == 8) // glass covered sees only  175
      {
        fs = 0;
        ft = 0;
        fm = fm * SiPMQE150;
        frad = frad * SiPMQE150;
        fx = fx * SiPMQE175;
      }
      else if (ic == 12) // PMT
      {
        fs = 0;
        ft = 0;
        fm = fm * PMTQE150;
        frad = frad * SiPMQE150;
        fx = fx * PMTQE175;
      }
      else
      { // all other sipms
        // multiply by efficiencies
        fs = fs * SiPMQE128Ham;
        ft = ft * SiPMQE128Ham;
        fm = fm * SiPMQE150;
        frad = frad * SiPMQE150;
        fx = fx * SiPMQE175;
      }

      // total light for channel and bin ic
      // double mval = fs + ft + fx + fm + lateBkg[ic];
      double mval = fs + ft + fx + fm;
      // printf("line501  model value chan %i sample %i x %E fs %E ft %E fx %E fm %E frad %E mval %E  \n", ic, j, x, fs, ft, fx, fm, frad, mval);
      //  mval =] fs + ft;
      //   for plottting components
      fitComp[ic][SINGLETCOMP][j] = fs;
      fitComp[ic][TRIPLETCOMP][j] = ft;
      fitComp[ic][XENONCOMP][j] = fx;
      fitComp[ic][RADIATIVECOMP][j] = frad;
      fitComp[ic][VISIBLECOMP][j] = fvis;
      // if (ic == 8 && j == 1500)
      //   printf("....line418 sample %i fx %E \n", j, fx);
      fitComp[ic][MIXEDCOMP][j] = fm;
      // fitComp[ic][BKGCOMP][j] = lateBkg[j];
      //  output fitted function
      fitWave[ic][j] = mval;

      // printf("line614 modelFitAll chan %i j %i --  fs %E ft %E fx %E fm %E frad %E mval %E \n", ic, j, fs, ft, fx, fm, frad, mval);

      /*******/
      if (mval <= 0)
      {
        // printf("line563 NEGATIVE model value chan %i sample %i xTrigger %i x %E fs %E ft %E fx %E fm %E bkg %E set to 0 \n", ic, j, j - int(xTrigger / 2.0), x, fs, ft, fx, fm, mval);
        continue;
      } /* for bug in not fixed TRIGSTART
       else if (j > int(xTrigger) / 2 && j < int(xTrigger) / 2 + 100) // example condition to check for a specific bin
         printf("line567 POSITIVE model value chan %i sample %i xT
         */

      // build the likelihood f to minimize
      double y = buff[ic][j]; // observe)d
      // we do an NLL
      double yterm = 0.0; // in this case Prob=1 so log=0
      if (y <= 0)
      {
        // skip zero bins
        continue;
      }

      yterm = y - y * log(y);
      if (isnan(yterm))
      {
        printf("line577 ibin YTERM is NAN chan %i sample  %i f=%E  y = %E fs %E ft %E fx %E fm %E bkg %E \n", ic, j, f, y, fs, ft, fx, fm, lateBkg[j]);
      }
      else
      { // sum up NLL as defined ROOT documentation
        f += mval - y * log(mval) - yterm;
        // printf("line555 modelFitAll chan %i j %i y %E fx %E frad %E fm %E mval %E nLL = %E\n", ic, j, y, fx, fm, frad, mval, f);
      }

      // ntScan->Fill(par[PPM], fx, f);
      //   leave warnning printout
      if (isnan(f))
      {
        printf("line265  ibin F is NAN chan %i sample %i f=%E mval = %E x = %E y = %E yterm %E rad3 %E krad23 %E fs %E ft %E fx %E fm %E frad %E \n", ic, j, f, mval, x, y, yterm, rad3, krad23, fs, ft, fx, fm, frad);
      }
    }

    // printf("channel modelFitAll channelic, ic %i nLL = %E\n", ic, f);
  } // loop over channels
  // printf("return modelFitAll nLL = %E\n", f);
}
