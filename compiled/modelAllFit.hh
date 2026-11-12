/*****model fitting to all channels at once revised Oct 2026 * version using new absorption model
                                                                  *and new light model with diffusion limited reaction rate and collision de -
    excitation quenching rate
            *Nov 3 2026 *
        ***/
// light components
#ifndef MODELALLFIT_H
#define MODELALLFIT_H
#include <cmath>
#include <vector>
#include "TMinuit.h"
#include "TSpline.h"
enum
{
  SINGLETCOMP,
  TRIPLETCOMP,
  MIXEDSINGLETCOMP,
  MIXEDTRIPLETCOMP,
  XENONCOMP,
  VISIBLECOMP,
  NUMCOMP
};
// fit parameters
enum moldelFitParameters
{
  NORM = 0,
  PPM,
  TAUT,
  TAUQUENCH,
  TAUMIXS,
  TAUMIXT,
  TAUMIXR,
  BRAD,
  TAUQUECHMIX,
  THELEVEL,
  NPARS
};

enum
{
  NCHAN = 13,
  NCHANPMT = 12,
  MAXSAMPLE = 7500,
  NSIPMLEVELS = 4,
  NLEVEL = 7
  /******
    levels 0,1,2,3 are SIPMS
    level 4 is PMT, level 5 is SIPM 0 and level 6 is SIPM 8
  ***/
};

TMinuit *gMinuit;
// fcn prototype
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
// return fitted function
static bool fixPar[NPARS] = {true, true, true, true, true, true, true, true, true, true};
static double buff[NLEVEL][MAXSAMPLE]; // buffer to store light curve data
static double fitWave[NLEVEL][MAXSAMPLE];
static double fitComp[NLEVEL][NUMCOMP][MAXSAMPLE];
static double lateBkg[NLEVEL]; // no longer used
vector<vector<int>> sensorLevels;

static const double inf = std::numeric_limits<double>::infinity();
std::vector<double> startIdx;
std::vector<double> startVal;
static int nominalTrigger = 700;                         // nominal trigger bin Nov 6
double RsingeltTriplet = 0.25;                           // fraction of singlet/triplet light that is prompt (1/4 for LAr, 1/3 for LXe)
double Sfrac = RsingeltTriplet / (1. + RsingeltTriplet); //

/*** define singlet integral region Sept 28 2026 */
static int earlyBin = 1400 / 2; //
static int lateBin = 1600 / 2;

/***** units are nanoseconds ****/
static double shift = 0.;
double tResolution = 11.43; // 11.43; // base on fit to singlet level 1 Nov 11;
// static double tResolution = 20.0;
/* old consants used in old light model calcLightComponents */
static double tSinglet0 = 7.0;
static double tMix0 = 4700.;
static double tXe0 = 20.0;
static double tTriplet0 = 944.0; // from fit in range document in fitting 1600.0; // 2100.0;

// from paper
static double kUnit0 = 1.0E-4;       // unit to convert to inverse nanoseconds
static double kqZero = 1.3 * kUnit0; // kq in the paper collisiona de-excitation quenching rate
static double kxZero = 2.9 * kUnit0; // kx in the paper diffusion limited reaction rate /[PPM]
static double LY = 41.;              // Doke, April 2009 https://arxiv.org/abs/0910.4956v1
// LEGEND value 25.6;                              //  photons/kev LEGEND number , ref see Doke
static double nPhotons = 60. * LY;                  // 60 keV gamma
static double visibleYield = 0.2;                   // photons/keV
static double nPhotonsVisible = 60. * visibleYield; //

/* from Bondar arXiv:2202.09 radiative component
static double r2ConstDefault = 3.76E-5;
static double r3ConstDefault = 1.95E-4;
static double C1ConstDefault = 0.09;
*/

/****absprption model Doug Oct. 2 2026 with PPM_MOLAR
static double lambda1ConstDefault = 16.59;   // cm at 0.1 PPM_MOLAR
static double lambda2ConstDefault = 75.81;   // cm at 0.1 PPM_MOLAR
static double lambda3ConstDefault = 2995.82; // cm at 0.1 PPM_MOLAR
static double AConstantDefault = 0.4417;
static double CConstantDefault = 0.1936;
*****/

/**** absprption model Doug Oct. 6 2026 with PPM_MOLAR
    static double lambda1ConstDefault = 9.20;    // cm at 0.1 PPM_MOLAR
    static double lambda2ConstDefault = 57.96;   // cm at 0.1 PPM_MOLAR
    static double lambda3ConstDefault = 4391.16; // cm at 0.1 PPM_MOLAR
    static double AConstantDefault = 0.3720;
    static double CConstantDefault = 0.1653;
******/

/****
       absprption model Doug Oct. 20 2026 with PPM_MOLAR
******/
static double lambda1ConstDefault = 4.49;    // cm at 0.1 PPM_MOLAR
static double lambda2ConstDefault = 38.55;   // cm at 0.1 PPM_MOLAR
static double lambda3ConstDefault = 5694.54; // cm at 0.1 PPM_MOLAR
static double AConstantDefault = 0.4010;
static double CConstantDefault = 0.1399;

static double backgroundConstDefault = 0.04153;
//
static int iTrigger = 729;

// values in samples
static int iLowFitRange = 600; //
                               // int iHighFitRange = 7500;
static int iHighFitRange = 15000 / 2;

// NOT USING THIS
/* Ion-beam excitation of liquid argon M. Hofmann et al.  Eur. Phys. J. C (2013) 73:2618 */
/*
Electron transport and electron–ion recombination in liquid argon simulation based on the Cohen–Lekner theory
  Mariusz Wojcik a, Tachiya b doi.org/10.1016/S0009-2614(02)01177-6
*/
// static double trecomb = 2.976081E-03; // ns
// NOT USING THIS

/* new model constants */
/*taus*/
static double tausConstant = 4;
/*taut*/
static double tautConstant = 1184;
/*taudx_S*/
static double taudx_SConstant = 50.0;
/*taudx_T*/
static double taudx_TConstant = 3000.0;
/*taudx_X*/
static double taudx_XConstant = 800.0;
/*taudn*/
static double taudnConstant = inf;
/*tauquench*/
static double tauquenchConstant = 2250;
/*tauquenchmix*/
static double tauquenchmixConstant = 3200;
/*taumix_S*/
static double taumix_SConstant = 35.0;
/*taumix_T*/
static double taumix_TConstant = 4700.0;
/*taumix_R*/
static double taumix_RConstant = 250.0;
/*taux*/
static double tauxConstant = 20;
/*k_rec*/
static double k_recConstant = 1.0 / 36.0;
/* Absorb defaults to 0.3 */
/*b_rad*/
static double b_radConstant = 0.70;
// used to calculate taudx_S
static double RArArsingletConstant = 52.0e-8; // Foerster Resonant Energy Transfer radius, only happens for the singlet (dipole) state (just guess for now) {VARY}?
double fitNorm;                               // not "norm" -- collides/ambiguous with std::norm (<complex>) once pulled into scope by ROOT headers
double tTriplet;
double tauQuench;
double tauMixS;
double tauMixT;
double tauMixR;
double bRad;
double tauQuenchMix;
double k_rec = 1.0 / 36.0;

static double lpar[NPARS];       // pass parameters to light model
static Double_t step[NPARS];     ///< Step sizes for Minuit parameter exploration
static TString lparNames[NPARS]; // parameter names
static TString compNames[NPARS]; // parameter names
static bool badChannel[NCHAN];
static bool isCovered[NCHAN];

// effiecienies
static double SiPMQE128Ham = 0.15; // 0.15;
static double SiPMQE150 = 0.238;
static double SiPMQE175 = 0.238;
static double SiPMQuartz = 0.93;
// static double PMTQE175 = 0.38;
static double PMTQE150 = 0.01;
static double PMTQE175 = 0.38;
static double PMTQE400 = 0.35;

// SiPM (VUV4) photon detection efficiency vs wavelength [nm], ported from
// sipm_detector_qe() in model_with_absorption_and_convolution7.py -- same
// wavelength/PDE arrays and cubic interpolation as SipmDetectorQE() in
// lar_analytic.C (there via TSpline3, matching scipy's interp1d(kind='cubic'))
static const double sipmQEWavelength[] = {120, 128, 150, 175, 200, 250, 300, 350, 400, 450,
                                          500, 550, 600, 650, 700, 750, 800, 850, 900};
static const double sipmQEPde[] = {0.12, 0.15, 0.20, 0.24, 0.28, 0.33, 0.38, 0.44, 0.48, 0.50,
                                   0.46, 0.38, 0.28, 0.20, 0.14, 0.09, 0.05, 0.03, 0.01};
static const int sipmQENpts = sizeof(sipmQEWavelength) / sizeof(sipmQEWavelength[0]);

// evaluate the SiPM QE curve at a single wavelength [nm]; 0 outside the
// tabulated 120-900 nm range
static double SipmDetectorQE(double wavelength)
{
  static TSpline3 sipmQESpline("sipmQE", const_cast<double *>(sipmQEWavelength),
                               const_cast<double *>(sipmQEPde), sipmQENpts);
  if (wavelength < sipmQEWavelength[0] || wavelength > sipmQEWavelength[sipmQENpts - 1])
    return 0.0;
  return sipmQESpline.Eval(wavelength);
}

// vector form, matching the python/lar_analytic.C call signature
static std::vector<double> SipmDetectorQE(const std::vector<double> &wavelengths)
{
  std::vector<double> out(wavelengths.size());
  for (size_t i = 0; i < wavelengths.size(); ++i)
    out[i] = SipmDetectorQE(wavelengths[i]);
  return out;
}

// piecewise-linear interpolation, matching scipy's
// interp1d(kind='linear', bounds_error=False, fill_value=0.0): returns 0
// outside [x[0], x[n-1]] instead of extrapolating
static double LinearInterpFillZero(const double *x, const double *y, int n, double xv)
{
  if (xv < x[0] || xv > x[n - 1])
    return 0.0;
  for (int i = 0; i < n - 1; ++i)
  {
    if (xv >= x[i] && xv <= x[i + 1])
    {
      double frac = (x[i + 1] == x[i]) ? 0.0 : (xv - x[i]) / (x[i + 1] - x[i]);
      return y[i] + frac * (y[i + 1] - y[i]);
    }
  }
  return 0.0;
}

// polycarbonate window transmission vs wavelength [nm], ported from
// polycarbonate_trans() in model_with_absorption_and_convolution7.py -- same
// wavelength/transmission arrays and linear interpolation as
// PolycarbonateTrans() in lar_analytic.C
static const double polycarbTransWavelength[] = {120, 128, 150, 175, 200, 250, 300, 380, 410, 450,
                                                 500, 550, 600, 650, 700, 750, 800, 850, 900};
static const double polycarbTransValue[] = {0.00, 0.00, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.85, 0.85,
                                            0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85};
static const int polycarbTransNpts = sizeof(polycarbTransWavelength) / sizeof(polycarbTransWavelength[0]);

// evaluate the polycarbonate transmission curve at a single wavelength [nm];
// 0 outside the tabulated 120-900 nm range, clamped to [0, 1]
static double PolycarbonateTrans(double wavelength)
{
  double v = LinearInterpFillZero(polycarbTransWavelength, polycarbTransValue, polycarbTransNpts, wavelength);
  return std::clamp(v, 0.0, 1.0);
}

// vector form, matching the python/lar_analytic.C call signature
static std::vector<double> PolycarbonateTrans(const std::vector<double> &wavelengths)
{
  std::vector<double> out(wavelengths.size());
  for (size_t i = 0; i < wavelengths.size(); ++i)
    out[i] = PolycarbonateTrans(wavelengths[i]);
  return out;
}

// quartz window transmission vs wavelength [nm], ported from
// quartz_trans() in model_with_absorption_and_convolution7.py -- same
// wavelength/transmission arrays and linear interpolation as QuartzTrans()
// in lar_analytic.C
static const double quartzTransWavelength[] = {120, 130, 140, 150, 155, 160, 165, 170, 175, 180,
                                               200, 300, 400, 500, 600, 700, 800, 850, 900};
static const double quartzTransValue[] = {0.0036, 0.0022, 0.0014, 0.0035, .1344, .3840, .7008, .8605, .8862,
                                          .8946, .9084, .9247, .9289, .9307, .9318, .9324, .9329, .9331, .9331};
static const int quartzTransNpts = sizeof(quartzTransWavelength) / sizeof(quartzTransWavelength[0]);

// evaluate the quartz transmission curve at a single wavelength [nm]; 0
// outside the tabulated 120-900 nm range, clamped to [0, 1]
static double QuartzTrans(double wavelength)
{
  double v = LinearInterpFillZero(quartzTransWavelength, quartzTransValue, quartzTransNpts, wavelength);
  return std::clamp(v, 0.0, 1.0);
}

// vector form, matching the python/lar_analytic.C call signature
static std::vector<double> QuartzTrans(const std::vector<double> &wavelengths)
{
  std::vector<double> out(wavelengths.size());
  for (size_t i = 0; i < wavelengths.size(); ++i)
    out[i] = QuartzTrans(wavelengths[i]);
  return out;
}

static double Gaussian(double xv, double amp, double mu, double sigma)
{
  return amp * std::exp(-((xv - mu) * (xv - mu)) / (2 * sigma * sigma));
}

static double Trapezoid(const std::vector<double> &y, const std::vector<double> &x)
{
  double sum = 0;
  for (size_t i = 1; i < x.size(); ++i)
    sum += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]);
  return sum;
}

// 128nm emission line shape and Rayleigh baseline, ported from
// generate_emission_spectrum() in model_with_absorption_and_convolution7.py
static const double emissionA_E = 8132.0, emissionMu_E = 128.2, emissionSigma_E = 3.006;
static const double emissionRayleighBaseline = 1.00;

// 128nm absorption-cross-section double-Gaussian fit parameters, ported from
// the same python script's fit to Neumeier_absorption.csv / 1ppmabsorption.csv
// / 10ppmabsorption.csv / 1000ppmabsorption.csv (see FitAbsorptionCrossSection()
// in lar_analytic.C, which this was copied from); re-running that fit here
// (root -l -b -q lar_analytic.C -e 'FitAbsorptionCrossSection()') reproduced:
//   A1_cs = 5.23, E1_cs = 9.80, sigma1_cs = 0.05
//   A2_cs = 1.95, E2_cs = 10.14, sigma2_cs = 0.11
// kept as plain (non-const) globals, like in lar_analytic.C, so they can be
// overwritten with a fresh fit result if the Neumeier data is ever refit.
static double emissionA1_cs = 5.23, emissionE1_cs = 9.80, emissionSigma1_cs = 0.05;
static double emissionA2_cs = 1.95, emissionE2_cs = 10.14, emissionSigma2_cs = 0.11;

// emission spectrum vs wavelength [nm] for one time sample: y128/y150/y175/y500
// are the already-time-integrated singlet+triplet, mixed, xenon(175nm), and
// visible/NBrS(500nm) light yields; x is distance [cm], c is dopant
// concentration [ppm] -- both only affect the 128nm line's self-absorption.
// Ported from generate_emission_spectrum() in
// model_with_absorption_and_convolution7.py, matching GenerateEmissionSpectrum()
// in lar_analytic.C.
static std::vector<double> GenerateEmissionSpectrum(const std::vector<double> &wavelengths, double y128,
                                                    double y150, double y175, double y500, double x, double c)
{
  const size_t n = wavelengths.size();
  std::vector<double> shape128(n, 0.0);
  for (size_t i = 0; i < n; ++i)
  {
    const double wl = wavelengths[i];
    if (wl >= 110 && wl <= 900)
    {
      const double em128 = Gaussian(wl, emissionA_E, emissionMu_E, emissionSigma_E);
      const double E = 1239.84 / wl;
      const double cs1 = emissionA1_cs * std::exp(-((E - emissionE1_cs) * (E - emissionE1_cs)) / (2 * emissionSigma1_cs * emissionSigma1_cs));
      const double cs2 = emissionA2_cs * std::exp(-((E - emissionE2_cs) * (E - emissionE2_cs)) / (2 * emissionSigma2_cs * emissionSigma2_cs));
      const double tScaled = std::exp(-(cs1 + cs2) * c * x) * emissionRayleighBaseline;
      shape128[i] = em128 * tScaled;
    }
  }
  const double area128 = Trapezoid(shape128, wavelengths);
  if (area128 > 0)
    for (auto &v : shape128)
      v = v / area128 * y128;

  const double sig150 = 2.8, sig175 = 4.5;
  const double amp150 = y150 / (sig150 * std::sqrt(2 * M_PI));
  const double amp175 = y175 / (sig175 * std::sqrt(2 * M_PI));

  std::vector<double> shapeVis(n);
  for (size_t i = 0; i < n; ++i)
  {
    const double wl = wavelengths[i];
    const double nbrsContinuum = (wl >= 250) ? std::exp(-std::pow((wl - 560.0) / 220.0, 4)) : 0.0;
    const double nirHump = 0.6 * Gaussian(wl, 1.0, 760.0, 55.0);
    shapeVis[i] = nbrsContinuum + nirHump;
  }
  const double areaVis = Trapezoid(shapeVis, wavelengths);
  if (areaVis > 0)
    for (auto &v : shapeVis)
      v = v / areaVis * y500;

  std::vector<double> result(n);
  for (size_t i = 0; i < n; ++i)
    result[i] = shape128[i] + Gaussian(wavelengths[i], amp150, 149, sig150) +
                Gaussian(wavelengths[i], amp175, 175, sig175) + shapeVis[i];
  return result;
}

bool geoVersionOld = false;

static double distanceLevel[NLEVEL];
static void setDistanceLevels()
{
  if (!geoVersionOld) // new geometry if geoVersionOld is false
  {
    distanceLevel[0] = 1.062; // trigger
    distanceLevel[1] = 11.48;
    distanceLevel[2] = 21.41;
    distanceLevel[3] = 31.35;
    distanceLevel[4] = 42.70; // PMT
    distanceLevel[5] = 31.35; // sipm 0
    distanceLevel[6] = 11.48; // sipm 8
  }
  else
  {
    distanceLevel[0] = 1.251; // trigger
    distanceLevel[1] = 11.789;
    distanceLevel[2] = 21.723;
    distanceLevel[3] = 31.668;
    distanceLevel[4] = 42.017; // PMT
    distanceLevel[5] = 31.668;
    distanceLevel[6] = 11.789;
  }
}

static void setParNames() //
{
  lparNames[NORM] = TString("norm");
  lparNames[PPM] = TString("ppm");
  lparNames[TAUT] = TString("tauT");
  lparNames[TAUQUENCH] = TString("tauQuench");
  lparNames[TAUMIXS] = TString("tauMixS");
  lparNames[TAUMIXT] = TString("tauMixT");
  lparNames[TAUMIXR] = TString("tauMixR");
  lparNames[BRAD] = TString("bRad");
  lparNames[TAUQUECHMIX] = TString("tauQenchMix");
  lparNames[THELEVEL] = TString("theLevel");
}
static void setCompNames() //
{
  compNames[SINGLETCOMP] = TString("singletComp");
  compNames[TRIPLETCOMP] = TString("tripletCpmp");
  compNames[XENONCOMP] = TString("xenonComp");
  compNames[MIXEDSINGLETCOMP] = TString("mixedSingletComp");
  compNames[MIXEDTRIPLETCOMP] = TString("mixedTripletComp");
  compNames[VISIBLECOMP] = TString("visibleComp");
}
// populate lpar[NPARS] with default starting values for each fit parameter,
// using the model constants declared above and indexed by moldelFitParameters
static void fillDefaultParameters(double ppm = 0.01, int theLevel = -1)
{
  lpar[NORM] = nPhotons;
  lpar[PPM] = ppm;
  lpar[TAUT] = tautConstant;
  lpar[TAUQUENCH] = tauquenchConstant;
  lpar[TAUMIXS] = taumix_SConstant;
  lpar[TAUMIXT] = taumix_TConstant;
  lpar[TAUMIXR] = taumix_RConstant;
  lpar[BRAD] = b_radConstant;
  lpar[TAUQUECHMIX] = tauquenchmixConstant;
  lpar[THELEVEL] = double(theLevel); // fit all levels by default

  // push the defaults into gMinuit too, if it has been created yet (setupMinuit()
  // runs after setupModelAllFit(), so gMinuit is still null the very first time
  // fillDefaultParameters() is called)
  if (gMinuit)
  {
    int ierflg = 0;
    for (int ip = 0; ip < NPARS; ++ip)
    {
      if (ip == THELEVEL)
      {
        // THELEVEL is a level selector (-1 = fit all, else a level index), always
        // held FIX by the caller right after being set -- not a continuously-varying
        // physical parameter, so it can't use the 0.1x/10x-of-value bound formula
        // below (which inverts for its negative default and breaks for 0).
        step[ip] = 0;
        gMinuit->mnparm(ip, lparNames[ip].Data(), lpar[ip], step[ip], -1.0, double(NLEVEL - 1), ierflg);
      }
      else
      {
        step[ip] = 1.E-6 * lpar[ip];
        double_t arglist[10];
        arglist[0] = ip + 1; // par
        gMinuit->mnparm(ip, lparNames[ip].Data(), lpar[ip], step[ip], 0.1 * lpar[ip], 100. * lpar[ip], ierflg);
        if (fixPar[ip])
          gMinuit->mnexcm("FIX", arglist, 1, ierflg);
      }
    }
  }
}

static void setupModelAllFit()
{
  /* set geomegtry version for modelAllFit,hh */
  // ntScan = new TNtuple("ntScan", "ntScan", "ppm:x:f");
  sensorLevels.resize(NLEVEL);
  sensorLevels[0] = {9, 10, 11}; // trigger
  sensorLevels[1] = {6, 7};      // sipm 6
  sensorLevels[2] = {3, 4, 5};   // sipm 3
  sensorLevels[3] = {1, 2};      // sipm 1
  sensorLevels[4] = {12};        // PMT
  sensorLevels[5] = {0};         // sipm 0
  sensorLevels[6] = {8};         // sipm 8

  geoVersionOld = false;
  setDistanceLevels();
  printf("setParNames and setCompNames\n");
  setParNames();
  setCompNames();

  // covered channels
  for (unsigned ic = 0; ic < NCHAN; ++ic)
    isCovered[ic] = false;

  // covered set to true
  isCovered[0] = true;
  isCovered[8] = true;
  // isCovered[3] = true;
}

// == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
//  MINUIT OPTIMIZER INITIALIZATION
// ============================================================================
static void setupMinuit()
{
  // Create Minuit minimizer with NPARS parameters
  gMinuit = new TMinuit(NPARS);
  gMinuit->SetFCN(fcn); // Set likelihood function pointer
  double_t arglist[10];
  int ierflg = 0;
  arglist[0] = 0.5; // for likelihood up from minimum for 1 sigma errors
  gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
  // setupModelAllFit() (in particular setParNames()) must run before
  // fillDefaultParameter(), since fillDefaultParameters() now also pushes
  // the defaults into gMinuit via mnparm(), which needs lparNames[] filled in
  setupModelAllFit();
  // fillDefaultParameters();
}

// level
static int getLevel(int ichan)
{
  int ilevel = -1; // triggger 9,10,11
  if (ichan == 9 || ichan == 10 || ichan == 11)
    ilevel = 0;
  if (ichan == 6 || ichan == 7)
    ilevel = 1;
  if (ichan == 3 || ichan == 4 || ichan == 5)
    ilevel = 2;
  if (ichan == 1 || ichan == 2)
    ilevel = 3;
  if (ichan == 12)
    ilevel = 4;
  if (ichan == 0)
    ilevel = 5;
  if (ichan == 8)
    ilevel = 6;
  return ilevel;
}

static void sumCurvesToLevels(std::vector<TH1D *> hchan, std::vector<TH1D *> &hsum, int theFileNumber)
{
  // make clones
  hsum.resize(NLEVEL);
  for (int ilevel = 0; ilevel < NLEVEL; ++ilevel)
  {
    hsum[ilevel] = (TH1D *)hchan[0]->Clone(Form("levelCurve%iFile%i", ilevel, theFileNumber));
    hsum[ilevel]->Reset("ICESM"); // reset
  }

  for (int ilevel = 0; ilevel < NLEVEL; ++ilevel)
  {
    /* sum over channels */
    for (int ichan = sensorLevels[ilevel][0]; ichan <= sensorLevels[ilevel].back(); ++ichan)
    {
      for (int ibin = 0; ibin < hchan[ichan]->GetNbinsX(); ++ibin)
      {
        hsum[ilevel]->SetBinContent(ibin, hsum[ilevel]->GetBinContent(ibin) + hchan[ichan]->GetBinContent(ibin));
        double error = sqrt(pow(hsum[ilevel]->GetBinError(ibin), 2.) + pow(hchan[ichan]->GetBinError(ibin), 2.));
        hsum[ilevel]->SetBinError(ibin, error);
      }
      // aftter summing average divide by nlevels to normalize to one sipm
      // printf("line 244 after summing average ilevel %i numberOfSensors %lu histEntries %f histEntriesScaled %f\n", ilevel, sensorLevels[ilevel].size(), hsum[ilevel]->Integral(), hsum[ilevel]->Integral() / double(sensorLevels[ilevel].size()));
    }
  }

  // normalize to single SIPM
  for (int ilevel = 0; ilevel < NLEVEL; ++ilevel)
  {
    hsum[ilevel]->Scale(1.0 / double(sensorLevels[ilevel].size()));
    // printf("line 229 after scaling ilevel %i numberOfSensors %lu integral %f \n", ilevel, sensorLevels[ilevel].size(), hsum[ilevel]->Integral());
  }
}

/**
 * @brief Populate histogram with fitted waveform samples from optimization result
 * @details Transfers fitWave array (computed by fcn() during minimization) into a TH1D
 *          histogram. Data is stored in modelAllFit.hh as fitWave[NCHAN][MAXSAMPLE]
 * @param level level index
 * @param hist Pointer to histogram to fill with fitted values
 */
static void fillFitWave(int ilevel, TH1D *hist)
{
  hist->Reset("ICES");
  for (int ib = 1; ib < hist->GetNbinsX(); ++ib)
  {
    /* debug
    if (ib > int(theTrigStart / 2.) && ib < int(theTrigStart / 2.) + 100) // example condition to check for a specific bin
      printf("line 255 in fillFitWave ichan %i sample %i val %E  \n", ichan, ib, fitWave[ichan][ib]);
      */
    double val = max(fitWave[ilevel][ib], 1.E-9);
    hist->SetBinContent(ib, val);
    hist->SetBinError(ib, 0.);
    hist->GetYaxis()->SetTitle("yield");
    hist->GetXaxis()->SetTitle("time [ns]");
  }
}

/**
 * @brief Populate histogram with fitted component waveform for given channel and component type
 * @details Extracts fitted component contributes (singlet, triplet, mixed) stored in
 *          fitComp[NCHAN][NUMCOMP][MAXSAMPLE] and fills histogram
 * @param ilevel level index
 * @param icomp Component index (e.g., singlet=0, triplet=1, mixed=2)
 * @param hist Pointer to histogram to fill
 */
static void fillCompWave(int ilevel, int icomp, TH1D *hist)
{
  hist->Reset("ICES");
  for (int ib = 1; ib < hist->GetNbinsX(); ++ib)
  {
    double val = max(fitComp[ilevel][icomp][ib], 1.E-9);
    hist->SetBinContent(ib, val);
    // if (icomp == XENONCOMP && ib > 600 && ib < 1000)
    //   printf("!!!! icomp %i chan %i sample %i val %E hist %E \n", icomp, ichan, ib, val, hist->GetBinContent(ib));
    hist->SetBinError(ib, 0);
    hist->GetYaxis()->SetTitle("yield");
    hist->GetXaxis()->SetTitle("time [ns]");
  }
}

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
  // Taken from fits to Neumeier data at 0.1 PPM_MOLAR and scaled;
  ppm = max(1.0E-9, ppm);
  double lambda1 = lambda1ConstDefault * 0.1 / ppm;
  double lambda2 = lambda2ConstDefault * 0.1 / ppm;
  double lambda3 = lambda3ConstDefault * 0.1 / ppm;
  double Tr128 = AConstantDefault * exp(-dist / lambda1) + CConstantDefault * exp(-dist / lambda2) + (1 - AConstantDefault - CConstantDefault) * exp(-dist / lambda3);
  // printf("Absorption: ppm %.3f dist %.3f A %.3E C %.3E lambda1 %.3f lambda2 %.3f lambda3 %.3f Tr128 %.3E \n", ppm, dist, Aconstant, Cconstant , lambda1, lambda2, lambda3, Tr128);
  return 1. - Tr128;
}
/* add labda1Const as argument*/
static double AbsorptionByValue(double ppm, double dist, double lambda1Const = lambda1ConstDefault)
{
  // Calculate absorption as a function of distance and xenon concentration.%
  // Taken from fits to Neumeier data at 0.1 PPM and scaled;
  ppm = max(1.0E-9, ppm);
  double lambda1 = lambda1Const * 0.1 / ppm;
  double lambda2 = lambda2ConstDefault * 0.1 / ppm;
  double lambda3 = lambda3ConstDefault * 0.1 / ppm;
  double Tr128 = AConstantDefault * exp(-dist / lambda1) + CConstantDefault * exp(-dist / lambda2) + (1 - AConstantDefault - CConstantDefault) * exp(-dist / lambda3);
  // printf("Absorption: ppm %.3f dist %.3f A %.3E C %.3E lambda1 %.3f lambda2 %.3f lambda3 %.3f Tr128 %.3E \n", ppm, dist, Aconstant, Cconstant , lambda1, lambda2, lambda3, Tr128);
  return 1. - Tr128;
}
/** exponential convolution with  gaussian time resolution ***/
static double expGaus(double x, double tau)
{
  x -= shift; // compensate for shift in mean due to smearing
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

/*****************************************************************************
 * closed-form solution of the linear LAr/Xe kinetics cascade, ported from
 * lar_analytic.C, matching the lar_system ODE in
 * model_with_absorption_and_convolution7.py
 *
 *   dS  = -lS*S
 *   dT  = -lT*T
 *   dMS = aSS*S + aTS*T - lMS*MS
 *   dMT = aST*S + aTT*T + b*MS - lMT*MT
 *   dX  = g*(MS + MT) - lX*X
 *   dR  = -k_rec*R
 *
 * with
 *   lS   = 1/taus + 1/taudx_S + 1/taudn + 1/tauquench
 *   lT   = 1/taut + 1/taudx_T + 1/taudn + 1/tauquench
 *   lMS  = 1/taumix_S + 1/taumix_R + 1/taudx_X + 1/taudn + 1/tauquenchmix
 *   lMT  = 1/taumix_T + 1/taudx_X + 1/taudn + 1/tauquenchmix
 *   lX   = 1/taux + 1/taudn + 1/tauquench
 *   b    = 1/taumix_R
 *   g    = 1/taudx_X
 *
 * S and T each feed BOTH MS and MT, split by the radiative-cascade branching
 * fraction b_rad (fraction of absorbed 128 nm photons -> MS vs -> MT), via
 * the shared absorbed-photon rate P_abs = Absorb*((1/taus)*S + (1/taut)*T);
 * S and T also each have a separate direct diffusion feed (taudx_S -> MS,
 * taudx_T -> MT) that is NOT split by b_rad:
 *   aSS = b_rad*Absorb/taus + 1/taudx_S        (S -> MS)
 *   aTS = b_rad*Absorb/taut                    (T -> MS)
 *   aST = (1-b_rad)*Absorb/taus                (S -> MT)
 *   aTT = (1-b_rad)*Absorb/taut + 1/taudx_T     (T -> MT)
 *
 * lS and lT carry NO absorption factor: Absorb only reduces the light that
 * escapes detection, it does not change the population decay rate, since
 * photons that get reabsorbed still remove population via the same
 * radiative channel, just redirected into MS/MT via P_abs above.
 *
 * Every component is a sum of chain terms
 *   Chain(l1,...,ln; t) = sum_i exp(-li t) / prod_{j!=i} (lj - li)
 * (the time convolution of exp(-l1 t), ..., exp(-ln t)), with
 *   integral_0^inf Chain dt = 1 / (l1*...*ln).
 *****************************************************************************/
struct LArRates
{
  double lS, lT, lMS, lMT, lX; // total loss rates (1/ns)
  double aSS, aTS, aST, aTT;   // feed rates S->MS, T->MS, S->MT, T->MT (1/ns)
  double b, g;                 // feed rates MS->MT, (MS,MT)->X (1/ns)
  double k;                    // recombination rate (1/ns)
  // characteristic emission rates (1/ns): the fraction of each total loss
  // rate that is the radiative (light-emitting) channel specifically, as
  // opposed to diffusion/quenching losses that do not emit a photon here
  double kEmitS, kEmitT, kEmitMS, kEmitMT, kEmitX; // = 1/taus, 1/taut, 1/taumix_S, 1/taumix_T, 1/taux
};

struct LArState
{
  double S, T, MS, MT, X, R;
};

void printLArState(LArState &state)
{
  printf("LArState: S %E T %E MS %E MT %E X %E R %E \n", state.S, state.T, state.MS, state.MT, state.X, state.R);
}

// 1/x, treating an infinite lifetime as a zero rate
static double LArInv(double x) { return std::isinf(x) ? 0.0 : 1.0 / x; }

static LArRates RatesFromTaus(double taus, double taut, double taudx_S, double taudx_T,
                              double taudx_X, double taudn, double tauquench, double tauquenchmix,
                              double taumix_S, double taumix_T, double taumix_R, double taux,
                              double k_rec, double Absorb, double b_rad)
{
  LArRates r;
  r.lS = LArInv(taus) + LArInv(taudx_S) + LArInv(taudn) + LArInv(tauquench);
  r.lT = LArInv(taut) + LArInv(taudx_T) + LArInv(taudn) + LArInv(tauquench);
  r.lMS = LArInv(taumix_S) + LArInv(taumix_R) + LArInv(taudx_X) + LArInv(taudn) + LArInv(tauquenchmix);
  r.lMT = LArInv(taumix_T) + LArInv(taudx_X) + LArInv(taudn) + LArInv(tauquenchmix);
  r.lX = LArInv(taux) + LArInv(taudn) + LArInv(tauquench);
  r.aSS = b_rad * Absorb * LArInv(taus) + LArInv(taudx_S);
  r.aTS = b_rad * Absorb * LArInv(taut);
  r.aST = (1.0 - b_rad) * Absorb * LArInv(taus);
  r.aTT = (1.0 - b_rad) * Absorb * LArInv(taut) + LArInv(taudx_T);
  r.b = LArInv(taumix_R);
  r.g = LArInv(taudx_X);
  r.k = k_rec;
  r.kEmitS = LArInv(taus);
  r.kEmitT = LArInv(taut);
  r.kEmitMS = LArInv(taumix_S);
  r.kEmitMT = LArInv(taumix_T);
  r.kEmitX = LArInv(taux);
  return r;
}

// Converts populations at a given time into detected light per unit time.
// transFactor is the fraction of S/T's own radiative decay that escapes
// detection (1-Absorb in the main sweep, Tr_detector in the Neumeier
// section); tSinceStart = t - t_start, used to gate the prompt NBrS term
// (zero before the excitation time).
LArState LightComponents(const LArState &y, double transFactor, double tSinceStart)
{

  double N_nbrs_prompt = fitNorm * 0.0005;
  double tau_nbrs = 0.5;
  double f_vis_slow = 0.0025;
  LArState light;
  light.S = transFactor * (1 / tSinglet0) * y.S;
  light.T = transFactor * (1 / tTriplet) * y.T;
  light.MS = (1 / tauMixS) * y.MS;
  light.MT = (1 / tauMixT) * y.MT;
  light.X = (1 / tauxConstant) * y.X;
  const double promptNbrs = (tSinceStart >= 0) ? (N_nbrs_prompt / tau_nbrs) * std::exp(-tSinceStart / tau_nbrs) : 0.0;
  light.R = promptNbrs + k_rec * y.R + f_vis_slow * (1 / tTriplet) * y.T;
  return light;
}

// Convolution of exponentials exp(-l t) for l in rates (rates must be distinct)
// Chain() is the single hottest function in the fit: fcn() calls
// AnalyticSolution() -> Chain() ~15 times per bin, per level, per MIGRAD
// call (roughly 7 levels * 6850 bins * ~1000 MIGRAD calls ~= 70 million
// AnalyticSolution calls, ~1 billion Chain calls). Every call site passes a
// brace-init list, e.g. Chain({r.lS, r.lMS}, t) -- taking that directly as
// std::initializer_list<double> (a view over a small stack-built array) lets
// the compiler skip building it; the old std::vector<double> parameter forced
// a heap allocation/free on every single call just to copy that same data,
// which dominated the runtime.
static double Chain(std::initializer_list<double> rates, double t)
{
  double out = 0.0;
  const double *r = rates.begin();
  const size_t n = rates.size();
  for (size_t i = 0; i < n; ++i)
  {
    double denom = 1.0;
    for (size_t j = 0; j < n; ++j)
      if (j != i)
        denom *= (r[j] - r[i]);
    out += std::exp(-r[i] * t) / denom;
  }

  return out;
}

// Populations at time t, measured from the start of integration
static LArState AnalyticSolution(double t, const LArState &y0, const LArRates &r)
{
  const double S0 = y0.S, T0 = y0.T, MS0 = y0.MS, MT0 = y0.MT, X0 = y0.X, R0 = y0.R;
  LArState y;
  // singlet population convolved with the timing resolution (Chain({r.lS}, t)
  // reduces to exp(-r.lS * t) for a single rate -- expGaus is that same term
  // smeared by the Gaussian resolution, matching calcLightComponents' fs term)
  const double singletShiftNs = -5. * 2.; // -40 samples, 2 ns/sample (bw)
  // y.S = S0 * exp(-r.lS * t);
  y.S = S0 * expGaus(t, 1.0 / r.lS);
  if (t < 0)
  {
    y.T = 0;
    y.R = 0;
    y.MS = 0;
    y.MT = 0;
    y.X = 0;
    return y;
  }
  // if (y.S > 1.E-1)
  //  printf("time %f y.S %E \n", t, y.S);
  y.T = T0 * Chain({r.lT}, t);
  y.R = R0 * Chain({r.k}, t);
  y.MS = MS0 * Chain({r.lMS}, t) + r.aSS * S0 * Chain({r.lS, r.lMS}, t) + r.aTS * T0 * Chain({r.lT, r.lMS}, t);
  y.MT = MT0 * Chain({r.lMT}, t) + r.aST * S0 * Chain({r.lS, r.lMT}, t) + r.aTT * T0 * Chain({r.lT, r.lMT}, t) +
         r.b * MS0 * Chain({r.lMS, r.lMT}, t) + r.b * r.aSS * S0 * Chain({r.lS, r.lMS, r.lMT}, t) +
         r.b * r.aTS * T0 * Chain({r.lT, r.lMS, r.lMT}, t);
  y.X = X0 * Chain({r.lX}, t) +
        r.g * (MS0 * Chain({r.lMS, r.lX}, t) + r.aSS * S0 * Chain({r.lS, r.lMS, r.lX}, t) +
               r.aTS * T0 * Chain({r.lT, r.lMS, r.lX}, t)) +
        r.g * (MT0 * Chain({r.lMT, r.lX}, t) + r.aST * S0 * Chain({r.lS, r.lMT, r.lX}, t) +
               r.aTT * T0 * Chain({r.lT, r.lMT, r.lX}, t) + r.b * MS0 * Chain({r.lMS, r.lMT, r.lX}, t) +
               r.b * r.aSS * S0 * Chain({r.lS, r.lMS, r.lMT, r.lX}, t) +
               r.b * r.aTS * T0 * Chain({r.lT, r.lMS, r.lMT, r.lX}, t));
  return y;
}

// Integrals from 0 to infinity of each population (multiply by kEmit for light yield)
static LArState TimeIntegrals(const LArState &y0, const LArRates &r)
{

  LArState I;
  I.S = y0.S / r.lS;
  I.T = y0.T / r.lT;
  I.MS = (y0.MS + r.aSS * I.S + r.aTS * I.T) / r.lMS;
  I.MT = (y0.MT + r.aST * I.S + r.aTT * I.T + r.b * I.MS) / r.lMT;
  I.X = (y0.X + r.g * (I.MS + I.MT)) / r.lX;
  I.R = y0.R / r.k;
  return I;
}

// Light emitted per unit time from each population: a state decaying with
// total rate l (summed over every loss channel: radiative emission plus any
// diffusion/quenching) emits light only through its own radiative channel,
// at rate kEmit <= l, so light(t) = kEmit * population(t) -- NOT population(t)
// divided by the total lifetime 1/l, which would overcount by including the
// non-radiative loss channels as if they emitted light too.
static LArState LightYield(const LArState &y, const LArRates &r)
{
  LArState light;
  light.S = r.kEmitS * y.S;
  light.T = r.kEmitT * y.T;
  light.MS = r.kEmitMS * y.MS;
  light.MT = r.kEmitMT * y.MT;
  light.X = r.kEmitX * y.X;
  light.R = r.k * y.R; // R decays only via recombination light emission
  return light;
}

/*
fcn is required by Minuit to have exactly these argements
returns likelihood value for some set of parameters
*/
static void printModel()
{

  printf(" \n\n >>> printModel lpars parameters <<< \n");
  for (int ii = 0; ii < NPARS; ++ii)
  {
    printf("\t  param %i %s %.4E  \n", ii, lparNames[ii].Data(), lpar[ii]);
  }

  if (gMinuit)
  {
    printf(" \n >>> printModel gMinuit parameters <<< \n");
    for (int ii = 0; ii < NPARS; ++ii)
    {
      double value = 0, error = 0;
      gMinuit->GetParameter(ii, value, error);
      printf("\t  gMinuit param %i %s value %.4E +/- %.4E \n", ii, lparNames[ii].Data(), value, error);
    }
  }

  // Rates derived from the current lpar[] values, at a representative
  // distance (10 cm) -- same inputs fcn() builds per-level, just evaluated
  // here once for inspection rather than inside the fit loop.
  {
    double ppm = max(1.0E-9, lpar[PPM]);
    double dist = 10.0; // representative distance (cm)
    double ab = Absorption(ppm, dist);
    LArRates r = RatesFromTaus(/*taus*/ tausConstant, /*taut*/ lpar[TAUT], /*taudx_S*/ taudx_SConstant, /*taudx_T*/ taudx_TConstant, /*taudx_X*/ taudx_XConstant, /*taudn*/ taudnConstant, /*tauquench*/ lpar[TAUQUENCH], /*tauquenchmix*/ lpar[TAUQUECHMIX], /*taumix_S*/ lpar[TAUMIXS], /*taumix_T*/ lpar[TAUMIXT], /*taumix_R*/ lpar[TAUMIXR], tauxConstant, /*k_rec*/ k_recConstant, /*Absorb*/ ab, /*b_rad*/ lpar[BRAD]);
    printf(" \n >>> printModel RatesFromTaus (ppm %.4f, dist %.1f cm, Absorb %.4E) <<< \n", ppm, dist, ab);
    printf("\t  lS %.4E  lT %.4E  lMS %.4E  lMT %.4E  lX %.4E\n", r.lS, r.lT, r.lMS, r.lMT, r.lX);
    printf("\t  aSS %.4E  aTS %.4E  aST %.4E  aTT %.4E\n", r.aSS, r.aTS, r.aST, r.aTT);
    printf("\t  b %.4E  g %.4E  k %.4E\n", r.b, r.g, r.k);
  }
}
/**
 * @brief Calculate the singlet/triplet/xenon/mixed/radiative light components at time x
 * @details Computes the model emission components (equation 6) for a given level/channel
 *          at a single time bin, then applies the per-level SiPM/PMT QE scaling (fs, ft,
 *          fx, fm, frad are modified in place -- this mirrors the per-channel branching
 *          used in fcn: level 5 = poly carbon covered, level 6 = glass covered, level 4 = PMT,
 *          all others = ordinary SiPM).
 */
static void calcLightComponents(int ilevel, double x, double ab, double alpha1, double alpha3,
                                double kx, double kxPrime, double tkxPrime,
                                double l1, double l3, double lX, double t1, double t3,
                                double c1, double c3, double tTriplet, double tMix,
                                double norm, double rad2, double rad3, double c1rad,
                                double &fs, double &ft, double &fx, double &fm, double &frad)
{
  // model emission components terms in equation 6
  // convoute with resolution using expGaus
  fs = (1. - ab) * alpha1 / tSinglet0 * expGaus(x, t1); // singlet
  ft = (1. - ab) * alpha3 / tTriplet * expGaus(x, t3);  // triplet

  // radiative component
  double krad3 = rad3 + 1. / tMix;
  double krad23 = rad2 - krad3;
  // here time x is nanoseconds
  double tnano = 2 * x;
  double pRadiative = rad2 * krad3 / krad23 * TMath::Exp(-krad3 * tnano) * (1. - TMath::Exp(-krad23 * tnano));
  pRadiative = max(pRadiative, 1.E-9);
  frad = c1rad * norm * abs(pRadiative);

  // xenon emission x_i terms in paper
  double xterm1 = c1 * kx * alpha1 / (l1 - kxPrime) * ((expGaus(x, tkxPrime) - expGaus(x, tXe0)) / (lX - kxPrime) - (expGaus(x, t1) - expGaus(x, tXe0)) / (lX - l1));
  double xterm3 = c3 * kx * alpha3 / (l3 - kxPrime) * ((expGaus(x, tkxPrime) - expGaus(x, tXe0)) / (lX - kxPrime) - (expGaus(x, t3) - expGaus(x, tXe0)) / (lX - l3));
  // divide by Xenon lifetime from equation 6
  fx = (xterm1 + xterm3) / tXe0;
  //  mixed component
  double mterm1 = alpha1 * c1 / (l1 - kxPrime) * (expGaus(x, tkxPrime) - expGaus(x, t1));
  double mterm3 = alpha3 * c3 / (l3 - kxPrime) * (expGaus(x, tkxPrime) - expGaus(x, t3));
  fm = (mterm1 + mterm3) / tMix; // mixed

  // additional factors depend on SIPM channel
  if (ilevel == 5) // poly carbon covered sees only above 400 nm
  {
    fs = 0;
    ft = 0;
    fm = 0;
    fx = fx * SiPMQE175;
    frad = frad * SiPMQE150;
  }
  else if (ilevel == 6) // glass covered sees only  175
  {
    fs = 0;
    ft = 0;
    fm = fm * SiPMQE150;
    frad = frad * SiPMQE150;
    fx = fx * SiPMQE175;
  }
  else if (ilevel == 4) // PMT
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
}

/*****************************************************************************
 **********************    fcn function fit by minuit  *******************
 ********************** Minuit specifies function arguments
 **********************r eturned f is NLL value    ************************
 **************************************************************************
 */
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // pack parameters into static array onto lightModel
  // for (int k = 0; 8 < npar; ++k)
  //  lpar[k] = par[k];
  f = 0;          // return value
  double bw = 2.; // ns

  // double bkg = par[BKGCONST];

  /******************************** loop over levels ***************************/
  // double chanList[3] = {8, 5, 0};
  for (int ilevel = 0; ilevel < NLEVEL; ++ilevel)
  {
    if (par[THELEVEL] >= 0 && ilevel != par[THELEVEL])
      continue; // only fit the specified level
    /* skip bad channels
    if (isBadChannel(ic))
    {
      printf("fcn skipping bad channel %i \n", ic);
      continue;
    }
    */

    // printf("line 687fcn level %i \n", ilevel);

    /** dist */
    double dist = distanceLevel[ilevel]; // from header distanceLevels.hh
    // SiPMQ128
    double SiPMQ128 = SiPMQE128Ham;

    /****
     * absorption as a function of distance and xenon concentration.%
    double lambda1 = par[LAMBDA1CONST] * 0.1 / ppm;
    double lambda2 = lambda2ConstDefault * 0.1 / ppm;
    double lambda3 = lambda3ConstDefault * 0.1 / ppm;
    double Tr128 = AConstantDefault * exp(-dist / lambda1) + CConstantDefault * exp(-dist / lambda2) + (1 - AConstantDefault - CConstantDefault) * exp(-dist / lambda3);
    double ab = 1. - Tr128;
     ****/
    // after calculating absorption get rates
    double ppm = max(1.0E-9, par[PPM]);
    double ab = Absorption(ppm, dist);
    // printf("line485*********** ilevel %i ppm %.3f ab %.3E *************\n", ilevel, ppm, ab);

    fitNorm = par[NORM];
    tTriplet = par[TAUT];
    tauQuench = par[TAUQUENCH];
    tauMixS = par[TAUMIXS];
    tauMixT = par[TAUMIXT];
    tauMixR = par[TAUMIXR];
    bRad = par[BRAD];
    tauQuenchMix = par[TAUQUECHMIX];

    LArRates larRates = RatesFromTaus(/*taus*/ tausConstant, /*taut*/ tTriplet, /*taudx_S*/ taudx_SConstant, /*taudx_T*/ taudx_TConstant, /*taudx_X*/ taudx_XConstant, /*taudn*/ taudnConstant, /*tauquench*/ tauQuench, /*tauquenchmix*/ tauQuenchMix, /*taumix_S*/ tauMixS, /*taumix_T*/ tauMixT, /*taumix_R*/ tauMixR, tauxConstant, /*k_rec*/ k_recConstant, /*Absorb*/ ab, /*b_rad*/ bRad);

    /* geometric efficiencies
    double fourPi = 2. * TMath::TwoPi();
    double effGeo = pow(0.6, 2.) / fourPi / pow(dist, 2.);
    double aPmt = TMath::Pi() / 4.0 * pow(6.40, 2); // R11410-20  Effective area : 64 mm dia units here are cm
    if (ilevel == 4)                                // PMT
      effGeo = aPmt / fourPi / pow(dist, 2.);
    double effGeoTrig = pow(0.6, 2.) / fourPi / pow(distanceLevel[0], 2.);
    */
    /*
        set to 1 using geo normalized data
    */
    double effGeo = 1.;

    //  MAXSAMPLE; // singlet MAXSAMPLE;
    //   singlet region
    //   iHighFitRange = 1500;
    /****
     ****   loop over bins to fit
     */
    for (int j = iLowFitRange; j < iHighFitRange; ++j) // 7500 is total samples
    {
      /* skip dip region for trigger sipms */
      bool dip = j > 1400 / 2 && j < 1700 / 2;
      // if (dip && ilevel == 0) //
      //   continue;
      /* set sample starting based on TRIGSTART*/
      // nominalTrigger and j are both bin indices (NOT ns); the previous
      // "double(j - int(xTrigger/2.0))" halved xTrigger instead of j-xTrigger,
      // and never multiplied the bin offset by the bin width (bw, declared
      // above but unused) to get ns, so by the real trigger bin the model
      // already thought ~388 ns had elapsed -- many e-foldings past the ~4 ns
      // singlet lifetime, making every light.S/T component numerically zero
      // across the whole fit range.
      double xTrigger = double(nominalTrigger);                                       // in samples/
      const double lightShiftNs = -42. * bw;                                          // shift entire light curve 40 samples in the negative direction
      double x = (double(j) - xTrigger + 0.5) * bw - 7. * tResolution - lightShiftNs; // ns since trigger

      // model emission components (singlet, triplet, mixedSinglet, mixedTriplet, xenon, visible/recomb),
      // computed from the closed-form LAr/Xe kinetics so ution (ported from lar_analytic.C)
      LArState larState0 = {Sfrac * fitNorm, (1.0 - Sfrac) * fitNorm, 0.0, 0.0, 0.0, 36.0};

      /*
      if (x < 0)
      {
        // before the trigger, nothing has happened yet -- AnalyticSolution's
        // Chain() terms diverge for negative t (exp(-rate*t) with rate>0,
        // t<0 blows up), so the pre-trigger state is just all-zero light,
        // not an extrapolation backward through the model.
        for (int icomp = 0; icomp < NUMCOMP; ++icomp)
          fitComp[ilevel][icomp][j] = 0.0;
        fitWave[ilevel][j] = 0.0;
        continue;
      }
        */
      // printLArState(larState0);
      LArState larPop = AnalyticSolution(x, larState0, larRates);
      // printf("x %.3f ns  S %.3E T %.3E MS %.3E MT %.3E X %.3E R %.3E \n", x, larPop.S, larPop.T, larPop.MS, larPop.MT, larPop.X, larPop.R);
      //  printLArState(larPop);

      LArState light = LightComponents(larPop, 1. - ab, x);
      // LArState light = TimeIntegrals(larPop, larRates)s/;

      /**** add efficiencies and QE factors to each component, depending on the level/channel ***/
      light.S *= SiPMQE128Ham;
      light.T *= SiPMQE128Ham;
      light.MS *= SiPMQE150;
      light.MT *= SiPMQE150;

      // additional factors depend on SIPM channel
      if (isCovered[ilevel]) // chan 0 is polycarbonat chan 8 8 glass covered sees only  175
      // multiply by efficiencies
      {
        light.S = 0;
        light.T = 0;
        light.MS = 0;
        light.MT = 0;
        light.R *= SiPMQuartz;
      }
      else if (ilevel == 4) // PMT
      {
        light.X *= PMTQE175;
        light.MS *= PMTQE150;
        light.MT *= PMTQE150;
        light.S = 0;
        light.T = 0;
      }
      /*
      else
      { // all other sipms
        light.S *= SiPMQE175;
      }
      */

      fitComp[ilevel][SINGLETCOMP][j] = light.S;
      fitComp[ilevel][TRIPLETCOMP][j] = light.T;
      fitComp[ilevel][MIXEDSINGLETCOMP][j] = light.MS;
      fitComp[ilevel][MIXEDTRIPLETCOMP][j] = light.MT;
      fitComp[ilevel][XENONCOMP][j] = light.X;
      fitComp[ilevel][VISIBLECOMP][j] = light.R; // R = recombination/visible light channel

      //  output fitted function
      double mval = light.S + light.T + light.MS + light.MT + light.X + light.R;
      fitWave[ilevel][j] = mval;

      // printf("line614 modelFitAll chan %i j %i --  mval %E \n", ilevel, j, mval);

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
      double y = buff[ilevel][j]; // observe)d
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
        printf("line577 ibin YTERM is NAN level %i sample  %i f=%E  y = %E mval %E bkg %E \n", ilevel, j, f, y, mval, lateBkg[j]);
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
        printf("line265  ibin F is NAN level %i sample %i f=%E mval = %E x = %E y = %E yterm %E \n", ilevel, j, f, mval, x, y, yterm);
      }
    }

    // printf("channel modelFitAll channelic, ic %i nLL = %E\n", ic, f);
  } // loop over channels
    // printf("return modelFitAll nLL = %E\n", f);
}
#endif // MODELALLFIT_HH