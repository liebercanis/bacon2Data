// April 28 2025
// distance levels for 13 channels July 16 2026
#include <iostream>
#include <fstream>
#include <numeric>
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TGraph.h"
#include "Math/Vector3D.h"
#include "TBRawRun.hxx"
#include "TBSimRun.hxx"
#include "triggerPeakFit.hh"
#include "TReadGains.hxx"
#include "TMinuit.h"
#include "modelAllFit.hh"
#include "failCodes.hh"

std::string sdate;
// time is in microseconds
using namespace TMath;
TFile *fout;
TRandom3 *ran;
TString geoName;
bool isFid;
double theDopant;

enum
{
  NPOINTS = 1000
};

TReadGains *readGains;

// for writing raw data
TBRawEvent *rawEvent;
TBRawRun *rawRun;
TBSimRun *simRun;
std::vector<uint16_t> wave;
TDirectory *scanDir;
TDirectory *modelDir;

bool writeRawData = true;
bool useMap = false;
int reportInterval = 1000;
double zZero = 0.3; // source position
TDirectory *eventDir;
TDirectory *histDir;
std::vector<TH1D *> hffit;            ///< Model histograms per channel (reserved for future use)
std::vector<TH1D *> hfitModel;        ///< Final fitted model histograms per channel
std::vector<vector<TH1D *>> hfitComp; ///< Final fitted model histograms per channel
std::vector<double> timeComp;
std::vector<std::vector<double>> compIntegral;
std::vector<double> compIntegralSum;
std::vector<int> numCompPhotons;

Double_t vstart[NPARS];      ///< Starting parameter values for Minuit minimization
static Double_t step[NPARS]; ///< Step sizes for Minuit parameter exploration
int theFitChannel = -1;      ///< Channel index to fit. Use -1 to simultaneously fit all 12 PMT channels.

double qsumNominalFromBtb = 3.11E4;

int trigCount9;
TNtuple *ntOrigin;
TNtuple *ntTrigCh;
TNtuple *ntTrig;
TNtuple *ntTern;
TNtuple *ntMean;
TNtuple *ntFit;
TNtuple *ntGammaPeak;
TNtuple *ntNorm;
TNtuple *ntScan;
/* geant maps */
TH3D *originPDF;     // pdf of event origins
TH3D *fluxMapChan9;  // geo efficiency values
TH3D *fluxMapChan10; ///
TH3D *fluxMapChan11; ///
TH2D *hXYMap;
TH3D *hXYZMap;
TH1D *hRadiusMap;
TH1D *hRhoMap;
TH1D *hPhiMap;
TH1D *hZMap;
TH2D *hRhoZMap;
TH3D *hRhoPhiZMap;

TH2D *hXYMapFid;
TH3D *hXYZMapFid;
TH1D *hRadiusMapFid;
TH1D *hRhoMapFid;
TH1D *hPhiMapFid;
TH1D *hZMapFid; ///
TH2D *hRhoZMapFid;
TH3D *hRhoPhiZMapFid;

TH1D *hEffGeo;
TH1D *hEffGeo9;
TH1D *hEffGeo10;
TH1D *hEffGeo11;
TH1D *hPoisson;
TH1D *hPhotonTrig[3];
TH1D *hPhotonTrigSum;

TH1D *hSignalPhotons[3];
TH1D *hNumberSPE[3];
TH1D *hSignalPhotonsEvent[3];
TH1D *hGeoEff[NCHAN];
TH1D *hPhoton[NCHAN];
TH1D *hPhotonTime[NCHAN];
TH1D *hPhotonSum[NCHAN];
TH1D *hSinglet[NCHAN];
TH1D *hConvolve[NCHAN];
TH1D *hSignalNb[NCHAN]; // no baseline
TH1D *hSignal[NCHAN];
TH1D *hSignalSumNoBaseline[NCHAN];
TH1D *hSignalSum[NCHAN];
TH1D *hSignalNorm[NCHAN];
TH1D *hSignalEff[NCHAN];

// for convolve test
TDirectory *testDir;
TNtuple *ntConvolve;
TNtuple *ntConvolveCheck;
TH1D *hSignalSumTest9;
TH1D *hPhotonTest9;
TH1D *hPhotonSumTest9;
TH1D *hConvolveTest9;

Long64_t totalPhotons;
Long64_t ncount[NCHAN];
Long64_t ncountSinglet[NCHAN];
double chanEff[NCHAN];
TH1D *hEventPass;
TH1D *hCount;
TH1D *hCountSinglet;
TH1D *hResponse;
TH1D *hTime;
TH1D *hTrigDiffTime;
TH1D *hTrigDiffTime30;
TH1D *hTrigDiffTimeCut;
TH1D *hGammaPeak;
TH2D *hTriangle;
TH2D *hTriangleCut;
TH2D *hTriangleMean;
// uint16_t maxAdc = pow(2, 14);
double sigmaNoise;
// Peak of the normalized Landau used as the SPE amplitude scale in convolve().
// The locally computed version in btb() is for reporting only; convolve() always uses this global.
double landauMax = 0.018063;
double nominalGeo;
// 2*14         // ns
/* parameters quoted in talk  "A new optical model for LEGEND-200
with remage" Manuel Huber <ge38nap@mytum.de>, Luigi Pertoldi
LEGEND collaboration meeting · March 25, 2025 */
double singletFrac = 0.14; // ; //;0.23;   // Segretto PHYSICAL REVIEW D 103, 043001 (2021)
//  LY is in modelFitAll.hh
double numPhotons = LY * 60.;
//*50. / 34.; // scale to gamma peak in data  60 keV gamma
double fillFactor = 1.0;
double reflection = 1.; //.8; // guess.. angular dependance?
int binWidth = 2;
double noiseToSignal = 0.04;
double baseline = 1100.;    // 1100; // ADC
double meanFreePath = 1.53; // from table in cm3frmom rtabtable in cm3frmom rtabtable in cm
double totalEventEffiency;
double triggerTimes[3];
// Minimum cosθ to avoid shadowing by the source holder: z=0.851 cm (SiPM plane), r=0.4 cm (holder radius).
double cosMin = 0.851 / sqrt(pow(0.4, 2) + pow(0.851, 2));
double maxTriggerTimeDifference = 24.0; // Aug 9
unsigned timeOffset = 13;               // changed from 17 may 13, 2024
double trigTimeShift[3];

double peakQsum[3];

// Note that the actual coordinate origin is centered on the source.
ROOT::Math::XYZVector eventOrigin(0, 0, 0);
// ROOT::Math::XYZVector eventOriginOffset(0, 0, zZero);

// z is positive into array!
// ROOT::Math::XYZVector positionSipm9(1.052, -0.608, 0.851);
// ROOT::Math::XYZVector positionSipm10(-1.052, -0.608, 0.851);
// ROOT::Math::XYZVector positionSipm11(0.000, 1.216, 0.851);

// z is positive into array
ROOT::Math::XYZVector positionSipm9;
ROOT::Math::XYZVector positionSipm10;
ROOT::Math::XYZVector positionSipm11;

TMinuit *gMinuit;
Double_t arglist[1];
int ierflg;
TGraph *gplot1 = NULL;
bool show = false;
bool doMinuit = false;

enum
{
  MAXSCANPLOTS = 10 // interaction point + total photons
};
/*
efficiencies  PMTQE175 = 0.38;
static double QEff128(double ppm, double dist)
*/
int nChannel[NCHAN];
double eff[NCHAN];
int triggerStart = binWidth * 730; // 730; sipm rise time convert to ns
double speMPV = double(triggerStart);
double speSigma = 20.; // ns from single PI data fit
TF1 *speLandau;

/// Plotting colors for each detector channel visualization
int colors[NCHAN] = {kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kSpring, kTeal, kAzure, kViolet, kPink, kGray};

void setupMinuit()
{
  for (int ich = 0; ich < NCHAN; ++ich)
    for (int isamp = 0; isamp < NCHAN; ++isamp)
      buff[ich][isamp] = 0;
  printf("setupMinuit with theDopant %.3f NPARS %i\n", theDopant, NPARS);
  setParNames();
  setCompNames();
  gMinuit = new TMinuit(NPARS);
  gMinuit->SetFCN(fcn); // Set likelihood function pointer
  arglist[0] = 0.5;     // for likelihood up from minimum for 1 sigma errors
  gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

  /* total photons per event LY defined in modelAllFit.hh */
  double startNorm = 60. * LY;
  Double_t arglist[10];
  int ierflg = 0;
  arglist[0] = 0.5; // for likelihood up from minimum for 1 sigma errors
  gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

  // ============================================================================
  //  FIT PARAMETER INITIALIZATION
  // ============================================================================
  // Set initial guesses for all fit parameters
  vstart[NORM] = startNorm;     ///< Photon yield per event
  vstart[TRIGSTART] = 2. * 700; ///< Trigger timing offset
  vstart[SFRAC] = 0.14;         ///< Singlet fraction (ref: Segretto 2021)
  vstart[PPM] = theDopant;      ///< Dopant concentration
  vstart[TAU3] = tTriplet0;     ///< Triplet decay time
  vstart[TAUM] = 4700.0;        ///< Mixed component decay time
  vstart[BKGCONST] = 4.0E-6;    ///< Constant background rate
  vstart[KXCONST] = 1.0;
  vstart[THECHANNEL] = theFitChannel; ///< Channel selection flag
  /*
  printf("starting parameter values \n");
  for (int ip = 0; ip < NPARS; ++ip)
    printf(" par %i %s start val %.3f \n", ip, lparNames[ip].Data(), vstart[ip]);
    */

  // copy into Minuit
  /* have to put some errors here otherwise it will be constant*/
  for (unsigned j = 0; j < NPARS; ++j)
  {
    step[j] = 1.E-6 * vstart[j];
    gMinuit->mnparm(j, lparNames[j].Data(), vstart[j], step[j], 0.1 * vstart[j], 10. * vstart[j], ierflg);
    lpar[j] = vstart[j];
  }
  for (int ip = 0; ip < NPARS; ++ip)
    printf(" par %i %s start val %.3f \n", ip, lparNames[ip].Data(), vstart[ip]);

  // == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
  //  PARAMETER CONSTRAINTS AND BOUNDARIES
  // ============================================================================
  // Fix parameters that are held constant during minimization
  // Note: Minuit uses 1-based indexing for parameters (adds 1 to C++ indices)

  // fix channel
  arglist[0] = THECHANNEL + 1; // channel
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  arglist[0] = TRIGSTART + 1; // trigger
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);
  // arglist[1] = 0.01 * vstart[TRIGSTART]; // low
  // arglist[2] = 10. * vstart[TRIGSTART];  // high
  // gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);

  arglist[0] = BKGCONST + 1; // par
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  arglist[0] = KXCONST + 1; // par
  // gMinuit->mnexcm("FIX", arglist, 1, ierflg);
  arglist[1] = 0.1 * vstart[KXCONST]; // low
  arglist[2] = 10. * vstart[KXCONST]; // high
  gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);

  // arglist[0] = SFRAC + 1; // kp
  // gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  // Set bounds for variable parameters to restrict optimization domain
  arglist[0] = NORM + 1;            // par
  arglist[1] = 0.01 * vstart[NORM]; // low
  arglist[2] = 10. * vstart[NORM];  // high
  gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);
  // gMinuit->mnexcm("FIX", arglist, 3, ierflg);

  arglist[0] = SFRAC + 1;     // par
  arglist[1] = vstart[SFRAC]; // low
  arglist[2] = vstart[SFRAC];
  // gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);
  gMinuit->mnexcm("FIX", arglist, 3, ierflg);

  // set limits ... here par starts with 1 so add 1
  arglist[0] = TAU3 + 1;         // par
  arglist[1] = 0.01 * tTriplet0; // low
  arglist[2] = 2.0 * tTriplet0;  // high
  // gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  // set limits ... here par starts with 1 so add 1
  arglist[0] = PPM + 1; // par
  arglist[1] = 0.0;     // low
  arglist[2] = 50.0;    // high
  gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);
  // gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  arglist[0] = TAUM + 1;     // par tau mixed
  arglist[1] = 0.01 * tMix0; // low
  arglist[2] = 10. * tMix0;  // high
  // gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);
  gMinuit->mnexcm("FIX", arglist, 3, ierflg);

  printf("\n...  call mnprin \n");
  double amin;
  gMinuit->mnprin(1, amin);
  // Evaluate likelihood at starting point to verify initialization
  double fval = 0;
  double gin[NPARS];
  int npar = NPARS;
  int llist = NPARS; ///< Number of parameters
  printModel(700, vstart);
  fcn(llist, gin, fval, vstart, ierflg);
  printf(" starting value >>>>   fval %E \n", fval);
  double fvalStart = fval;
  if (isnan(fval))
  {
    printf("gMinuit returns NAN\n");
    return;
  }
}
/**
 * @brief Populate histogram with fitted waveform samples from optimization result
 * @details Transfers fitWave array (computed by fcn() during minimization) into a TH1D
 *          histogram. Data is stored in modelAllFit.hh as fitWave[NCHAN][MAXSAMPLE]
 * @param ichan Channel index
 * @param hist Pointer to histogram to fill with fitted values
 */
void fillFitWave(int ichan, TH1D *hist)
{
  // std::cout << " fillFitWave " << ichan << "  " << hist->GetName() << std::endl;
  hist->Reset("ICES");
  for (int ib = 1; ib < hist->GetNbinsX(); ++ib)
  {
    double val = max(fitWave[ichan][ib], 1.E-9);
    hist->SetBinContent(ib, val);
    hist->SetBinError(ib, sqrt(val) / 10.);
    hist->GetYaxis()->SetTitle("yield");
    hist->GetXaxis()->SetTitle("time [ns]");
    hffit[ichan] = hist;
  }
}

/**
 * @brief Populate histogram with fitted component waveform for given channel and component type
 * @details Extracts fitted component contributes (singlet, triplet, mixed) stored in
 *          fitComp[NCHAN][NUMCOMP][MAXSAMPLE] and fills histogram
 * @param ichan Channel index
 * @param icomp Component index (e.g., singlet=0, triplet=1, mixed=2)
 * @param hist Pointer to histogram to fill
 */
void fillCompWave(int ichan, int icomp, TH1D *hist)
{
  // std::cout << " fillCompWave " << ichan << " comp  " << icomp << " " << hist->GetName() << std::endl;
  hist->Reset("ICES");
  for (int ib = 1; ib < hist->GetNbinsX(); ++ib)
  {
    double val = max(fitComp[ichan][icomp][ib], 1.E-9);
    hist->SetBinContent(ib, val);
    // if (icomp == XENONCOMP && ib > 600 && ib < 1000)
    //   printf("!!!! icomp %i chan %i sample %i val %E hist %E \n", icomp, ichan, ib, val, hist->GetBinContent(ib));
    hist->SetBinError(ib, 0);
    hist->GetYaxis()->SetTitle("yield");
    hist->GetXaxis()->SetTitle("time [ns]");
  }
  std::cout << std::endl;
}
void fillModel()
{
  modelDir->cd();
  // Initialize histogram vectors for all channels
  hffit.resize(NCHAN);
  hfitModel.resize(NCHAN);
  hfitComp.resize(NCHAN);
  for (unsigned ic = 0; ic < NCHANPMT; ++ic)
  {
    // printf("fill fitWaveFitChan%i", ic);
    TH1D *hFit = new TH1D(Form("fitWaveFitChan%i", ic), Form("fitWaveFitChan%i", ic), MAXSAMPLE, 0, 2 * MAXSAMPLE);
    hFit->Reset("ICES");
    hFit->SetTitle((Form("fitWaveFitChan%i %.3fPPM", ic, theDopant)));
    hFit->SetMarkerColor(colors[ic]);
    hFit->SetLineColor(colors[ic]);
    hfitModel[ic] = hFit;
    fillFitWave(ic, hFit);
  }

  // drawing
  // Store individual component contributions (singlet, triplet, mixed) for each channel
  for (unsigned ic = 0; ic < NCHAN; ++ic)
  {
    for (int icomp = 0; icomp < NUMCOMP; ++icomp)
    {
      /*
       TString tprint;
      tprint.Form("fillCompWaveFitChan%iComp%i name %s", ic, icomp, compNames[icomp].Data());
      cout << tprint << endl;
      */
      TH1D *hFit = new TH1D(Form("fitComp%sChan%i", compNames[icomp].Data(), ic), Form("fit%sChan%i", compNames[icomp].Data(), ic), MAXSAMPLE, 0, 2 * MAXSAMPLE);
      hFit->Reset("ICES");
      hFit->SetTitle((Form("fit%sChan%i %.3fPPM", compNames[icomp].Data(), ic, theDopant)));
      hFit->SetLineColor(colors[ic]);
      fillCompWave(ic, icomp, hFit); // only need one of these
      hfitComp[ic].push_back(hFit);
    }
  }

  compIntegral.resize(NCHAN);
  compIntegralSum.clear();
  compIntegralSum.resize(NCHAN);
  // calculate comp sizes
  for (unsigned ic = 0; ic < NCHAN; ++ic)
  {
    for (int icomp = 0; icomp < hfitComp[ic].size(); ++icomp)
    {
      compIntegral[ic].push_back(hfitComp[ic][icomp]->Integral());
      compIntegralSum[ic] += compIntegral[ic].back();
    }
  }

  if (show)
  {
    printf("MESSAGE ******** compIntegrals ******* \n");
    for (unsigned ic = 0; ic < NCHAN; ++ic)
    {
      for (int icomp = 0; icomp < hfitComp[ic].size(); ++icomp)
      {
        printf("chan %i comp %i compIntegral %.3E\n", ic, icomp, compIntegral[ic][icomp]);
      }
      printf("chan %i compIntegralSum %.3E\n\n", ic, compIntegralSum[ic]);
    }
  }
  // end of function
  fout->cd();
}

/* get times for sipm channel */
void getTime(int ic, int icomp, int nPhotons)
{
  timeComp.clear();
  for (int i = 0; i < nPhotons; ++i)
  {
    timeComp.push_back(hfitComp[ic][icomp]->GetRandom());
  }
  return;
}

double gainFunc(int ich)
{
  double g = readGains->sipmPeakGain[ich];
  return g;
}

// get geant4 map in dir bobj
bool getMap()
{
  bool rc = false;
  originPDF = NULL;
  fluxMapChan9 = NULL;
  fluxMapChan10 = NULL;
  fluxMapChan11 = NULL;
  TString mapFileName = TString(getenv("BOBJ")) + TString("/geantSim-2025-07-16-15-22-0.root");
  TFile *fmap = new TFile(mapFileName, "readonly");
  fmap->GetObject("OriginMap", originPDF);
  fmap->GetObject("FluxMapChan9", fluxMapChan9);
  fmap->GetObject("FluxMapChan10", fluxMapChan10);
  fmap->GetObject("FluxMapChan11", fluxMapChan11);
  if (originPDF && fluxMapChan9 && fluxMapChan10 && fluxMapChan11)
  {
    printf("getMap found histograms %s %s %s %s \n", originPDF->GetName(), fluxMapChan9->GetName(), fluxMapChan10->GetName(), fluxMapChan11->GetName());
    rc = true;
  }

  return rc;
}

// trigger condition return maximum time between 3 SIPM first photons
double eventTrigger(int iev)
{
  // printf("line104 %0.f %0.f %0.f \n", hPhoton[9]->GetEntries(), hPhoton[10]->GetEntries(), hPhoton[11]->GetEntries());
  // number of samples is 7500 each bin is 2 ns
  // hPhoton x-axis is in ns

  double tdiff = double(2 * MAXSAMPLE); // better to put in overflow
  // printf("line204 eventTrigger photon integrals 9) %.3f %.3f 10) %.3f %.3f 11) %.3f  %.3f\n", hPhoton[9]->GetEntries(), hSignalPhotonsEvent[0]->Integral(),
  //        hPhoton[10]->GetEntries(), hSignalPhotonsEvent[1]->Integral(), hPhoton[11]->GetEntries(), hSignalPhotonsEvent[2]->Integral());
  //  all must have at least 1 photon
  /*
  if (hSignalPhotonsEvent[0]->GetEntries() < 1)
    return tdiff;
  if (hSignalPhotonsEvent[1]->GetEntries() < 1)
    return tdiff;
  if (hSignalPhotonsEvent[2]->GetEntries() < 1)
    return tdiff;
  */

  // TH1D *hclone9 = (TH1D *)hSignalPhotonsEvent[0]->Clone();

  // all must have at least 1 photon
  double prompt9 = hPhoton[9]->Integral(730, 745) / gainFunc(9);
  double prompt10 = hPhoton[10]->Integral(730, 745) / gainFunc(10);
  double prompt11 = hPhoton[11]->Integral(730, 745) / gainFunc(11);

  hNumberSPE[0]->Fill(hPhoton[9]->GetEntries());
  hNumberSPE[1]->Fill(hPhoton[10]->GetEntries());
  hNumberSPE[2]->Fill(hPhoton[11]->GetEntries());

  bool canTrig = true;
  if (prompt9 < 1)
    canTrig = false;
  if (prompt10 < 1)
    canTrig = false;
  if (prompt11 < 1)
    canTrig = false;

  // if (!canTrig)
  // printf("line236 event %i %.0f %.0f %.0f\n", iev, hPhoton[9]->GetEntries(), hPhoton[10]->GetEntries(), hPhoton[11]->GetEntries());

  if (!canTrig)
  {
    if (show)
      printf("eventTrigger fails event %i %f %f %f \n", iev, prompt9, prompt10, prompt11);
    // printf("event %i cannot trig \n", iev);
    return tdiff;
  }

  /*
  if (hPhoton[9]->GetEntries() < 1)
    return tdiff;
  if (hPhoton[10]->GetEntries() < 1)
    return tdiff;
  if (hPhoton[11]->GetEntries() < 1)
    return tdiff;
    */

  ++trigCount9;

  // collect first times hPhoton x-axis each bin is 2ns
  std::vector<double> ftimes;

  for (int isipm = 9; isipm < 12; ++isipm)
  {
    /* using convolution waveform
    for (int ibin = 1; ibin < hSignalPhotonsEvent[isipm]->GetNbinsX(); ++ibin)
    {
      if (hSignalPhotonsEvent[isipm]->GetBinContent(ibin) > 1)
        printf("chan %i bin %i val %f \n", isipm, ibin, hSignalPhotonsEvent[isipm]->GetBinContent(ibin));
    }
       */

    /* based on photon arrival */
    for (int ibin = 730; ibin < 746; ++ibin)
    {
      if (hPhoton[isipm]->GetBinContent(ibin) > 0)
      {
        // ftimes.push_back(hPhoton[isipm]->GetBinLowEdge(ibin));
        ftimes.push_back(hPhoton[isipm]->GetBinCenter(ibin));
        break;
      }
    }
    /* alternative logic
    int nph = 0;
    for (int ibin = 1; ibin < hPhoton[isipm]->GetNbinsX(); ++ibin)
    {
      if (hPhoton[isipm]->GetBinContent(ibin) > 0)
        ++nph;
      if (nph > 0) // number of photons to trigger
      {
        ftimes.push_back(hPhoton[isipm]->GetBinCenter(ibin)); // the time this photon is detected
        break;
      }
    }
    */
  }
  if (ftimes.size() < 3)
  {
    printf("event %i ftimes %ld  trig \n", iev, ftimes.size());
    return tdiff;
  }
  // this event can trigger
  hEventPass->SetBinContent(3, hEventPass->GetBinContent(3) + 1);

  // Sort in ascending order (default)
  std::sort(ftimes.begin(), ftimes.end());
  tdiff = ftimes[2] - ftimes[0];

  // printf("line189 %i (%.0f %.0f %.0f)  tdiff %.0f \n", iev, ftimes[0], ftimes[1], ftimes[2], tdiff);
  //   convert to ns
  return tdiff;
}

//// https://mathworld.wolfram.com/TernaryDiagram.html
void makeTernary(double a, double b, double c, double &x, double &y)
{
  double s = a + b + c;
  x = 0.5 * (a + 2. * b) / s;
  y = sqrt(3.) / 2. * a / s;
}

TGraph *myScan(int thePar, double xlow, double xhigh)
{
  int maxPoints = 100;
  std::vector<double> xval;
  std::vector<double> yval;

  // get min parameters
  double fitVal[NPARS];
  double fitErr[NPARS];
  for (int ipar = 0; ipar < NPARS; ++ipar)
  {
    gMinuit->GetParameter(ipar, fitVal[ipar], fitErr[ipar]);
    // printf("line141 ipar %i par %f err %f \n", ipar, fitVal[ipar], fitErr[ipar]);
  }

  double *fGin;
  double nLL;
  for (int i = 0; i < maxPoints; ++i)
  {
    double x = xlow + double(i) * (xhigh - xlow) / double(maxPoints);
    fitVal[thePar] = x;
    gMinuit->Eval(thePar, fGin, nLL, &fitVal[0], 4);
    xval.push_back(x);
    if (nLL >= 100)
      nLL = 0.0; // for plotting
    yval.push_back(nLL);
    // printf("line53 i %i par nph %f r %f theta %f  phi %f \n", i, fitVal[0], fitVal[1], fitVal[2], fitVal[3]);
    ntScan->Fill(nLL, peakMeanQsum[0], peakMeanQsum[1], peakMeanQsum[2], peakFitQsum[0], peakFitQsum[1], peakFitQsum[2], fitVal[1], fitVal[2], fitVal[3]);
    // printf("mySCAN %i %f %f \n", i, x, nLL);
  }
  // make and return graph
  return new TGraph(maxPoints, &xval[0], &yval[0]);
}

// spe landau shape
static double myLandau(Double_t *xx, Double_t *par)
{
  double x = xx[0];
  return par[2] * TMath::Landau(x, par[0], par[1], true); // normalized
}

void convolve(TH1D *hist, double time, double gain) // time is when photon arrives
{
  int offsetBin = 728; // max bin of hResponse read offf of Response histogram
  int startBin = hist->FindBin(time);
  // hist->SetBinContent(startBin, hist->GetBinContent(startBin) + gain);
  for (int ib = startBin; ib < hist->GetNbinsX(); ++ib)
    hist->SetBinContent(ib, hist->GetBinContent(ib) + gain / landauMax * hResponse->GetBinContent(ib - startBin + offsetBin));
  // printf("convolve: gain %f landauMax %f hist %s time %f startBin %i offsetBin %i response %.0f integral %.0f \n", gain, landauMax, hist->GetName(), time, startBin, offsetBin, hResponse->Integral(), hist->Integral());
}

void convolveTest(int ntries)
{
  // loop over number of photons per event
  for (int i = 0; i < 20; ++i)
  {
    int nsinglet = i + 1;
    // reset histograms for this nsinglet number
    double timeShift = timeOffset + trigTimeShift[0];
    double singletTime;
    double singletTimeAve;
    // generate singlet
    // loop over trials
    printf("convolveTest nsinglet %i \n", nsinglet);
    for (int k = 0; k < ntries; ++k)
    {
      singletTimeAve = 0;
      // reset event histograms
      hPhotonTest9->Reset("ICESM");
      hPhotonSumTest9->Reset("ICESM");
      hConvolveTest9->Reset("ICESM");
      hSignalSumTest9->Reset("ICESM");
      double gain9 = readGains->sipmPeakGain[9];
      for (int j = 0; j < nsinglet; ++j)
      {
        singletTime = ran->Exp(tSinglet0);
        singletTimeAve += singletTime;
        double time = timeShift + triggerStart + singletTime;
        hPhotonTest9->Fill(time);
        hPhotonSumTest9->Fill(time);
        convolve(hConvolveTest9, time, gain9);
      } // singlet loop
      singletTimeAve /= double(ntries);
      // add noise and sum convolution
      // printf("nsinglet %i photon sum %.0f signal integral %.0f \n", i, hPhotonSumTest9->Integral(), hConvolveTest9->Integral());
      for (int ibin = 1; ibin <= hSignalSumTest9->GetNbinsX(); ++ibin)
      {
        double binNoise = ran->Gaus(0.0, sigmaNoise);
        hSignalSumTest9->SetBinContent(ibin, binNoise + hConvolveTest9->GetBinContent(ibin) + hSignalSumTest9->GetBinContent(ibin));
      } // photon generation loop
      // printf("nsinglet %i trial %i nphotons %.0f signal integral %.0f \n", i, k, hPhotonSumTest9->Integral(), hSignalSumTest9->Integral());
      ntConvolve->Fill(singletTimeAve, hPhotonSumTest9->GetEntries(), hConvolveTest9->Integral() / qsumNominalFromBtb, hSignalSumTest9->Integral() / qsumNominalFromBtb);
    } // tries
  } // singlet value
}

ROOT::Math::XYZVector getXYZVector(double r, double theta, double phi) // angles in radians
{
  // construct XYZ coordinates
  double cost = cos(theta);
  double sint = sin(theta);
  double cosp = cos(phi);
  double sinp = sin(phi);
  double x = r * sint * cosp;
  double y = r * sint * sinp;
  double z = r * cost;
  ROOT::Math::XYZVector pos(x, y, z);
  return pos;
}

double effGeoSim(int ichan) // uses PositionVector3D eventOrigin;
{
  bool isTrig = ichan == 9 || ichan == 10 || ichan == 11;
  if (!isTrig)
    return effGeoFunc(ichan);

  int ilevel = getLevel(ichan);
  // double e = 1.0;
  // if (ilevel != 0)
  //  return e;

  double e = effGeoFunc(ichan);
  if (eventOrigin.R() == 0.)
    return e;

  /* correct for gamma interaction position */
  /* get from map */
  if (useMap)
  {
    int xbin = originPDF->GetXaxis()->FindBin(eventOrigin.R());
    int ybin = originPDF->GetYaxis()->FindBin(cos(eventOrigin.Theta()));
    int zbin = originPDF->GetZaxis()->FindBin(eventOrigin.Phi());
    // int globalBin = originPDF->GetBin(eventOrigin.R(), eventOrigin.Theta(), eventOrigin.Phi());
    if (ichan == 9)
      e = fluxMapChan9->GetBinContent(xbin, ybin, zbin);
    else if (ichan == 10)
      e = fluxMapChan10->GetBinContent(xbin, ybin, zbin);
    else if (ichan == 11)
      e = fluxMapChan11->GetBinContent(xbin, ybin, zbin);
    // printf("Bins (%i,%i,%i) at %f %f %f  effGeo %E \n", xbin, ybin, zbin, eventOrigin.R(), eventOrigin.Theta(), eventOrigin.Phi(), e);
  }
  else
  {
    /* this was old method */
    /* det positions Georgia May 2025 */
    double trigRadius = positionSipm9.R();
    // double trigRadius = 1.;
    //   convert to radians
    double trigTheta = positionSipm9.Theta(); // 11,10,9
    // phi are different
    double trigPhi[3];
    trigPhi[0] = positionSipm9.Phi();  // 9
    trigPhi[1] = positionSipm10.Phi(); // 9
    trigPhi[2] = positionSipm11.Phi(); // 9

    ROOT::Math::XYZVector rSipm = getXYZVector(trigRadius, trigTheta, trigPhi[ichan - 9]);
    ROOT::Math::XYZVector relative = rSipm - eventOrigin;

    double distance2 = relative.Mag2();
    // Area of SiPMs is 6.0mm x 6.0mm
    //    Channels 6, 7, and 8 are at 11.6 cm
    //    from the source Channels 3, 4, and 5 are at 23.2 cm
    //    from the source Channels 0, 1, and 2 are at 34.8 cm from the source Channel 12 is at 36 cm from the source.
    double aPmt = TMath::Pi() / 4.0 * pow(6.4, 2); // R11410-20  Effective area : 64 mm dia
    double a = pow(0.6, 2.);                       // SIPM area

    // correct solid angle
    ROOT::Math::XYZVector runit = rSipm.Unit();
    ROOT::Math::XYZVector ounit = relative.Unit();
    double cos = runit.Dot(ounit);
    // printf(" chan %i cos %f \n", ichan, cos);
    if (cos < 0.)
      cos = 0; // origin behind sipm
    a *= cos;

    e = a / distance2 / (4.0 * TMath::Pi());
    // shift phi for printing.
    double localPhi = relative.Phi() * 360. / TMath::TwoPi();
    if (localPhi < 0)
      localPhi += 360.;
    if (0)
      printf("effGeoSim ichan  %i level %i  cos %f dist %f R,Theta,Phi (%f,%f,%f) area %f geo eff %f \n", ichan, ilevel, cos, sqrt(distance2), relative.R(), relative.Theta() * 360. / TMath::TwoPi(), localPhi, a, e);
  }
  return e;
}

void btb(int ngen = 10000000, double thePPM = 30.)
{
  /** use nominal gains  **/
  readGains = new TReadGains(false);
  for (int il = 0; il < NCHAN; ++il)
    nChannel[il] = 0;
  // zerp trig coount
  trigCount9 = 0;

  // trigger time shifts — only the last assignment (+5 ns) takes effect; earlier lines are dead
  trigTimeShift[0] = 0.;
  trigTimeShift[0] = -5.;
  trigTimeShift[0] = +5.;

  geoVersionOld = false;
  setDistanceLevels(geoVersionOld);
  if (geoVersionOld)
    geoName = TString("OLD");
  else
    geoName = TString("NEW");

  if (geoVersionOld)
  {
    // z is positive into array set XYZ coordinates
    positionSipm9.SetCoordinates(0.795, -0.459, 0.851);
    positionSipm10.SetCoordinates(-0.795, -0.459, 0.851);
    positionSipm11.SetCoordinates(0.000, 0.918, 0.851);
    printf(" btb sim OLD geometry NOMAP generate ngen =  %i LY %.1f photons/kev * 60 = %.1f nominalGain %f nominalTrigGain %f \n", ngen, LY, numPhotons, readGains->sipmPeakGain[8], readGains->sipmPeakGain[9]);
  }
  else
  {
    // z is positive into array
    /* chan 9 1.062 0.795 -0.459 -0.534
      chan 10 1.062 -0.795 -0.459 -0.534
      chan 11 1.062 0.000 0.918 -0.534
    */
    positionSipm9.SetCoordinates(0.795, -0.459, 0.534);
    positionSipm10.SetCoordinates(-0.795, -0.459, 0.534);
    positionSipm11.SetCoordinates(0.000, 0.918, 0.534);
    printf(" btb sim NEW geometry NOMAP generate ngen =  %i LY %.1f photons/kev * 60 = %.1f nominalGain %f nominalTrigGain %f \n", ngen, LY, numPhotons, readGains->sipmPeakGain[8], readGains->sipmPeakGain[9]);
  }

  printf("level distances 0 = %.3f 1= %.3f 2= %.3f 3 %.3f 4 %.3f \n", distanceLevel[0], distanceLevel[1], distanceLevel[2], distanceLevel[3], distanceLevel[4]);

  printf(" trigger sipm positions :  \n");
  printf(" \t sipm 9 : x %.3f y %.3f z %.3f R %.3f rho %.3f phi %.3f  \n",
         positionSipm9.X(), positionSipm9.Y(), positionSipm9.Z(),
         positionSipm9.R(), positionSipm9.Rho(), positionSipm9.Phi() * 180 / TMath::Pi());

  printf(" \t sipm 10 : x %.3f y %.3f z %.3f R %.3f rho %.3f phi %.3f  \n",
         positionSipm10.X(), positionSipm10.Y(), positionSipm10.Z(),
         positionSipm10.R(), positionSipm10.Rho(), positionSipm10.Phi() * 180 / TMath::Pi());

  printf(" \t sipm 11: x %.3f y %.3f z %.3f R %.3f rho %.3f phi %.3f \n",
         positionSipm11.X(), positionSipm11.Y(), positionSipm11.Z(),
         positionSipm11.R(), positionSipm11.Rho(), positionSipm11.Phi() * 180 / TMath::Pi());

  if (useMap)
  {
    if (!getMap())
    {
      printf("no geant4 map\n");
      exit(0);
    }
  }

  /* nominal geo with SIPM at with R=0 origin detector origin */
  double trigDistanceR = positionSipm9.R();
  nominalGeo = pow(0.6, 2.) / pow(trigDistanceR, 2.) / (4.0 * TMath::Pi());
  double nPhotonsNominal = nPhotons * pow(nominalGeo * SiPMQE128Ham, 3.0);
  printf("\n ******************* nominal photon %.0f trig R %.3f nominal goem eff %.3f sum on 3 sipms %.2f ********************\n",
         nPhotons, trigDistanceR, nominalGeo, nPhotonsNominal);

  /* channel efficiences */
  double nominalGeo9;
  for (int i = 0; i < NCHAN - 1; ++i)
  {
    double effGeoSimi = effGeoSim(i);
    double trigRadius = distanceLevel[getLevel(i)];
    nominalGeo = pow(0.6, 2.) / pow(trigRadius, 2.) / (4.0 * TMath::Pi());
    if (i == 9)
      nominalGeo9 = nominalGeo;
    printf("chan %i eventOrgin R %.3f nominal effGeoSim %E (%E)  ratio %E \n", i, eventOrigin.R(), effGeoFunc(i), nominalGeo, nominalGeo9 / nominalGeo);
    eff[i] = effGeoSimi * SiPMQE128Ham * fillFactor;
    // double dist = distanceLevel[level(i)];
  }
  eff[NCHAN - 1] = 0.0; // pmt sees light > 175 nm

  printf(" efficiences at 128 nm PPM=%f\n", thePPM);
  for (int i = 0; i < NCHAN; ++i)
    printf("\t chan %i eff %f \n", i, eff[i]);

  ran = new TRandom3();
  // open raw output file/raw/
  rawRun = NULL;

  time_t rawtime;
  struct tm *timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);

  int day = timeinfo->tm_mday;         // 1-31, as-is
  int month = timeinfo->tm_mon + 1;    // tm_mon is 0-11, so add 1
  int year = timeinfo->tm_year + 1900; // tm_year is years since 1900
  // correct format is month day year
  // char output[30];
  // strftime(output, 30, "%Y-%m-%d-%H-%M", timeinfo);
  // TString tdateTag = TString(output);
  TString tdateTag;
  tdateTag.Form("%02i_%02i_%i", month, day, year);
  TString fullname;
  if (geoVersionOld)
    fullname = (Form("btbSimOLD-%s-%i.root", tdateTag.Data(), ngen));
  else
    fullname = (Form("btbSimNEW-run-%s-%i-PPM-%05i.root", tdateTag.Data(), ngen, int(thePPM * 1000.)));
  fout = new TFile(fullname, "recreate"); // DEF made to update rather than recreate so that it doesn't write over a file already made.
  modelDir = fout->mkdir("modelDir");
  fout->cd();
  theDopant = thePPM;
  printf("opened output file %s date %s dopant %.3f\n", fout->GetName(), tdateTag.Data(), theDopant);
  printf("absorbtion factor %.3f PPM \n", theDopant);

  setupMinuit();
  fillModel();

  /*
  for (unsigned ic = 0; ic < NCHANPMT; ++ic)
    for (unsigned icomp = 0; icomp < hfitComp[ic].size(); ++icomp)
      printf("ic %i icomp %i %s \n", ic, icomp, hfitComp[ic][icomp]->GetName());
      */

  cout << "tdateTag = " << tdateTag << endl;

  std::vector<double> distance;
  std::vector<double> abdist;
  distance.resize(NPOINTS);
  abdist.resize(NPOINTS);
  for (int i = 0; i < NPOINTS; ++i)
  {
    distance[i] = double(i) * 0.05;
    abdist[i] = Absorbtion(theDopant, distance[i]);
    // printf("%i ppm %f dist %f A %.3E\n", i, theDopant, distance[i], abdist[i]);
  }

  TGraph *gAbDist = new TGraph(NPOINTS, &distance[0], &abdist[0]);
  gAbDist->SetMarkerStyle(21);
  gAbDist->SetMarkerSize(.7);
  gAbDist->SetName("absorbtion-distance");
  gAbDist->SetTitle(Form("absorbtion factor at %.3f PPM", theDopant));
  gAbDist->GetXaxis()->SetTitle("distance [cm]");
  gAbDist->GetYaxis()->SetTitle("absorbtion factor");
  TCanvas *cabDist = new TCanvas("absorbtionDist", "absorbtionDist");
  cabDist->SetGrid();
  cabDist->SetLogx();
  gAbDist->Draw("ap");
  fout->Append(gAbDist);

  scanDir = fout->mkdir("scanDir");

  // make output tree
  simRun = new TBSimRun("sim0");
  simRun->clear();

  if (writeRawData)
  {
    rawRun = new TBRawRun(tdateTag);
    rawRun->updateTime(rawtime);
    rawRun->btree->SetTitle("simulation");
    // rawRun->print();a
  }
  ntConvolveCheck = new TNtuple("ntConvolveCheck", "ntConvolveCheck", "chan:nsinglet:ntriplet:nph:qconv:qsum");
  ntNorm = new TNtuple("ntNorm", "norm and efficieincy", "ch:eff:nph:signal:signalEff:signalNorm");
  ntOrigin = new TNtuple("ntOrigin", " event origin ", "ev:r:cos:theta:phi:x:y:z");
  ntTrig = new TNtuple("ntTrig", " trigger info by event  ", "ev:trigSum:nph9:nph10:nph11:r:rho:phi:x:y:z:tdiff:pass");
  ntTrigCh = new TNtuple("ntTrigCh", " trigger info by channel ", "ev:ch:qsum:psum:nph:rho:phi:x:y:z:effgeo");
  ntTern = new TNtuple("ntTern", " trigger sipm", "eventRho:eventPhi:eventZ:nph9:nph10:nph11:qsum9:qsum10:qsum11:mean9:mean10:mean11:xq:yq");
  ntMean = new TNtuple("ntMean", "trigger means ", "ev:eventRho:eventPhi:eventZ:qsum9:qsum10:qsum11:mean9:mean10:mean11:xternq:yternq");
  ntFit = new TNtuple("ntFit", "trigger peak fit", "ev:numPhotons:eventR:eventCos:eventPhi:qsum9:qsum10:qsum11:fitR:fitCos:fitPhi:errR:errTheta:ierr");
  ntScan = new TNtuple("ntScan", "scan", "nll:mean9:mean10:mean11:qsum9:qsum10:qsum11:r:theta:phi");
  ntGammaPeak = new TNtuple("ntGammaPeak", "nt gamma peak", "ev:nph:ph9:ph10:ph11:nph9:nph10:nph11:sum");

  hPoisson = new TH1D("Poisson", " total photons in event  ", 200, 0, 200.);
  hPhotonTrig[0] = new TH1D("PhotonTrig9", " total sipm 9 trigger photons in event  ", MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
  hPhotonTrig[1] = new TH1D("PhotonTrig10", " total sipm 10 trigger photons in event  ", MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
  hPhotonTrig[2] = new TH1D("PhotonTrig11", " total sipm 11 trigger photons in event  ", MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
  hSignalPhotons[0] = new TH1D("SignalPhotons9", "photons/event channel 9", MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
  hSignalPhotons[1] = new TH1D("SignalPhotons10", "photons/event channel 10", MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
  hSignalPhotons[2] = new TH1D("SignalPhotons11", "photons/event channel 11", MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
  hNumberSPE[0] = new TH1D("NumberSPE9", "photons/event channel 9", 50, 0, 50);
  hNumberSPE[1] = new TH1D("NumberSPE10", "photons/event channel 10", 50, 0, 50);
  hNumberSPE[2] = new TH1D("NumberSPE11", "photons/event channel 11", 50, 0, 50);

  hSignalPhotonsEvent[0] = new TH1D("SignalPhotonsEvent9", "photons in event  channel 9", MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
  hSignalPhotonsEvent[1] = new TH1D("SignalPhotonsEvent10", "photons in event channel 10", MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
  hSignalPhotonsEvent[2] = new TH1D("SignalPhotonsEvent11", "photons in event channel 11", MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
  hSignalPhotonsEvent[0]->SetDirectory(nullptr);
  hSignalPhotonsEvent[1]->SetDirectory(nullptr);
  hSignalPhotonsEvent[2]->SetDirectory(nullptr);

  hEffGeo = new TH1D("EffGeo", " geometric efficiency / nominal ", 150, 0, 1.5);
  hEffGeo9 = new TH1D("EffGeo9", " ch 9 geometric efficiency / nominal ", 150, 0, 1.5);
  hEffGeo10 = new TH1D("EffGeo10", " ch 10 geometric efficiency / nominal ", 150, 0, 1.5);
  hEffGeo11 = new TH1D("EffGeo11", " ch 11 geometric efficiency / nominal ", 150, 0, 1.5);

  /* define ntuples amd histograms here */
  hXYMap = new TH2D("XYMap", "event Y versus X  [cm] ", 100, -4., 4., 100, -4., 4.);
  hXYMap->GetXaxis()->SetTitle("event X [cm]");
  hXYMap->GetYaxis()->SetTitle("event Y [cm]");

  hXYZMap = new TH3D("XYZMap", "event x y z  [cm] ", 100, -4., 4., 100, -4., 4., 50, 0, 4.);
  hXYZMap->GetXaxis()->SetTitle("event X [cm]");
  hXYZMap->GetYaxis()->SetTitle("event Y [cm]");
  hXYZMap->GetZaxis()->SetTitle("event Z [cm]");
  //
  hRadiusMap = new TH1D("RadiusMap", "event radius [cm] ", 100, 0., 10.);
  hRadiusMap->GetXaxis()->SetTitle("event radius R [cm]");
  //
  hRhoMap = new TH1D("RhoMap", "event cylindrical rho [cm] ", 100, 0., 4.);
  hRhoMap->GetXaxis()->SetTitle("cylindrical rho [cm]");
  //
  hZMap = new TH1D("ZMap", "event cylindrical Z [cm] ", 2000, -10., 10.);
  hZMap->GetXaxis()->SetTitle("event Z [cm]");
  //
  hPhiMap = new TH1D("PhiMap", "event phi", 100, -TMath::Pi(), TMath::Pi());
  hRhoZMap = new TH2D("RhoZMap", "cylindrical rho z  map ", 100, 0., 2., 100, 0., 4.);
  hRhoZMap->GetXaxis()->SetTitle("cylindrical rho [cm]");
  hRhoZMap->GetYaxis()->SetTitle("Z");
  //
  hRhoPhiZMap = new TH3D("RhoZPhiMap", "cylindrical rho phi z  map ", 100, 0., 2., 100, -TMath::Pi(), TMath::Pi(), 100, 0., 4.);
  hRhoPhiZMap->GetXaxis()->SetTitle("cylindrical rho [cm]");
  hRhoPhiZMap->GetYaxis()->SetTitle("phi");
  hRhoPhiZMap->GetZaxis()->SetTitle("Z");

  /** fid  */
  hXYMapFid = new TH2D("XYMapFid", "event Y versus X  [cm] ", 100, -4., 4., 100, -4., 4.);
  hXYMapFid->GetXaxis()->SetTitle("event X [cm]");
  hXYMapFid->GetYaxis()->SetTitle("event Y [cm]");

  hXYZMapFid = new TH3D("XYZMapFid", "event x y z  [cm] ", 100, -4., 4., 100, -4., 4., 50, 0, 4.);
  hXYZMapFid->GetXaxis()->SetTitle("event X [cm]");
  hXYZMapFid->GetYaxis()->SetTitle("event Y [cm]");
  hXYZMapFid->GetZaxis()->SetTitle("event Z [cm]");
  //
  hRadiusMapFid = new TH1D("RadiusMapFid", "event radius [cm] ", 100, 0., 10.);
  hRadiusMapFid->GetXaxis()->SetTitle("event radius R [cm]");
  //
  hRhoMapFid = new TH1D("RhoMapFid", "event cylindrical rho [cm] ", 100, 0., 4.);
  hRhoMapFid->GetXaxis()->SetTitle("cylindrical rho [cm]");
  //
  hZMapFid = new TH1D("ZMapFid", "event cylindrical Z [cm] ", 2000, -10., 10.);
  hZMapFid->GetXaxis()->SetTitle("event Z [cm]");
  //
  hPhiMapFid = new TH1D("PhiMapFid", "event phi", 100, -TMath::Pi(), TMath::Pi());
  hRhoZMapFid = new TH2D("RhoZMapFid", "cylindrical rho z  map ", 100, 0., 2., 100, 0., 4.);
  hRhoZMapFid->GetXaxis()->SetTitle("cylindrical rho [cm]");
  hRhoZMapFid->GetYaxis()->SetTitle("Z");
  //
  hRhoPhiZMapFid = new TH3D("RhoZPhiMapFid", "cylindrical rho phi z  map ", 100, 0., 2., 100, -TMath::Pi(), TMath::Pi(), 100, 0., 4.);
  hRhoPhiZMapFid->GetXaxis()->SetTitle("cylindrical rho [cm]");
  hRhoPhiZMapFid->GetYaxis()->SetTitle("phi");
  hRhoPhiZMapFid->GetZaxis()->SetTitle("Z");

  //
  hEventPass = new TH1D("hEventPass", "event pass", 4, 0, 4);
  hTriangle = new TH2D("Triangle", "ytern vs xtern", 100, 0., 1., 100, 0., 1.);
  hTriangle->GetXaxis()->SetTitle("xtern");
  hTriangle->GetYaxis()->SetTitle("ytern");
  hTriangleCut = new TH2D("TriangleCut", "ytern vs xtern", 100, 0., 1., 100, 0., 1.);
  hTriangleCut->GetXaxis()->SetTitle("xtern");
  hTriangleCut->GetYaxis()->SetTitle("ytern");
  hTriangleMean = new TH2D("TriangleMean", "ytern vs xtern", 100, 0., 1., 100, 0., 1.);
  hTriangleMean->GetXaxis()->SetTitle("xtern");
  hTriangleMean->GetYaxis()->SetTitle("ytern");

  hCount = new TH1D("Count", "hit count", 13, 0, 13);
  hCountSinglet = new TH1D("CountSinglet", "hit count singlet", 13, 0, 13);
  hTime = new TH1D("Time", "photon time ", 7500, 0, 2 * 7500);
  hTrigDiffTime = new TH1D("TrigDiffTime", " time difference ", MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
  hTrigDiffTime->GetXaxis()->SetTitle("max time diff [ns]");
  hTrigDiffTime30 = new TH1D("TrigDiffTime30", " time difference <30 ns", MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
  hTrigDiffTime30->GetXaxis()->SetTitle("max time diff [ns]");
  hTrigDiffTimeCut = new TH1D("TrigDiffTimeCut", Form(" time difference <%.0f ns", maxTriggerTimeDifference), MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
  hTrigDiffTimeCut->GetXaxis()->SetTitle("max time diff [ns]");
  hGammaPeak = new TH1D("GammaPeak", "sum of trigger sipms", 200, 0, 200);
  hGammaPeak->GetXaxis()->SetTitle("trigger sipm sum [SPE]");
  hGammaPeak->GetYaxis()->SetTitle("events per bin");

  // landau response function
  speLandau = new TF1("myLandau", myLandau, 0, MAXSAMPLE * binWidth, 3);
  // set SPE response parameters
  speLandau->SetParName(1, "MPV");
  speLandau->SetParameter(0, speMPV);
  speLandau->SetParName(1, "sigma");
  speLandau->SetParameter(1, speSigma);
  speLandau->SetParName(2, "norm");
  speLandau->SetParameter(2, 1); // single SPE
  // modelFit::modelFit(int theFit, int ichan, double ppm)
  hResponse = new TH1D("Response", "sipm response", MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
  hResponse->GetXaxis()->SetTitle("time [ns]");
  hResponse->GetYaxis()->SetTitle("photons/2ns");

  // fill response
  for (int ib = 1; ib < hResponse->GetNbinsX(); ++ib)
    hResponse->SetBinContent(ib, speLandau->Eval(hResponse->GetBinCenter(ib)) * double(binWidth));
  double landauMax = hResponse->GetBinContent(hResponse->GetMaximumBin());
  printf(" landau response integral %E  gain %f  landauMax %f SPE %E \n", hResponse->Integral(), gainFunc(9), landauMax, gainFunc(9) / landauMax);

  sigmaNoise = gainFunc(9) * noiseToSignal;
  TH1D *hPhotonSum9 = new TH1D("PhotonSum9", "photon sipm 9 sum/event", 50, 0., 50.);
  TH1D *hSingletSum9 = new TH1D("SingletSum9", "photon sipm 9 sum/event", 50, 0., 50.);
  TH1D *hPhotonAll = new TH1D("PhotonAll", "all photons ", 200, 0.5 * numPhotons, 1.5 * numPhotons);
  TH1D *hPhotonTrigSum = new TH1D("PhotonTrigSum", "photon trig sum ", 200, 0., 40.);
  TH1D *hPhotonSumCut = new TH1D("PhotonSumCut", "photon sum cut", 200, 0., 40.);

  // make individual light curves
  /*
  for (int ih = 0; ih < NCHAN; ++ih)
    models[ih] = new modelFit(MODELALL, ih, 0);
  */

  for (int ih = 0; ih < NCHAN; ++ih)
  {
    int ilevel = getLevel(ih);
    hGeoEff[ih] = new TH1D(Form("GeoEff%i", ih), Form("Photon%i-level%i", ih, ilevel), 5000, 0, 100);
  }

  eventDir = fout->mkdir("eventDir");
  histDir = fout->mkdir("histDir");
  histDir->cd();

  TH1D *hNoise = new TH1D("Noise", "Noise", 200, -10 * sigmaNoise, 10 * sigmaNoise);
  histDir->cd();

  for (int ih = 0; ih < NCHAN; ++ih)
  {
    int ilevel = getLevel(ih);
    hPhotonTime[ih] = new TH1D(Form("PhotonTime%i", ih), Form("PhotonTime%i-level%i", ih, ilevel), MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
    hPhotonTime[ih]->GetXaxis()->SetTitle("time [ns]");
    hPhotonTime[ih]->GetYaxis()->SetTitle("photons/2ns");
    // hPhotonTime[ih]->SetDirectory(nullptr);
  }

  for (int ih = 0; ih < NCHAN; ++ih)
  {
    // modelFit::modelFit(int theFit, int ichan, double ppm)
    int ilevel = getLevel(ih);
    hPhotonSum[ih] = new TH1D(Form("PhotonSum%i", ih), Form("PhotonSum%i-level%i", ih, ilevel), MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
    hPhotonSum[ih]->GetXaxis()->SetTitle("time [ns]");
    hPhotonSum[ih]->GetYaxis()->SetTitle("photons/2ns");
    //
  }
  for (int ih = 0; ih < NCHAN; ++ih)
  {
    int ilevel = getLevel(ih);
    // modelFit::modelFit(int theFit, int ichan, double ppm)
    hPhoton[ih] = new TH1D(Form("Photon%i", ih), Form("Photon%i-level%i", ih, ilevel), MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
    hPhoton[ih]->GetXaxis()->SetTitle("time [ns]");
    hPhoton[ih]->GetYaxis()->SetTitle("photons/2ns");
    hPhoton[ih]->SetDirectory(nullptr);

    // modelFit::modelFit(int theFit, int ichan, double ppm)
    hSinglet[ih] = new TH1D(Form("Singlet%i", ih), Form("Singlet%i-level%i", ih, ilevel), MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
    hSinglet[ih]->GetXaxis()->SetTitle("time [ns]");
    hSinglet[ih]->GetYaxis()->SetTitle("photons/2ns");
    hSinglet[ih]->SetDirectory(nullptr);
    //
    hConvolve[ih] = new TH1D(Form("Convolve%i", ih), Form("Convolve%i-level%i", ih, ilevel), MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
    hConvolve[ih]->GetXaxis()->SetTitle("time [ns]");
    hConvolve[ih]->GetYaxis()->SetTitle("photons/2ns");
    hConvolve[ih]->SetDirectory(nullptr);

    //
    hSignalNb[ih] = new TH1D(Form("SignalNb%i", ih), Form("SignalNb%i-level%i", ih, ilevel), MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
    hSignalNb[ih]->GetXaxis()->SetTitle("time [ns]");
    hSignalNb[ih]->GetYaxis()->SetTitle("photons/2ns");
    hSignalNb[ih]->SetDirectory(nullptr);
    //
    hSignal[ih] = new TH1D(Form("Signal%i", ih), Form("Signal%i-level%i", ih, ilevel), MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
    hSignal[ih]->GetXaxis()->SetTitle("time [ns]");
    hSignal[ih]->GetYaxis()->SetTitle("photons/2ns");
    hSignal[ih]->SetDirectory(nullptr);

  } //
  for (int ih = 0; ih < NCHAN; ++ih)
  {
    int ilevel = getLevel(ih);

    hSignalSumNoBaseline[ih] = new TH1D(Form("SignalSumNoBaseline%i", ih), Form("SignalSumNoBaseline%i-level%i", ih, ilevel), MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
    hSignalSumNoBaseline[ih]->GetXaxis()->SetTitle("time [ns]");
    hSignalSumNoBaseline[ih]->GetYaxis()->SetTitle("photons/2ns");

    hSignalSum[ih] = new TH1D(Form("SignalSum%i", ih), Form("SignalSum%i-level%i", ih, ilevel), MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
    hSignalSum[ih]->GetXaxis()->SetTitle("time [ns]");
    hSignalSum[ih]->GetYaxis()->SetTitle("photons/2ns");
    hSignalNorm[ih] = new TH1D(Form("SignalNorm%i", ih), Form("SignalNorm%i-level%i", ih, ilevel), MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
    hSignalNorm[ih]->GetXaxis()->SetTitle("time [ns]");
    hSignalNorm[ih]->GetYaxis()->SetTitle("photons/2ns");

    hSignalEff[ih] = new TH1D(Form("SignalEff%i", ih), Form("SignalEff%i-level%i", ih, ilevel), MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
    hSignalEff[ih]->GetXaxis()->SetTitle("time [ns]");
    hSignalEff[ih]->GetYaxis()->SetTitle("photons/2ns");
  }

  // test histograms
  testDir = fout->mkdir("testDir");
  testDir->cd();
  ntConvolve = new TNtuple("ntConvolve", "ntConvolve", "time:nph:qconv:qsum");

  hPhotonTest9 = new TH1D("PhotonTest9", "ConvolveTest9", MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
  hPhotonTest9->GetXaxis()->SetTitle("time [ns]");
  hPhotonTest9->GetYaxis()->SetTitle("photons/2ns");

  hPhotonSumTest9 = new TH1D("PhotonSumTest9", "ConvolveTest9", MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
  hPhotonSumTest9->GetXaxis()->SetTitle("time [ns]");
  hPhotonSumTest9->GetYaxis()->SetTitle("photons/2ns");

  hConvolveTest9 = new TH1D("ConvolveTest9", "ConvolveTest9", MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
  hConvolveTest9->GetXaxis()->SetTitle("time [ns]");
  hConvolveTest9->GetYaxis()->SetTitle("photons/2ns");

  hSignalSumTest9 = new TH1D("SignalSumTest9", "SignalSumTest9", MAXSAMPLE, 0, MAXSAMPLE * (binWidth));
  hSignalSumTest9->GetXaxis()->SetTitle("time [ns]");
  hSignalSumTest9->GetYaxis()->SetTitle("photons/2ns");

  /* end of test dir*/
  /* end of define ntuples amd histograms here */
  // print info for channel
  printf("\n******* efficienies ****\n ");
  for (int ich = 8; ich >= 0; --ich)
  {
    double eff = effGeoFunc(ich);
    int ilevel = getLevel(ich);
    printf(" chan %i level %i distance %f eff %E total eff %E\n", ich, ilevel, distanceLevel[ilevel], eff, eff * SiPMQE128Ham * fillFactor);
  }

  // trigger sipms
  double effTrigger = effGeoSim(9);
  for (int ich = 9; ich < NCHANPMT; ++ich)
  {
    int ilevel = distanceLevel[ich];
    printf(" chan %i level %i origin(%f,%f,%f) distance %f eff %E total eff %E\n", ich, ilevel,
           eventOrigin.X(), eventOrigin.Y(), eventOrigin.Z(), distanceLevel[ilevel], effTrigger, effTrigger * SiPMQE128Ham * fillFactor);
  }

  printf("\t\t nominal yield nphotons %.0f  3 SIPM sum %f\n,", numPhotons, 3. * numPhotons * effTrigger * SiPMQE128Ham * fillFactor);
  // print info for pmt
  double effPmt = effGeoFunc(NCHANPMT);
  int ilevel = distanceLevel[4];
  printf(" chan %i level %i distance %f eff %E \n", 12, ilevel, distanceLevel[ilevel], effPmt);
  printf("***********\n\n\n ");

  /* static double tTriplet0 = 1600.0; ns
    static double tSinglet0 = 7.0; ns
  */

  /*******
   *            loop over events
   *******/
  totalPhotons = 0;
  int nsinglet;
  int ntriplet;
  int nTrigger = 0;
  for (int iev = 0; iev < ngen; ++iev) // start of event loop
  {
    if (iev / 1000 * 1000 == iev)
    {
      printf("...btbSim event %i passed %lld \n", iev, ntFit->GetEntries());
      fflush(stdout);
    }
    hSignalPhotonsEvent[0]->Reset("ICESM");
    hSignalPhotonsEvent[1]->Reset("ICESM");
    hSignalPhotonsEvent[2]->Reset("ICESM");
    for (int ich = 0; ich < NCHAN; ++ich)
    {
      // histogram reset
      hPhoton[ich]->Reset("ICESM");
      hSinglet[ich]->Reset("ICESM");
      hConvolve[ich]->Reset("ICESM");
      hSignalNb[ich]->Reset("ICESM");
      hSignal[ich]->Reset("ICESM");
      // hSignalSumNoBaseline[ich]->Reset("ICESM");
      // hSignalSumNoBaseline[ich]->Reset("ICESM");
    }

    hEventPass->SetBinContent(1, hEventPass->GetBinContent(1) + 1); // generated over 4 PI

    // zero trigger times array
    for (int i = 0; i < 3; ++i)
      triggerTimes[i] = 0.;

    int nPhotonsEvent = (int)ran->Gaus(numPhotons, sqrt(numPhotons));
    totalPhotons += nPhotonsEvent;
    hPhotonAll->Fill(nPhotonsEvent);

    /********  generate gamma position  *********/
    // double gammaCosTheta = 2. * ran->Rndm() - 1.; // cos flat from 1 to -1
    /* cut out the events blocked by the source holder
    if (gammaCosTheta < cosMin)
      continue;
    */

    // generate r,cosTheta from geant4 map;

    double gammaR = 0;
    double gammaCosTheta = -1;
    double gammaPhi = 0;
    /********  generate gamma position  from geant4 map  *********/
    if (useMap)
    {
      originPDF->GetRandom3(gammaR, gammaCosTheta, gammaPhi, ran);
    }
    /* generate from Rndm */
    else
    {
      gammaR = abs(ran->Exp(meanFreePath));             // See TMath Double_t TRandom::Exp	(	Double_t	tau	)	returns exp( -t/tau )
      gammaCosTheta = 2.0 * ran->Rndm() - 1.;           // cos range ig -1 to +1
      gammaPhi = (2. * ran->Rndm() - 1.) * TMath::Pi(); // -pi to pi
    }
    eventOrigin = getXYZVector(gammaR, acos(gammaCosTheta), gammaPhi);
    /*  add z offset zZero if we are using new geometry */
    // if (!geoVersionOld)
    //   eventOrigin = eventOrigin;
    //   +eventOriginOffset;
    double localPhi = eventOrigin.Phi() * 360. / TMath::TwoPi();
    if (localPhi < 0)
      localPhi += 360.;
    hZMap->Fill(eventOrigin.Z());
    hRadiusMap->Fill(eventOrigin.R());
    hRhoMap->Fill(eventOrigin.Rho());
    hRhoZMap->Fill(eventOrigin.Rho(), eventOrigin.Z());
    hRhoPhiZMap->Fill(eventOrigin.Rho(), eventOrigin.Phi(), eventOrigin.Z());
    hPhiMap->Fill(eventOrigin.Phi());
    hXYMap->Fill(eventOrigin.X(), eventOrigin.Y());
    hXYZMap->Fill(eventOrigin.X(), eventOrigin.Y(), eventOrigin.Z());

    // cut -Z (up going in btb) events
    /* this is how events were generated */
    /*
    if (gammaCosTheta < 0)
    {
      if (show)
        printf(" skip event %i %f \n", iev, gammaCosTheta);
      continue;
    }
    */

    /** convolve test  ***/
    if (iev == -1)
    {
      printf("******************* run convolution test  ********************\n");
      convolveTest(1000);
    }

    // fiducial cut
    double effGeoSim9 = effGeoSim(9);
    double effGeoSim10 = effGeoSim(10);
    double effGeoSim11 = effGeoSim(11);
    // Z > 0 means the interaction point is on the detector side of the source (array at +Z)
    isFid = false;
    if (eventOrigin.Z() > 0.0 && effGeoSim9 > 0 && effGeoSim10 > 0 && effGeoSim11 > 0)
      isFid = true;

    if (!isFid)
    {
      if (show)
        printf(" skip event %i XYZ %.3f %.3f %.3f \n", iev, eventOrigin.X(), eventOrigin.Y(), eventOrigin.Z());
      continue;
    }
    hZMapFid->Fill(eventOrigin.Z());
    hRadiusMapFid->Fill(eventOrigin.R());
    hRhoMapFid->Fill(eventOrigin.Rho());
    hRhoZMapFid->Fill(eventOrigin.Rho(), eventOrigin.Z());
    hRhoPhiZMapFid->Fill(eventOrigin.Rho(), eventOrigin.Phi(), eventOrigin.Z());
    hPhiMapFid->Fill(eventOrigin.Phi());
    hXYMapFid->Fill(eventOrigin.X(), eventOrigin.Y());
    hXYZMapFid->Fill(eventOrigin.X(), eventOrigin.Y(), eventOrigin.Z());

    hEventPass->SetBinContent(2, hEventPass->GetBinContent(2) + 1); // for fiducial events
    // eventOrigin.SetZ(abs(eventOrigin.Z()));

    ntOrigin->Fill(iev, eventOrigin.R(), cos(eventOrigin.Theta()), eventOrigin.Theta() * 360. / TMath::TwoPi(), localPhi, eventOrigin.X(), eventOrigin.Y(), eventOrigin.Z());
    // printf(" event origin x %f y %f z %f \n", eventOrigin.X(), eventOrigin.Y(), eventOrigin.Z());

    if (iev / reportInterval * reportInterval == iev)
      printf("... event %i total photon %0.f (rho,z,phi) = (%f, %f, %f) (r,theta,Phi) = (%f , %f ,%f ) \n", iev, double(totalPhotons), eventOrigin.Rho(), eventOrigin.Z(), eventOrigin.Phi() * 360. / TMath::TwoPi(), eventOrigin.R(), eventOrigin.Theta() * 360. / TMath::TwoPi(), eventOrigin.Phi() * 360. / TMath::TwoPi());

    // loop over channels
    for (int ich = 0; ich < NCHAN; ++ich)
    {
      // set the nominal gain from file modelFitGamma.hh
      double timeShift = 0;
      if (ich > 8 && ich < 12)
      {
        timeShift = timeOffset + trigTimeShift[ich - 9]; // trig amp delay and gain shift effect on trigger threshold
      }
      //
      sigmaNoise = gainFunc(ich) * noiseToSignal;
      bool invert = ich > 8; // trigger SiPMs (9-11) have inverted ADC polarity: stored as 2^14 - signal

      if (rawRun)
      {
        rawEvent = rawRun->getDet(ich); // If channel branch doesn't exist getDet calls addDet
        rawEvent->clear();
        rawEvent->channel = unsigned(ich);
        rawEvent->trigger = unsigned(iev);
      }
      TDet *det = simRun->getDet(ich); // If channel branch doesn't exist getDet calls addDet
      det->clear();
      det->event = iev;
      det->trigger = iev;
      det->nspe = nsinglet + ntriplet;

      // simRun->btree->GetListOfBranches()->ls();
      // printf(" simRun ev %i  channel %i \n", iev, ich);

      // rawEvent->time = EventInfo->TriggerTimeTag;
      bool isTrig = false;
      isTrig = ich == 9 || ich == 10 || ich == 11;
      double effGeoSimi = effGeoSim(ich);
      hGeoEff[ich]->Fill(nominalGeo9 / effGeoSimi);
      // printf("xxxxv event %i chan %i (r,cosTheta,phi) (%f,%f,%f) effGeo %E  \n", iev, ich, eventOrigin.R(), eventOrigin.Theta() * 360. / TMath::TwoPi(), localPhi, effGeoSimi);
      if (isTrig)
        hEffGeo->Fill(effGeoSimi / nominalGeo);
      if (ich == 9)
        hEffGeo9->Fill(effGeoSimi / nominalGeo);
      if (ich == 10)
        hEffGeo10->Fill(effGeoSimi / nominalGeo);
      if (ich == 11)
        hEffGeo11->Fill(effGeoSimi / nominalGeo);

      // compInt includes QE
      double eff = effGeoSimi * fillFactor * reflection;
      if (show)
      {
        printf("MESSAGE ******** compIntegrals ******* \n");
        for (int icomp = 0; icomp < NUMCOMP; ++icomp)
        {
          printf("chan %i comp %i compIntegral %.3E\n", ich, icomp, compIntegral[ich][icomp]);
        }
        printf("chan %i compIntegralSum %.3E\n\n", ich, compIntegralSum[ich]);
      }

      // get commponent photons
      numCompPhotons.clear();
      numCompPhotons.resize(NUMCOMP);
      for (int icomp = 0; icomp < NUMCOMP; ++icomp)
      {
        numCompPhotons[icomp] = int(nPhotonsEvent * eff * compIntegral[ich][icomp] / compIntegralSum[ich]);
        // printf("line 1457 iev %i ich %i icomp %i nphotons %i comp fraction %E nphotons %i\n", iev, ich, icomp, nPhotonsEvent, eff * compIntegral[ich][icomp] / compIntegralSum[ich], numCompPhotons[icomp]);
      }

      double nsmean = double(nPhotonsEvent) * eff * singletFrac;
      double ntmean = double(nPhotonsEvent) * eff - nsmean;
      nsinglet = ran->Poisson(nsmean);
      ntriplet = ran->Poisson(ntmean);
      hPoisson->Fill(nsinglet + ntriplet);

      ncount[ich] += nsinglet + ntriplet;
      ncountSinglet[ich] += nsinglet;
      totalPhotons += nsinglet + ntriplet;

      /*
          big loop over time components
      */
      for (int icomp = 0; icomp < NUMCOMP; ++icomp)
      {
        // printf("call getTime ev %d chan %i \n", iev, ich);
        getTime(ich, icomp, numCompPhotons[icomp]);

        // printf("line 1492 iev %i ich %i icomp %i nphotons %i size %lu \n", iev, ich, icomp, numCompPhotons[icomp], timeComp.size());

        // fill timeComp photons
        for (unsigned iphoton = 0; iphoton < timeComp.size(); ++iphoton)
        {
          hPhotonTime[ich]->Fill(timeComp[iphoton]);

          double time = timeShift + triggerStart + timeComp[iphoton];
          double gain = gainFunc(ich);
          hPhoton[ich]->Fill(time, gain);
          hSinglet[ich]->Fill(time, gain);
          // hPhotonSum[ich]->Fill(time, gain);
          hPhotonSum[ich]->Fill(time);
          nChannel[ich] = nChannel[ich] + 1;
          convolve(hConvolve[ich], time, gainFunc(ich));
          TH1D *hist = hConvolve[ich];

          if (ich == 9)
            hPhotonTrig[0]->Fill(time);
          if (ich == 10)
            hPhotonTrig[1]->Fill(time);
          if (ich == 11)
            hPhotonTrig[2]->Fill(time);

          // printf("event %i chan %i  max value %E\n", iev, ich, hist->GetBinContent(hist->GetMaximumBin()));
          hTime->Fill(time);
          // make a TDetHit for photon
          TDetHit hit;
          hit.startTime = double(hTime->FindBin(time)); // convert to samples
          // printf("line579 time %f %f bin %i  \n", time, hit.startTime, hPhoton[ich]->FindBin(time));
          hit.qpeak = gainFunc(ich);
          det->hits.push_back(hit);
        }
      }

      // if (ich == 9)
      //   printf("line960 event %i photons 9 entries %f photon trig0 %f \n", iev, hPhoton[9]->GetEntries(), hPhotonTrig[0]->GetEntries());

      // add baseline and noise
      for (int ibin = 1; ibin <= hSignal[ich]->GetNbinsX(); ++ibin)
      {
        double binNoise = ran->Gaus(0.0, sigmaNoise);
        hNoise->Fill(binNoise);
        hSignalNb[ich]->SetBinContent(ibin, binNoise + hConvolve[ich]->GetBinContent(ibin));
        hSignal[ich]->SetBinContent(ibin, baseline + binNoise + hConvolve[ich]->GetBinContent(ibin));
        hSignalSumNoBaseline[ich]->SetBinContent(ibin, binNoise +
                                                           hConvolve[ich]->GetBinContent(ibin) + hSignalSumNoBaseline[ich]->GetBinContent(ibin));
        hSignalSum[ich]->SetBinContent(ibin, baseline + binNoise +
                                                 hConvolve[ich]->GetBinContent(ibin) + hSignalSum[ich]->GetBinContent(ibin));
      }
      /*
      if (ich == 9)
        printf(">>>> chan %i nsinglet %i ntriplet %i nph %.0f convolve %0.3f qsum %.0f\n", ich, nsinglet, ntriplet,
               hPhoton[ich]->Integral(), hConvolve[ich]->Integral(), hSignalSumNoBaseline[ich]->Integral());
               */
      ntConvolveCheck->Fill(float(ich), float(nsinglet), float(ntriplet), hPhoton[ich]->Integral(), hConvolve[ich]->Integral(), hSignalNb[ich]->Integral());

      /* event histograms */
      TString histName;
      if (eventDir->GetList()->GetEntries() < 100 && hPhoton[ich]->GetEntries() > 0)
      {
        eventDir->cd();
        histName.Form("hPhotonCh%iEv%i", ich, iev);
        TH1D *hPhotonEvent = (TH1D *)hPhoton[ich]->Clone(histName);
        hPhotonEvent->SetTitle(histName);
        //
        histName.Form("hConvolCh%iEv%i", ich, iev);
        TH1D *hConvolveEvent = (TH1D *)hConvolve[ich]->Clone(histName);
        hConvolveEvent->SetTitle(histName);

        histName.Form("hSignalNoBaselineCh%iEv%i", ich, iev);
        TH1D *hNoBaseline = (TH1D *)hSignalNb[ich]->Clone(histName);
        hNoBaseline->SetTitle(histName);
        //
        histName.Form("hSignalCh%iEv%i", ich, iev);
        TH1D *hSignalEvent = (TH1D *)hSignal[ich]->Clone(histName);
        hSignalEvent->SetTitle(histName);
      }

      // file wave for this channel
      if (rawRun)
      {
        wave.clear();
        if (invert) // trigger sipm and ADC
        {
          for (int ibin = 1; ibin <= hSignal[ich]->GetNbinsX(); ++ibin)
          {
            uint16_t adc = pow(2, 14) - hSignal[ich]->GetBinContent(ibin);
            /*
            uint16_t unAdc = -1. * (adc - pow(2, 14));
            if (ibin == 1)
              printf("line831 iev %i channnel %i bin %i adc %u %u \n", iev, ich, ibin, adc, unAdc);
              */
            wave.push_back(adc);
          }
        }
        else
        {
          for (int ibin = 1; ibin <= hSignal[ich]->GetNbinsX(); ++ibin)
            wave.push_back((uint16_t)hSignal[ich]->GetBinContent(ibin));
          // for (int ibin = 1; ibin <= hSignal[ich]->GetNbinsX(); ++ibin)
          //  printf("iev %i ich %i ibin %i short %u float %f \n", iev, ich, ibin, (uint16_t)hSignal[ich]->GetBinContent(ibin), hSignal[ich]->GetBinContent(ibin));
        }
        rawEvent->rdigi = wave;
      }
      double qsum = 0;
      for (int ibin = 1; ibin <= hSignal[ich]->GetNbinsX(); ++ibin)
        qsum += (hSignal[ich]->GetBinContent(ibin) - baseline) / gainFunc(ich) * landauMax;

      double psum = 0;
      for (int ibin = 1; ibin <= hPhoton[ich]->GetNbinsX(); ++ibin)
        psum += hPhoton[ich]->GetBinContent(ibin) / gainFunc(ich);

      ntTrigCh->Fill(iev, ich, qsum, psum, hPhoton[ich]->GetEntries(), eventOrigin.Rho(), localPhi, eventOrigin.X(), eventOrigin.Y(), eventOrigin.Z(), effGeoSimi / nominalGeo);
      // printf("line508 ch %i nPhotonsEvent %i eff %E nPhotonsEvent*eff %.0f nhotons %i %i \n", ich, nPhotonsEvent, eff, nPhotonsEvent * eff, nsinglet + ntriplet, int(hPhoton[ich]->GetEntries()));
      // if (iev / 1 * 1 == iev && ich < 12 && ich > 8)

      // printf("line993 ich %i %f \n", ich, hPhoton[ich]->GetEntries());
    } // end channel loop

    /* summed convolved waveforms  */
    for (int ibin = 0; ibin < hSignalPhotons[0]->GetNbinsX(); ++ibin)
    {
      hSignalPhotons[0]->SetBinContent(ibin, hSignalPhotons[0]->GetBinContent(ibin) + hConvolve[9]->GetBinContent(ibin) / readGains->sipmSumGain[9]);
      hSignalPhotons[1]->SetBinContent(ibin, hSignalPhotons[1]->GetBinContent(ibin) + hConvolve[10]->GetBinContent(ibin) / readGains->sipmSumGain[10]);
      hSignalPhotons[2]->SetBinContent(ibin, hSignalPhotons[2]->GetBinContent(ibin) + hConvolve[11]->GetBinContent(ibin) / readGains->sipmSumGain[11]);
      hSignalPhotonsEvent[0]->SetBinContent(ibin, hSignalPhotonsEvent[0]->GetBinContent(ibin) + hConvolve[9]->GetBinContent(ibin) / readGains->sipmSumGain[9]);
      hSignalPhotonsEvent[1]->SetBinContent(ibin, hSignalPhotonsEvent[1]->GetBinContent(ibin) + hConvolve[10]->GetBinContent(ibin) / readGains->sipmSumGain[10]);
      hSignalPhotonsEvent[2]->SetBinContent(ibin, hSignalPhotonsEvent[2]->GetBinContent(ibin) + hConvolve[11]->GetBinContent(ibin) / readGains->sipmSumGain[11]);
    }
    // printf("line991 event %i  photon integrals 9 %.3f %.3f 10 %.3f %.3f 11 %.3f  %.3f\n", iev, hPhoton[9]->Integral(), hSignalPhotonsEvent[9]->Integral(),
    //        hPhoton[10]->Integral(), hSignalPhotonsEvent[10]->Integral(), hPhoton[11]->Integral(), hSignalPhotonsEvent[11]->Integral());
    // printf("line993 event %i  photon integrals 9) %.3f %.3f 10) %.3f %.3f 11) %.3f  %.3f\n", iev, hPhoton[9]->GetEntries(), hSignalPhotons[0]->Integral(),
    //       hPhoton[10]->GetEntries(), hSignalPhotons[1]->Integral(), hPhoton[11]->GetEntries(), hSignalPhotons[2]->Integral());

    // ensure the event triggers
    double maxTriggerDiff = eventTrigger(iev);
    hTrigDiffTime->Fill(maxTriggerDiff);
    // printf("line1276....%i %i %f \n", iev, isFid, maxTriggerDiff);
    //  if (maxTriggerDiff < maxTriggerTimeDifference)
    //  hTrigDiffTime10->Fill(maxTriggerDiff);
    if (maxTriggerDiff < 30)
      hTrigDiffTime30->Fill(maxTriggerDiff);

    // event passes trigger
    bool trigPass = false;
    if (maxTriggerDiff < maxTriggerTimeDifference)
    {
      hEventPass->SetBinContent(4, hEventPass->GetBinContent(4) + 1);
      trigPass = true;
    }
    else if (maxTriggerDiff < double(2 * MAXSAMPLE))
      printf("fails trigger time cut event %i tdiff %f \n", iev, maxTriggerDiff);

    // as in real data qsum
    peakQsum[0] = hSignalNb[9]->Integral() / qsumNominalFromBtb;
    peakQsum[1] = hSignalNb[10]->Integral() / qsumNominalFromBtb;
    peakQsum[2] = hSignalNb[11]->Integral() / qsumNominalFromBtb;

    // printf("line990 photons %f signal %f normed %f \n", hPhoton[9]->GetEntries(), hSignalNb[9]->Integral(), peakQsum[0]);

    double gammaPeakSum = peakQsum[0] + peakQsum[1] + peakQsum[2];

    /* trigger pass cut*/
    /*
    double prompt9 = hPhoton[9]->Integral(0, int(maxTriggerDiff / 2.)) / gain;
    double prompt10 = hPhoton[10]->Integral(0, int(maxTriggerDiff / 2.)) / gain;
    double prompt11 = hPhoton[11]->Integral(0, int(maxTriggerDiff / 2.)) / gain;
    double prompt9 = hPhoton[9]->Integral(730, 745) / gain;
    double prompt10 = hPhoton[10]->Integral(730, 745) / gain;
    double prompt11 = hPhoton[11]->Integral(730, 745) / gain;
    // printf("event %i %.0f, %.0f, %.0f \n", iev, prompt9, prompt10, prompt11);
    if (prompt9 > 0 && prompt10 > 0 && prompt11 > 0)
      hEventPass->SetBinContent(3, hEventPass->GetBinContent(3) + 1);
    */

    if (!trigPass)
    {
      if (show)
        printf(" event %i does not trigger %f \n", iev, maxTriggerDiff);
      continue;
    }
    ++nTrigger;
    // printf("triggered event % i prompt sums %.0f, %.0f, %.0f \n", iev, prompt9, prompt10, prompt11);
    hTrigDiffTimeCut->Fill(maxTriggerDiff);
    double photonSum = hPhoton[9]->GetEntries() + hPhoton[10]->GetEntries() + hPhoton[11]->GetEntries();
    hPhotonTrigSum->Fill(photonSum);

    ntTrig->Fill(iev, photonSum, hPhoton[9]->GetEntries(), hPhoton[10]->GetEntries(), hPhoton[11]->GetEntries(), eventOrigin.R(), eventOrigin.Rho(), localPhi, eventOrigin.X(), eventOrigin.Y(), eventOrigin.Z(), maxTriggerDiff, trigPass);

    if (show)
      printf("xxx event %i nph %.0f %.0f %.0f\n", iev, hPhoton[9]->GetEntries(), hPhoton[10]->GetEntries(), hPhoton[11]->GetEntries());

    // cut on fitted radius
    // printf("line577 fitted radius %f\n", fitVal[1]);
    // tell triggerPeakFit the qsum normalize to the total photons
    peakFitQsum[0] = hPhoton[9]->GetEntries();
    peakFitQsum[1] = hPhoton[10]->GetEntries();
    peakFitQsum[2] = hPhoton[11]->GetEntries();

    // make fraction cut
    /* try a cut like TUM */
    double qFraction[3];
    double trigQSum = peakQsum[0] + peakQsum[1] + peakQsum[2];
    qFraction[0] = peakQsum[0] / trigQSum;
    qFraction[1] = peakQsum[1] / trigQSum;
    qFraction[2] = peakQsum[2] / trigQSum;

    // Now ready for minimization step with MIGRAD
    // set starting param values
    double fitVal[NPARS];
    double fitErr[NPARS];

    for (unsigned i = 0; i < 3; ++i)
    {
      fitVal[i] = 0;
      fitErr[i] = 0;
    }
    double amin = 0;

    double step = 0.0001;

    // printf("event %i fill fit ntuple\n", iev);
    //  fill fit ntuple
    if (nTrigger / reportInterval * reportInterval == nTrigger)
    {
      printf(".x.x.x report event %i nLL %f nphotons %i  singlet %i triplet %i tot  %i  photons (%.0f, %.0f, %.0f sum %.0f )  qsum(%f,%f,%f) mean(%f,%f,%f) \n", iev, amin, nPhotonsEvent, nsinglet, ntriplet, nsinglet + ntriplet, hPhoton[9]->GetEntries(), hPhoton[10]->GetEntries(), hPhoton[11]->GetEntries(), photonSum, peakFitQsum[0], peakFitQsum[1], peakFitQsum[2], peakMeanQsum[0], peakMeanQsum[1], peakMeanQsum[2]);
    }

    double xternMean, yternMean;
    makeTernary(peakMeanQsum[0], peakMeanQsum[1], peakMeanQsum[2], xternMean, yternMean);
    hTriangleMean->Fill(xternMean, yternMean);

    // use photon  number
    double xternPh, yternPh;
    makeTernary(peakFitQsum[0], peakFitQsum[1], peakFitQsum[2], xternPh, yternPh);

    // use hSignal integral
    double xternQ, yternQ;
    makeTernary(peakQsum[0], peakQsum[1], peakQsum[2], xternQ, yternQ);
    hTriangle->Fill(xternQ, yternQ);

    bool passTriangle = false;
    if ((xternQ > 0.3 && xternQ < 0.65) && (xternQ > 0.1 && yternQ < 0.5))
      passTriangle = true;

    if (passTriangle)
    {
      hTriangleCut->Fill(xternQ, yternQ);
      hPhotonSumCut->Fill(photonSum);
    }
    // printf("line1197 %f\n", gammaPeakSum);
    hGammaPeak->Fill(gammaPeakSum);
    ntGammaPeak->Fill(iev, nPhotonsEvent, hPhoton[9]->GetEntries(), hPhoton[10]->GetEntries(), hPhoton[11]->GetEntries(), peakQsum[0], peakQsum[1], peakQsum[2], gammaPeakSum);

    // if (!passTriangle)
    //   printf("line959 event %i !passTriangle %f %f \n", iev, xternQ, yternQ);

    // printf("line892 nph  (%.0f  %.0f  %.0f)  qsum (%.3f   %.3f  %.3f) xtern (%.3f %.3f)  ytern (%.3f %.3f)  \n", peakFitQsum[0], peakFitQsum[1], peakFitQsum[2], peakQsum[0], peakQsum[1], peakQsum[2], xternPh, xternQ, yternPh, yternQ);

    ntTern->Fill(eventOrigin.Rho(), eventOrigin.Phi(), eventOrigin.Z(), peakFitQsum[0], peakFitQsum[1], peakFitQsum[2], qFraction[0], qFraction[1], qFraction[2], peakMeanQsum[0], peakMeanQsum[1], peakMeanQsum[2], xternQ, yternQ);

    // printf(" nTrigger %i %.0f %.0f %.0f (%f %f %f)  xtern %f ytern %f \n", nTrigger, hPhoton[9]->GetEntries(), hPhoton[10]->GetEntries(), hPhoton[11]->GetEntries(), peakFitQsum[0], peakFitQsum[1], peakFitQsum[2], xternQ, yternQ);

    ntMean->Fill(iev, photonSum, eventOrigin.Rho(), eventOrigin.Phi(), eventOrigin.Z(),
                 peakFitQsum[0], peakFitQsum[1], peakFitQsum[2], peakMeanQsum[0], peakMeanQsum[1], peakMeanQsum[2], xternQ, yternQ);
    ntFit->Fill(iev, fitVal[0], eventOrigin.Rho(), eventOrigin.Phi(), eventOrigin.Z(),
                peakFitQsum[0], peakFitQsum[1], peakFitQsum[2], fitVal[1], cos(fitVal[2]), fitVal[3], fitErr[1], fitErr[2], ierflg);

    if (rawRun)
    {
      printf("fill rawRun %i\n", iev);
      rawRun->fill();
    }

    printf("fill simRun %i\n", iev);
    simRun->fill();
    hPhotonSum9->Fill(hPhoton[9]->GetEntries());
    hSingletSum9->Fill(hSinglet[9]->GetEntries());
    totalEventEffiency = hPhotonSumCut->Integral() / hPhotonAll->Integral();

  } // end of event loop

  // normalize
  double totalPass = hEventPass->GetBinContent(1);

  for (int ich = 0; ich < NCHAN; ++ich)
  {
    double eff = effGeoFunc(ich) * SiPMQE128Ham;
    printf("*** normalize hSignal norm to total pass %.0f geo %.E gain %.E \n", totalPass, eff, gainFunc(ich));
    for (int ibin = 1; ibin <= hSignal[ich]->GetNbinsX(); ++ibin)
    {
      // calculate channel efficiency
      /* have to normalize by gain */
      chanEff[ich] = (nChannel[ich]) / double(totalPhotons);
      hSignalNorm[ich]->SetBinContent(ibin, hSignalSumNoBaseline[ich]->GetBinContent(ibin) / double(totalPass) / gainFunc(ich));
      hSignalEff[ich]->SetBinContent(ibin, hSignalSumNoBaseline[ich]->GetBinContent(ibin) / double(totalPass) / chanEff[ich] / gainFunc(ich));
    }
  }

  for (int ich = 0; ich < NCHAN; ++ich)
  {
    hCount->SetBinContent(ich + 1, double(ncount[ich]) / double(ngen));
    hCountSinglet->SetBinContent(ich + 1, double(ncountSinglet[ich]) / double(ngen));
  }

  // summary
  printf("****** generated %i events.\nphoton count: \n", ngen);
  for (int ih = 0; ih < NCHAN; ++ih)
  {
    if (nTrigger > 1)
      printf(" chan %i (singlet,total) photons (%i, %i)  (singlet,total) photons/triggered (%.3f, %.3f)  \n",
             ih, (int)ncountSinglet[ih], (int)ncount[ih], double(ncountSinglet[ih]) / double(nTrigger), double(ncount[ih]) / double(nTrigger));
    else
      printf("nTrigger = 0\n");
  }
  double sipmPhotonSum = 0;
  printf("******** Photon count by level: *******\n");
  for (int ich = 0; ich < NCHAN; ++ich)
  {
    double eff = effGeoFunc(ich) * SiPMQE128Ham;
    printf("channel %i nphotons %i nphotons/total %.2E eff %.2E photon sum %.2E no baseline Sum %.2E  norm %.2f signalEff %.2F \n", ich, nChannel[ich], chanEff[ich], eff, hPhotonSum[ich]->Integral(), hSignalSumNoBaseline[ich]->Integral(), hSignalNorm[ich]->Integral(), hSignalEff[ich]->Integral());

    ntNorm->Fill(ich, eff, hPhotonSum[ich]->Integral(), hSignalSumNoBaseline[ich]->Integral(), hSignalEff[ich]->Integral(), hSignalNorm[ich]->Integral());
  }
  // fout->ls();
  // hEventPass->Print("all");
  std::vector<TString> passLabel;
  passLabel.resize(hEventPass->GetNbinsX());
  passLabel[0] = TString("all");
  passLabel[1] = TString("fid");
  passLabel[2] = TString("canTrig");
  passLabel[3] = TString("triggered");
  for (int ibin = 1; ibin <= hEventPass->GetNbinsX(); ++ibin)
    printf(" cut %i  %s  %.0f \n", ibin, passLabel[ibin - 1].Data(), hEventPass->GetBinContent(ibin));

  double trigRate = 0;
  double gammaRate = 105450. * 0.36;
  if (hEventPass->GetBinContent(1) > 0)
    trigRate = gammaRate * hEventPass->GetBinContent(4) / hEventPass->GetBinContent(1);

  // printf(" hTrigDiffTime10 entries %f \n", hTrigDiffTime10->GetEntries());

  printf("trigger photons summary: \n");
  for (int itrig = 0; itrig < 3; ++itrig)
    printf("chan %i sum photons %.3E integral in SPE %.3E gain %f\n", itrig + 9, hPhotonTrig[itrig]->Integral(),
           hSignalPhotons[itrig]->Integral(), hSignalPhotons[itrig]->Integral() / hPhotonTrig[itrig]->Integral());

  // multiply single sipm rate by gamma rate from source
  printf("trig sippm count %i rate %.3E Hz trig rate %.3E \n", trigCount9, double(trigCount9) / double(ngen) * gammaRate, trigRate);

  double photonTrigSum = 0;
  for (int i = 0; i < 3; ++i)
  {
    printf("trig NSPE count %i mean %.2f   \n", i + 9, hNumberSPE[i]->GetMean());
    photonTrigSum += hNumberSPE[i]->GetMean();
  }

  printf("\n\n ******* gamma peak photon sum %.3f numPhotons %.0f ratio %.3f  mean from histo %.3f ******* \n\n", photonTrigSum, numPhotons, photonTrigSum / double(numPhotons), hGammaPeak->GetMean());

  for (int ich = 0; ich < 13; ++ich)
    printf(" effScaleFactor[%i] = %.3E ;\n", ich, hGeoEff[ich]->GetMean());

  printf("efficiencies:\n");
  for (int ich = 0; ich < NCHAN; ++ich)
  {
    double eff = effGeoFunc(ich) * SiPMQE128Ham;
    printf("\t channel %i chanEff %.2E nominal eff %.2E ratio %.3E \n", ich, chanEff[ich], eff, chanEff[ich] / eff);
  }

  printf("********* end of btb with ngen %i triggers %i rawRun %i LY %.2f coincidence %.3f rate %.3f Hz ********\n", ngen, nTrigger, int(rawRun->btree->GetEntries()), LY, maxTriggerTimeDifference, trigRate);
  // simRun->print();
}

// static TBRun *theTBRun;
int main(int argc, char *argv[])
{
  /** geoVersion  **/
  geoVersionOld = true;

  int ngen = 100000;

  std::cout << "  usage: btbSim <ngen> default 1000000  " << argv[0] << std::endl;
  printf("\n ");
  if (argc > 1)
  {
    ngen = atoi(argv[1]);
  }
  double thePPM = 30;
  if (argc > 2)
  {
    thePPM = atof(argv[2]);
  }
  /* dont need this prints all the time
  double par[NPAR] = {numPhotons, 0, 0, 0};
  peakFitQsum[0] = 1.;
  peakFitQsum[1] = 2.;
  peakFitQsum[2] = 3.;
  triggerPeakFitShow = true;
  // set total photon yield
  peakFit(par);
  triggerPeakFitShow = false;
  */

  printf("***** START of btb geometry ngen = %.0E PPM %.03f *****\n", double(ngen), thePPM);
  btb(ngen, thePPM);
  printf("... end %s version %s ngen %i passed %lld total efficiency %.3f  file %s exit\n", argv[0], geoName.Data(), ngen, ntFit->GetEntries(), totalEventEffiency, fout->GetName());
  //  eventDir->ls();
  // printf("hGammaPeak entries %f \n", hGammaPeak->GetEntries());
  histDir->ls();
  fout->Write();
  fout->Close();
  exit(0);
}
