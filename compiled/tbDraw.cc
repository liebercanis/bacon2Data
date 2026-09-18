/**
 * @file tbDraw.cc
 * @date June 23 2026
 */
#include <iostream>
#include <fstream>
#include "TGraph.h"
#include "TMinuit.h"
#include "distanceLevels.hh"
#include "modelAllFit.hh"
#include "failCodes.hh"

// ============================================================================
//  GLOBAL CONSTANTS AND DEFINITIONS
// ============================================================================
// Time resolution: time is measured in microseconds
using namespace TMath;
// ============================================================================
//  GLOBAL OBJECTS AND DATA STRUCTURES
// ============================================================================
TFile *fin;            ///< Input ROOT file handle
TFile *fout;           ///< Output ROOT file handle for results
TNtuple *ntParScan;    ///< N-tuple storing parameter scan results
TMultiGraph *mgAbDist; ///< Multi-graph for absorption vs distance
TMultiGraph *mgAbPPM;  ///< Multi-graph for absorption vs PPM

std::vector<TH1D *> hnorm;               ///< Normalized detector response histograms per channel
std::vector<TH1D *> hcurve;              ///< Raw curve histograms per channel
std::vector<TH1D *> hffit;               ///< Fitted waveforms per channel
std::vector<TH1D *> hmodel;              ///< Model histograms per channel (reserved for future use)
std::vector<TH1D *> hfitModel;           ///< Final fitted model histograms per channel
std::vector<double> abValueByConstant1;  /// absorption by constant values
std::vector<double> abValueByConstant10; /// absorption by constant values
std::vector<double> singletByAbsorb;     /// singlet by absorption values
std::vector<double> absorbConstValue;    /// singlet by absorption values

// ============================================================================
//  ANALYSIS PARAMETERS
// ============================================================================
enum
{
  NPOINTS = 10000
};
int nominalTrigger = 729;  ///< Expected trigger timing bin
double singletStart = 700; ///< Singlet scintillation region start [ns]
double singletEnd = 750;   ///< Singlet scintillation region end [ns]

static Double_t vstart[NPARS]; ///< Starting parameter values for Minuit minimization
static Double_t step[NPARS];   ///< Step sizes for Minuit parameter exploration

/// Plotting colors for each detector channel visualization
int colors[NCHAN] = {kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kSpring, kTeal, kAzure, kViolet, kPink, kGray};
/// Legible palette for overlaying curves across absorbConst values
int absorbColors[] = {kRed, kGreen + 2, kBlue, kMagenta, kCyan + 2, kOrange + 7,
                      kSpring, kTeal + 2, kAzure + 2, kViolet, kPink + 9, kBlack};

static double AbsorptionByValue(double ppm, double dist, double absorb1Const = absorb1ConstDefault)
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
  absorb1Const = max(1.0E-9, absorb1Const); // avoid 0/0 -> NaN when dist is also 0
  double lambda1 = absorb1Const * 0.1 / ppm;
  double lambda2 = 3.38 * 0.1 / ppm;
  double lambda3 = 50.5 * 0.1 / ppm;
  double Tr128 = Aconstant * exp(-dist / lambda1) + Cconstant * exp(-dist / lambda2) + (1 - Aconstant - Cconstant) * exp(-dist / lambda3);
  // printf("Absorption: ppm %.3f dist %.3f A %.3E C %.3E lambda1 %.3f lambda2 %.3f lambda3 %.3f Tr128 %.3E \n", ppm, dist, Aconstant, Cconstant, lambda1, lambda2, lambda3, Tr128);
  return 1. - Tr128;
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
  std::cout << " fillFitWave " << ichan << "  " << hist->GetName() << std::endl;
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

double AbsorptionForConstant(double ppm, double dist, double absorb1Const = absorb1ConstDefault)
{
  /* new fit from Doug August 19*/
  ppm = max(1.0E-9, ppm);
  absorb1Const = max(1.0E-9, absorb1Const); // avoid 0/0 -> NaN when dist is also 0
  double lambda1 = absorb1Const * 0.1 / ppm;
  double lambda2 = 3.38 * 0.1 / ppm;
  double lambda3 = 50.5 * 0.1 / ppm;
  double Tr128 = Aconstant * exp(-dist / lambda1) + Cconstant * exp(-dist / lambda2) + (1 - Aconstant - Cconstant) * exp(-dist / lambda3);
  // printf("Absorption: ppm %.3f dist %.3f A %.3E C %.3E lambda1 %.3f lambda2 %.3f lambda3 %.3f Tr128 %.3E \n", ppm, dist, Aconstant, Cconstant, lambda1, lambda2, lambda3, Tr128);
  return 1. - Tr128;
}

// ============================================================================
//  MAIN FITTING ROUTINE
// ============================================================================

/**
 * @brief Main fitting routine: optimize model parameters to match detector data
 * @details Performs likelihood minimization using ROOT Minuit optimizer with chi-squared
 *          goodness-of-fit metric. Outputs fitted waveforms, component decomposition,
 *          parameter scan plots, and comparison canvases.
 * @param theFitChannel Channel index to fit. Use -1 to simultaneously fit all 12 PMT channels.
 *                     This enables global optimization of parameters shared across detectors.
 */
void tbDrawAbsorb(int theFitChannel = 7, double dopant = 10., double absorb1Const = absorb1ConstDefault, int icolor = 0)
{

  std::vector<double> distance(NPOINTS); //
  std::vector<double> abdist(NPOINTS);   //
  std::vector<double> ppm(NPOINTS);      //
  std::vector<double> ab(NPOINTS);       //
  double theDopant = dopant;
  for (int i = 0; i < NPOINTS; ++i)
  {
    distance[i] = double(i) * .005;
    abdist[i] = AbsorptionForConstant(theDopant, distance[i], absorb1Const);
    // printf("ppm %f dist %f A %.3E\n", theDopant, distance[i], ab[i]);
  }

  // draw absorption
  double dist = 1.; // cm
  for (int i = 0; i < NPOINTS; ++i)
  {
    ppm[i] = double(i) * .001;
    ab[i] = AbsorptionByValue(ppm[i], dist, absorb1Const);
  }
  abValueByConstant1.push_back(AbsorptionByValue(theDopant, 1., absorb1Const));
  abValueByConstant10.push_back(AbsorptionByValue(theDopant, 10., absorb1Const));

  TGraph *gAbDist = new TGraph(NPOINTS, &distance[0], &abdist[0]);
  gAbDist->SetMarkerStyle(21);
  gAbDist->SetMarkerSize(.7);
  gAbDist->SetName(Form("absorbtionVsDistancePPM%.3fAbsConst%.0E", theDopant, absorb1Const));
  gAbDist->SetTitle(Form("absorbtionVsDistancePPM%.3fAbsConst%.0E", theDopant, absorb1Const));
  gAbDist->GetXaxis()->SetTitle(Form("distance [cm] absorb1Const: %.3f", absorb1Const));
  gAbDist->SetLineColor(absorbColors[icolor % 12]);
  gAbDist->SetMarkerColor(absorbColors[icolor % 12]);
  TCanvas *cabDist = new TCanvas("absorbtionDist", "absorbtionDist");
  cabDist->SetGrid();
  cabDist->SetLogx();
  gAbDist->Draw("ap");
  mgAbDist->Add(gAbDist);
  fout->Add(gAbDist);

  TGraph *gAbPPM = new TGraph(NPOINTS, &ppm[0], &ab[0]);
  gAbPPM->SetMarkerStyle(21);
  gAbPPM->SetMarkerSize(.7);
  gAbPPM->SetName(Form("absorbtionVsPPMAt1cmAbsConst%.0E", absorb1Const));
  gAbPPM->SetTitle(Form("absorbtionVsPPMAt1cmAbsConst%.0E", absorb1Const));
  gAbPPM->GetXaxis()->SetTitle("dopant [PPM]");
  gAbPPM->GetYaxis()->SetTitle("absorbtion factor");
  gAbPPM->SetLineColor(absorbColors[icolor % 12]);
  gAbPPM->SetMarkerColor(absorbColors[icolor % 12]);
  TCanvas *cabPPM = new TCanvas("absorbtionPPM", "absorbtionPPM");
  cabPPM->SetGrid();
  cabPPM->SetLogx();
  gAbPPM->Draw("ap");
  mgAbPPM->Add(gAbPPM);
  fout->Add(gAbPPM);

  // ============================================================================
  //  ANALYSIS CONFIGURATION
  // ============================================================================

  // Initialize model parameters and optical properties from modelAllFit.hh
  setupModelAllFit();

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
  vstart[THECHANNEL] = theFitChannel;  ///< Channel selection flag
  vstart[ABSORB1CONST] = absorb1Const; ///< Absorption constant for xenon
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
  // for (int ip = 0; ip < NPARS; ++ip)
  //   printf(" par %i %s start val %.3f \n", ip, lparNames[ip].Data(), vstart[ip]);

  // copy into Minuit
  /* have to put some errors here otherwise it will be constant*/
  for (unsigned j = 0; j < NPARS; ++j)
  {
    step[j] = 1.E-6 * vstart[j];
    gMinuit->mnparm(j, lparNames[j].Data(), vstart[j], step[j], 0.1 * vstart[j], 10. * vstart[j], ierflg);
    lpar[j] = vstart[j];
  }

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
  fcn(llist, gin, fval, lpar, ierflg);
  printf(" after call to fcn starting value >>>>   fval %E \n", fval);
  double fvalStart = fval;
  if (isnan(fval))
  {
    printf("gMinuit returns NAN\n");
    return;
  }

  // fill fit wave

  // Initialize histogram vectors for all channels
  hffit.resize(NCHAN);
  hfitModel.resize(NCHAN);
  if (theFitChannel < 0)
  {
    for (unsigned ic = 0; ic < NCHANPMT; ++ic)
    {
      printf("fill fitWaveFitChan%i", ic);
      TH1D *hFit = new TH1D(Form("fitWaveFitChan%i", ic), Form("fitWaveFitChan%i", ic), MAXSAMPLE, 0, 2 * MAXSAMPLE);
      hFit->Reset("ICES");
      hFit->SetTitle((Form("fitWaveFitChan%i %.3fPPM", ic, dopant)));
      hFit->SetMarkerColor(colors[ic]);
      hFit->SetLineColor(colors[ic]);
      hfitModel[ic] = hFit;
      fillFitWave(ic, hFit);
    }
  }
  else // only 1 channel
  {
    printf("fill fitWaveFitChan%i", theFitChannel);
    TH1D *hFit = new TH1D(Form("fitWaveFitChan%iAbs%.0E", theFitChannel, absorb1Const), Form("fitWaveFitChan%iAbs%.0E", theFitChannel, absorb1Const), MAXSAMPLE, 0, 2 * MAXSAMPLE);
    hFit->Reset("ICES");
    hFit->SetTitle((Form("fitWaveFitChan%i %.3fPPM", theFitChannel, dopant)));
    hFit->SetMarkerColor(colors[theFitChannel]);
    hFit->SetLineColor(colors[theFitChannel]);
    hfitModel[theFitChannel] = hFit;
    fillFitWave(theFitChannel, hFit);
    int earlyBin = 1350 / 2;
    int lateBin = 1460 / 2;
    singletByAbsorb.push_back(hFit->Integral(earlyBin, lateBin));
  }

  // drawing
  // Store individual component contributions (singlet, triplet, mixed) for each channel
  if (theFitChannel < 0)
  {
    for (unsigned ic = 0; ic < NCHANPMT; ++ic)
    {
      for (int icomp = 0; icomp < NUMCOMP; ++icomp)
      {
        TH1D *hFit = new TH1D(Form("fit%sChan%iAbs%.0E", compNames[icomp].Data(), ic, absorb1Const), Form("fit%sChan%iAbs%.0E", compNames[icomp].Data(), ic, absorb1Const), MAXSAMPLE, 0, 2 * MAXSAMPLE);
        hFit->Reset("ICES");
        hFit->SetTitle((Form("fit%sChan%i %.3fPPM", compNames[icomp].Data(), ic, theDopant)));
        hFit->SetLineColor(colors[ic]);
        fillCompWave(ic, icomp, hFit); // only need one of these
      }
    }
  }
  else
  {
    for (int icomp = 0; icomp < NUMCOMP; ++icomp)
    {
      TH1D *hFit = new TH1D(Form("fit%sChan%iAbs%.0E", compNames[icomp].Data(), theFitChannel, absorb1Const), Form("fit%sChan%iAbs%.0E", compNames[icomp].Data(), theFitChannel, absorb1Const), MAXSAMPLE, 0, 2 * MAXSAMPLE);
      hFit->Reset("ICES");
      hFit->SetTitle((Form("fit%sChan%i %.3fPPM", compNames[icomp].Data(), theFitChannel, theDopant)));
      hFit->SetLineColor(colors[theFitChannel]);
      fillCompWave(theFitChannel, icomp, hFit);
    }
  }

  printf("\n...  finished tbDraw absorb %.0E \n", absorb1Const);
}

void tbDraw(int theFitChannel = 7, double dopant = 10)
{
  fout = new TFile(Form("tbDrawPPM%.2f.root", dopant), "recreate");
  // == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
  //  MINUIT OPTIMIZER INITIALIZATION
  // ============================================================================
  // Create Minuit minimizer with NPARS parameters
  TMinuit *gMinuit = new TMinuit(NPARS);
  gMinuit->SetFCN(fcn); // Set likelihood function pointer

  mgAbDist = new TMultiGraph("mgAbDist", "absorption vs distance, all absorbConst");
  mgAbPPM = new TMultiGraph("mgAbPPM", "absorption vs PPM, all absorbConst");

  int nAbsorbConst = 20; /// number of absorbConst values
  for (int i = 0; i < nAbsorbConst; ++i)
  {
    double ab = pow(2., i + 1);
    //  *absorb1ConstDefault;
    absorbConstValue.push_back(ab);
    printf("\n\n\t\t******* tbDraw: calling tbDrawAbsorb for absorbConst %.3E******", ab);
    tbDrawAbsorb(theFitChannel, dopant, ab, i);
  }

  // Draw multi-graphs for absorption vs distance and absorption vs PPM
  printf("Number of graphs in mgAbDist: %d\n", mgAbDist->GetListOfGraphs()->GetSize());
  TCanvas *canAbDist = new TCanvas("canAbDist", "absorption vs distance, all absorbConst");
  canAbDist->SetGrid();
  canAbDist->SetLogx();
  mgAbDist->GetXaxis()->SetTitle("distance [cm]");
  mgAbDist->GetYaxis()->SetTitle("absorption factor");
  mgAbDist->Draw("a");
  canAbDist->BuildLegend();
  canAbDist->Print("mgAbDist.pdf");
  fout->Add(mgAbDist);
  fout->Append(canAbDist);

  printf("Number of graphs in mgAbPPM: %d\n", mgAbPPM->GetListOfGraphs()->GetSize());
  TCanvas *canAbPPM = new TCanvas("canAbPPM", "absorption vs PPM, all absorbConst");
  canAbPPM->SetGrid();
  canAbPPM->SetLogx();
  mgAbPPM->GetXaxis()->SetTitle("dopant [PPM]");
  mgAbPPM->GetYaxis()->SetTitle("absorption factor");
  mgAbPPM->Draw("a");
  canAbPPM->BuildLegend();
  canAbPPM->Print("mgAbPPM.pdf");
  fout->Add(mgAbPPM);
  fout->Append(canAbPPM);

  if (singletByAbsorb.size() != absorbConstValue.size())
  {
    printf("tbDraw: size mismatch absorbConstValue %lu singletByAbsorb %lu, skipping graph\n",
           absorbConstValue.size(), singletByAbsorb.size());
    return;
  }
  TGraph *gSingletByAbsorb = new TGraph(absorbConstValue.size(), &absorbConstValue[0], &singletByAbsorb[0]);
  gSingletByAbsorb->SetName(Form("singletByAbsorbChan%iPPM%.3f", theFitChannel, dopant));
  gSingletByAbsorb->SetTitle(Form("singletByAbsorbChan%iPPM%.3f", theFitChannel, dopant));
  gSingletByAbsorb->SetMarkerStyle(21);
  gSingletByAbsorb->GetXaxis()->SetTitle("absorbConst");
  gSingletByAbsorb->GetYaxis()->SetTitle("singlet integral");
  TCanvas *canSingletByAbsorb = new TCanvas(Form("singletByAbsorbChan%iPPM%.3f", theFitChannel, dopant), "singlet integral vs absorbConst");
  canSingletByAbsorb->SetGrid();
  canSingletByAbsorb->SetLogx();
  gSingletByAbsorb->Draw("ap");
  canSingletByAbsorb->Print("singletByAbsorb.pdf");

  if (abValueByConstant1.size() != absorbConstValue.size())
  {
    printf("tbDraw: size mismatch absorbConstValue %lu abValueByConstant %lu, skipping graph\n",
           absorbConstValue.size(), abValueByConstant1.size());
    fout->Write();
    return;
  }

  // at dist 1 cm
  TGraph *gAbValueByConstant1 = new TGraph(absorbConstValue.size(), &absorbConstValue[0], &abValueByConstant1[0]);
  gAbValueByConstant1->SetName(Form("abValueByConstantAt1cmPPM%.3f", dopant));
  gAbValueByConstant1->SetTitle(Form("abValueByConstantAt1cmPPM%.3f", dopant));
  gAbValueByConstant1->SetMarkerStyle(21);
  gAbValueByConstant1->GetXaxis()->SetTitle("absorbConst");
  gAbValueByConstant1->GetYaxis()->SetTitle(Form("absorption factor at %.3f PPM @ 1 cm", dopant));
  TCanvas *canAbValueByConstant1 = new TCanvas("canAbValueByConstant1", "absorption factor vs absorbConst");
  canAbValueByConstant1->SetGrid();
  canAbValueByConstant1->SetLogx();
  gAbValueByConstant1->Draw("ap");
  canAbValueByConstant1->Print("abValueByConstant1.pdf");
  fout->Add(gAbValueByConstant1);
  fout->Append(canAbValueByConstant1);

  // at dist 10 cm
  TGraph *gAbValueByConstant10 = new TGraph(absorbConstValue.size(), &absorbConstValue[0], &abValueByConstant10[0]);
  gAbValueByConstant10->SetName(Form("abValueByConstantAt10cmPPM%.3f", dopant));
  gAbValueByConstant10->SetTitle(Form("abValueByConstantAt10cmPPM%.3f", dopant));
  gAbValueByConstant10->SetMarkerStyle(21);
  gAbValueByConstant10->GetXaxis()->SetTitle("absorbConst");
  gAbValueByConstant10->GetYaxis()->SetTitle(Form("absorption factor at %.3f PPM @ 10 cm", dopant));
  TCanvas *canAbValueByConstant10 = new TCanvas("canAbValueByConstant10", "absorption factor vs absorbConst");
  canAbValueByConstant10->SetGrid();
  canAbValueByConstant10->SetLogx();
  gAbValueByConstant10->Draw("ap");
  canAbValueByConstant10->Print("abValueByConstant10.pdf");
  fout->Add(gAbValueByConstant10);
  fout->Append(canAbValueByConstant10);

  fout->Write();

  printf("\n...  finished tbDraw \n");
}
