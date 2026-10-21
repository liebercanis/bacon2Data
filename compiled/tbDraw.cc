/**
 * @file tbDraw.cc
 * @date June 23 2026
 */
#include <iostream>
#include <fstream>
#include "TGraph.h"
#include "TMinuit.h"
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

std::vector<TH1D *> hnorm;                                 ///< Normalized detector response histograms per channel
std::vector<TH1D *> hlevel;                                ///< Normalized detector response histograms per channel
std::vector<TH1D *> hcurve;                                ///< Raw curve histograms per channel
std::vector<TH1D *> hffit;                                 ///< Fitted waveforms per channel
std::vector<TH1D *> hmodel;                                ///< Model histograms per channel (reserved for future use)
std::vector<TH1D *> hfitModel;                             ///< Final fitted model histograms per channel
std::vector<double> abValueByConstant1;                    /// absorption by constant values
std::vector<double> abValueByConstant10;                   /// absorption by constant values
std::vector<double> singletByAbsorb;                       /// singlet by absorption values
std::vector<double> absorbConstValue;                      /// singlet by absorption values
std::vector<std::vector<double>> singletByDopantAndAbsorb; /// singlet by dopant, absorption values
std::vector<double> ppmFile;                               /// ppm values for files

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

int colors[NCHAN] = {kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kSpring, kTeal, kAzure, kViolet, kPink, kGray};
/// Legible palette for overlaying curves across absorbConst values
int absorbColors[] = {kRed, kGreen + 2, kBlue, kMagenta, kCyan + 2, kOrange + 7,
                      kSpring, kTeal + 2, kAzure + 2, kViolet, kPink + 9, kBlack};

/**
 * @brief Populate histogram with fitted waveform samples from optimization result
 * @details Transfers fitWave array (computed by fcn() during minimization) into a TH1D
 *          histogram. Data is stored in modelAllFit.hh as fitWave[NCHAN][MAXSAMPLE]
 * @param ichan Levelnel index
 * @param hist Pointer to histogram to fill with fitted values
 */

// ============================================================================
//  MAIN FITTING ROUTINE
// ============================================================================

/**
 * @brief Main fitting routine: optimize model parameters to match detector data
 * @details Performs likelihood minimization using ROOT Minuit optimizer with chi-squared
 *          goodness-of-fit metric. Outputs fitted waveforms, component decomposition,
 *          parameter scan plots, and comparison canvases.
 * @param  theFitLevelLevel index to fit. Use -1 to simultaneously fit all 12 PMT channels.
 *                     This enables global optimization of parameters shared across detectors.
 */
void tbDrawAbsorb(int theFitLevel = 0, double dopant = 10., double lambda1Const = lambda1ConstDefault, int jdop = 0, int iabs = 0)
{

  std::vector<double> distance(NPOINTS); //
  std::vector<double> abdist(NPOINTS);   //
  std::vector<double> ppm(NPOINTS);      //
  std::vector<double> ab(NPOINTS);       //
  double theDopant = dopant;
  for (int i = 0; i < NPOINTS; ++i)
  {
    distance[i] = double(i) * .005;
    abdist[i] = AbsorptionByValue(theDopant, distance[i], lambda1Const);
    // printf("ppm %f dist %f A %.3E\n", theDopant, distance[i], ab[i]);
  }

  // draw absorption
  double dist = 1.; // cm
  for (int i = 0; i < NPOINTS; ++i)
  {
    ppm[i] = double(i) * .001;
    ab[i] = AbsorptionByValue(ppm[i], dist, lambda1Const);
  }
  abValueByConstant1.push_back(AbsorptionByValue(theDopant, 1., lambda1Const));
  abValueByConstant10.push_back(AbsorptionByValue(theDopant, 10., lambda1Const));

  TGraph *gAbDist = new TGraph(NPOINTS, &distance[0], &abdist[0]);
  gAbDist->SetMarkerStyle(21);
  gAbDist->SetMarkerSize(.7);
  gAbDist->SetName(Form("absorbtionVsDistancePPM%.3fAbsConst%.0E", theDopant, lambda1Const));
  gAbDist->SetTitle(Form("absorbtionVsDistancePPM%.3fAbsConst%.0E", theDopant, lambda1Const));
  gAbDist->GetXaxis()->SetTitle(Form("distance [cm] lambda1Const: %.3f", lambda1Const));
  gAbDist->SetLineColor(absorbColors[iabs % 12]);
  gAbDist->SetMarkerColor(absorbColors[iabs % 12]);
  TCanvas *cabDist = new TCanvas(Form("cabAbsorbtionVsDistanceJdop%iIabs%i", jdop, iabs), Form("canAsorbtionVsDistancePPM%.3fAbsConst%.0E", theDopant, lambda1Const));
  cabDist->SetGrid();
  cabDist->SetLogx();
  gAbDist->Draw("ap");
  mgAbDist->Add(gAbDist);
  fout->Add(gAbDist);

  TGraph *gAbPPM = new TGraph(NPOINTS, &ppm[0], &ab[0]);
  gAbPPM->SetMarkerStyle(21);
  gAbPPM->SetMarkerSize(.7);
  gAbPPM->SetName(Form("absorbtionVsPPMAt1cmAbsConst%.0E", lambda1Const));
  gAbPPM->SetTitle(Form("absorbtionVsPPMAt1cmAbsConst%.0E", lambda1Const));
  gAbPPM->GetXaxis()->SetTitle("dopant [PPM]");
  gAbPPM->GetYaxis()->SetTitle("absorbtion factor");
  gAbPPM->SetLineColor(absorbColors[iabs % 12]);
  gAbPPM->SetMarkerColor(absorbColors[iabs % 12]);
  TCanvas *cabPPM = new TCanvas(Form("canAsorbtionVsPPMAt1cmAbsConstJdop%iIabs%i", jdop, iabs), Form("canAbsorbtionVsPPMAt1cmAbsConst%.0E", lambda1Const));
  cabPPM->SetGrid();
  cabPPM->SetLogx();
  gAbPPM->Draw("ap");
  mgAbPPM->Add(gAbPPM);
  fout->Add(gAbPPM);

  // ============================================================================
  //  ANALYSIS CONFIGURATION
  // ============================================================================

  /* total photons per event LY defined in modelAllFit.hh */
  double startNorm = 60. * LY;
  Double_t arglist[10];
  int ierflg = 0;

  // ============================================================================
  //  FIT PARAMETER INITIALIZATION
  // ============================================================================
  // Set initial guesses for all fit parameters
  lpar[NORM] = startNorm;     ///< Photon yield per event
  lpar[TRIGSTART] = 2. * 700; ///< Trigger timing offset
  lpar[SFRAC] = 0.14;         ///< Singlet fraction (ref: Segretto 2021)
  lpar[PPM] = theDopant;      ///< Dopant concentration
  lpar[TAU3] = tTriplet0;     ///< Triplet decay time
  lpar[TAUM] = 4700.0;        ///< Mixed component decay time
  lpar[BKGCONST] = 4.0E-6;    ///< Constant background rate
  lpar[KXCONST] = 1.0;
  lpar[THELEVEL] = theFitLevel;      ///< Levelnel selection flag
  lpar[LAMBDA1CONST] = lambda1Const; ///< Absorption constant for xenon

  // Initialize model parameters and optical properties from modelAllFit.hh
  setupModelAllFit();

  // setup minuit from modelAllFit.hh
  setupMinuit();
  /*
  printf("starting parameter values \n");
  for (int ip = 0; ip < NPARS; ++ip)
    printf(" par %i %s start val %.3f \n", ip, lparNames[ip].Data(), lpar[ip]);
    */
  arglist[0] = KXCONST + 1; // par
  // gMinuit->mnexcm("FIX", arglist, 1, ierflg);
  arglist[1] = 0.1 * lpar[KXCONST]; // low
  arglist[2] = 10. * lpar[KXCONST]; // high
  gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);

  // arglist[0] = SFRAC + 1; // kp
  // gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  // Set bounds for variable parameters to restrict optimization domain
  arglist[0] = NORM + 1;          // par
  arglist[1] = 0.01 * lpar[NORM]; // low
  arglist[2] = 10. * lpar[NORM];  // high
  gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);
  // gMinuit->mnexcm("FIX", arglist, 3, ierflg);

  arglist[0] = SFRAC + 1;   // par
  arglist[1] = lpar[SFRAC]; // low
  arglist[2] = lpar[SFRAC];
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

  arglist[0] = LAMBDA1CONST + 1;    // par tau mixed
  arglist[1] = 0.01 * lambda1Const; // low
  arglist[2] = 10. * lambda1Const;  // high
  // gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);
  printf("\t fix LAMBDA1CONST %.0f \n", arglist[0]);
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);

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
  if (theFitLevel < 0)
  {
    for (unsigned ilevel = 0; ilevel < NLEVEL; ++ilevel)
    {
      printf("fill fitWaveFitLevel%i", ilevel);
      TH1D *hFit = new TH1D(Form("fitWaveFitLevel%ijdop%iiabs%i", theFitLevel, jdop, iabs), Form("fitWaveFitLevel%ilambda1%.3E", theFitLevel, lambda1Const), MAXSAMPLE, 0, 2 * MAXSAMPLE);

      hFit->Reset("ICES");
      hFit->SetTitle((Form("fitWaveFitLevel%i %.3fPPM", ilevel, dopant)));
      hFit->SetMarkerColor(colors[ilevel]);
      hFit->SetLineColor(colors[ilevel]);
      hfitModel[ilevel] = hFit;
      fillFitWave(ilevel, hFit);
      fout->Add(hFit);
    }
  }
  else // only 1 channel
  {
    printf("fill fitWaveFitLevel%i", theFitLevel);
    TH1D *hFit = new TH1D(Form("fitWaveFitLevel%ijdop%iiabs%i", theFitLevel, jdop, iabs), Form("fitWaveFitLevel%iAbs%.3E", theFitLevel, lambda1Const), MAXSAMPLE, 0, 2 * MAXSAMPLE);
    hFit->Reset("ICES");
    hFit->SetTitle((Form("fitWaveFitLevel%i %.3fPPM", theFitLevel, dopant)));
    hFit->SetMarkerColor(colors[theFitLevel]);
    hFit->SetLineColor(colors[theFitLevel]);
    hfitModel[theFitLevel] = hFit;
    fillFitWave(theFitLevel, hFit);
    singletByDopantAndAbsorb[jdop][iabs] = hFit->Integral(earlyBin, lateBin);

    fout->Add(hFit);
  }

  // drawing
  // Store individual component contributions (singlet, triplet, mixed) for each channel
  if (theFitLevel < 0)
  {
    for (unsigned ic = 0; ic < NLEVEL; ++ic)
    {
      for (int icomp = 0; icomp < NUMCOMP; ++icomp)
      {
        TH1D *hFit = new TH1D(Form("fit%sLevel%ijsop%iiabs%iAbs%.0E", compNames[icomp].Data(), ic, jdop, iabs, lambda1Const), Form("fit%sLevel%iAbs%.0E", compNames[icomp].Data(), ic, lambda1Const), MAXSAMPLE, 0, 2 * MAXSAMPLE);
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
      TH1D *hFit = new TH1D(Form("fit%sLevel%ijsop%iiabs%iAbs%.0E", compNames[icomp].Data(), theFitLevel, jdop, iabs, lambda1Const), Form("fit%sLevel%iAbs%.0E", compNames[icomp].Data(), theFitLevel, lambda1Const), MAXSAMPLE, 0, 2 * MAXSAMPLE);
      hFit->Reset("ICES");
      hFit->SetTitle((Form("fit%sLevel%i %.3fPPM", compNames[icomp].Data(), theFitLevel, theDopant)));
      hFit->SetLineColor(colors[theFitLevel]);
      fillCompWave(theFitLevel, icomp, hFit);
      // if (icomp == 0)
      //   singletByDopantAndAbsorb[jdop][iabs] = hFit->Integral(earlyBin, lateBin);
    }
  }

  printf("\n...  finished tbDraw absorb %.0E jdop %i iabs %i \n", lambda1Const, jdop, iabs);
}

void tbDraw(int theFitLevel = 0, double dopant = 10)
{
  fout = new TFile(Form("tbDrawPPM%.2f.root", dopant), "recreate");
  // == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
  //  MINUIT OPTIMIZER INITIALIZATION
  // ============================================================================
  // Create Minuit minimizer with NPARS parameters

  mgAbDist = new TMultiGraph("mgAbDist", "absorption vs distance, all absorbConst");
  mgAbPPM = new TMultiGraph("mgAbPPM", "absorption vs PPM, all absorbConst");

  ppmFile.push_back(0.001);
  ppmFile.push_back(0.01);
  ppmFile.push_back(0.03);
  ppmFile.push_back(0.05);
  ppmFile.push_back(0.1);
  ppmFile.push_back(0.3);
  ppmFile.push_back(0.5);
  ppmFile.push_back(1.);
  ppmFile.push_back(2.);
  ppmFile.push_back(5.);
  ppmFile.push_back(10.);
  ppmFile.push_back(15.);
  ppmFile.push_back(30.);

  std::vector<double> dopantVals;
  int nAbsorbConst = 10;        /// number of absorbConst values
  int nDopant = ppmFile.size(); /// number of dopant values
  singletByDopantAndAbsorb.resize(nDopant, std::vector<double>(nAbsorbConst, 0.));
  /* loop over dopants */
  for (int jdop = 0; jdop < nDopant; ++jdop)
  {
    double dopantj = ppmFile[jdop]; // dopant values from the ppmFile vector
    dopantVals.push_back(dopantj);

    /* loop over absorption  */
    double lambda1ConstMin = 0.1 * lambda1ConstDefault;
    double lambda1ConstMax = 10.0 * lambda1ConstDefault;
    double absStep = (lambda1ConstMax - lambda1ConstMin) / double(nAbsorbConst);
    for (int iabs = 0; iabs < nAbsorbConst; ++iabs)
    {
      double ab = lambda1ConstMin + double(iabs) * absStep;
      //  *absorb1ConstDefault;
      absorbConstValue.push_back(ab);
      printf("\n\n\t\t******* tbDraw: calling tbDrawAbsorb for iabs = %d, absorbConst %.3E dopant %.3f ******", iabs, ab, dopantj);
      tbDrawAbsorb(theFitLevel, dopantj, ab, jdop, iabs);
    }
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

  // singlet integral vs dopant, one curve per absorbConst value
  TMultiGraph *mgSingletByDopant = new TMultiGraph("mgSingletByDopant", "singlet integral vs dopant");

  std::vector<double> singletVals;
  for (int iabs = 0; iabs < nAbsorbConst; ++iabs)
  {
    // fill singlet vals by PPM for this iabs value
    singletVals.clear();
    for (unsigned jdop = 0; jdop < dopantVals.size(); ++jdop)
      singletVals.push_back(singletByDopantAndAbsorb[jdop][iabs]);
    //
    TGraph *gSingletByDopant = new TGraph(dopantVals.size(), &dopantVals[0], &singletVals[0]);
    gSingletByDopant->SetName(Form("gSingletByDopantIabs%i", iabs));
    gSingletByDopant->SetTitle(Form("lambda1=%.0E", absorbConstValue[iabs]));
    gSingletByDopant->SetLineColor(absorbColors[iabs % 12]);
    gSingletByDopant->SetMarkerColor(absorbColors[iabs % 12]);
    gSingletByDopant->SetMarkerStyle(21);
    fout->Add(gSingletByDopant);
    mgSingletByDopant->Add(gSingletByDopant, "lp");
  }

  TCanvas *canSingletByDopant = new TCanvas("canSingletByDopant", "singlet integral vs dopant");
  canSingletByDopant->SetGrid();
  mgSingletByDopant->GetXaxis()->SetTitle("dopant [PPM]");
  mgSingletByDopant->GetYaxis()->SetTitle(Form("singlet integral level %i", theFitLevel));
  canSingletByDopant->SetLogx();
  mgSingletByDopant->Draw("a");
  canSingletByDopant->BuildLegend();
  canSingletByDopant->Print("singletByDopant.pdf");
  fout->Add(mgSingletByDopant);
  fout->Append(canSingletByDopant);

  // fout->ls();
  fout->Write();

  printf("\n...  finished tbDraw \n");
}
