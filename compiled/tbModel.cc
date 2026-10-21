/**
 * @file tbModel.cc
 * @date June 23 2026
 */
#include <iostream>
#include <fstream>
#include "TGraph.h"
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
TFile *fin;         ///< Input ROOT file handle
TFile *fout;        ///< Output ROOT file handle for results
TNtuple *ntParScan; ///< N-tuple storing parameter scan results

std::vector<TH1D *> hnorm;     ///< Normalized detector response histograms per channel
std::vector<TH1D *> hcurve;    ///< Raw curve histograms per channel
std::vector<TH1D *> hffit;     ///< Fitted waveforms per channel
std::vector<TH1D *> hmodel;    ///< Model histograms per channel (reserved for future use)
std::vector<TH1D *> hfitModel; ///< Final fitted model histograms per channel
std::vector<double> ppmFile;
std::vector<std::vector<std::vector<TH1D *>>> hCompWaves; ///< Component waveforms per channel and component type
std::vector<std::vector<double>> singletByFile;           /// full waveforms per channel

double theLambda1Constant = lambda1ConstDefault; ///< Absorption constant for current analysis

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
double theTrigStart = 1400;

/// Plotting colors for each detector channel visualization
int colors[NCHAN] = {kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kSpring, kTeal, kAzure, kViolet, kPink, kGray};

void plotComponents(int ifile, int ilevel)
{

  TCanvas *can = new TCanvas(Form("compLevel%iFile%i", ilevel, ifile), Form("compLevel%iFile%i", ilevel, ifile));
  for (int icomp = 0; icomp < NUMCOMP; ++icomp)
  {
    hCompWaves[ifile][ilevel][icomp]->GetYaxis()->SetRangeUser(1.e-5, 10.);
    hCompWaves[ifile][ilevel][icomp]->SetLineColor(colors[icomp]);
    hCompWaves[ifile][ilevel][icomp]->Draw(icomp == 0 ? "hist" : "hist same");
  }
  can->SetLogy();
  can->BuildLegend();
  can->Print(".pdf");
}

// ============================================================================
//  MAIN FITTING ROUTINE
// ============================================================================

/**
 * @brief Main fitting routine: optimize model parameters to match detector data
 * @details Performs likelihood minimization using ROOT Minuit optimizer with chi-squared
 *          goodness-of-fit metric. Outputs fitted waveforms, component decomposition,
 *          parameter scan plots, and comparison canvases.
 * @param theLevel level index to fit. Use -1 to simultaneously fit all levels.
 *                     This enables global optimization of parameters shared across detectors.
 */
void tbDraw(int theLevel = 7, int ifile = 0)
{ //
  // draw absorption
  double dopant = ppmFile[ifile];
  double dist = 10.; // cm
  double ppm[NPOINTS];
  double ab[NPOINTS];
  for (int i = 0; i < NPOINTS; ++i)
  {
    ppm[i] = double(i + 1) * .001;
    ab[i] = Absorption(ppm[i], dist);
    // printf("ppm %f dist %f A %.3E\n", ppm[i], dist, ab[i]);
  }

  TGraph *gAb = new TGraph(NPOINTS, &ppm[0], &ab[0]);
  gAb->SetMarkerStyle(21);
  gAb->SetMarkerSize(.7);
  gAb->SetName("absorbtion");
  gAb->SetTitle("absorbtion model");
  gAb->GetXaxis()->SetTitle("PPM");
  gAb->GetYaxis()->SetTitle(Form("absorbtion factor at  %.0f cm", dist));
  TCanvas *cab = new TCanvas("absorbtion", "absorbtion");
  cab->SetGrid();
  gAb->Draw("ap");

  double distance[NPOINTS]; // cm
  double abdist[NPOINTS];
  double theDopant = dopant;
  for (int i = 0; i < NPOINTS; ++i)
  {
    distance[i] = double(i) * .005;
    abdist[i] = Absorption(theDopant, distance[i]);
    // printf("ppm %f dist %f A %.3E\n", theDopant, distance[i], ab[i]);
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

  // ============================================================================
  //  ANALYSIS CONFIGURATION
  // ============================================================================

  fout = new TFile(Form("tbModelPPM%i.root", ifile), "recreate");

  std::vector<unsigned> badList;
  setBadChannels(badList);

  /* total photons per event LY defined in modelAllFit.hh */
  // ============================================================================
  //  FIT PARAMETER INITIALIZATION
  // ============================================================================
  // Set initial guesses for all fit parameters
  double startNorm = 60. * LY;
  printf("trigstart %f\n", theTrigStart);
  /**** Set initial guesses for all fit parameters ****/
  lpar[NORM] = startNorm;         ///< Photon yield per event
  lpar[TRIGSTART] = theTrigStart; ///< Trigger timing offset
  // for Gamma I_s/I_t=0.3
  lpar[SFRAC] = 0.23;     ///< Singlet fraction (ref: Segretto 2021)
  lpar[PPM] = dopant;     ///< Dopant concentration
  lpar[TAU3] = tTriplet0; ///< Triplet decay time
  lpar[TAUM] = 4700.0;    ////< Mixed component decay time
  lpar[BKGCONST] = 0.0;   ///< Constant background rate
  lpar[KXCONST] = 1.0;    ///< Rate of transfer to mixed state
  // radiative
  lpar[R2CONST] = r2ConstDefault;
  lpar[R3CONST] = r3ConstDefault;
  lpar[C1CONST] = C1ConstDefault;
  lpar[LAMBDA1CONST] = theLambda1Constant;

  // setup minuit from modelAllFit.hh
  setupMinuit();

  // == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
  //  PARAMETER CONSTRAINTS AND BOUNDARIES
  // ============================================================================
  // Fix parameters that are held constant during minimization
  // Note: Minuit uses 1-based indexing for parameters (adds 1 to C++ indices)

  // fix level
  /*
  arglist[0] = THELEVEL + 1; // level
  printf("fix level %.0f \n", arglist[0]);
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);
  */
  Double_t arglist[10];
  int ierflg = 0;
  arglist[0] = TRIGSTART + 1; // trigger
  printf("fix TRIGSTART %.0f \n", arglist[0]);
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);
  // arglist[1] = 0.01 * lpar[TRIGSTART]; // low
  // arglist[2] = 10. * lpar[TRIGSTART];  // high
  // gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);

  arglist[0] = BKGCONST + 1; // par
  printf("fix BKGCONST %.0f \n", arglist[0]);
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  arglist[0] = KXCONST + 1; // par
  printf("\t KXCONST %.0f \n", arglist[0]);
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
  printf("\t limit NORM %.0f \n", arglist[0]);
  gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);
  // gMinuit->mnexcm("FIX", arglist, 3, ierflg);

  arglist[0] = SFRAC + 1;          // par
  arglist[1] = lpar[SFRAC];        // low
  arglist[2] = lpar[SFRAC] * 1.01; // high
  printf("\t FIX SFRAC %.0f \n", arglist[0]);
  // gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  // set limits ... here par starts with 1 so add 1
  arglist[0] = TAU3 + 1;         // par
  arglist[1] = 0.01 * tTriplet0; // low
  arglist[2] = 2.0 * tTriplet0;  // high
  // gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);
  printf("\t fix TAU3 %.0f \n", arglist[0]);
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  // set limits ... here par starts with 1 so add 1
  arglist[0] = PPM + 1; // par
  arglist[1] = 0.0;     // low
  arglist[2] = 50.0;    // high
  printf("\t set PPM %.0f \n", arglist[0]);
  gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);
  // gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  arglist[0] = TAUM + 1;     // par tau mixed
  arglist[1] = 0.01 * tMix0; // low
  arglist[2] = 10. * tMix0;  // high
  // gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);
  printf("\t fix TAUM %.0f \n", arglist[0]);
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  arglist[0] = LAMBDA1CONST + 1;           // par tau mixed
  arglist[1] = 0.01 * lambda1ConstDefault; // low
  arglist[2] = 10. * lambda1ConstDefault;  // high
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
  printf(" starting value >>>>   fval %E \n", fval);
  double fvalStart = fval;
  if (isnan(fval))
  {
    printf("gMinuit returns NAN\n");
    return;
  }

  // fill fit wave and make singlet graph
  vector<double> singletValues;
  vector<double> levelNumber;

  // Initialize histogram vectors for all levels
  printf("\t\t start filling histograms \n");
  hffit.resize(NLEVEL);
  hfitModel.resize(NLEVEL);

  if (theLevel < 0)
  {
    for (unsigned ic = 0; ic < NLEVEL; ++ic)
    {
      // printf("fill fitWaveFitLevel%i", ic);
      TH1D *hFit = new TH1D(Form("fitWaveFitLevel%i", ic), Form("fitWaveFitLevel%i", ic), MAXSAMPLE, 0, 2 * MAXSAMPLE);
      hFit->Reset("ICES");
      hFit->SetTitle((Form("fitWaveFitLevel%i %.3fPPM", ic, dopant)));
      hFit->SetMarkerColor(colors[ic]);
      hFit->SetLineColor(colors[ic]);
      hfitModel[ic] = hFit;
      // printf("hfitModel name %s \n", hFit->GetName());
      fillFitWave(ic, hFit);
      levelNumber.push_back(double(ic));
      singletValues.push_back(hFit->Integral(earlyBin, lateBin));
      singletByFile[ic].push_back(singletValues.back());
    }
    // make singlet graph
    TGraph *singletGraph = new TGraph(singletValues.size(), &levelNumber[0], &singletValues[0]);
    singletGraph->SetName(Form("singletGraphPPM%i", ifile));
    singletGraph->SetTitle(Form("singletGraphPPM%i", ifile));
    fout->Add(singletGraph);
  }
  else // only 1 channel
  {
    // printf("fill fitWaveFitLevel%i", theLevel);
    TH1D *hFit = new TH1D(Form("fitWaveFitLevel%i", theLevel), Form("fitWaveFitLevel%i", theLevel), MAXSAMPLE, 0, 2 * MAXSAMPLE);
    hFit->Reset("ICES");
    hFit->SetTitle((Form("fitWaveFitLevel%i %.3fPPM", theLevel, dopant)));
    hFit->SetMarkerColor(colors[theLevel]);
    hFit->SetLineColor(colors[theLevel]);
    hfitModel[theLevel] = hFit;
    // printf("hfitModel name %s \n", hFit->GetName());
    fillFitWave(theLevel, hFit);
  }

  // fill individual component contributions (singlet, triplet, mixed) for each level

  if (theLevel < 0)
  {
    for (unsigned ic = 0; ic < NLEVEL; ++ic)
    {
      for (int icomp = 0; icomp < NUMCOMP; ++icomp)
      {
        TH1D *hFit = new TH1D(Form("fit%sLevel%i", compNames[icomp].Data(), ic), Form("fit%sLevel%i", compNames[icomp].Data(), ic), MAXSAMPLE, 0, 2 * MAXSAMPLE);
        hFit->Reset("ICES");
        hFit->SetTitle((Form("fit%sLevel %i %.3fPPM", compNames[icomp].Data(), ic, theDopant)));
        hFit->SetLineColor(colors[ic]);
        fillCompWave(ic, icomp, hFit); // only need one of these
        hCompWaves[ifile][ic].push_back(hFit);
        /*
        if (icomp == 0)
        {
          singletValues.push_back(hFit->Integral(earlyBin, lateBin));
          singletByFile[ic].push_back(singletValues.back());
        }
        */
      }
    }
  }
  else
  {
    for (int icomp = 0; icomp < NUMCOMP; ++icomp)
    {
      TH1D *hFit = new TH1D(Form("fit%sLevel%i", compNames[icomp].Data(), theLevel), Form("fit%sLevel%i", compNames[icomp].Data(), theLevel), MAXSAMPLE, 0, 2 * MAXSAMPLE);
      hFit->Reset("ICES");
      hFit->SetTitle((Form("fit%sLevel%i %.3fPPM", compNames[icomp].Data(), theLevel, theDopant)));
      hFit->SetLineColor(colors[theLevel]);
      fillCompWave(theLevel, icomp, hFit);
      hCompWaves[ifile][theLevel].push_back(hFit);
      // printf("hFit name %s \n", hFit->GetName());
      //  printf("hFit name %s \n", hCompWaves[ifile][theLevel].back()->GetName());
    }
  }

  fout->Write();

  printf("\n...  finished tbDraw  file %i \n", ifile);
}
void tbModel()
{

  // Initialize model parameters and optical properties from modelAllFit.hh
  setupModelAllFit();

  // ppm of file

  ppmFile.push_back(1.E-3);
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

  std::vector<double> abValues;

  for (int ifile = 0; ifile < ppmFile.size(); ++ifile)
  {
    abValues.push_back(Absorption(ppmFile[ifile], 10.)); // dist = 10 cm, matches tbDraw()'s reference distance
    for (int ilevel = 0; ilevel < NLEVEL; ++ilevel)
    {
      double ab = Absorption(ppmFile[ifile], distanceLevel[ilevel]);
      // printf("ppm %f dist %f A %.3E\n", ppmFile[ifile], distanceLevel[ilevel], ab);
    }
  }

  hCompWaves.resize(ppmFile.size());
  for (size_t ifile = 0; ifile < ppmFile.size(); ++ifile)
  {
    hCompWaves[ifile].resize(NCHAN);
  }

  singletByFile.resize(NCHAN);

  printf("\t number of file %li lambda1Constant = %.3E\n", ppmFile.size(), lpar[LAMBDA1CONST]);
  gStyle->SetOptStat(0);
  for (unsigned ifile = 0; ifile < ppmFile.size(); ++ifile)
  {
    printf("tbDraw ifile %i \n", ifile);
    tbDraw(-1, ifile);
    for (int ichan = 0; ichan < NLEVEL; ++ichan)
    {
      /*for (int icomp = 0; icomp < hCompWaves[ifile][ichan].size(); ++icomp)
      {
        if (hCompWaves[ifile][ichan][icomp])
        {
          printf("tbDraw file %i chan %i comp %i name %s \n", ifile, ichan, icomp, hCompWaves[ifile][ichan][icomp]->GetName());
        }
        else
          printf(" nocompwave ifile %i chan %i comp %i \n", ifile, ichan, icomp);
      }*/
      plotComponents(ifile, ichan); // ifile,ichan
    }
  }

  printf("\t number of file %li lambda1Constant = %.3E\n", ppmFile.size(), lpar[LAMBDA1CONST]);

  // singlet integral vs file (PPM), one graph per level
  TMultiGraph *mgSingletByFile = new TMultiGraph("mgSingletByFile", "singlet integral vs PPM, all levels");
  for (unsigned ic = 0; ic < NLEVEL; ++ic)
  {
    if (singletByFile[ic].size() != ppmFile.size())
    {
      printf("mgSingletByFile: level %u size mismatch ppmFile %lu singletByFile %lu, skipping\n",
             ic, ppmFile.size(), singletByFile[ic].size());
      continue;
    }
    TGraph *gSingletChan = new TGraph(ppmFile.size(), &ppmFile[0], &singletByFile[ic][0]);
    gSingletChan->SetName(Form("gSingletByFileLevel%u", ic));
    gSingletChan->SetTitle(Form("level %u", ic));
    gSingletChan->SetLineColor(colors[ic]);
    gSingletChan->SetMarkerColor(colors[ic]);
    gSingletChan->SetMarkerStyle(21);
    mgSingletByFile->Add(gSingletChan, "lp");
    fout->Add(gSingletChan);
  }
  TCanvas *canSingletByFile = new TCanvas("canSingletByFile", "singlet integral vs PPM, all channels");
  canSingletByFile->SetGrid();
  canSingletByFile->SetLogx();
  mgSingletByFile->GetXaxis()->SetTitle("PPM");
  mgSingletByFile->GetYaxis()->SetTitle("singlet integral");
  mgSingletByFile->Draw("a");
  canSingletByFile->BuildLegend();
  canSingletByFile->Print("singletByFile.pdf");
  fout->Add(mgSingletByFile);
  fout->Append(canSingletByFile);

  // absorption value vs PPM
  if (abValues.size() != ppmFile.size())
  {
    printf("abValues: size mismatch ppmFile %lu abValues %lu, skipping graph\n", ppmFile.size(), abValues.size());
    return;
  }
  TGraph *gAbValues = new TGraph(ppmFile.size(), &ppmFile[0], &abValues[0]);
  gAbValues->SetName("gAbValues");
  gAbValues->SetTitle("absorption factor vs PPM (10 cm)");
  gAbValues->SetMarkerStyle(21);
  gAbValues->GetXaxis()->SetTitle("PPM");
  gAbValues->GetYaxis()->SetTitle("absorption factor at 10 cm");
  TCanvas *canAbValues = new TCanvas("canAbValues", "absorption factor vs PPM");
  canAbValues->SetGrid();
  canAbValues->SetLogx();
  gAbValues->Draw("ap");
  canAbValues->Print("abValues.pdf");
  fout->Add(gAbValues);
  fout->Append(canAbValues);
}
