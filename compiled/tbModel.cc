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
TFile *fin;         ///< Input ROOT file handle
TFile *fout;        ///< Output ROOT file handle for results
TNtuple *ntParScan; ///< N-tuple storing parameter scan results

std::vector<TH1D *> hnorm;     ///< Normalized detector response histograms per channel
std::vector<TH1D *> hcurve;    ///< Raw curve histograms per channel
std::vector<TH1D *> hffit;     ///< Fitted waveforms per channel
std::vector<TH1D *> hmodel;    ///< Model histograms per channel (reserved for future use)
std::vector<TH1D *> hfitModel; ///< Final fitted model histograms per channel
std::vector<double> ppmFile;

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
void tbDraw(int theFitChannel = 7, double dopant = 30.)
{
  // draw absorption
  double dist = 10.; // cm
  double ppm[NPOINTS];
  double ab[NPOINTS];
  for (int i = 0; i < NPOINTS; ++i)
  {
    ppm[i] = double(i) * .001;
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

  fout = new TFile(Form("tbModelPPM%.2f.root", theDopant), "recreate");

  // == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
  //  MINUIT OPTIMIZER INITIALIZATION
  // ============================================================================
  // Create Minuit minimizer with NPARS parameters
  TMinuit *gMinuit = new TMinuit(NPARS);
  gMinuit->SetFCN(fcn); // Set likelihood function pointer

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
  printf(" starting value >>>>   fval %E \n", fval);
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
    TH1D *hFit = new TH1D(Form("fitWaveFitChan%i", theFitChannel), Form("fitWaveFitChan%i", theFitChannel), MAXSAMPLE, 0, 2 * MAXSAMPLE);
    hFit->Reset("ICES");
    hFit->SetTitle((Form("fitWaveFitChan%i %.3fPPM", theFitChannel, dopant)));
    hFit->SetMarkerColor(colors[theFitChannel]);
    hFit->SetLineColor(colors[theFitChannel]);
    hfitModel[theFitChannel] = hFit;
    fillFitWave(theFitChannel, hFit);
  }

  // drawing
  // Store individual component contributions (singlet, triplet, mixed) for each channel
  if (theFitChannel < 0)
  {
    for (unsigned ic = 0; ic < NCHANPMT; ++ic)
    {
      for (int icomp = 0; icomp < NUMCOMP; ++icomp)
      {
        TH1D *hFit = new TH1D(Form("fit%sChan%i", compNames[icomp].Data(), ic), Form("fit%sChan%i", compNames[icomp].Data(), ic), MAXSAMPLE, 0, 2 * MAXSAMPLE);
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
      TH1D *hFit = new TH1D(Form("fit%sChan%i", compNames[icomp].Data(), theFitChannel), Form("fit%sChan%i", compNames[icomp].Data(), theFitChannel), MAXSAMPLE, 0, 2 * MAXSAMPLE);
      hFit->Reset("ICES");
      hFit->SetTitle((Form("fit%sChan%i %.3fPPM", compNames[icomp].Data(), theFitChannel, theDopant)));
      hFit->SetLineColor(colors[theFitChannel]);
      fillCompWave(theFitChannel, icomp, hFit);
    }
  }

  fout->Write();

  printf("\n...  finished tbDraw \n");
}
void tbModel()
{
  // Initialize model parameters and optical properties from modelAllFit.hh
  setupModelAllFit();

  // ppm of file

  ppmFile.push_back(1.E-2);
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

  for (int ifile = 0; ifile < ppmFile.size(); ++ifile)
  {
    {
      for (int ilevel = 0; ilevel < NLEVELS; ++ilevel)
      {
        double ab = Absorption(ppmFile[ifile], distanceLevel[ilevel]);
        printf("ppm %f dist %f A %.3E\n", ppmFile[ifile], distanceLevel[ilevel], ab);
      }
    }
  }

  for (unsigned ifile = 0; ifile < ppmFile.size(); ++ifile)
    tbDraw(-2, ppmFile[ifile]);
}
