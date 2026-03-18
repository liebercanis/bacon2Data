/**
 * @file tbFit.cc
 * @brief Time-based waveform fitting using Minuit optimizer for BACoN detector analysis
 * @details Performs chi-squared minimization fitting to detector curves using the ROOT Minuit library.
 *          Supports fitting individual channels or all channels simultaneously.
 *          Output includes fitted waveforms, component decomposition, and parameter scans.
 * @author M.Gold
 * @date February 17, 2026
 * @note Supersedes tbFitAll.cc (deprecated). Uses modelAllFit.hh for model definitions.
 */
#include <iostream>
#include <fstream>
#include "TGraph.h"
#include "TMinuit.h"
#include "modelAllFit.hh"

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

// ============================================================================
//  ANALYSIS PARAMETERS
// ============================================================================
int nominalTrigger = 729;  ///< Expected trigger timing bin
double singletStart = 700; ///< Singlet scintillation region start [ns]
double singletEnd = 750;   ///< Singlet scintillation region end [ns]

static Double_t vstart[NPARS]; ///< Starting parameter values for Minuit minimization
static Double_t step[NPARS];   ///< Step sizes for Minuit parameter exploration

/// Plotting colors for each detector channel visualization
int colors[NCHAN] = {kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kSpring, kTeal, kAzure, kViolet, kPink, kGray};

// ============================================================================
//  UTILITY FUNCTIONS FOR FITTING AND ANALYSIS
// ============================================================================

/**
 * @brief Calculate reduced chi-squared between data and model for a given channel
 * @param ic Channel index
 * @return Chi-squared per degree of freedom (divided by number of bins)
 */
double doChsiq(int ic)
{
  TH1D *hdata = hnorm[ic];
  TH1D *hmodel = hffit[ic];

  double xval = 0;
  double mean = 0;
  double chisq = 0;
  int nbins = hdata->GetNbinsX();
  for (int ibin = 0; ibin < nbins; ++ibin)
  {
    xval = hdata->GetBinContent(ibin);
    mean = hmodel->GetBinContent(ibin);
    if (mean > 0)
      chisq += pow((xval - mean), 2) / mean; // assuming error is sqrt(mean)
    // printf("ic %i bin %i x %f mean %f chsq %f\n", ic, ibin, xval, mean, chisq);
  }

  return chisq / double(nbins); /// per degree of freedom
}

/**
 * @brief Create comparison canvas displaying data vs fitted model for channel range
 * @param i1 Starting channel index
 * @param i2 Ending channel index
 * @param canName Canvas name and title
 * @return Pointer to created TCanvas with 2x2 panel layout
 */
TCanvas *makeCanFit(int i1, int i2, TString canName)
{
  printf(" makeCanFit %s from %i to %i size %lu \n", canName.Data(), i1, 12, hnorm.size());
  bool firstPlot = true;
  TCanvas *can = new TCanvas(canName, canName);
  can->Divide(2, 2);
  int ipanel = 0;
  for (int i = i1; i <= i2; ++i)
  {
    ++ipanel;
    if (isBadChannel(i))
      continue;
    hnorm[i]->GetXaxis()->SetRangeUser(1000, 15000);
    hfitModel[i]->GetXaxis()->SetRangeUser(1000, 15000);
    hfitModel[i]->SetLineWidth(2);
    can->cd(ipanel);
    gPad->SetLogy();
    hnorm[i]->Draw("");
    hfitModel[i]->Draw("HISTSAME");
  }
  // can->BuildLegend();
  can->SetLogy();
  can->Print(".pdf");
  return can;
}

/**
 * @brief Extract and store late-time background levels for each detector channel
 * @details Uses two methods: (1) integral method over 6000-7500 ns range, and
 *          (2) polynomial fit to 12000-15000 ns tail region. Late background is stored
 *          in the global lateBkg array (defined in modelAllFit.hh).
 */
void fillLateBkg()
{

  for (unsigned ih = 0; ih < hnorm.size(); ++ih)
  {
    // integrate bin ranges
    lateBkg[ih] = hnorm[ih]->Integral(6000, 7500) / double(1500);
    TFitResultPtr fitptr = hnorm[ih]->Fit("pol1", "QS", "", 12000, 15000);
    Int_t fitStatus = fitptr;
    // printf("fitptr %i \n", fitStatus);
    //  hnorm[ih]->GetListOfFunctions()->ls();
    //  gROOT->GetListOfFunctions()->ls();
    TF1 *fit = (TF1 *)hnorm[ih]->GetListOfFunctions()->FindObject("pol1");
    if (fit)
    {
      printf("%u %s late int inte %E per bin fit %E slope %E \n", ih, hnorm[ih]->GetName(), lateBkg[ih], fit->GetParameter(0), fit->GetParameter(1));
      // or just use average
      // lateBkg[ih] = fit->GetParameter(0);
    }
  }
  for (unsigned ih = 0; ih < hnorm.size(); ++ih)
  {
    // printf("%u %s inte %E \n", ih, hnorm[ih]->GetName(), hnorm[ih]->Integral());
    hnorm[ih]->GetListOfFunctions()->Clear();
  }
  for (unsigned ih = 0; ih < hnorm.size(); ++ih)
  {
    // printf("%u %s inte %E \n", ih, hnorm[ih]->GetName(), hnorm[ih]->Integral());
    hnorm[ih]->GetListOfFunctions()->Clear();
  }
}

/**
 * @brief Perform a likelihood scan over a specified parameter range
 * @param thePar Parameter index to scan (e.g., NORM, TAU3, PPM)
 * @param xlow Lower bound for scan range
 * @param xhigh Upper bound for scan range
 * @return TGraph containing parameter values vs negative log-likelihood
 */
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
    yval.push_back(nLL);
    // printf("line53 i %i par nph %f r %f theta %f  phi %f \n", i, fitVal[0], fitVal[1], fitVal[2], fitVal[3]);
    ntScan->Fill(nLL, fitVal[1], fitVal[2], fitVal[3]);
    // printf("mySCAN par %i x= %f nLL %E \n", i, x, nLL);
  }
  // make and return graph
  return new TGraph(maxPoints, &xval[0], &yval[0]);
}

/**
 * @brief Execute parameter scan and save results to output file
 * @param thePar Parameter index to scan
 * @return TGraph with scan results, configured with titles and saved to output
 */
TGraph *parameterScan(int thePar)
{
  TGraph *graph = myScan(thePar, 0.001 * lpar[thePar], 2. * lpar[thePar]);
  graph->SetName(Form("ScanPar%i", thePar));
  graph->SetTitle(Form("ScanPar%i", thePar));
  graph->GetYaxis()->SetTitle("FCN likelihood value");
  graph->GetXaxis()->SetTitle(Form("parameter %s", lparNames[thePar].Data()));
  fout->Add(graph);
  return graph;
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
  // std::cout << " fillFitWave " << ichan << "  " << hist->GetName() << std::endl;
  hist->Reset("ICES");
  for (int ib = 1; ib < hist->GetNbinsX(); ++ib)
  {
    double val = max(fitComp[ichan][icomp][ib], 1.E-9);
    // if (ichan == 8 && ib == 1500)
    //   printf("!!!! chan %i sample %i val %E \n", ichan, ib, val);
    hist->SetBinContent(ib, val);
    hist->SetBinError(ib, 0);
    hist->GetYaxis()->SetTitle("yield");
    hist->GetXaxis()->SetTitle("time [ns]");
  }
  std::cout << std::endl;
}

/**
 * @brief Load all detector waveform histograms from input ROOT file
 * @details Scans input file for histograms containing "CurveChan" (raw curves) or
 *          "NormChan" (normalized curves). Both types are stored in global vectors.
 */
void getCurves()
{

  TIter next(fin->GetListOfKeys());
  TKey *key;
  int ifile = 0;
  while (TKey *key = (TKey *)next())
  {
    TClass *cl = gROOT->GetClass(key->GetClassName());

    if (!cl->InheritsFrom("TH1D"))
      continue;
    TH1D *h = (TH1D *)key->ReadObj();

    if (TString(h->GetName()).Contains("CurveChan"))
      hcurve.push_back(h);

    if (TString(h->GetName()).Contains("NormChan"))
    {
      hnorm.push_back(h);
      fout->Append(h);
    }
  }
}

/**
 * @brief Perform constant background level fit to high-time tail region of histogram
 * @details Fits histogram to a constant (pol0) function in the range 12000-15000 ns.
 *          Used to estimate background level for subtraction. Returns fitted constant.
 * @param hist Pointer to histogram to fit
 * @return Fitted background level (constant parameter value)
 */
double fitBack(TH1D *hist)
{
  double lowCut = 12000;
  double highCut = 15000;
  TF1 *gfit = NULL;
  double ave = 0;
  int lowBins = hist->FindBin(lowCut);
  int highBins = hist->FindBin(highCut);
  auto fitBack = new TF1("fitBack", "pol0", lowCut, highCut);
  fitBack->SetParameter(0, 1E-1);
  fitBack->SetParLimits(0, 1.E-9, 1E2);

  hist->Fit("fitBack", "LF", " ", lowCut, highCut);
  gfit = (TF1 *)hist->GetListOfFunctions()->FindObject("fitBack");
  if (gfit)
  {
    ave = gfit->GetParameter(0);
  }
  else
    printf("P1 Fit to hist fails \n");
  double aveb = hist->Integral(lowCut, highCut) / double(highBins - lowBins);
  printf(" \t\t ave %E aveb %E \n\n", ave, aveb);
  return ave;
}

/**
 * @brief Open and validate input ROOT data file
 * @param fileName Path to input ROOT file for analysis
 * @return True if file successfully opened and assigned to global fin; false otherwise
 */
bool openFile(TString fileName)
{
  // open input file and make some histograms
  printf(" looking for file %s\n", fileName.Data());

  bool exists = false;
  FILE *aFile;
  aFile = fopen(fileName.Data(), "r");
  if (aFile)
  {
    fclose(aFile);
    exists = true;
  }
  if (!exists)
  {
    printf(" couldnt open file %s\n", fileName.Data());
    return false;
  }

  fin = new TFile(fileName, "readonly");
  printf(" opened file %s\n", fileName.Data());
  return true;
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
void tbFit(int theFitChannel = -1)
{

  // Initialize histogram vectors for all channels
  hffit.resize(NCHAN);
  hmodel.resize(NCHAN);

  // ============================================================================
  //  INPUT FILE SELECTION AND VALIDATION
  // ============================================================================
  TString inputFile = TString("post-anaCRun-btbSimNEW-2026-02-13-100000-7857.root");
  inputFile = TString("post-11_19_2025-11_19_2025-1371746.root"); // Override with data file
  // inputFile = TString("anaCRun-btbSimNEW-2026-02-23-10-18-100000-0.root");

  if (!openFile(inputFile))
    return;

  // ============================================================================
  //  ANALYSIS CONFIGURATION
  // ============================================================================
  double dopant = 1.E-2; ///< Dopant concentration [PPM] for simulations

  // Initialize model parameters and optical properties from modelAllFit.hh
  setupModelAllFit();

  // Determine whether input is simulation or experimental data
  bool isSim = false;
  if (inputFile.Contains("btb")) // "btb" in filename indicates batch simulation
  {
    isSim = true;
  }

  // Create output file with naming convention reflecting data type and dopant level
  if (isSim)
    fout = new TFile(Form("tbFitSimPPM%.2f.root", dopant), "recreate");
  else
    fout = new TFile(Form("tbFitPPM%.2f.root", dopant), "recreate");

  if (geoVersionOld)
    printf("OLD level distances 0 = %.3f 1= %.3f 2= %.3f 3 %.3f 4 %.3f \n", distanceLevel[0], distanceLevel[1], distanceLevel[2], distanceLevel[3], distanceLevel[4]);
  else
    printf("NEW level distances 0 = %.3f 1= %.3f 2= %.3f 3 %.3f 4 %.3f \n", distanceLevel[0], distanceLevel[1], distanceLevel[2], distanceLevel[3], distanceLevel[4]);

  // ============================================================================
  //  DATA LOADING AND PREPROCESSING
  // ============================================================================
  // Load all detector waveforms from input file hcurve and hnorm === hcurve is used for fit
  getCurves();

  // Extract and characterize late-time background for each channel
  // either from average or from fit to poly1
  fillLateBkg();

  // Transfer histogram data to buffer arrays for fitting algorithm
  printf("fill buff \n");
  for (unsigned ichan = 0; ichan < NCHANPMT; ++ichan)
  {
    printf(".... fill buffer for channel %i  hist %s maximum bin %i \n", ichan, hnorm[ichan]->GetName(), hnorm[ichan]->GetMaximumBin());
    // fill data buffer
    for (int isample = 1; isample < MAXSAMPLE; ++isample)
    {
      buff[ichan][isample - 1] = hnorm[ichan]->GetBinContent(isample); // C starts array from zero
    }
  }

  // ============================================================================
  //  MINUIT OPTIMIZER INITIALIZATION
  // ============================================================================
  // Create Minuit minimizer with NPARS parameters
  TMinuit *gMinuit = new TMinuit(NPARS);
  gMinuit->SetFCN(fcn); // Set likelihood function pointer

  double currentValue;
  double currentError;

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
  vstart[NORM] = startNorm;                                              ///< Photon yield per event
  vstart[TRIGSTART] = 2. * hnorm[9]->GetMaximumBin() - 10 * tResolution; ///< Trigger timing offset
  vstart[SFRAC] = 0.14;                                                  ///< Singlet fraction (ref: Segretto 2021)
  vstart[PPM] = dopant;                                                  ///< Dopant concentration
  vstart[TAU3] = tTriplet0;                                              ///< Triplet decay time
  vstart[TAUM] = 4700.0;                                                 ///< Mixed component decay time
  vstart[BKGCONST] = 4.0E-6;                                             ///< Constant background rate
  vstart[BKGTAU] = 5000.;                                                ///< Background decay timescale
  vstart[THECHANNEL] = theFitChannel;                                    ///< Channel selection flag
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

  // ============================================================================
  //  PARAMETER CONSTRAINTS AND BOUNDARIES
  // ============================================================================
  // Fix parameters that are held constant during minimization
  // Note: Minuit uses 1-based indexing for parameters (adds 1 to C++ indices)

  arglist[0] = TRIGSTART + 1; // trigger
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  arglist[0] = BKGCONST + 1; // par
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  arglist[0] = BKGTAU + 1; // par
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  // arglist[0] = SFRAC + 1; // kp
  // gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  // Set bounds for variable parameters to restrict optimization domain
  arglist[0] = NORM + 1;         // par
  arglist[1] = 0.01 * startNorm; // low
  arglist[2] = 10. * startNorm;  // high
  gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);
  // gMinuit->mnexcm("FIX", arglist, 3, ierflg);

  arglist[0] = SFRAC + 1; // par
  arglist[1] = 0.01;      // low
  arglist[2] = 1.0;
  gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);
  // gMinuit->mnexcm("FIX", arglist, 3, ierflg);

  // set limits ... here par starts with 1 so add 1
  arglist[0] = TAU3 + 1;         // par
  arglist[1] = 0.01 * tTriplet0; // low
  arglist[2] = 2.0 * tTriplet0;  // high
  gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);

  // set limits ... here par starts with 1 so add 1
  arglist[0] = PPM + 1; // par
  arglist[1] = 0.0;     // low
  arglist[2] = 100.;    // high
  // gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  arglist[0] = TAUM + 1;     // par tau mixed
  arglist[1] = 0.01 * tMix0; // low
  arglist[2] = 10. * tMix0;  // high
  // gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);
  gMinuit->mnexcm("FIX", arglist, 3, ierflg);

  // ============================================================================
  //  VERIFICATION AND MINIMIZATION
  // ============================================================================
  double amin;
  printf("call mnprin starting values \n");
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

  // ============================================================================
  //  LIKELIHOOD MINIMIZATION (MIGRAD ALGORITHM)
  // ============================================================================
  // Configure MIGRAD minimizer parameters
  arglist[0] = 1000000; ///< Maximum function calls permitted
  arglist[1] = 1.E-5;   ///< Convergence tolerance on FCN value

  // Execute MIGRAD minimization
  gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  printf("\n...  after fit, call to  MIGRAD returns ierflg %i \n", ierflg);

  // ============================================================================
  //  RESULT EXTRACTION AND REPORTING
  // ============================================================================
  // Retrieve final minimization statistics
  Double_t edm, errdef;
  Int_t nvpar, nparx, icstat;
  printf("...  call mnstat \n");
  gMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
  /*
  Prints the values of the parameters at the time of the call.
  According to the value of IKODE, the printout is: IKODE=INKODE=

  0 only info about function value
  1 parameter values, errors, limits
  2 values, errors, step sizes, internal values
  3 values, errors, step sizes, first derivs.
  4 values, parabolic errors, MINOS errors
  when INKODE=5, MNPRIN chooses IKODE=1,2, or 3, according to fISW[1]
  */
  printf("\n...  call mnprin \n");
  gMinuit->mnprin(1, amin);

  // ============================================================================
  //  FITTED WAVEFORM STORAGE AND OUTPUT
  // ============================================================================
  // Allocate histograms for fitted waveforms across all PMT channels
  hfitModel.resize(NCHANPMT);

  for (unsigned ic = 0; ic < NCHANPMT; ++ic)
  {
    TH1D *hFit = (TH1D *)hnorm[ic]->Clone(Form("fitWaveFitChan%i", ic));
    hFit->Reset("ICES");
    hFit->SetTitle((Form("fitWaveFitChan%i", ic)));
    hFit->SetLineColor(colors[ic]);
    hfitModel[ic] = hFit;
    fillFitWave(ic, hFit);
  }

  // Create subdirectory in output file for component decomposition histograms
  TDirectory *compDir = fout->mkdir("components");
  compDir->cd();

  // Store individual component contributions (singlet, triplet, mixed) for each channel
  if (theFitChannel == -1)
  {
    for (unsigned ic = 0; ic < NCHANPMT; ++ic)
    {
      for (int icomp = 0; icomp < NUMCOMP; ++icomp)
      {
        TH1D *hFit = (TH1D *)hnorm[ic]->Clone(Form("fit%sChan%i", compNames[icomp].Data(), ic));
        hFit->Reset("ICES");
        hFit->SetTitle((Form("fit%sChan%i", compNames[icomp].Data(), ic)));
        hFit->SetLineColor(colors[ic]);
        fillCompWave(ic, icomp, hFit); // only need one of these
      }
    }
  }
  else
  {
    for (int icomp = 0; icomp < NUMCOMP; ++icomp)
    {
      TH1D *hFit = (TH1D *)hnorm[theFitChannel]->Clone(Form("fit%sChan%i", compNames[icomp].Data(), theFitChannel));
      hFit->Reset("ICES");
      hFit->SetTitle((Form("fit%sChan%i", compNames[icomp].Data(), theFitChannel)));
      hFit->SetLineColor(colors[theFitChannel]);
      fillCompWave(theFitChannel, icomp, hFit);
    }
  }

  // ============================================================================
  //  PARAMETER SCAN AND VISUALIZATION
  // ============================================================================
  // Create n-tuple to store parameter scan results
  ntParScan = new TNtuple("ntParScan", "parameter scan", "nll:fitVal1:fitVal2:fitVal3");

  // Scan TAU3 (triplet decay time) parameter across physically motivated range
  int thePar = TAU3;
  printf("scan parameter %i %s from %f to %f \n", thePar, lparNames[thePar].Data(), 0.001 * lpar[thePar], 2. * lpar[thePar]);
  parameterScan(thePar);

  thePar = NORM;
  printf("scan parameter %i %s from %f to %f \n", thePar, lparNames[thePar].Data(), 0.001 * lpar[thePar], 2. * lpar[thePar]);
  parameterScan(thePar);

  if (theFitChannel == -1)
  {
    makeCanFit(0, 2, TString("canFitLevel0"));
    makeCanFit(3, 5, TString("canFitLevel1"));
    makeCanFit(6, 8, TString("canFitLevel2"));
    makeCanFit(9, 11, TString("canFitTrig"));
  }

  // ============================================================================
  //  GOODNESS-OF-FIT EVALUATION
  // ============================================================================
  // Compute and report reduced chi-squared for each non-flagged channel
  for (unsigned ic = 0; ic < NCHANPMT; ++ic)
  {
    if (!isBadChannel(ic))
      printf("channel %i chisq/dof %.3E \n", ic, doChsiq(ic));
  }

  printf(" shift %.2f peak curve chan 9 %d peak fit %d ns \n", shift, 2 * hcurve[9]->GetMaximumBin(), 2 * hffit[9]->GetMaximumBin());

  printf("\n...  finished tbFit \n");
}
