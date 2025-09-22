#include <iostream>
#include <fstream>
#include "TMinuit.h"
#include "TFile.h"
#include "TString.h"
#include "TH1D.h"
#include "TGraph.h"
#include "modelAllFit.hh"
// time is in microseconds
using namespace TMath;
TFile *fin;
TFile *fout;
TDirectory *runSumDir;
std::vector<TH1D *> hwave;
std::vector<TH1D *> hffit;
std::vector<TH1D *> hmodel;
std::vector<TH1D *> hffitPmt;
std::vector<TH1D *> hffitChan;
TString histSet;
double dopant[2];
TString summaryFile[2];
double ylow = 1300;
double yhigh = 4000;

/*
kWhite  = 0,   kBlack  = 1,   kGray    = 920,  kRed    = 632,  kGreen  = 416,
kBlue   = 600, kYellow = 400, kMagenta = 616,  kCyan   = 432,  kOrange = 800,
kSpring = 820, kTeal   = 840, kAzure   =  860, kViolet = 880,  kPink   = 900
*/
int colors[NCHAN] = {kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kSpring, kTeal, kAzure, kViolet, kPink, kGray};

void fillHistogram(TF1 *func, TH1D *hist)
{
  if (!func)
  {
    printf("!!! fillHistogram called with null func !!!!! \n");
    return;
  }
  if (hist == NULL)
    return;
  for (int ib = 1; ib < hist->GetNbinsX(); ++ib)
  {
    double xbin = hist->GetBinCenter(ib);
    double fbin = func->Eval(xbin);
    hist->SetBinContent(ib, fbin);
    hist->SetBinError(ib, 0);
    hist->GetYaxis()->SetTitle("yield");
    hist->GetXaxis()->SetTitle("time [ns]");
  }
}

void fillFitWave(int ichan, TH1D *hist)
{
  // std::cout << " fillFitWave " << ichan << "  " << hist->GetName() << std::endl;
  hist->Reset("ICES");
  for (int ib = 1; ib < hist->GetNbinsX(); ++ib)
  {
    double val = max(fitWave[ichan][ib], 1.E-9);
    hist->SetBinContent(ib, val);
    hist->SetBinError(ib, 0);
    hist->GetYaxis()->SetTitle("yield");
    hist->GetXaxis()->SetTitle("time [ns]");
    // if (ichan == 9 && (ib > 1000 && ib < 1100))
    //   printf("chan %i sample %i buff  %E  fitWave %E ... ", ichan, ib, buff[ichan][ib], hist->GetBinContent(ib));
  }
  std::cout << std::endl;
}

double fitBack(TH1D *hist)
{
  double lowCut = 8000;
  double highCut = 15000;
  TF1 *gfit = NULL;
  double ave = 0;
  int lowBins = hist->FindBin(lowCut);
  int highBins = hist->FindBin(highCut);
  auto fitBack = new TF1("fitBack", "pol0", lowCut, highCut);
  fitBack->SetParameter(0, 1E-5);
  fitBack->SetParLimits(0, 1.E-12, 1.E2);

  hist->Fit("fitBack", "LF", " ", lowCut, highCut);
  gfit = (TF1 *)hist->GetListOfFunctions()->FindObject("fitBack");
  if (gfit)
  {
    ave = gfit->GetParameter(0);
  }
  else
    printf("P1 Fit to hist fails \n");
  double aveb = hist->Integral(lowCut, highCut) / double(highBins - lowBins);
  printf(" BBBB background fit ave %E aveb %E \n\n", ave, aveb);
  return ave;
}

int openFile(int fileNum = 0)
{
  TString fileName = summaryFile[fileNum];
  // open input filex
  hffit.resize(NCHAN);
  hmodel.resize(NCHAN);
  hffitPmt.resize(NCHAN);

  printf(" looking for summary file %s\n", fileName.Data());

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
  // file exists, so open with TFile
  fin = new TFile(fileName, "readonly");
  printf(" opened file %s\n", fileName.Data());
  runSumDir = NULL;
  // get subdirectory pointer
  fin->GetObject("runSumDir", runSumDir);
  if (!runSumDir)
  {
    printf(" no runSumDir in file %s\n", fileName.Data());
    return false;
  }

  /* get sum histos from file0 */
  TIter next(runSumDir->GetListOfKeys());
  TKey *key;
  int iGot = 0;
  while (TKey *key = (TKey *)next())
  {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1D"))
      continue;
    TH1D *h = (TH1D *)key->ReadObj();
    TString hname(h->GetName());
    if (hname.Contains("RunPeakWave"))
    {
      h->SetName(Form("%s", h->GetName()));
      h->SetTitle(Form("%s", h->GetName()));
      fout->Add(h);
      hwave.push_back(h);
      ++iGot;
    }
  }
  return iGot;
}

void tbFitAll(int fileNum = 0)
{
  summaryFile[0] = TString("summary-05_19_2025-05_19_2025-nfiles-14-created-2025-09-09-15-58.root");
  summaryFile[1] = TString("summary-05_27_2025-05_27_2025-nfiles-22-created-2025-09-09-15-56.root");
  dopant[0] = 0.05;
  dopant[1] = 0.00;

  fout = new TFile(Form("tbFitAllPPM%.2f.root", dopant[fileNum]), "recreate");
  printf(" opened output file %s dopant %f \n", fout->GetName(), dopant[fileNum]);
  // get data histograms from file
  int iGot = openFile(fileNum);
  printf(" got %i data histograms %lu \n", iGot, hwave.size());
  if (iGot < 1)
    return;

  printf("max bins \n");
  for (int ichan = 0; ichan < hwave.size(); ++ichan)
    printf("%i %s max bin %i \n", ichan, hwave[ichan]->GetName(), hwave[ichan]->GetMaximumBin());

  /* fill buffer */
  printf("fill buff \n");
  for (unsigned ichan = 0; ichan < NCHAN; ++ichan)
  {
    fout->Add(hwave[ichan]);
    std::cout << ".... for channel %i " << ichan << "  " << hwave[ichan]->GetName() << " maximum bin " << hwave[ichan]->GetMaximumBin() << std::endl;
    // fill data buffer
    for (int isample = 1; isample < MAXSAMPLE; ++isample)
    {
      buff[ichan][isample] = hwave[ichan]->GetBinContent(isample);
      // if (buff[ichan][isample] == 0)
      //   printf("sample %i val %E  ... ", isample, buff[ichan][isample]);
    }
  }

  // clones to store fit histos
  /*
  for (unsigned ih = 0; ih < hwave.size(); ++ih)
  {
    if (hwave[ih] == NULL)
      continue;
    hwave[ih]->GetListOfFunctions()->Clear();
    hffit[ih] = (TH1D *)hwave[ih]->Clone(Form("hFitCh%i", ih));
    fout->Add(hffit[ih]);
    hmodel[ih] = (TH1D *)hwave[ih]->Clone(Form("hModelCh%i", ih));
    fout->Add(hmodel[ih]);
  }
    */

  // clones to store fit components histos
  /*
  for (unsigned ih = 0; ih < NCOMP; ++ih)
  {
    TString histString;
    histString.Form("ModelPmt%s", compName[ih].Data());
    hffitPmt[ih] = (TH1D *)hwave[0]->Clone(histString);
    hffitPmt[ih]->SetTitle(histString);
    fout->Add(hffitPmt[ih]);
    hffitChan[ih] = (TH1D *)hwave[0]->Clone(Form("ModelCh7%s", compName[ih].Data()));
    hffitChan[ih]->SetTitle(Form("ModelCh7%s", compName[ih].Data()));
    fout->Add(hffitChan[ih]);
  }
    */

  printf("setParNames\n");
  setParNames();

  // Set starting values and step sizes for parameters
  static Double_t vstart[NPARS];
  static Double_t step[NPARS];
  /*
    lparNames[NORM] = TString("norm");
    lparNames[SFRAC] = TString("sfrac");
    lparNames[PPM] = TString("ppm");
    lparNames[TAU3] = TString("tau3");
    lparNames[TAUM] = TString("taumix");
    lparNames[BKGCONST] = TString("bkgconst");
    lparNames[BKGTAU] = TString("bkgtau");
  */

  // fit starting values
  vstart[NORM] = 2.06322e+03;
  vstart[SFRAC] = 0.25; //
  vstart[PPM] = dopant[fileNum];
  vstart[TAU3] = 1600.0;
  vstart[TAUM] = 4700.0;
  vstart[BKGCONST] = 0.0;
  vstart[BKGTAU] = 5000.;

  printf("starting parameter values \n");
  for (int ip = 0; ip < NPARS; ++ip)
    printf(" par %i %s start val %f \n", ip, lparNames[ip].Data(), vstart[ip]);

  TMinuit *gMinuit = new TMinuit(NPARS); // initialize TMinuit with a maximum of 5 params
  gMinuit->SetFCN(fcn);

  Double_t arglist[10];

  int ierflg = 0;
  arglist[0] = 0.5; // for likelihood
  gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

  for (unsigned j = 0; j < NPARS; ++j)
  {
    step[j] = 1.E-6 * vstart[j];
    gMinuit->mnparm(j, lparNames[j].Data(), vstart[j], step[j], 0, 0, ierflg);
    lpar[j] = vstart[j];
  }

  double ftChan[NCHAN];
  double fsChan[NCHAN];
  double xChan[NCHAN];
  for (int ichan = 0; ichan < NCHAN; ++ichan)
    xChan[ichan] = double(ichan);

  int ibin = 705;
  printModel(iTrigger, lpar, fsChan, ftChan);

  TGraph *fsGraph = new TGraph(NCHAN, xChan, fsChan);
  fsGraph->SetName(Form("fsAt%i", ibin));
  fsGraph->SetTitle(Form("fast component at sample %i", ibin));
  fsGraph->SetMarkerStyle(21);
  fsGraph->SetMarkerColor(kBlue);
  fout->Append(fsGraph);

  double currentValue;
  double currentError;

  /******************************/
  // fix parameters minuit is fortran!
  /******************************/
  arglist[0] = TAUM + 1; // par tau3
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  arglist[0] = SFRAC + 1; // kp
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  arglist[0] = BKGCONST + 1; //
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  arglist[0] = BKGTAU + 1; //
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  /*
      set limits ... here par starts with 1 so add 1
   */
  // set limits ... here par starts with 1 so add 1
  arglist[0] = TAU3 + 1; // par
  arglist[1] = 100.;     // low
  arglist[2] = 10000.;   // high
  gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);

  // set limits ... here par starts with 1 so add 1
  arglist[0] = PPM + 1; // par
  arglist[1] = 0.0;     // low
  arglist[2] = 100.;    // high
  gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);

  printf(" \n\n >>> modelFit start parameters fit ppm %f \n", vstart[1]);
  for (int ii = 0; ii < NPARS; ++ii)
  {
    gMinuit->GetParameter(ii, currentValue, currentError);
    printf("\t  param %i %s %.4E  \n", ii, lparNames[ii].Data(), currentValue);
  }

  double amin;
  gMinuit->mnprin(1, amin);

  double fval = 0;
  double gin[NPARS];
  int npar = NPARS;

  // Call the function once and get the return value.
  int llist = NPARS; // Number of parameters
  fcn(llist, gin, fval, lpar, ierflg);
  printf(" starting value >>>>   fval %E \n", fval);
  if (isnan(fval))
  {
    printf("gMinuit returns NAN\n");
    return;
  }
  // gMinuit->mnexcm("CALL FCN", vstart, llist, ierflg);

  // fill fit function histogram
  fout->cd(); // add to output file
  for (int ichan = 0; ichan < NCHAN; ++ichan)
  {
    TH1D *hFit = (TH1D *)hwave[ichan]->Clone(Form("fitWaveDefaultChan%i", ichan));
    hFit->SetLineColor(colors[ichan]);
    fillFitWave(ichan, hFit);
  }

  // fout->ls();
  fout->Write();

  // minimize with MIGRAD
  // Now ready for minimization step
  arglist[0] = 100000; // maxcalls
  arglist[1] = 1.E-2;  // tolerance

  /* MIGrad[maxcalls][tolerance]*/
  gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  // Print results
  Double_t edm, errdef;
  Int_t nvpar, nparx, icstat;
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
  gMinuit->mnprin(1, amin);

  // put parameters into model par array
  for (int k = 0; k < NPARS; ++k)
  {
    gMinuit->GetParameter(k, currentValue, currentError);
    lpar[k] = currentValue;
    printf("\t copy %s new value %f \n", lparNames[k].Data(), lpar[k]);
  }

  for (int ichan = 0; ichan < NCHAN; ++ichan)
  {
    TH1D *hFit = (TH1D *)hwave[ichan]->Clone(Form("fitWaveFitChan%i", ichan));
    hFit->SetLineColor(colors[ichan]);
    fillFitWave(ichan, hFit);
  }

  fout->Write();
}