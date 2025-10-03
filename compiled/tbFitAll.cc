#include <iostream>
#include <fstream>
#include "TMinuit.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TGraph.h"
#include "modelAllFit.hh"
// time is in microseconds
using namespace TMath;
TFile *fin;
TFile *fout;
TDirectory *histoDir;
std::vector<TH1D *> hwave;
std::vector<TH1D *> hffit;
std::vector<TH1D *> hmodel;
std::vector<TH1D *> hffitPmt;
std::vector<TH1D *> hffitChan;
TString histSet;
double dopant[3];
TString summaryFile[3];
double ylow = 1300;
double yhigh = 4000;
bool isSim;
double xTrigger;

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
  }
}

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
  TString fileString;
  histoDir = nullptr;
  if (isSim)
  {
    isSim = true;
    fin->GetObject("sumDir", histoDir);
    fileString = TString("sumPeakWave");
    printf("line127 file %s IS SIMULATION\n", fileName.Data());
  }
  else
  {
    // get subdirectory pointer
    fin->GetObject("runSumDir", histoDir);
    fileString = TString("RunPeakWave");
    printf("line139 file %s IS BTB Data\n", fileName.Data());
  }

  if (!histoDir)
  {
    printf(" no histoDir in file %s\n", fileName.Data());
    return false;
  }

  /* get sum histos from file0 */
  TIter next(histoDir->GetListOfKeys());
  TKey *key;
  int iGot = 0;
  while ((key = (TKey *)next()))
  {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1D"))
      continue;
    TH1D *h = (TH1D *)key->ReadObj();
    TString hname(h->GetName());
    if (hname.Contains(fileString))
    {
      h->SetName(Form("%s", h->GetName()));
      h->SetTitle(Form("%s", h->GetName()));
      // see if we already have this cycle histo
      bool addIt = true;
      for (unsigned ihist = 0; ihist < hwave.size(); ++ihist)
      {
        if (h->GetName() == hwave[ihist]->GetName())
        {
          addIt = false;
          printf("skip %s cycle %i\n", h->GetName(), int(key->GetCycle()));
        }
      }
      if (addIt)
      {
        printf("add %s cycle %i\n", h->GetName(), int(key->GetCycle()));
        TString oldName(h->GetName());
        // rename so sim a data have same name. first get the chan number
        TString tchan = oldName(oldName.Last('e') + 1, 1);
        TString newName(Form("RunPeakWave%s", tchan.Data()));
        h->SetName(newName);
        fout->Add(h);
        hwave.push_back(h);
        ++iGot;
      }
    }
  }
  // fin->Close(); // cannot close because histograms are on input file
  return iGot;
}

void tbFitAll(int fileNum = 2)
{

    summaryFile[0] = TString("summary-05_19_2025-05_19_2025-nfiles-14-created-2025-09-09-15-58.root");
  summaryFile[1] = TString("summary-05_27_2025-05_27_2025-nfiles-22-created-2025-09-09-15-56.root");
  summaryFile[2] = TString("caenData/anaCRun-btbSimOffset-2025-09-09-16-09-1000000-0.root"); // new geometry
  dopant[0] = 0.05;
  dopant[1] = 0.00;
  dopant[2] = 0.00;

  geoVersionOld = true;
  setDistanceLevels(geoVersionOld);

  // is this simulation?
  isSim = false;
  if (summaryFile[fileNum].Contains("btb"))
  {
    isSim = true;
  }
  if (isSim)
    fout = new TFile(Form("tbFitAllSimPPM%.2f.root", dopant[fileNum]), "recreate");
  else
    fout = new TFile(Form("tbFitAllPPM%.2f.root", dopant[fileNum]), "recreate");

  if (geoVersionOld)
    printf("OLD level distances 0 = %.3f 1= %.3f 2= %.3f 3 %.3f 4 %.3f \n", distanceLevel[0], distanceLevel[1], distanceLevel[2], distanceLevel[3], distanceLevel[4]);
  else
    printf("NEW level distances 0 = %.3f 1= %.3f 2= %.3f 3 %.3f 4 %.3f \n", distanceLevel[0], distanceLevel[1], distanceLevel[2], distanceLevel[3], distanceLevel[4]);

  /* channel efficiences */
  for (int ichan = 0; ichan < NCHAN - 1; ++ichan)
  {
    printf("chan %i nominal effGeo  %E   \n", ichan, effGeoFunc(ichan));
  }

  printf(" opened output file %s dopant %f \n", fout->GetName(), dopant[fileNum]);
  // get data histograms from file
  int iGot = openFile(fileNum);
  printf(" got %i data histograms %lu \n", iGot, hwave.size());
  if (iGot < 1)
    return;

  // fout->ls();

  /* fill buffer */
  printf("fill buff \n");
  for (unsigned ichan = 0; ichan < NCHAN; ++ichan)
  {
    printf(".... fill buffer for channel %i  hist %s maximum bin %i \n", ichan, hwave[ichan]->GetName(), hwave[ichan]->GetMaximumBin());
    // fill data buffer
    for (int isample = 1; isample < MAXSAMPLE; ++isample)
    {
      buff[ichan][isample - 1] = hwave[ichan]->GetBinContent(isample); // C starts array from zero
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

  printf("setParNames and setCompNames\n");
  setParNames();
  setCompNames();

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

  if (isSim)
    xTrigger = 698;
  else
    xTrigger = 695;

  // fit starting values
  vstart[NORM] = 2.06322e+03;
  vstart[TRIGSTART] = xTrigger;
  vstart[SFRAC] = 0.20; // btbSim value
  vstart[PPM] = dopant[fileNum];
  vstart[TAU3] = 1600.0;
  vstart[TAUM] = 4700.0;
  vstart[BKGCONST] = 0.0;
  vstart[BKGTAU] = 5000.;
  // fout->ls();

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

  int ibin = 1500;
  // int(xTrigger);
  printModel(ibin, lpar, fsChan, ftChan);

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

  arglist[0] = TRIGSTART + 1; // trigger
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);

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
  arglist[0] = NORM + 1; // par
  arglist[1] = 1.;       // low
  arglist[2] = 1.0E12;   // high
  gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);

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

  arglist[0] = PPM + 1; //
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);

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
  double fvalStart = fval;
  if (isnan(fval))
  {
    printf("gMinuit returns NAN\n");
    return;
  }
  // gMinuit->mnexcm("CALL FCN", vstart, llist, ierflg);

  // fill fit function histogram
  fout->cd(); // add to output file
  // fout->ls();
  for (int ichan = 0; ichan < NCHAN; ++ichan)
  {
    TH1D *hFit = (TH1D *)hwave[ichan]->Clone(Form("fitWaveDefaultChan%i", ichan));
    hFit->Reset("ICES");
    hFit->SetTitle((Form("fitWaveDefaultChan%i", ichan)));
    hFit->SetLineColor(colors[ichan]);
    fillFitWave(ichan, hFit);
  }

  // minimize with MIGRAD
  // Now ready for minimization step
  arglist[0] = 1000000; // maxcalls
  arglist[1] = 1.E-5;   // tolerance

  /* MIGrad[maxcalls][tolerance]*/
  gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  printf("\n...  call to  MIGRAD returns ierflg %i \n", ierflg);

  // Print results
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

  // put parameters into model par array
  for (int k = 0; k < NPARS; ++k)
  {
    gMinuit->GetParameter(k, currentValue, currentError);
    lpar[k] = currentValue;
    printf("\t param %s new value %.4E  \n", lparNames[k].Data(), lpar[k]);
  }
  fcn(llist, gin, fval, lpar, ierflg);
  printf("  >>>> starting fcn %E ending fcn  fval %E \n", fvalStart, fval);

  if (isnan(fval))
  {
    printf("gMinuit returns NAN\n");
    return;
  }
  ibin = 1500;
  // int(xTrigger);
  printModel(ibin, lpar, fsChan, ftChan);
  // printf("... add fit waves to ouput file \n");
  fout->cd(); // add to output file
  // fout->ls();
  for (int ichan = 0; ichan < NCHAN; ++ichan)
  {
    TH1D *hFit = (TH1D *)hwave[ichan]->Clone(Form("fitWaveFitChan%i", ichan));
    hFit->Reset("ICES");
    hFit->SetTitle((Form("fitWaveFitChan%i", ichan)));
    hFit->SetLineColor(colors[ichan]);
    fillFitWave(ichan, hFit);
  }

  TDirectory *compDir = fout->mkdir("components");
  compDir->cd();

  // plot by channel first
  for (int ichan = 0; ichan < NCHAN; ++ichan)
  {
    for (int icomp = 0; icomp < NUMCOMP; ++icomp)
    {
      TH1D *hFit = (TH1D *)hwave[ichan]->Clone(Form("fit%sChan%i", compNames[icomp].Data(), ichan));
      hFit->Reset("ICES");
      hFit->SetTitle((Form("fit%sChan%i", compNames[icomp].Data(), ichan)));
      hFit->SetLineColor(colors[ichan]);
      fillCompWave(ichan, icomp, hFit);
    }
  }

  fout->cd();
  // graph single peaks
  std::vector<double> dataPeak;
  std::vector<double> fitPeak;
  std::vector<double> fchan;
  std::vector<double> dataCorrPeak;
  std::vector<double> fitCorrPeak;
  printf("max bins \n");
  for (int ichan = 0; ichan < 13; ++ichan)
  {
    int peakBin = hwave[ichan]->GetMaximumBin();
    if (fitWave[ichan][peakBin] == 0)
      continue;
    fchan.push_back(ichan);
    double eff = effGeoFunc(ichan);
    dataPeak.push_back(buff[ichan][peakBin]);
    fitPeak.push_back(fitWave[ichan][peakBin]);
    dataCorrPeak.push_back(buff[ichan][peakBin] / eff);
    fitCorrPeak.push_back(fitWave[ichan][peakBin] / eff);
    printf("ichan %i peak bin %i  xTrigger buff %.3E model %.3E effGeo %.3E data ratio %.3E fit ratio %.3E\n",
           ichan, peakBin, buff[ichan][peakBin], fitWave[ichan][peakBin],
           eff, buff[ichan][peakBin] / eff, fitWave[ichan][peakBin] / eff);
  }

  TGraph *grData = new TGraph(fchan.size(), &fchan[0], &dataPeak[0]);
  grData->SetName("dataPeak");
  grData->SetTitle("dataPeak");
  grData->GetHistogram()->GetXaxis()->SetTitle("chan");
  grData->GetHistogram()->GetYaxis()->SetTitle("value");
  grData->SetMarkerStyle(21);
  grData->SetMarkerColor(kBlue);

  TGraph *grFit = new TGraph(fchan.size(), &fchan[0], &fitPeak[0]);
  grFit->SetName("fitPeak");
  grFit->SetTitle("fitPeak");
  grFit->SetMarkerStyle(22);
  grFit->SetMarkerColor(kRed);

  TCanvas *canPeak = new TCanvas("singletPeak", "singlet peak");
  grData->Draw("ap");
  grFit->Draw("psame");
  gPad->SetLogy();
  canPeak->BuildLegend();
  canPeak->SetGrid();

  TGraph *grCorrData = new TGraph(fchan.size(), &fchan[0], &dataCorrPeak[0]);
  grCorrData->SetName("dataCorrPeak");
  grCorrData->SetTitle("dataCorrPeak");
  grCorrData->GetHistogram()->GetXaxis()->SetTitle("chan");
  grCorrData->GetHistogram()->GetYaxis()->SetTitle("efficiency corrected value");
  grCorrData->SetMarkerStyle(21);
  grCorrData->SetMarkerColor(kBlue);

  TGraph *grCorrFit = new TGraph(fchan.size(), &fchan[0], &fitCorrPeak[0]);
  grCorrFit->SetName("fitCorrPeak");
  grCorrFit->SetTitle("fitCorrPeak");
  grCorrFit->SetMarkerStyle(22);
  grFit->SetMarkerColor(kRed);

  TCanvas *canCorrPeak = new TCanvas("singletCorrPeak", "singlet peak eff corrected");
  grCorrData->Draw("ap");
  grCorrFit->Draw("psame");
  gPad->SetLogy();
  canCorrPeak->BuildLegend();
  canCorrPeak->SetGrid();

  fout->Append(grData);
  fout->Append(grFit);

  // grData->Print("all");
  // grFit->Print("all");

  // eff geo corrected

  // delete gMinuit;
  // fout->Purge(1);
  // fout->ls();
  fout->Write();
  fin->Close(); // cannot close before write because histos are on input file
}