/*
root macro to fit data with model
using modelAllFit.hh
Nov 13 2025
******* deprecated as of Feb 2026
use tvbFit.cc instead
***********
*/
#include <iostream>
#include <string>
#include <fstream>
#include "TMinuit.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "modelAllFit.hh"
// time is in microseconds
using namespace TMath;
TFile *fin;
TFile *fout;
std::vector<TH1D *> hnorm;
std::vector<TH1D *> hcurve;
std::vector<TH1D *> hffit;
std::vector<TH1D *> hmodel;

TString histSet;
double dopant[4];
TString summaryFile[4];
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
      hnorm.push_back(h);
  }
}

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
    hffit[ichan] = hist;
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

bool openFile(int fileNum = 0)
{
  TString fileName = summaryFile[fileNum];
  // open input filex
  hffit.resize(NCHAN);
  hmodel.resize(NCHAN);

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
  return exists;
}

void tbFitAll(int fileNum = 0)
{

  summaryFile[0] = TString("post-anaCRun-btbSimNEW-2026-02-13-100000-7857.root");
  dopant[0] = 0.00;

  geoVersionOld = false;
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
  if (!openFile(fileNum))
    return;

  getCurves();
  printf(" got data histograms %lu \n", hcurve.size());
  double startNorm = 0;
  for (unsigned ich = 0; ich < hcurve.size(); ++ich)
  {
    printf(" chan %i integral s %f \n", ich, hcurve[ich]->Integral());
    if (ich < NCHANPMT)
      startNorm += hcurve[ich]->Integral();
  }
  printf("starting norm %E \n", startNorm);
  // fout->ls();

  /* fill buffer */
  printf("fill buff \n");
  for (unsigned ichan = 0; ichan < NCHANPMT; ++ichan)
  {
    printf(".... fill buffer for channel %i  hist %s maximum bin %i \n", ichan, hcurve[ichan]->GetName(), hcurve[ichan]->GetMaximumBin());
    // fill data buffer
    for (int isample = 1; isample < MAXSAMPLE; ++isample)
    {
      buff[ichan][isample - 1] = hcurve[ichan]->GetBinContent(isample); // C starts array from zero
      // if (buff[ichan][isample] == 0)
      //   printf("sample %i val %E  ... ", isample, buff[ichan][isample]);
    }
  }

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
    xTrigger = 705;
  else
    xTrigger = 705;

  // fit starting values
  vstart[NORM] = startNorm;
  vstart[TRIGSTART] = xTrigger;
  vstart[SFRAC] = 0.14; //;0.23;   // Segretto PHYSICAL REVIEW D 103, 043001 (2021)
  vstart[PPM] = dopant[fileNum];
  vstart[TAU3] = 1600.0;
  vstart[TAUM] = 4700.0;
  vstart[BKGCONST] = 0.0;
  vstart[BKGTAU] = 5000.;
  // fout->ls();

  printf("starting parameter values \n");
  for (int ip = 0; ip < NPARS; ++ip)
    printf(" par %i %s start val %f \n", ip, lparNames[ip].Data(), vstart[ip]);

  // make TMinuit class instance and set minimization function
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

  /** look at model prior to fitting  */
  // plot by channel first
  for (int ichan = 0; ichan < hcurve.size(); ++ichan)
  {
    if (ichan != 8)
      continue;
    for (int icomp = 0; icomp < NUMCOMP; ++icomp)
    {
      TH1D *hFit = (TH1D *)hcurve[ichan]->Clone(Form("fit%sChan%i", compNames[icomp].Data(), ichan));
      hFit->Reset("ICES");
      hFit->SetTitle((Form("fit%sChan%i", compNames[icomp].Data(), ichan)));
      hFit->SetLineColor(colors[ichan]);
      fillCompWave(ichan, icomp, hFit);
    }
  }

  // minimize with MIGRADfill
  // Now ready for minimization step
  arglist[0] = 1000000; // maxcalls
  arglist[1] = 1.E-5;   // tolerance

  /* MIGrad[maxcalls][tolerance]*/
  gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  printf("\n...  after fit, call to  MIGRAD returns ierflg %i \n", ierflg);

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
    TH1D *hFit = (TH1D *)hcurve[ichan]->Clone(Form("fitWaveFitChan%i", ichan));
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
      TH1D *hFit = (TH1D *)hcurve[ichan]->Clone(Form("fit%sChan%i", compNames[icomp].Data(), ichan));
      hFit->Reset("ICES");
      hFit->SetTitle((Form("fit%sChan%i", compNames[icomp].Data(), ichan)));
      hFit->SetLineColor(colors[ichan]);
      fillCompWave(ichan, icomp, hFit);
    }
  }

  fout->cd();
  return;

  // graph single peaks
  std::vector<double> dataPeak;
  std::vector<double> fitPeak;
  std::vector<double> fchan;
  std::vector<double> dataCorrPeak;
  std::vector<double> fitCorrPeak;
  std::vector<double> dataPeakBin;
  printf("peak fit results: \n");
  for (int ichan = 0; ichan < 13; ++ichan)
  {
    int peakBin = hcurve[ichan]->GetMaximumBin();
    dataPeakBin.push_back(hcurve[ichan]->GetBinContent(peakBin));
    if (fitWave[ichan][peakBin] == 0)
      continue;
    fchan.push_back(ichan);
    double eff = effGeoFunc(ichan);
    dataPeak.push_back(buff[ichan][peakBin]);
    fitPeak.push_back(fitWave[ichan][peakBin]);
    dataCorrPeak.push_back(buff[ichan][peakBin] / eff);
    fitCorrPeak.push_back(fitWave[ichan][peakBin] / eff);
    printf("ichan %i peak bin %i  xTrigger buff %.3E model %.3E effGeo %.3E data eff corr %.3E fit eff corr %.3E\n",
           ichan, peakBin, buff[ichan][peakBin], fitWave[ichan][peakBin],
           eff, buff[ichan][peakBin] / eff, fitWave[ichan][peakBin] / eff);
  }

  // sort peak bins
  std::sort(dataPeakBin.begin(), dataPeakBin.end());
  printf(" sorted peak bins %f %f \n", dataPeakBin[0], dataPeakBin[dataPeakBin.size() - 1]);

  // make multigraph
  printf("number of peaks is %lu \n", fchan.size());

  TGraph *grData = new TGraph(fchan.size(), &fchan[0], &dataPeak[0]);
  grData->SetName("dataPeak");
  grData->SetTitle("dataPeak");
  grData->SetMarkerStyle(21);
  grData->SetMarkerColor(kBlue);

  TGraph *grFit = new TGraph(fchan.size(), &fchan[0], &fitPeak[0]);
  grFit->SetName("fitPeak");
  grFit->SetTitle("fitPeak");
  grFit->SetMarkerStyle(22);
  grFit->SetMarkerColor(kRed);

  TMultiGraph *gmultPeak = new TMultiGraph();
  gmultPeak->Add(grData);
  gmultPeak->Add(grFit);
  TCanvas *canPeak = new TCanvas("singletPeak", "singlet peak");
  gmultPeak->GetXaxis()->SetTitle("channel");
  gmultPeak->GetYaxis()->SetTitle("singlet peak value");
  gmultPeak->Draw("apm");
  // gPad->SetLogy();
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
  grCorrFit->SetMarkerColor(kRed);

  TMultiGraph *gmultCorr = new TMultiGraph();
  gmultCorr->Add(grCorrData);
  gmultCorr->Add(grCorrFit);
  TCanvas *canCorr = new TCanvas("singletCorrPeak", "singlet corrected peak");
  gmultCorr->GetXaxis()->SetTitle("channel");
  gmultCorr->GetYaxis()->SetTitle("singlet peak value");
  gmultCorr->Draw("apm");
  // gPad->SetLogy();
  canCorr->BuildLegend();
  canCorr->SetGrid();

  fout->Append(grData);
  fout->Append(grFit);

  /* make comparison plots */
  TCanvas *canChan;
  bool firstPlot = true;
  for (int i = 0; i < 13; ++i)
  {
    if (fitWave[i][hcurve[i]->GetMaximumBin()] == 0)
      continue;
    // int ichan = atoi(string(fsname.substr(fsname.find_last_of("n") + 1, 1)).c_str());
    printf("%s %s \n", hcurve[i]->GetName(), hffit[i]->GetName());
    // hcurve[i]->Draw("");
    canChan = new TCanvas(Form("lightCurveChan%i", i), Form("lightCurveChan%i", i));
    printf("ffit %i max %f \n", i, hffit[i]->GetBinContent(hffit[i]->GetMaximumBin()));
    printf("chan %i data bin %i val %f peak bin fit %i val %f \n", i,
           hcurve[i]->GetMaximumBin(), hcurve[i]->GetBinContent(hcurve[i]->GetMaximumBin()),
           hffit[i]->GetMaximumBin(), hffit[i]->GetBinContent(hffit[i]->GetMaximumBin()));
    hffit[i]->GetYaxis()->SetRangeUser(10, 1.1 * hffit[i]->GetBinContent(hffit[i]->GetMaximumBin()));
    hffit[i]->GetXaxis()->SetRangeUser(650., 5000.);
    hcurve[i]->GetYaxis()->SetRangeUser(10, 1.1 * hffit[i]->GetBinContent(hffit[i]->GetMaximumBin()));
    hcurve[i]->GetXaxis()->SetRangeUser(650., 5000.);
    hcurve[i]->Draw("");
    hffit[i]->Draw("sames");
    canChan->BuildLegend();
    canChan->SetLogy();
  }

  // grData->Print("all");
  // grFit->Print("all");

  // eff geo corrected

  // delete gMinuit;
  // fout->Purge(1);
  // fout->ls();
  fout->Write();
  // fin->Close(); // cannot close before write because histos are on input file
}