/*
   read and fit gains from summary
   updated sept 16 2024
   updated June 14 2025!
   new version jan 12 2026
*/
#include <ctime>
#include <iostream>
#include "TDirectory.h"
#include "TGraphErrors.h"

TFile *fin;
TFile *fout;
TString tag;
// double nominalGain = 227.4; // average
double theNominalGain;
/**************** define nominal gains ***************/
double nominalGain = 170.;     // was 160.0; set Jue 13 2025
double nominalTrigGain = 700.; //
double nominalQsumGain = 7050;
double nominalQsumTrigGain = 2.9E4;
double nominalPmtGain = 502.;
double nominalQsumPmtGain = 1713;
/******/
std::vector<TH1D *> peakList;
std::vector<TH1D *> sumList;

TGraphErrors *gSavedGain;
TGraphErrors *gNewGain;

std::vector<double> firstPeak;
std::vector<int> firstPeakChan;

std::vector<double> sipmSavedGain;
std::vector<double> sipmSavedGainError;
/* for gain graph */
std::vector<double> sipmGain;
std::vector<double> sipmGainError;
std::vector<double> sipmNumber;
std::vector<double> sipmNumberError;
std::vector<double> fSpeNumber;
std::vector<double> fSpeNumberError;
std::vector<double> fFitAdc;
std::vector<double> fFitAdcError;

std::string sdate;

enum
{
  CHANNELS = 14,
  NONSUMCHANNELS = CHANNELS - 1
};

enum
{
  MAXPOINTS = 4
};
int colors[11] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 30, 40};

/*  start of code */
string currentDate()
{
  time_t rawtime;
  struct tm *timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  char output[30];
  strftime(output, 30, "%Y-%m-%d-%H-%M", timeinfo);
  return string(output);
}

int getChan(TH1D *h)
{
  TString hname(h->GetName());
  int ichan = TString(hname(hname.Last('n') + 1, hname.Length())).Atoi();
  return ichan;
}

// Define a linear fit
double fline(double *x, double *par)
{
  return par[0] + x[0] * par[1];
}

int findPeakBin(TH1D *h, double fitStart, double fitEnd, double &integral)
{
  int ilow = h->FindBin(fitStart);
  int ihigh = h->FindBin(fitEnd);
  int ipeakBin = ilow;
  integral = 0;
  double ymax = h->GetBinContent(ilow);
  for (int i = ilow; i < ihigh; ++i)
  {
    integral += h->GetBinContent(i);
    if (h->GetBinContent(i) > ymax)
    {
      ymax = h->GetBinContent(i);
      ipeakBin = i;
    }
  }
  return ipeakBin;
}

bool readGains(TString fileName)
{
  TFile *gfin = new TFile(fileName, "readonly");
  if (gfin->IsZombie())
  {
    cout << "Error opening file" << fileName << endl;
    return false;
  }
  cout << " opened sipm gain file " << fileName << endl;
  gSavedGain = NULL;
  gfin->GetObject("gGain", gSavedGain);
  if (gSavedGain == NULL)
  {
    cout << "no gGain in file " << endl;
    return false;
  }

  printf("number saved gain points %i \n", gSavedGain->GetN());
  sipmSavedGain.clear();
  sipmSavedGainError.clear();
  sipmSavedGain.resize(NONSUMCHANNELS);
  sipmSavedGainError.resize(NONSUMCHANNELS);

  for (int i = 0; i < gSavedGain->GetN(); ++i)
  {
    int index = int(gSavedGain->GetPointX(i));
    sipmSavedGain[index] = gSavedGain->GetPointY(i);
    sipmSavedGainError[index] = gSavedGain->GetErrorY(i);
  }

  printf("\t\t\t stored gains %lu \n", sipmSavedGain.size());
  for (unsigned long j = 0; j < sipmSavedGain.size(); ++j)
  {
    printf(" %lu  saved gain %.4f error %.4f   \n", j, sipmSavedGain[j], sipmSavedGainError[j]);
  }
  return true;
}
void getHistosFromFile()
{
  // get histos from file
  TDirectory *gainSumDir = nullptr;
  fin->GetObject("gainDir", gainSumDir);

  //
  TIter next(gainSumDir->GetListOfKeys());
  TKey *key;
  while (TKey *key = (TKey *)next())
  {
    TClass *cl = gROOT->GetClass(key->GetClassName());

    if (!cl->InheritsFrom("TH1D"))
      continue;
    TH1D *h = (TH1D *)key->ReadObj();
    fout->Append(h);

    if (TString(h->GetName()).Contains("Peak"))
      peakList.push_back(h);

    if (TString(h->GetName()).Contains("Sum"))
      sumList.push_back(h);
  }
}

void doGains(TString type)
{
  firstPeak.clear();
  firstPeakChan.clear();
  // find first peaks
  cout << " do gains for type " << type << endl;
  if (type.Contains("Peak")) // for peak hits
  {
    for (unsigned i = 0; i < NONSUMCHANNELS; ++i)
    {
      firstPeak.push_back(peakList[i]->GetBinCenter(peakList[i]->GetMaximumBin()));
      firstPeakChan.push_back(getChan(peakList[i]));
    }
  }
  else // for sum hits
  {
    for (unsigned i = 0; i < NONSUMCHANNELS; ++i)
    {
      firstPeak.push_back(sumList[i]->GetBinCenter(sumList[i]->GetMaximumBin()));
      firstPeakChan.push_back(getChan(sumList[i]));
    }
  }

  for (unsigned i = 0; i < firstPeak.size(); ++i)
    printf(" type %s chan %i peak at %.2f\n", type.Data(), firstPeakChan[i], firstPeak[i]);

  // just do for chan 9 for starters
  // collect the points

  fSpeNumber.clear();
  fSpeNumberError.clear();
  fFitAdc.clear();
  fFitAdcError.clear();
  fSpeNumber.resize(MAXPOINTS);
  fSpeNumberError.resize(MAXPOINTS);
  fFitAdc.resize(MAXPOINTS);
  fFitAdcError.resize(MAXPOINTS);

  for (int i = 0; i < NONSUMCHANNELS; ++i)
  {
    // collect the points//
    printf("for channel %i :\n", i);
    unsigned nToFit = 0;
    // histogram to fit to from the appropriate list
    TH1D *hFit = peakList[i];
    if (type.Contains("Sum"))
      hFit = sumList[i];
    //
    for (int ipeak = 0; ipeak < MAXPOINTS; ++ipeak)
    {
      fSpeNumber[ipeak] = ipeak + 1;
      fSpeNumberError[ipeak] = 0;
      double width = 0;
      if (type.Contains("Sum"))
        width = 400;
      else
        width = 10;
      double fitStart = firstPeak[i] * double(ipeak + 1) - width;
      double fitEnd = firstPeak[i] * double(ipeak + 1) + width;
      double integral = 0;
      int ipeakBin = findPeakBin(hFit, fitStart, fitEnd, integral);
      if (integral < 10.)
        continue;
      ++nToFit;
      fFitAdc[ipeak] = hFit->GetBinLowEdge(ipeakBin);
      fFitAdcError[ipeak] = width / sqrt(integral); // gaussian error sigma/sqrt(N)
      // printf("point %i %f %f peak bin %i val %f integral %f \n", ipeak, fitStart, fitEnd, ipeakBin, peakList[i]->GetBinLowEdge(ipeakBin), integral);
    }

    // fit the line

    /* collect list of ipeakBins for integral > something and fill a TGraphErrors
     */
    printf("THEFIT type %s channel %i fit to %i points \n", type.Data(), i, nToFit);
    TGraphErrors *g = new TGraphErrors(nToFit, &fSpeNumber[0], &fFitAdc[0], &fSpeNumberError[0], &fFitAdcError[0]);
    g->SetName(Form("Graph%sChan%i", type.Data(), i));
    g->SetTitle(Form("Graph to fit for Chan%i", i));
    fout->Add(g);

    // fit line
    g->Fit("myLine", "Q");
    // g->GetHistogram()->GetListOfFunctions()->ls();
    TF1 *gFit = g->GetFunction("myLine");
    if (gFit != nullptr) // good fit
    {
      TCanvas *gcan = new TCanvas(Form("Gain%sMarkerChan%i", type.Data(), i), Form("%schan%i", type.Data(), i));
      gPad->SetLogy(0);
      gStyle->SetOptFit();
      // g->GetHistogram()->GetXaxis()->SetRangeUser(0, 4);
      // g->GetHistogram()->GetYaxis()->SetRangeUser(0, 1.5E5);
      gFit->SetLineStyle(5);
      gFit->SetLineWidth(1);
      g->SetMarkerStyle(20);
      g->SetMarkerSize(.5);
      g->Draw("APE1");
      gFit->Draw("same");
      gPad->SetGrid();
      gcan->Print(".pdf");
      fout->Add(g);
      fout->Add(gcan);

      printf("%s LINEFIT %i slope %f error %f \n", type.Data(), i, gFit->GetParameter(1), gFit->GetParError(1));
      sipmGain[i] = gFit->GetParameter(1);
      sipmGainError[i] = gFit->GetParError(1);
      sipmNumber[i] = double(i);
      sipmNumberError[i] = 0;
    }
    else // fit fails
    {
      printf("line269 !!!!!! fit to myLine fails for hist %i \n", i);
      for (unsigned ip = 0; ip < fSpeNumber.size(); ++ip)
        printf("peak %i %f %f \n", ip, fSpeNumber[ip], fFitAdc[ip]);
      sipmGain[i] = firstPeak[i];
      sipmGainError[i] = 0;
      sipmNumber[i] = double(i);
      sipmNumberError[i] = 0;
    }
  }

  // make final graph of gains
  TGraphErrors *absGain = new TGraphErrors(sipmGain.size(), &sipmNumber[0], &sipmGain[0], &sipmNumberError[0], &sipmGainError[0]);
  TString absoluteName;
  absoluteName = Form("gain%s", type.Data());
  TString absoluteTitle;
  absoluteTitle = Form("absolute %s gain  %s", type.Data(), tag.Data());
  absGain->SetName(absoluteName);
  absGain->SetTitle(absoluteTitle);
  absGain->GetHistogram()->GetYaxis()->SetTitle(absoluteTitle.Data());
  absGain->GetHistogram()->GetXaxis()->SetTitle("channel");
  absGain->SetTitle(absoluteTitle);
  absGain->SetMarkerStyle(23);
  absGain->SetMarkerColor(kRed);
  absGain->SetMarkerSize(1.3);
  fout->Append(absGain);
}

/***********************
 * macro main entry
 * ******************* */
void gainPost() // default all
{
  sipmGain.resize(13);
  sipmGainError.resize(13);
  sipmNumber.resize(13);
  sipmNumberError.resize(13);

  sdate = currentDate();
  printf(" making cans on %s \n", sdate.c_str());

  // put in explicit file name and get tag
  TString fileName("compiled/post-10_06_2025-10_06_2025-1951999.root");
  tag = TString(fileName(fileName.First("-") + 1, 21));
  cout << " gains from file " << fileName << " with date tag " << tag << endl;

  // open file
  fin = new TFile(fileName, "readonly");
  if (fin->IsZombie())
    return;

  // open output file
  std::string sdate = currentDate();
  fout = new TFile(Form("gains-%s-%s.root", tag.Data(), sdate.c_str()), "recreate");

  // function to fit line
  TF1 *line = new TF1("myLine", fline, 0, 2.E5, 2);

  // read old gain file for comparison
  TString savedGainTag = TString("gains-2024-02-15-17-26-alpha");
  TString gainFileName = TString(getenv("BOBJ")) + savedGainTag + TString(".root");
  printf("read gains from file %s \n", gainFileName.Data());
  if (!readGains(gainFileName))
  {
    printf("no gain file %s so exit \n", gainFileName.Data());
    exit(0);
  }
  gSavedGain->GetHistogram()->GetYaxis()->SetTitle("absolute saved gain ");
  gSavedGain->GetHistogram()->GetXaxis()->SetTitle("channel");
  gSavedGain->SetName(savedGainTag);
  gSavedGain->SetTitle(Form("saved gains %s ", savedGainTag.Data()));
  fout->Append(gSavedGain);

  getHistosFromFile();
  printf("line171 have peak gain hists %lu and sum gain hists %lu \n", peakList.size(), sumList.size());
  doGains("Peak");
  doGains("Sum");

  fout->ls();
  fout->Write();
}
