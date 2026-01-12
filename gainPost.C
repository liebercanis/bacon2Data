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
bool doFit = false;
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
std::vector<double> sipmSavedGain;
std::vector<double> sipmSavedGainError;
std::vector<double> sipmGain;
std::vector<double> sipmGainError;

std::vector<double> firstPeak;
std::vector<int> firstPeakChan;
std::vector<double> firstSum;
std::vector<int> firstSumChan;

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

    if (TString(h->GetName()).Contains("Peak"))
      peakList.push_back(h);

    if (TString(h->GetName()).Contains("Sum"))
      sumList.push_back(h);
  }
}

/***********************
 * macro main entry
 * ******************* */
void gainPost() // default all
{
  sdate = currentDate();
  printf(" making cans on %s \n", sdate.c_str());

  // put in explicit file name and get tag
  TString fileName("compiled/post-10_06_2025-10_06_2025-969976.root");
  TString tag = TString(fileName(fileName.First("-") + 1, 21));
  cout << " gains from file " << fileName << " with date tag " << tag << endl;

  // open file
  fin = new TFile(fileName, "readonly");
  if (fin->IsZombie())
    return;

  // open output file
  std::string sdate = currentDate();
  fout = new TFile(Form("gainPeak-%s-%s.root", tag.Data(), sdate.c_str()), "recreate");

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

  // find first peaks
  for (unsigned i = 0; i < NONSUMCHANNELS; ++i)
  {
    firstPeak.push_back(peakList[i]->GetBinCenter(peakList[i]->GetMaximumBin()));
    firstPeakChan.push_back(getChan(peakList[i]));
  }

  for (unsigned i = 0; i < firstPeak.size(); ++i)
    printf(" %i %s chan %i peak at %.2f\n", i, peakList[i]->GetName(), firstPeakChan[i], firstPeak[i]);

  // find first sum peaks
  for (unsigned i = 0; i < NONSUMCHANNELS; ++i)
  {
    firstSum.push_back(sumList[i]->GetBinCenter(sumList[i]->GetMaximumBin()));
    firstSumChan.push_back(getChan(sumList[i]));
  }

  for (unsigned i = 0; i < firstPeak.size(); ++i)
    printf("%i %s chan %i sum at %.2f\n", i, sumList[i]->GetName(), firstSumChan[i], firstSum[i]);

  // just do for chan 9 for starters
  // collect the points
  std::vector<double> fSpeNumber;
  std::vector<double> fSpeNumberError;
  std::vector<double> fFitAdc;
  std::vector<double> fFitAdcError;
  fSpeNumber.resize(MAXPOINTS);
  fSpeNumberError.resize(MAXPOINTS);
  fFitAdc.resize(MAXPOINTS);
  fFitAdcError.resize(MAXPOINTS);

  for (int i = 0; i < NONSUMCHANNELS; ++i)
  {
    // collect the points//
    printf("for channel %i :\n", i);
    unsigned nToFit = 0;
    for (int ipeak = 0; ipeak < 4; ++ipeak)
    {
      fSpeNumber[ipeak] = ipeak + 1;
      fSpeNumberError[ipeak] = 0;
      double width = firstPeak[i] / 10.;
      double fitStart = firstPeak[i] * double(ipeak + 1) - width;
      double fitEnd = firstPeak[i] * double(ipeak + 1) + width;
      double integral = 0;
      int ipeakBin = findPeakBin(peakList[i], fitStart, fitEnd, integral);
      if (integral < 10.)
        continue;
      ++nToFit;
      fFitAdc[ipeak] = peakList[i]->GetBinLowEdge(ipeakBin);
      fFitAdcError[ipeak] = width;
      printf("point %i %f %f peak bin %i val %f integral %f \n", ipeak, fitStart, fitEnd, ipeakBin, peakList[i]->GetBinLowEdge(ipeakBin), integral);
    }

    // fit the line

    /* collect list of ipeakBins for integral > something and fill a TGraphErrors
     */
    TGraphErrors *g = new TGraphErrors(nToFit, &fSpeNumber[0], &fFitAdc[0], &fSpeNumberError[0], &fFitAdcError[0]);
    g->SetName(Form("GraphChan%i", i));
    g->SetTitle(Form("Graph to fit for Chan%i", i));
    fout->Append(g);
    /* fit line and get slope and error if nToFit greater than 1
    collect into new graph of new gains
    */
  }

  fout->ls();
  fout->Write();
}
