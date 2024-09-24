#include <ctime>
#include <iostream>
#include <iterator>
#include <locale>
#include <TString.h>
#include <TROOT.h>
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TDatime.h"
#include "TDirectory.h"
#include "TSystemDirectory.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TList.h"
#include "TF1.h"
#include "TNtuple.h"
#include <TROOT.h>
#include <TKey.h>
#include <TBranch.h>
#include <TBranchElement.h>
#include <TVirtualFFT.h>
#include <TChain.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TBranchElement.h>
#include <TFile.h>
#include <Rtypes.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TLeaf.h>
#include <TFormula.h>
#include <TStyle.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
//
TNtuple *ntThreshAll;
TNtuple *ntThreshold;
Float_t event;
Float_t chan;
Float_t sampleLow;
Float_t ddigiLow;
Float_t sampleHigh;
Float_t ddigiHigh;
Float_t maxBin;
Float_t adcMax;

Float_t ev;
Float_t ch;
Float_t samp;
Float_t ddigi;

std::vector<TH1D *> vddigiLow;
std::vector<TH1D *> vddigiHigh;
std::vector<TH1D *> vadcMax;
TH2D *hThresh;
TH2D *hThreshTrig;
TH2D *hThreshPmt;
TH2D *hThreshLow;
TH2D *hThreshTrigLow;
TH2D *hThreshPmtLow;
enum
{
  NCHAN = 13
};

void Init()
{
  // Set branch addresses and branch pointers
  if (!ntThreshold)
    return;
  ntThreshold->SetBranchAddress("event", &event);
  ntThreshold->SetBranchAddress("chan", &chan);
  ntThreshold->SetBranchAddress("sampleLow", &sampleLow);
  ntThreshold->SetBranchAddress("ddigiLow", &ddigiLow);
  ntThreshold->SetBranchAddress("sampleHigh", &sampleHigh);
  ntThreshold->SetBranchAddress("ddigiHigh", &ddigiHigh);
  ntThreshold->SetBranchAddress("maxBin", &maxBin);
  ntThreshold->SetBranchAddress("adcMax", &adcMax);

  if (!ntThreshAll)
    return;
  ntThreshAll->SetBranchAddress("event", &ev);
  ntThreshAll->SetBranchAddress("chan", &ch);
  ntThreshAll->SetBranchAddress("sample", &samp);
  ntThreshAll->SetBranchAddress("ddigi", &ddigi);
}

void thresh(TString tag = "09_10_2024")
{
  ntThreshAll = NULL;
  ntThreshold = NULL;
  TString fileName;
  fileName.Form("caenData/anaCRun-run-%s-diffStep5.root", tag.Data());
  TFile *fin = new TFile(fileName, "readonly");

  if (fin->IsZombie())
    return;

  printf("opened file %s \n", fileName.Data());

  fin->GetObject("ntThresholdAll", ntThreshAll);
  if (!ntThreshAll)
  {
    printf("no ntuple ntThresholdAll \n");
    return;
  }

  fin->GetObject("ntThresold", ntThreshold);
  if (!ntThreshold)
  {
    printf("no ntuple ntThreshold \n");
    return;
  }
  // output file
  TFile *fout = new TFile(Form("ntThreshold%s.root", tag.Data()), "recreate");

  TString hname;
  TString htitle;
  for (int i = 0; i < NCHAN; ++i)
  {
    hname.Form("DigiLowChan%i", i);
    htitle.Form("DigiLowChan%i", i);
    vddigiLow.push_back(new TH1D(hname, htitle, 1000, 0, 1600));
    hname.Form("DigiHighChan%i", i);
    htitle.Form("DigiHighChan%i", i);
    vddigiHigh.push_back(new TH1D(hname, htitle, 1000, 0, 1600));

    hname.Form("AdcMaxChan%i", i);
    htitle.Form("AdcMaxChan%i", i);
    vadcMax.push_back(new TH1D(hname, htitle, 1000, 0, 1600));
  }

  hThreshLow = new TH2D("ThreshLow", "non trigger thresh", 100, 0, 1000, 100, 0, 1000);
  hThresh = new TH2D("Thresh", "non trigger thresh", 100, 0, 1000, 100, 0, 1000);
  hThreshLow->SetMarkerStyle(21);
  hThreshLow->SetMarkerColor(kRed);

  hThreshTrigLow = new TH2D("ThreshTrigLow", "trigger thresh", 1000, 0, 10000, 1000, 0, 10000);
  hThreshTrig = new TH2D("ThreshTrig", "trigger thresh", 1000, 0, 10000, 1000, 0, 10000);

  hThreshPmtLow = new TH2D("ThreshPmtLow", "trigger PMT", 100, 0, 100, 100, 0, 100);
  hThreshPmt = new TH2D("ThreshPmt", "trigger PMT", 100, 0, 100, 100, 0, 100);

  fout->ls();

  Init();
  Long64_t nEntries = ntThreshAll->GetEntries();
  printf("looping over ntThreshAll %lld \n", nEntries);

  nEntries = ntThreshold->GetEntries();
  printf("looping over ntThresh %lld \n", nEntries);

  // loop over entries
  for (Long64_t entry = 0; entry < nEntries; ++entry)
  {
    ntThreshold->GetEntry(entry);
    printf("ev %i chan %i %f %i %f  \n", int(event), int(chan), ddigiLow, int(maxBin), adcMax);
    vddigiLow[int(chan)]->Fill(-1. * ddigiLow); // flip the sign because these are negative
    vddigiHigh[int(chan)]->Fill(ddigiHigh);
    vadcMax[int(chan)]->Fill(adcMax);

    if (chan < 9 && adcMax < 75)
      hThreshLow->Fill(-ddigiLow, ddigiHigh);
    else if (chan < 9)
      hThresh->Fill(-ddigiLow, ddigiHigh);

    if (chan > 8 && adcMax < 75)
      hThreshTrigLow->Fill(-ddigiLow, ddigiHigh);
    else if (chan > 8)
      hThreshTrig->Fill(-ddigiLow, ddigiHigh);

    if (chan == 12 && adcMax < 20)
      hThreshPmtLow->Fill(-ddigiLow, ddigiHigh);
    else if (chan == 12)
      hThreshPmt->Fill(-ddigiLow, ddigiHigh);
  }
  fout->Write();
}
