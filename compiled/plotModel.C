#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSpline.h"
#include "modelFit.hh"

std::vector<TH1D *> hwave0;
std::vector<TH1D *> hfit0;
TFile *fout;
TSpline5 *spline5 = 0;
int colors[NCHAN] = {kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kSpring, kTeal, kAzure, kViolet, kPink, kGray};

enum
{
  MAXCHAN = 4
};
int chanList[MAXCHAN] = {9, 8, 5, 0};

// get smallest nonzero bin
double getMinBin(TH1D *h)
{
  double min = 1.E9;
  for (int ibin = 0; ibin < h->GetNbinsX(); ++ibin)
    if (h->GetBinContent(ibin) > 0 && h->GetBinContent(ibin) < min)
    {
      min = h->GetBinContent(ibin);
    }
  return min;
}

void plotModel(TString fileName)
{
  TFile *fin = new TFile(fileName, "readonly");
  printf(" file is %s \n", fileName.Data());

  /* get sum histos from file0
  TIter next(fin->GetListOfKeys());
  TKey *key;
  while (TKey *key = (TKey *)next())
  {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1D"))
      continue;
    TH1D *h = (TH1D *)key->ReadObj();
    TString hname(h->GetName());
    // cout << " hist " << hname << endl;
    if (hname.Contains("fitWave"))
    {
      h->SetName(Form("%s", h->GetName()));
      h->SetTitle(Form("%s", h->GetName()));
      hwave0.push_back(h);
    }
      */
  TH1D *h = nullptr;
  TString hname;
  for (int ichan = 0; ichan < MAXCHAN; ++ichan)
  {
    hname.Form("fitWaveFitChan%i", chanList[ichan]);
    fin->GetObject(hname, h);
    if (h)
      cout << "got " << hname << endl;
    hfit0.push_back(h);
    //
    hname.Form("RunPeakWave%i", chanList[ichan]);
    fin->GetObject(hname, h);
    if (h)
      cout << "got " << hname << endl;
    hwave0.push_back(h);
  }
  printf(" got %lu %lu \n", hwave0.size(), hfit0.size());

  if (hwave0.size() != MAXCHAN || hfit0.size() != MAXCHAN)
  {
    return;
  }

  for (int ichan = 0; ichan < MAXCHAN; ++ichan)
  {
    printf(" %i data %s fit %s \n", ichan, hwave0[ichan]->GetName(), hfit0[ichan]->GetName());
  }

  return;

  TString canName;
  canName.Form("model-%s", fileName.Data());
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas *can1 = new TCanvas(canName, canName);
  bool first = true;
  // for (int ichan = hwave0.size() - 1; ichan >= 0; --ichan)
  // find max bin
  double maxBin = 0;
  double minBin = 1E9;
  for (int ichan = 0; ichan < hwave0.size(); ++ichan)
  {
    if (hwave0[ichan]->GetMaximum() > maxBin)
      maxBin = hwave0[ichan]->GetMaximum();
    hwave0[ichan]->SetLineColor(colors[ichan]);
    minBin = getMinBin(hwave0[ichan]);
  }
  maxBin *= 1.1;
  printf("min bin set to %E max bin set to %E\n", minBin, maxBin);

  for (int ihist = 0; ihist < hwave0.size(); ++ihist) // skip PMT for now
  {
    hwave0[ihist]->GetYaxis()->SetRangeUser(minBin, 1.1 * maxBin);
    if (first)
    {
      hwave0[ihist]->Draw();
      first = false;
    }
    else
      hwave0[ihist]->Draw("sames");
  }
  can1->BuildLegend();
  gPad->SetLogy();

  for (int ihist = 0; ihist < hwave0.size(); ++ihist) // skip PMT for now
  {
    minBin = getMinBin(hwave0[ihist]);
    maxBin = hwave0[ihist]->GetMaximum();
    printf(" hist %i %s and %s %E to %E \n", ihist, hwave0[ihist]->GetName(), hfit0[ihist]->GetName(), minBin, maxBin);
    hwave0[ihist]->GetYaxis()->SetRangeUser(minBin, 1.1 * maxBin);
    hfit0[ihist]->GetYaxis()->SetRangeUser(minBin, 1.1 * maxBin);
    canName.Form("channel-%i", chanList[ihist]);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TCanvas *can1 = new TCanvas(canName, canName);
    hwave0[ihist]->Draw();
    hfit0[ihist]->Draw("sames");
    gPad->SetLogy();
  }
}
