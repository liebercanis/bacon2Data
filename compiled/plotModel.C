#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSpline.h"
#include "modelFit.hh"

std::vector<TH1D *> hwave0;
TFile *fout;
TSpline5 *spline5 = 0;

void plotModel(double PPM = 0.0)
{
  TString fileName = TString(Form("tbFitAllPPM%.2f.root", PPM));
  TFile *fin = new TFile(fileName, "readonly");
  printf(" file is %s \n", fileName.Data());

  /* get sum histos from file0 */
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
  }
  cout << " number waves " << fin->GetName() << " " << hwave0.size() << endl;
  if (hwave0.size() < 1)
  {
    return;
  }

  for (int ichan = 0; ichan < hwave0.size(); ++ichan)
    cout << hwave0[ichan]->GetName() << endl;

  TString canName;
  canName.Form("model-%s", fileName.Data());
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas *can1 = new TCanvas(canName, canName);
  bool first = true;
  // for (int ichan = hwave0.size() - 1; ichan >= 0; --ichan)
  // find max bin
  double maxBin = 0;
  for (int ichan = 0; ichan < hwave0.size(); ++ichan)
  {
    if (hwave0[ichan]->GetMaximum() > maxBin)
      maxBin = hwave0[ichan]->GetMaximum();
  }
  maxBin *= 1.1;

  double chanList[5] = {9, 8, 5, 0, 12};

  for (int ilevel = 0; ilevel < 5; ++ilevel) // skip PMT for now
  {
    int ichan = chanList[ilevel];
    // hwave0[ichan]->GetYaxis()->SetRangeUser(1.E-9, maxBin);
    if (first)
    {
      hwave0[ichan]->Draw();
      first = false;
    }
    else
      hwave0[ichan]->Draw("sames");
  }
  can1->BuildLegend();
  gPad->SetLogy();
}
