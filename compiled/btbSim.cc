//  simple sim of btb
// April 28 2025
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TNtuple.h"
#include "TFile.h"
#include "modelFitGamma.hh"
#include "TBRawEvent.hxx"
#include "TBRawRun.hxx"
std::string sdate;
// time is in microseconds
using namespace TMath;
TFile *fout;
TRandom3 *ran;

// for writing raw data
TBRawEvent *rawEvent;
TBRawRun *rawRun;
std::vector<uint16_t> wave;

bool writeRawData = true;

modelFit *models[NCHAN];
TH1D *hPhoton[NCHAN];
TH1D *hSignal[NCHAN];
Long64_t totalPhotons;
Long64_t ncount[NCHAN];
TH1D *hCount;
TH1D *hResponse;
double gain = 1.E2;
// 2*14         // ns

double LY = 25.6; //  photone/kev Doke
double numPhotons = 60 * LY;
double singletFrac = 0.20;
int binWidth = 2;

int triggerStart = 700;
double speMPV = double(triggerStart);
double speSigma = 20.; // ns from single PI data fit
TF1 *speLandau;

// spe landau shape
static double myLandau(Double_t *xx, Double_t *par)
{
  double x = xx[0];
  return par[2] * TMath::Landau(x, par[0], par[1], true); // normalized
}

void convolve(TH1D *hist, double time) // time is when photon arrives
{
  int offsetBin = hist->FindBin(540.); // read offf of Response histogram
  int startBin = hist->FindBin(time);
  // printf("convolve: %s time %f startBin %i offsetBin %i\n", hist->GetName(), time, startBin, offsetBin);
  for (int ib = startBin; ib < hist->GetNbinsX(); ++ib)
    hist->SetBinContent(ib, hist->GetBinContent(ib) + gain * hResponse->GetBinContent(ib - startBin + offsetBin));
}
void btb(int ngen = 10000000)
{
  printf(" btb sim generate ngen =  %i \n", ngen);

  ran = new TRandom3();
  // open raw output file
  rawRun = NULL;
  char dateTag[10];
  time_t rawtime;
  time(&rawtime);
  struct tm *timeinfo;
  timeinfo = localtime(&rawtime);
  strftime(&dateTag[0], 10, "%m_%d_%Y", timeinfo);
  TString tdateTag = TString(dateTag).Data();
  TString fullname = (Form("btbSim-%s-%i.root", tdateTag.Data(), ngen));
  fout = new TFile(fullname, "recreate"); // DEF made to update rather than recreate so that it doesn't write over a file already made.
  printf("opened output file %s date %s \n", fout->GetName(), tdateTag.Data());

  if (writeRawData)
  {
    rawRun = new TBRawRun(TString(dateTag));
    rawRun->updateTime(rawtime);
    rawRun->btree->SetTitle("simulation");
    // rawRun->print();
  }

  hCount = new TH1D("Count", "hit count", 13, 0, 13);
  // landau response function
  speLandau = new TF1("myLandau", myLandau, 0, totalBins * theBinWidth, 3);
  // set SPE response parameters
  speLandau->SetParName(1, "MPV");
  speLandau->SetParameter(0, speMPV);
  speLandau->SetParName(1, "sigma");
  speLandau->SetParameter(1, speSigma);
  speLandau->SetParName(2, "norm");
  speLandau->SetParameter(2, 1); // single SPE
  // modelFit::modelFit(int theFit, int ichan, double ppm)
  hResponse = new TH1D("Response", "sipm response", totalBins, 0, totalBins * (theBinWidth));
  hResponse->GetXaxis()->SetTitle("time [ns]");
  hResponse->GetYaxis()->SetTitle("photons/2ns");

  // fill response
  for (int ib = 1; ib < hResponse->GetNbinsX(); ++ib)
    hResponse->SetBinContent(ib, speLandau->Eval(hResponse->GetBinCenter(ib)) * double(binWidth));
  printf(" landau response integral %E \n", hResponse->Integral());

  // make individual light curves
  for (int ih = 0; ih < NCHAN; ++ih)
    models[ih] = new modelFit(MODELALL, ih, 0);

  TDirectory *histDir = fout->mkdir("histDir");
  histDir->cd();
  for (int ih = 0; ih < NCHAN; ++ih)
  {
    // modelFit::modelFit(int theFit, int ichan, double ppm)
    hPhoton[ih] = new TH1D(Form("Photon%i", ih), Form("Photon%i-level%i", ih, level(ih)), totalBins, 0, totalBins * (theBinWidth));
    hPhoton[ih]->GetXaxis()->SetTitle("time [ns]");
    hPhoton[ih]->GetYaxis()->SetTitle("photons/2ns");
    hPhoton[ih]->SetDirectory(nullptr);
    hSignal[ih] = new TH1D(Form("Signal%i", ih), Form("Signal%i-level%i", ih, level(ih)), totalBins, 0, totalBins * (theBinWidth));
    hSignal[ih]->GetXaxis()->SetTitle("time [ns]");
    hSignal[ih]->GetYaxis()->SetTitle("photons/2ns");
    hSignal[ih]->SetDirectory(nullptr);
  }

  // print info for channel
  for (int ich = NCHAN - 2; ich >= 0; --ich)
  {
    double eff = effGeoFunc(ich);
    int ilevel = level(ich);
    printf(" chan %i level %i distance %f eff %E \n", ich, ilevel, distanceLevel[ilevel], eff);
  }

  // print info for pmt
  double eff = effGeoFunc(12);
  int ilevel = level(12);
  printf(" chan %i level %i distance %f eff %E \n", 12, ilevel, distanceLevel[ilevel], eff);

  /* static double tTriplet0 = 1600.0; ns
    static double tSinglet0 = 7.0; ns
  */

  // loop over events
  totalPhotons = 0;
  for (int iev = 0; iev < ngen; ++iev)
  {

    int nPhotonsEvent = (int)ran->Gaus(numPhotons, sqrt(numPhotons));
    totalPhotons += nPhotonsEvent;
    if (iev / 100 * 100 == iev)
      printf("... event %i total photons %E \n", iev, double(totalPhotons));

    // loop over channels
    for (int ich = 0; ich < NCHAN; ++ich)
    {
      hPhoton[ich]->Reset("ICESM");
      hSignal[ich]->Reset("ICESM");
      if (rawRun)
      {
        rawEvent = rawRun->getDet(ich); // If channel branch doesn't exist getDet calls addDet
        rawEvent->clear();
        rawEvent->channel = ich;
        rawEvent->trigger = iev;
      }
      // rawEvent->time = EventInfo->TriggerTimeTag;

      double eff = effGeoFunc(ich);
      double nsmean = double(nPhotonsEvent) * eff * singletFrac;
      double ntmean = double(nPhotonsEvent) * eff - nsmean;
      int nsinglet = ran->Poisson(nsmean);
      int ntriplet = ran->Poisson(ntmean);

      ncount[ich] += nsinglet + ntriplet;

      if (iev / 100000 * 100000 == iev)
        printf("event %i ich %i eff %E nPhotons  %E \n", iev, ich, eff, double(ncount[ich]));

      // singlet times
      for (int it = 0; it < nsinglet; ++it)
      {
        double time = triggerStart + ran->Exp(tSinglet0);
        hPhoton[ich]->Fill(time);
        convolve(hSignal[ich], time);
      }

      // triplet times
      for (int it = 0; it < ntriplet; ++it)
      {
        double time = triggerStart + ran->Exp(tTriplet0);
        hPhoton[ich]->Fill(time);
        convolve(hSignal[ich], time);
      }
      // file wave for this channel
      if (rawRun)
      {
        wave.clear();
        for (int ibin = 1; ibin < hSignal[ich]->GetNbinsX(); ++ibin)
          wave.push_back((uint16_t)hSignal[ich]->GetBinContent(ibin));
        rawEvent->rdigi = wave;
      }
    } // end channel loop
    // end of event loop

    /* event histograms */
    TString histName;
    fout->cd();
    if (histDir->GetList()->GetEntries() < 100)
      for (int ih = 0; ih < NCHAN; ++ih)
      {
        histDir->cd();
        histName.Form("hPhotonCh%iEv%i", ih, iev);
        TH1D *hPhotonEvent = (TH1D *)hPhoton[ih]->Clone(histName);
        hPhotonEvent->SetTitle(histName);
        //
        histName.Form("hSignalCh%iEv%i", ih, iev);
        TH1D *hSignalEvent = (TH1D *)hSignal[ih]->Clone(histName);
        hSignalEvent->SetTitle(histName);
      }
  }
  if (rawRun)
    rawRun->fill();

  for (int ich = 0; ich < NCHAN; ++ich)
    hCount->SetBinContent(ich + 1, ncount[ich]);
  // summary
  printf("generated %i \n", ngen);
  for (int ih = 1; ih < NCHAN; ++ih)
  {
    printf(" chan %i photons %i\n", ih, (int)hCount->GetBinContent(ih));
  }
  fout->ls();
  fout->Write();
  fout->Close();

  printf("end of btbgen  ngen %i exit\n", ngen);
}

// static TBRun *theTBRun;
int main(int argc, char *argv[])
{
  int ngen = 1000000;

  std::cout << "  usage: btbSim <ngen> default 1000000  " << argv[0] << std::endl;
  printf("\n ");
  if (argc > 1)
  {
    ngen = atoi(argv[1]);
  }

  btb(ngen);
  printf("... %s ngen %i exit\n", argv[0], ngen);
  exit(0);
}
