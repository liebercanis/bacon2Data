//  simple sim of btb
// April 28 2025
#include <iostream>
#include <fstream>
#include <numeric>
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TNtuple.h"
#include "TFile.h"
#include "Math/Vector3D.h"
#include "modelFitGamma.hh"
#include "TBRawRun.hxx"
#include "TBSimRun.hxx"
std::string sdate;
// time is in microseconds
using namespace TMath;
TFile *fout;
TRandom3 *ran;

// for writing raw data
TBRawEvent *rawEvent;
TBRawRun *rawRun;
TBSimRun *simRun;
std::vector<uint16_t> wave;

bool writeRawData = true;

modelFit *models[NCHAN];
TNtuple *ntOrigin;
TNtuple *ntTrigCh;
TNtuple *ntTrig;
TH1D *hPhoton[NCHAN];
TH1D *hConvolve[NCHAN];
TH1D *hSignal[NCHAN];
Long64_t totalPhotons;
Long64_t ncount[NCHAN];
TH1D *hCount;
TH1D *hResponse;
TH1D *hTime;
uint16_t maxAdc = pow(2, 14);
double gain = nominalGain;
double landauMax = 0.018063;
// 2*14         // ns
/* parameters quoted in talk  "A new optical model for LEGEND-200
with remage" Manuel Huber <ge38nap@mytum.de>, Luigi Pertoldi
LEGEND collaboration meeting · March 25, 2025 */
double LY = 25.6; //  photone/kev Doke
double numPhotons = 60 * LY;
double singletFrac = 0.20;
int binWidth = 2;
double noiseToSignal = 0.04;
double baseline = 1100.; // 1100; // ADC
double thePPM = 0.0;
double meanFreePath = 1.0; // guess for 60kev gamma in cm

ROOT::Math::XYZVector eventOrigin(0, 0, 0);

/*
efficiencies  PMTQE175 = 0.38;
static double QEff128(double ppm, double dist)
*/
double eff[NCHAN];
int triggerStart = 2 * 700; // 730; sipm rise time
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
  int offsetBin = hist->FindBin(640.); // read offf of Response histogram
  int startBin = hist->FindBin(time);
  // printf("convolve: %s time %f startBin %i offsetBin %i\n", hist->GetName(), time, startBin, offsetBin);
  for (int ib = startBin; ib < hist->GetNbinsX(); ++ib)
    hist->SetBinContent(ib, hist->GetBinContent(ib) + gain / landauMax * hResponse->GetBinContent(ib - startBin + offsetBin));
}

ROOT::Math::XYZVector getXYZVector(double r, double theta, double phi) // angles in radians
{
  // construct XYZ coordinates
  double cost = cos(theta);
  double sint = sin(theta);
  double cosp = cos(phi);
  double sinp = sin(phi);
  double x = r * sint * cosp;
  double y = r * sint * sinp;
  double z = r * cost;
  ROOT::Math::XYZVector pos(x, y, z);
  return pos;
}

double effGeoSim(int ichan, bool show = false) // uses PositionVector3D eventOrigin;
{
  bool isTrig = ichan == 9 || ichan == 10 || ichan == 11;
  if (!isTrig)
    return effGeoFunc(ichan);

  int ilevel = level(ichan);
  double e = 1.0;
  if (ilevel != 0)
    return e;

  // det positions Georgia May 2025
  double trigRadius = 1.486;
  // convert to radians
  double trigTheta = 55.06 / 360. * TMath::TwoPi(); // 11,10,9
  double trigPhi[3];
  trigPhi[0] = 0.;                           // 9
  trigPhi[1] = 240. / 360. * TMath::TwoPi(); // 10
  trigPhi[2] = 120. / 360. * TMath::TwoPi(); // 11

  ROOT::Math::XYZVector rSipm = getXYZVector(trigRadius, trigTheta, trigPhi[ichan - 9]);
  ROOT::Math::XYZVector relative = rSipm - eventOrigin;

  double distance2 = relative.Mag2();
  /* Area of SiPMs is 6.0mm x 6.0mm
      Channels 6, 7, and 8 are at 11.6 cm
      from the source Channels 3, 4, and 5 are at 23.2 cm
      from the source Channels 0, 1, and 2 are at 34.8 cm from the source Channel 12 is at 36 cm from the source.
  */
  double aPmt = TMath::Pi() / 4.0 * pow(6.4, 2); // R11410-20  Effective area : 64 mm dia
  double a = pow(0.6, 2.);                       // SIPM area
  double b = 4.0 * TMath::Pi();

  // correct solid angle
  ROOT::Math::XYZVector runit = rSipm.Unit();
  ROOT::Math::XYZVector ounit = relative.Unit();
  double cos = runit.Dot(ounit);
  // printf(" chan %i cos %f \n", ichan, cos);
  if (cos < 0.)
    cos = 0; // origin behind sipm
  a *= cos;

  e = a / b / distance2;
  // shift phi for printing.
  double localPhi = relative.Phi() * 360. / TMath::TwoPi();
  if (localPhi < 0)
    localPhi += 360.;
  if (show)
    printf("effGeoSim ichan  %i level %i  cos %f dist %f R,Theta,Phi (%f,%f,%f) area %f eff %f \n", ichan, ilevel, cos, sqrt(distance2), relative.R(), relative.Theta() * 360. / TMath::TwoPi(), localPhi, a, e);
  return e;
}

void btb(int ngen = 10000000)
{
  printf(" btb sim generate ngen =  %i LY %.1f photons/kev * 60 = %.1f \n", ngen, LY, numPhotons);

  /* channel efficiences */
  for (int i = 0; i < NCHAN - 1; ++i)
  {
    double dist = distanceLevel[level(i)];
    eff[i] = QEff128(thePPM, dist);
  }
  eff[NCHAN - 1] = 0.0; // pmt sees light > 175 nm

  printf(" efficiences at 128 nm PPM=%f\n", thePPM);
  for (int i = 0; i < NCHAN; ++i)
    printf("\t chan %i eff %f \n", i, eff[i]);

  ran = new TRandom3();
  // open raw output file/raw/
  rawRun = NULL;

  time_t rawtime;
  struct tm *timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  char output[30];
  strftime(output, 30, "%Y-%m-%d-%H-%M", timeinfo);
  TString tdateTag = TString(output);
  TString fullname = (Form("btbSim-%s-%i.root", tdateTag.Data(), ngen));
  fout = new TFile(fullname, "recreate"); // DEF made to update rather than recreate so that it doesn't write over a file already made.
  printf("opened output file %s date %s \n", fout->GetName(), tdateTag.Data());
  cout << tdateTag << endl;

  // make output tree
  simRun = new TBSimRun("sim0");
  simRun->clear();

  if (writeRawData)
  {
    rawRun = new TBRawRun(tdateTag);
    rawRun->updateTime(rawtime);
    rawRun->btree->SetTitle("simulation");
    // rawRun->print();
  }

  ntOrigin = new TNtuple("ntOrigin", " event origin ", "ev:r:cos:theta:phi:x:y:z");
  ntTrigCh = new TNtuple("ntTrigCh", " trigger info by channel ", "ev:ch:qsum:psum:nph:r:theta:phi:x:y:z");
  ntTrig = new TNtuple("ntTrig", " trigger info by event  ", "ev:nph9:nph10:nph11:r:theta:phi:x:y:z");
  hCount = new TH1D("Count", "hit count", 13, 0, 13);
  hTime = new TH1D("Time", "photon time ", 7500, 0, 2 * 7500);
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
  double landauMax = hResponse->GetBinContent(hResponse->GetMaximumBin());
  printf(" landau response integral %E  gain %f  landauMax %f SPE %E \n", hResponse->Integral(), gain, landauMax, gain / landauMax);

  double sigmaNoise = gain * noiseToSignal;
  TH1D *hNoise = new TH1D("Noise", "Noise", 200, -10 * sigmaNoise, 10 * sigmaNoise);

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
    //
    hConvolve[ih] = new TH1D(Form("Convolve%i", ih), Form("Convolve%i-level%i", ih, level(ih)), totalBins, 0, totalBins * (theBinWidth));
    hConvolve[ih]->GetXaxis()->SetTitle("time [ns]");
    hConvolve[ih]->GetYaxis()->SetTitle("photons/2ns");
    hConvolve[ih]->SetDirectory(nullptr);
    //
    hSignal[ih] = new TH1D(Form("Signal%i", ih), Form("Signal%i-level%i", ih, level(ih)), totalBins, 0, totalBins * (theBinWidth));
    hSignal[ih]->GetXaxis()->SetTitle("time [ns]");
    hSignal[ih]->GetYaxis()->SetTitle("photons/2ns");
    hSignal[ih]->SetDirectory(nullptr);
  }

  // print info for channel
  printf("\n******* efficienies ****\n ");
  for (int ich = 8; ich >= 0; --ich)
  {
    double eff = effGeoFunc(ich);
    int ilevel = level(ich);
    printf(" chan %i level %i distance %f eff %E \n", ich, ilevel, distanceLevel[ilevel], eff);
  }

  // trigger sipms
  for (int ich = 9; ich < 12; ++ich)
  {
    double eff = effGeoSim(ich, true);
  }
  printf("***********\n ");

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

    /********  generate gamma position  *********/
    double gammaR = abs(ran->Exp(meanFreePath));
    double gammaCosTheta = 2. * ran->Rndm() - 1.;   // cos flat from 1 to -1
    double gammaPhi = TMath::TwoPi() * ran->Rndm(); // phi flat from 0 to 2pi
    eventOrigin = getXYZVector(gammaR, acos(gammaCosTheta), gammaPhi);
    double localPhi = eventOrigin.Phi() * 360. / TMath::TwoPi();
    if (localPhi < 0)
      localPhi += 360.;

    ntOrigin->Fill(iev, eventOrigin.R(), cos(eventOrigin.Theta()), eventOrigin.Theta() * 360. / TMath::TwoPi(), localPhi), eventOrigin.X(), eventOrigin.Y(), eventOrigin.Z();

    if (iev / 100 * 100 == iev)
      printf("... event %i total photon %0.f (r,cosTheta,phi) = (%f, %f, %f) (r,theta,Phi) = (%f , %f ,%f ) \n", iev, double(totalPhotons), gammaR, gammaCosTheta, gammaPhi, eventOrigin.R(), eventOrigin.Theta() * 360. / TMath::TwoPi(), localPhi);
    // get event position

    // loop over channels
    for (int ich = 0; ich < NCHAN; ++ich)
    {
      bool invert = ich > 8; // invert trigger 9,10,11 and PMT
      hPhoton[ich]->Reset("ICESM");
      hConvolve[ich]->Reset("ICESM");
      hSignal[ich]->Reset("ICESM");
      if (rawRun)
      {
        rawEvent = rawRun->getDet(ich); // If channel branch doesn't exist getDet calls addDet
        rawEvent->clear();
        rawEvent->channel = unsigned(ich);
        rawEvent->trigger = unsigned(iev);
      }
      TDet *det = simRun->getDet(ich); // If channel branch doesn't exist getDet calls addDet
      det->clear();
      det->event = iev;
      det->trigger = iev;
      // simRun->btree->GetListOfBranches()->ls();
      // printf(" simRun ev %i  channel %i \n", iev, ich);

      // rawEvent->time = EventInfo->TriggerTimeTag;
      double eff = effGeoSim(ich) * SiPMQE128Ham;
      double nsmean = double(nPhotonsEvent) * eff * singletFrac;
      double ntmean = double(nPhotonsEvent) * eff - nsmean;
      int nsinglet = ran->Poisson(nsmean);
      int ntriplet = ran->Poisson(ntmean);

      ncount[ich] += nsinglet + ntriplet;

      if (iev / 100 * 100 == iev)
        if (ich < 12 && ich > 8)
          printf("event %i nphotons %i ich %i eff %E singlet %i triplet %i tot  %i \n", iev, nPhotonsEvent, ich, eff, nsinglet, ntriplet, nsinglet + ntriplet);

      // singlet times
      for (int it = 0; it < nsinglet; ++it)
      {
        double time = triggerStart + ran->Exp(tSinglet0);
        hPhoton[ich]->Fill(time, gain);
        convolve(hConvolve[ich], time);
        TH1D *hist = hConvolve[ich];
        // printf("event %i chan %i  max value %E\n", iev, ich, hist->GetBinContent(hist->GetMaximumBin()));
        hTime->Fill(time);
        // make a TDetHit for photon
        TDetHit hit;
        hit.startTime = (UInt_t)hTime->FindBin(time);
        hit.qpeak = gain;
        det->hits.push_back(hit);
      }
      // triplet times
      for (int it = 0; it < ntriplet; ++it)
      {
        double time = triggerStart + ran->Exp(tTriplet0);
        hPhoton[ich]->Fill(time, gain);
        convolve(hConvolve[ich], time);
        hTime->Fill(time);
        TDetHit hit;
        hit.startTime = (UInt_t)hTime->FindBin(time);
        hit.qpeak = gain;
        det->hits.push_back(hit);
        // printf(" \t\t after triplets %iev %ch %lu \n", iev, ich, det->hits.size());
      }
      // add baseline and noise
      for (int ibin = 1; ibin <= hSignal[ich]->GetNbinsX(); ++ibin)
      {
        double binNoise = ran->Gaus(0.0, sigmaNoise);
        hNoise->Fill(binNoise);
        hSignal[ich]->SetBinContent(ibin, baseline + binNoise + hConvolve[ich]->GetBinContent(ibin));
      }

      // file wave for this channel
      if (rawRun)
      {
        wave.clear();
        if (invert) // trigger sipm and ADC
        {
          for (int ibin = 1; ibin <= hSignal[ich]->GetNbinsX(); ++ibin)
          {
            uint16_t adc = maxAdc - hSignal[ich]->GetBinContent(ibin);
            wave.push_back(adc);
          }
        }
        else
        {
          for (int ibin = 1; ibin <= hSignal[ich]->GetNbinsX(); ++ibin)
            wave.push_back((uint16_t)hSignal[ich]->GetBinContent(ibin));
          // for (int ibin = 1; ibin <= hSignal[ich]->GetNbinsX(); ++ibin)
          //  printf("iev %i ich %i ibin %i short %u float %f \n", iev, ich, ibin, (uint16_t)hSignal[ich]->GetBinContent(ibin), hSignal[ich]->GetBinContent(ibin));
        }
        rawEvent->rdigi = wave;
        // printf(" .... ich %i wave size %lu \n", ich, wave.size());
        double qsum = 0;
        for (int ibin = 1; ibin <= hSignal[ich]->GetNbinsX(); ++ibin)
          qsum += (hSignal[ich]->GetBinContent(ibin) - baseline) / gain * landauMax;

        double psum = 0;
        for (int ibin = 1; ibin <= hPhoton[ich]->GetNbinsX(); ++ibin)
          psum += hPhoton[ich]->GetBinContent(ibin) / gain;

        ntTrigCh->Fill(iev, ich, qsum, psum, hPhoton[ich]->GetEntries(), eventOrigin.R(), eventOrigin.Theta() * 360. / TMath::TwoPi(), localPhi, eventOrigin.X(), eventOrigin.Y(), eventOrigin.Z());
      }
    } // end channel loop
    ntTrig->Fill(iev, hPhoton[9]->GetEntries(), hPhoton[10]->GetEntries(), hPhoton[11]->GetEntries(), eventOrigin.R(), eventOrigin.Theta() * 360. / TMath::TwoPi(), localPhi, eventOrigin.X(), eventOrigin.Y(), eventOrigin.Z());

    /* event histograms */
    TString histName;
    fout->cd();
    if (histDir->GetList()->GetEntries() < 100)
      for (int ih = 0; ih < NCHAN; ++ih)
      {
        if (hPhoton[ih]->GetEntries() < 1)
          continue;
        histDir->cd();
        histName.Form("hPhotonCh%iEv%i", ih, iev);
        TH1D *hPhotonEvent = (TH1D *)hPhoton[ih]->Clone(histName);
        hPhotonEvent->SetTitle(histName);
        //
        histName.Form("hConvolCh%iEv%i", ih, iev);
        TH1D *hConvolveEvent = (TH1D *)hConvolve[ih]->Clone(histName);
        hConvolveEvent->SetTitle(histName);
        //
        histName.Form("hSignalCh%iEv%i", ih, iev);
        TH1D *hSignalEvent = (TH1D *)hSignal[ih]->Clone(histName);
        hSignalEvent->SetTitle(histName);
      }
    if (rawRun)
      rawRun->fill();
    simRun->fill();
  } // end of event loop

  for (int ich = 0; ich < NCHAN; ++ich)
    hCount->SetBinContent(ich + 1, ncount[ich]);
  // summary
  printf("generated %i \n", ngen);
  for (int ih = 1; ih < NCHAN; ++ih)
  {
    printf(" chan %i photons %i\n", ih, (int)hCount->GetBinContent(ih));
  }
  // fout->ls();
  printf("end of btbgen  ngen %i %s exit\n", ngen, fout->GetName());
  simRun->print();
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
  printf("... %s ngen %i file %s exit\n", argv[0], ngen, fout->GetName());
  // fout->ls();
  fout->Write();
  fout->Close();
  exit(0);
}
