// simple sim of btb
// April 28 2025
#include <iostream>
#include <fstream>
#include <numeric>
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TGraph.h"
#include "Math/Vector3D.h"
#include "modelFitGamma.hh"
#include "TBRawRun.hxx"
#include "TBSimRun.hxx"
#include "triggerPeakFit.hh"
#include "TMinuit.h"
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
TDirectory *scanDir;

bool writeRawData = true;
bool useMap = false;
int reportInterval = 100;

modelFit *models[NCHAN];
TNtuple *ntOrigin;
TNtuple *ntTrigCh;
TNtuple *ntTDiff;
TNtuple *ntTrig;
TNtuple *ntTern;
TNtuple *ntMean;
TNtuple *ntFit;
TNtuple *ntScan;
/* geant maps */
TH3D *originPDF;     // pdf of event origins
TH3D *fluxMapChan9;  // geo efficiency values
TH3D *fluxMapChan10; ///
TH3D *fluxMapChan11; ////
TH1D *hRadiusMap;
TH1D *hRhoMap;
TH1D *hPhiMap;
TH1D *hZMap;
TH3D *hRhoPhiZMap;
TH1D *hEffGeo;

TH1D *hPhoton[NCHAN];
TH1D *hConvolve[NCHAN];
TH1D *hSignal[NCHAN];
TH1D *hSignalSum[NCHAN];
Long64_t totalPhotons;
Long64_t ncount[NCHAN];
TH1D *hEventPass;
TH1D *hCount;
TH1D *hResponse;
TH1D *hTime;
TH1D *hTrigDiffTime;
TH2D *hTriangle;
uint16_t maxAdc = pow(2, 14);
double gain;
double sigmaNoise;
double landauMax = 0.018063;
double nominalGeo;
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
double meanFreePath = 1.53; // from table in cm3frmom rtabtable in cm3frmom rtabtable in cm
double totalEventEffiency;
double triggerTimes[3];
double cosMin = 0.851 / sqrt(pow(0.4, 2) + pow(0.851, 2));
double maxTriggerTimeDiffernce = 30.;
unsigned timeOffset = 13; // changed from 17 may 13, 2024

ROOT::Math::XYZVector eventOrigin(0, 0, 0);
ROOT::Math::XYZVector positionSipm9(1.052, -0.608, -0.851);
ROOT::Math::XYZVector positionSipm10(-1.052, -0.608, -0.851);
ROOT::Math::XYZVector positionSipm11(0.000, 1.216, -0.851);

TMinuit *gMinuit;
Double_t arglist[1];
int ierflg;
TGraph *gplot1 = NULL;
bool show = false;
bool doMinuit = false;

enum
{
  NPAR = 4 // interaction point + total photons
};

enum
{
  MAXSCANPLOTS = 10 // interaction point + total photons
};

/*
efficiencies  PMTQE175 = 0.38;
static double QEff128(double ppm, double dist)
*/
double eff[NCHAN];
int triggerStart = theBinWidth * 730; // 730; sipm rise time convert to ns
double speMPV = double(triggerStart);
double speSigma = 20.; // ns from single PI data fit
TF1 *speLandau;

// get geant4 map in dir bobj
bool getMap()
{
  bool rc = false;
  originPDF = NULL;
  fluxMapChan9 = NULL;
  fluxMapChan10 = NULL;
  fluxMapChan11 = NULL;
  TString mapFileName = TString(getenv("BOBJ")) + TString("/geantSim-2025-07-16-15-22-0.root");
  TFile *fmap = new TFile(mapFileName, "readonly");
  fmap->GetObject("OriginMap", originPDF);
  fmap->GetObject("FluxMapChan9", fluxMapChan9);
  fmap->GetObject("FluxMapChan10", fluxMapChan10);
  fmap->GetObject("FluxMapChan11", fluxMapChan11);
  if (originPDF && fluxMapChan9 && fluxMapChan10 && fluxMapChan11)
  {
    printf("getMap found histograms %s %s %s %s \n", originPDF->GetName(), fluxMapChan9->GetName(), fluxMapChan10->GetName(), fluxMapChan11->GetName());
    rc = true;
  }

  return rc;
}

// trigger condition return maximum time betweed trigger photons
double eventTrigger()
{
  // printf("line104 %0.f %0.f %0.f \n", hPhoton[9]->GetEntries(), hPhoton[10]->GetEntries(), hPhoton[11]->GetEntries());
  // number of samples is 7500 each bin is 2 ns
  double tdiff = double(2 * 7500);
  // all must have at least 1 photon
  if (hPhoton[9]->GetEntries() < 1)
    return tdiff;
  if (hPhoton[10]->GetEntries() < 1)
    return tdiff;
  if (hPhoton[11]->GetEntries() < 1)
    return tdiff;

  // min time between 9,10
  double tdiff910 = double(2 * 7500);
  for (int ibin9 = 1; ibin9 < hPhoton[9]->GetNbinsX(); ++ibin9)
  {
    if (hPhoton[9]->GetBinContent(ibin9) == 0)
      continue;
    double time9 = hPhoton[9]->GetBinCenter(ibin9);
    for (int ibin10 = 1; ibin10 < hPhoton[10]->GetNbinsX(); ++ibin10)
    {
      if (hPhoton[10]->GetBinContent(ibin10) == 0)
        continue;
      double time10 = hPhoton[10]->GetBinCenter(ibin10);
      if (abs(time9 - time10) < tdiff910)
        tdiff910 = abs(time9 - time10);
    }
  }

  // min time between 9,11
  double tdiff911 = double(2 * 7500);
  for (int ibin9 = 1; ibin9 < hPhoton[9]->GetNbinsX(); ++ibin9)
  {
    if (hPhoton[9]->GetBinContent(ibin9) == 0)
      continue;
    double time9 = hPhoton[9]->GetBinCenter(ibin9);

    for (int ibin11 = 1; ibin11 < hPhoton[11]->GetNbinsX(); ++ibin11)
    {
      if (hPhoton[11]->GetBinContent(ibin11) == 0)
        continue;
      double time11 = hPhoton[11]->GetBinCenter(ibin11);
      if (abs(time9 - time11) < tdiff911)
        tdiff911 = abs(time9 - time11);
    }
  }

  // min time between 10,11
  double tdiff1011 = double(2 * 7500);
  for (int ibin10 = 1; ibin10 < hPhoton[10]->GetNbinsX(); ++ibin10)
  {
    if (hPhoton[10]->GetBinContent(ibin10) == 0)
      continue;
    double time10 = hPhoton[10]->GetBinCenter(ibin10);
    for (int ibin11 = 1; ibin11 < hPhoton[11]->GetNbinsX(); ++ibin11)
    {
      if (hPhoton[11]->GetBinContent(ibin11) == 0)
        continue;
      double time11 = hPhoton[11]->GetBinCenter(ibin11);
      if (abs(time10 - time11) < tdiff1011)
        tdiff1011 = abs(time10 - time11);
    }
  }

  // return the maximum time
  tdiff = TMath::Max(tdiff910, tdiff911);
  tdiff = TMath::Max(tdiff, tdiff1011);

  ntTDiff->Fill(tdiff910, tdiff911, tdiff1011, tdiff);

  // printf("line169 tdiff 9 10 %f tdiff 9 11 %f tdiff 10 %f tdiff %f \n", tdiff910, tdiff911, tdiff1011, tdiff);

  return tdiff;
}

//// https://mathworld.wolfram.com/TernaryDiagram.html
void makeTernary(double a, double b, double c, double &x, double &y)
{
  double s = a + b + c;
  x = 0.5 * (a + 2. * b) / s;
  y = sqrt(3.) / 2. * a / s;
}

TGraph *myScan(int thePar, double xlow, double xhigh)
{
  int maxPoints = 100;
  std::vector<double> xval;
  std::vector<double> yval;

  // get min parameters
  double fitVal[NPAR];
  double fitErr[NPAR];
  for (int ipar = 0; ipar < NPAR; ++ipar)
  {
    gMinuit->GetParameter(ipar, fitVal[ipar], fitErr[ipar]);
    // printf("line141 ipar %i par %f err %f \n", ipar, fitVal[ipar], fitErr[ipar]);
  }

  double *fGin;
  double nLL;
  for (int i = 0; i < maxPoints; ++i)
  {
    double x = xlow + double(i) * (xhigh - xlow) / double(maxPoints);
    fitVal[thePar] = x;
    gMinuit->Eval(thePar, fGin, nLL, &fitVal[0], 4);
    xval.push_back(x);
    if (nLL >= 100)
      nLL = 0.0; // for plotting
    yval.push_back(nLL);
    // printf("line53 i %i par nph %f r %f theta %f  phi %f \n", i, fitVal[0], fitVal[1], fitVal[2], fitVal[3]);
    ntScan->Fill(nLL, peakMeanQsum[0], peakMeanQsum[1], peakMeanQsum[2], peakFitQsum[0], peakFitQsum[1], peakFitQsum[2], fitVal[1], fitVal[2], fitVal[3]);
    // printf("mySCAN %i %f %f \n", i, x, nLL);
  }
  // make and return graph
  return new TGraph(maxPoints, &xval[0], &yval[0]);
}

// spe landau shape
static double myLandau(Double_t *xx, Double_t *par)
{
  double x = xx[0];
  return par[2] * TMath::Landau(x, par[0], par[1], true); // normalized
}

void convolve(TH1D *hist, double time) // time is when photon arrives
{
  int offsetBin = 728; // max bin of hResponse read offf of Response histogram
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

double effGeoSim(int ichan) // uses PositionVector3D eventOrigin;
{
  bool isTrig = ichan == 9 || ichan == 10 || ichan == 11;
  if (!isTrig)
    return effGeoFunc(ichan);

  int ilevel = level(ichan);
  double e = 1.0;
  if (ilevel != 0)
    return e;

  e = nominalGeo;
  if (eventOrigin.R() == 0.)
    return e;

  /* get from map */
  if (useMap)
  {
    int xbin = originPDF->GetXaxis()->FindBin(eventOrigin.R());
    int ybin = originPDF->GetYaxis()->FindBin(cos(eventOrigin.Theta()));
    int zbin = originPDF->GetZaxis()->FindBin(eventOrigin.Phi());
    // int globalBin = originPDF->GetBin(eventOrigin.R(), eventOrigin.Theta(), eventOrigin.Phi());
    if (ichan == 9)
      e = fluxMapChan9->GetBinContent(xbin, ybin, zbin);
    else if (ichan == 10)
      e = fluxMapChan10->GetBinContent(xbin, ybin, zbin);
    else if (ichan == 11)
      e = fluxMapChan11->GetBinContent(xbin, ybin, zbin);
    // printf("Bins (%i,%i,%i) at %f %f %f  effGeo %E \n", xbin, ybin, zbin, eventOrigin.R(), eventOrigin.Theta(), eventOrigin.Phi(), e);
  }
  else
  {
    /* this was old method */
    /* det positions Georgia May 2025 */
    double trigRadius = positionSipm9.R();
    // double trigRadius = 1.;
    //   convert to radians
    double trigTheta = positionSipm9.Theta(); // 11,10,9
    // phi are different
    double trigPhi[3];
    trigPhi[0] = positionSipm9.Phi();  // 9
    trigPhi[1] = positionSipm10.Phi(); // 9
    trigPhi[2] = positionSipm11.Phi(); // 9

    ROOT::Math::XYZVector rSipm = getXYZVector(trigRadius, trigTheta, trigPhi[ichan - 9]);
    ROOT::Math::XYZVector relative = rSipm - eventOrigin;

    double distance2 = relative.Mag2();
    // Area of SiPMs is 6.0mm x 6.0mm
    //    Channels 6, 7, and 8 are at 11.6 cm
    //    from the source Channels 3, 4, and 5 are at 23.2 cm
    //    from the source Channels 0, 1, and 2 are at 34.8 cm from the source Channel 12 is at 36 cm from the source.
    double aPmt = TMath::Pi() / 4.0 * pow(6.4, 2); // R11410-20  Effective area : 64 mm dia
    double a = pow(0.6, 2.);                       // SIPM area

    // correct solid angle
    ROOT::Math::XYZVector runit = rSipm.Unit();
    ROOT::Math::XYZVector ounit = relative.Unit();
    double cos = runit.Dot(ounit);
    // printf(" chan %i cos %f \n", ichan, cos);
    if (cos < 0.)
      cos = 0; // origin behind sipm
    a *= cos;

    e = a / distance2 / (4.0 * TMath::Pi());
    // shift phi for printing.
    double localPhi = relative.Phi() * 360. / TMath::TwoPi();
    if (localPhi < 0)
      localPhi += 360.;
    if (0)
      printf("effGeoSim ichan  %i level %i  cos %f dist %f R,Theta,Phi (%f,%f,%f) area %f geo eff %f \n", ichan, ilevel, cos, sqrt(distance2), relative.R(), relative.Theta() * 360. / TMath::TwoPi(), localPhi, a, e);
  }
  return e;
}

void btb(int ngen = 10000000)
{
  printf(" btb sim NOMAP generate ngen =  %i LY %.1f photons/kev * 60 = %.1f nominalGain %f nominalTrigGain %f \n", ngen, LY, numPhotons, nominalGain, nominalTrigGain);

  printf(" trigger sipm positions :  \n");
  printf(" \t sipm 9 : R %f Theta %f phi %f  \n", positionSipm9.R(), positionSipm9.Theta(), positionSipm9.Phi());
  printf(" \t sipm 10 : R %f Theta %f phi %f \n", positionSipm10.R(), positionSipm10.Theta(), positionSipm10.Phi());
  printf(" \t sipm 11 : R %f Theta %f phi %f  \n", positionSipm11.R(), positionSipm11.Theta(), positionSipm11.Phi());

  if (useMap)
  {
    if (!getMap())
    {
      printf("no geant4 map\n");
      exit(0);
    }
  }

  /* nominal geo at detector origin */
  double trigRadius = 1.486;
  nominalGeo = pow(0.6, 2.) / pow(trigRadius, 2.) / (4.0 * TMath::Pi());

  for (int ichan = 0; ichan < NCHAN; ++ichan)
    models[ichan] = new modelFit(4, ichan, thePPM);

  /* channel efficiences */
  for (int i = 0; i < NCHAN - 1; ++i)
  {
    double effGeoSimi = effGeoSim(i);
    printf("chan %i nominal effGeoSim %E \n", i, effGeoSimi);
    eff[i] = effGeoSimi * SiPMQE128Ham * fillFactor;
    // double dist = distanceLevel[level(i)];
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
  scanDir = fout->mkdir("scanDir");

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
  hEffGeo = new TH1D("EffGeo", " geometric efficiency / nominal ", 150, 0, 1.5);
  /* define ntuples amd histograms here */
  hRadiusMap = new TH1D("RadiusMap", "event radius [cm] ", 100, 0., 10.);
  hRhoMap = new TH1D("RhoMap", "event cylindrical rho [cm] ", 100, 0., 4.);
  hZMap = new TH1D("ZMap", "event cylindrical Z [cm] ", 100, 0., 10.);
  hPhiMap = new TH1D("PhiMap", "event phi", 100, -TMath::Pi(), TMath::Pi());
  hRhoPhiZMap = new TH3D("RhoZPhiMap", "cylindrical rho phi z  map ", 100, 0., 2., 100, -TMath::Pi(), TMath::Pi(), 100, 0., 4.);
  hEventPass = new TH1D("hEventPass", "event pass", 3, 0, 3);
  hTriangle = new TH2D("Triangle", "ytern vs xtern", 100, 0., 1., 100, 0., 1.);
  ntOrigin = new TNtuple("ntOrigin", " event origin ", "ev:r:cos:theta:phi:x:y:z");
  ntTrigCh = new TNtuple("ntTrigCh", " trigger info by channel ", "ev:ch:qsum:psum:nph:r:theta:phi:x:y:z");
  ntTDiff = new TNtuple("ntTDiff", "trig time differences", "tdiff910:tdiff911:tdiff1011:tdiff");
  ntTrig = new TNtuple("ntTrig", " trigger info by event  ", "ev:nph9:nph10:nph11:r:theta:phi:x:y:z");
  ntTern = new TNtuple("ntTern", " trigger sipm", "eventR:eventCos:eventPhi:qsum9:qsum10:qsum11:mean9:mean10:mean11:xmean:ymean:xq:yq");
  ntMean = new TNtuple("ntMean", "trigger means ", "ev:numPhotons:eventR:eventCos:eventPhi:qsum9:qsum10:qsum11:mean9:mean10:mean11:xternq:yternq");
  ntFit = new TNtuple("ntFit", "trigger peak fit ", "ev:numPhotons:eventR:eventCos:eventPhi:qsum9:qsum10:qsum11:fitR:fitCos:fitPhi:errR:errTheta:ierr");
  ntScan = new TNtuple("ntScan", "scan", "nll:mean9:mean10:mean11:qsum9:qsum10:qsum11:r:theta:phi");
  hCount = new TH1D("Count", "hit count", 13, 0, 13);
  hTime = new TH1D("Time", "photon time ", 7500, 0, 2 * 7500);
  hTrigDiffTime = new TH1D("TrigDiffTime", " time difference ", 7500, 0, 2 * 7500);
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

  sigmaNoise = gain * noiseToSignal;
  TH1D *hNoise = new TH1D("Noise", "Noise", 200, -10 * sigmaNoise, 10 * sigmaNoise);
  TH1D *hPhotonAll = new TH1D("PhotonAll", "all photons ", 200, 0.5 * numPhotons, 1.5 * numPhotons);
  TH1D *hPhotonSum = new TH1D("PhotonSum", "photon sum", 200, 0., 40.);
  TH1D *hPhotonSumCut = new TH1D("PhotonSumCut", "photon sum cut", 200, 0., 40.);

  // make individual light curves
  /*
  for (int ih = 0; ih < NCHAN; ++ih)
    models[ih] = new modelFit(MODELALL, ih, 0);
  */

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

    //
    hSignalSum[ih] = new TH1D(Form("SignalSum%i", ih), Form("SignalSum%i-level%i", ih, level(ih)), totalBins, 0, totalBins * (theBinWidth));
    hSignalSum[ih]->GetXaxis()->SetTitle("time [ns]");
    hSignalSum[ih]->GetYaxis()->SetTitle("photons/2ns");
  }
  /* end of define ntuples amd histograms here */

  // print info for channel
  printf("\n******* efficienies ****\n ");
  for (int ich = 8; ich >= 0; --ich)
  {
    double eff = effGeoFunc(ich);
    int ilevel = level(ich);
    printf(" chan %i level %i distance %f eff %E total eff %E\n", ich, ilevel, distanceLevel[ilevel], eff, eff * SiPMQE128Ham * fillFactor);
  }

  // trigger sipms
  double effTrigger = effGeoSim(9);
  for (int ich = 9; ich < 12; ++ich)
  {
    int ilevel = level(ich);
    printf(" chan %i level %i origin(%f,%f,%f) distance %f eff %E total eff %E\n", ich, ilevel,
           eventOrigin.X(), eventOrigin.Y(), eventOrigin.Z(), distanceLevel[ilevel], effTrigger, effTrigger * SiPMQE128Ham * fillFactor);
  }

  printf("\t\t nominal yield nphotons %.0f  3 SIPM sum %f\n,", numPhotons, 3. * numPhotons * effTrigger * SiPMQE128Ham * fillFactor);
  // print info for pmt
  double effPmt = effGeoFunc(12);
  int ilevel = level(12);
  printf(" chan %i level %i distance %f eff %E \n", 12, ilevel, distanceLevel[ilevel], effPmt);
  printf("***********\n\n\n ");

  /* static double tTriplet0 = 1600.0; ns
    static double tSinglet0 = 7.0; ns
  */

  /*******
   *            loop over events
   *******/
  totalPhotons = 0;
  int nsinglet;
  int ntriplet;
  int nTrigger = 0;
  for (int iev = 0; iev < ngen; ++iev) // start of event loop
  {

    hEventPass->SetBinContent(1, hEventPass->GetBinContent(1) + 1);

    // zero trigger times array
    for (int i = 0; i < 3; ++i)
      triggerTimes[i] = 0.;

    int nPhotonsEvent = (int)ran->Gaus(numPhotons, sqrt(numPhotons));
    totalPhotons += nPhotonsEvent;
    hPhotonAll->Fill(nPhotonsEvent);

    /********  generate gamma position  *********/
    // double gammaCosTheta = 2. * ran->Rndm() - 1.; // cos flat from 1 to -1
    /* cut out the events blocked by the source holder
    if (gammaCosTheta < cosMin)
      continue;
    */

    // generate r,cosTheta from geant4 map;

    double gammaR = 0;
    double gammaCosTheta = -1;
    double gammaPhi = 0;
    /********  generate gamma position  from geant4 map  *********/
    if (useMap)
    {
      originPDF->GetRandom3(gammaR, gammaCosTheta, gammaPhi, ran);
    }
    /* generate from Rndm */
    else
    {
      gammaR = abs(ran->Exp(meanFreePath));
      // gammaCosTheta = 2 * ran->Rndm() - 1.;
      gammaCosTheta = ran->Rndm();                      // use only positive z
      gammaPhi = (2. * ran->Rndm() - 1.) * TMath::Pi(); // -pi to pi
    }
    eventOrigin = getXYZVector(gammaR, acos(gammaCosTheta), gammaPhi);
    double localPhi = eventOrigin.Phi() * 360. / TMath::TwoPi();
    if (localPhi < 0)
      localPhi += 360.;
    hZMap->Fill(eventOrigin.Z());
    hRadiusMap->Fill(eventOrigin.R());
    hRhoMap->Fill(eventOrigin.Rho());
    hRhoPhiZMap->Fill(eventOrigin.Rho(), eventOrigin.Phi(), eventOrigin.Z());
    hPhiMap->Fill(eventOrigin.Phi());

    // cut -Z (up going in btb) events
    if (gammaCosTheta < 0)
    {
      if (show)
        printf(" skip event %i %f \n", iev, gammaCosTheta);
      continue;
    }

    hEventPass->SetBinContent(2, hEventPass->GetBinContent(2) + 1);
    // eventOrigin.SetZ(abs(eventOrigin.Z()));

    ntOrigin->Fill(iev, eventOrigin.R(), cos(eventOrigin.Theta()), eventOrigin.Theta() * 360. / TMath::TwoPi(), localPhi, eventOrigin.X(), eventOrigin.Y(), eventOrigin.Z());
    // printf(" event origin x %f y %f z %f \n", eventOrigin.X(), eventOrigin.Y(), eventOrigin.Z());

    if (iev / reportInterval * reportInterval == iev)
      printf("... event %i total photon %0.f (r,cosTheta,phi) = (%f, %f, %f) (r,theta,Phi) = (%f , %f ,%f ) \n", iev, double(totalPhotons), gammaR, gammaCosTheta, gammaPhi, eventOrigin.R(), eventOrigin.Theta() * 360. / TMath::TwoPi(), localPhi);
    // get event position

    // loop over channels
    for (int ich = 0; ich < NCHAN; ++ich)
    {
      // set the nominal gain from file modelFitGamma.hh
      gain = nominalGain;
      double timeShift = 0;
      if (ich > 8 && ich < 12)
      {
        gain = nominalTrigGain;
        timeShift = timeOffset; // trig amp delay
      }
      //
      sigmaNoise = gain * noiseToSignal;
      bool invert = ich > 8; // invert trigger 9,10,11 and PMT
      // histogram reset
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
      det->nspe = nsinglet + ntriplet;
      // simRun->btree->GetListOfBranches()->ls();
      // printf(" simRun ev %i  channel %i \n", iev, ich);

      // rawEvent->time = EventInfo->TriggerTimeTag;
      bool isTrig = false;
      isTrig = ich == 9 || ich == 10 || ich == 11;
      double effGeoSimi = effGeoSim(ich);
      // printf("xxxxv event %i chan %i (r,cosTheta,phi) (%f,%f,%f) effGeo %E  \n", iev, ich, eventOrigin.R(), eventOrigin.Theta() * 360. / TMath::TwoPi(), localPhi, effGeoSimi);
      if (isTrig)
        hEffGeo->Fill(effGeoSimi / nominalGeo);

      double eff = effGeoSimi * SiPMQE128Ham * fillFactor;
      double nsmean = double(nPhotonsEvent) * eff * singletFrac;
      double ntmean = double(nPhotonsEvent) * eff - nsmean;
      nsinglet = ran->Poisson(nsmean);
      ntriplet = ran->Poisson(ntmean);

      ncount[ich] += nsinglet + ntriplet;

      // singlet times
      for (int it = 0; it < nsinglet; ++it)
      {
        double time = timeShift + triggerStart + ran->Exp(tSinglet0);
        hPhoton[ich]->Fill(time, gain);
        convolve(hConvolve[ich], time);
        TH1D *hist = hConvolve[ich];
        // printf("event %i chan %i  max value %E\n", iev, ich, hist->GetBinContent(hist->GetMaximumBin()));
        hTime->Fill(time);
        // make a TDetHit for photon
        TDetHit hit;
        hit.startTime = double(hTime->FindBin(time)); // convert to samples
        // printf("line579 time %f %f bin %i  \n", time, hit.startTime, hPhoton[ich]->FindBin(time));
        hit.qpeak = gain;
        det->hits.push_back(hit);
      }
      // triplet times
      for (int it = 0; it < ntriplet; ++it)
      {
        double time = timeShift + triggerStart + ran->Exp(tTriplet0);
        hPhoton[ich]->Fill(time, gain);
        convolve(hConvolve[ich], time);
        hTime->Fill(time);
        TDetHit hit;
        hit.startTime = double(hTime->FindBin(time)); // convert to samples
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
        hSignalSum[ich]->SetBinContent(ibin, baseline + binNoise +
                                                 hConvolve[ich]->GetBinContent(ibin) + hSignalSum[ich]->GetBinContent(ibin));
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

        // printf("line508 ch %i nPhotonsEvent %i eff %E nPhotonsEvent*eff %.0f nhotons %i %i \n", ich, nPhotonsEvent, eff, nPhotonsEvent * eff, nsinglet + ntriplet, int(hPhoton[ich]->GetEntries()));
      }
      // if (iev / 1 * 1 == iev && ich < 12 && ich > 8)

    } // end channel loop

    // ensure the event triggers
    double maxTriggerDiff = eventTrigger();
    hTrigDiffTime->Fill(maxTriggerDiff);
    if (maxTriggerDiff > maxTriggerTimeDiffernce)
    {
      if (show)
        printf(" event %i does not trigger %f \n", iev, maxTriggerDiff);
      continue;
    }
    ++nTrigger;
    hEventPass->SetBinContent(3, hEventPass->GetBinContent(3) + 1);

    if (show)
      printf("xxx event %i nph %.0f %.0f %.0f\n", iev, hPhoton[9]->GetEntries(), hPhoton[10]->GetEntries(), hPhoton[11]->GetEntries());

    double photonSum = hPhoton[9]->GetEntries() + hPhoton[10]->GetEntries() + hPhoton[11]->GetEntries();
    hPhotonSum->Fill(photonSum);
    // cut on fitted radius
    // printf("line577 fitted radius %f\n", fitVal[1]);
    // tell triggerPeakFit the qsum normalize to the total photons
    peakFitQsum[0] = hPhoton[9]->GetEntries() / photonSum;
    peakFitQsum[1] = hPhoton[10]->GetEntries() / photonSum;
    peakFitQsum[2] = hPhoton[11]->GetEntries() / photonSum;

    ntTrig->Fill(iev, hPhoton[9]->GetEntries(), hPhoton[10]->GetEntries(), hPhoton[11]->GetEntries(), eventOrigin.R(), eventOrigin.Theta() * 360. / TMath::TwoPi(), localPhi, eventOrigin.X(), eventOrigin.Y(), eventOrigin.Z());

    hPhotonSumCut->Fill(photonSum);
    // Now ready for minimization step with MIGRAD
    // set starting param values
    double fitVal[NPAR];
    double fitErr[NPAR];

    for (unsigned i = 0; i < 3; ++i)
    {
      fitVal[i] = 0;
      fitErr[i] = 0;
    }
    double amin = 0;

    double step = 0.0001;
    if (doMinuit)
    {
      // Set starting values and step sizes for parameters
      gMinuit->mnparm(0, "yield", numPhotons, step, 0., 10. * numPhotons, ierflg);
      gMinuit->mnparm(1, "fitR", 0.1, step, 0., 2., ierflg);
      gMinuit->mnparm(2, "fitTheta", 0, step, 0., TMath::Pi(), ierflg);
      gMinuit->mnparm(3, "fitPhi", 0, step, -TMath::Pi(), TMath::Pi(), ierflg);
      // gMinuit->FixParameter(0);
      // minimize
      gMinuit->mnexcm("MIGRAD", arglist, 0, ierflg);
      if (ierflg != 0)
        printf("\t\t ***** MIGRAD event %i error code %i ******\n", iev, ierflg);

      /* get results */
      double edm, errdef;
      int nvpar, nparx, icstat;
      gMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);

      for (int ipar = 0; ipar < NPAR; ++ipar)
      {
        gMinuit->GetParameter(ipar, fitVal[ipar], fitErr[ipar]);
        if (show)
          printf("\t\t       event %i par %i fit %E err %E \n", iev, ipar, fitVal[ipar], fitErr[ipar]);
      }

      if (scanDir->GetList()->GetEntries() < MAXSCANPLOTS && ierflg == 0)
      {
        printf("scan perameter event %i\n", iev);
        gMinuit->SetGraphicsMode(kTRUE);
        gplot1 = myScan(1, 0, 2.);
        // gMinuit->mncomd("scan 1", ierflg);
        // gplot1 = (TGraph *)gMinuit->GetPlot();
        TGraph *gsave = (TGraph *)gplot1->Clone(Form("scan1Ev%i", iev));
        // gplot1->SetPoint(0, gplot1->GetPointX(0), gplot1->GetPointY(1)); // first point is NAN
        gsave->SetTitle(Form("scan of parameter 1 ev %i min nLL %.3E rmin = %.3f", iev, amin, fitVal[1]));
        // gplot1->Draw("al");
        scanDir->Add(gsave);
      }
    }
    else // calculate peakMeanQsum for this eventaOrigin
    {
      // set paramters
      gMinuit->mnparm(0, "yield", numPhotons, step, 0., 10. * numPhotons, ierflg);
      // eventOrigin
      gMinuit->mnparm(1, "fitR", eventOrigin.R(), step, 0., 2., ierflg);
      gMinuit->mnparm(2, "fitTheta", eventOrigin.Theta(), step, 0., TMath::Pi(), ierflg);
      gMinuit->mnparm(3, "fitPhi", eventOrigin.Phi(), step, -TMath::Pi(), TMath::Pi(), ierflg);
      // fill parameter array
      for (int ipar = 0; ipar < NPAR; ++ipar)
        gMinuit->GetParameter(ipar, fitVal[ipar], fitErr[ipar]);

      // calculate peakMeanQsum
      peakFit(fitVal);

      if (show)
      {
        printf("peakFit paramters: ");
        for (int ipar = 0; ipar < NPAR; ++ipar)
          printf("\t\t       event %i par %i fit %E err %E \n", iev, ipar, fitVal[ipar], fitErr[ipar]);
        printf(" peakFitQsum %f %f %f \n", peakFitQsum[0], peakFitQsum[1], peakFitQsum[2]);
      }
    }

    // printf("event %i fill fit ntuple\n", iev);
    //  fill fit ntuple
    if (nTrigger / reportInterval * reportInterval == nTrigger)
    {
      if (doMinuit)
        printf(".x.x.x report event %i nLL %f nphotons %i  singlet %i triplet %i tot  %i qsum(%f,%f,%f) mean(%f,%f,%f) fit(%f,%f,%f)\n", iev, amin, nPhotonsEvent, nsinglet, ntriplet, nsinglet + ntriplet, peakFitQsum[0], peakFitQsum[1], peakFitQsum[2], peakMeanQsum[0], peakMeanQsum[1], peakMeanQsum[2], fitVal[0], fitVal[1], fitVal[2]);
      else
        printf(".x.x.x report event %i nLL %f nphotons %i  singlet %i triplet %i tot  %i  photons (%.0f, %.0f, %.0f sum %.0f )  qsum(%f,%f,%f) mean(%f,%f,%f) \n", iev, amin, nPhotonsEvent, nsinglet, ntriplet, nsinglet + ntriplet, hPhoton[9]->GetEntries(), hPhoton[10]->GetEntries(), hPhoton[11]->GetEntries(), photonSum, peakFitQsum[0], peakFitQsum[1], peakFitQsum[2], peakMeanQsum[0], peakMeanQsum[1], peakMeanQsum[2]);
    }

    // printf("btbsim::  9 %f 10 %f 11 %f \n", peakMeanQsum[0], peakMeanQsum[1], peakMeanQsum[2]);
    double xternMean, yternMean;
    makeTernary(peakMeanQsum[0], peakMeanQsum[1], peakMeanQsum[2], xternMean, yternMean);

    double xternQ, yternQ;
    makeTernary(peakFitQsum[0], peakFitQsum[1], peakFitQsum[2], xternQ, yternQ);

    ntTern->Fill(eventOrigin.R(), cos(eventOrigin.Theta()), eventOrigin.Phi(), peakFitQsum[0], peakFitQsum[1], peakFitQsum[2], peakMeanQsum[0], peakMeanQsum[1], peakMeanQsum[2], xternMean, yternMean, xternQ, yternQ);

    hTriangle->Fill(xternQ, yternQ);

    ntMean->Fill(iev, photonSum, eventOrigin.R(), cos(eventOrigin.Theta()), eventOrigin.Phi(),
                 peakFitQsum[0], peakFitQsum[1], peakFitQsum[2], peakMeanQsum[0], peakMeanQsum[1], peakMeanQsum[2], xternQ, yternQ);
    ntFit->Fill(iev, fitVal[0], eventOrigin.R(), cos(eventOrigin.Theta()), eventOrigin.Phi(),
                peakFitQsum[0], peakFitQsum[1], peakFitQsum[2], fitVal[1], cos(fitVal[2]), fitVal[3], fitErr[1], fitErr[2], ierflg);

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
    totalEventEffiency = hPhotonSumCut->Integral() / hPhotonAll->Integral();
  } // end of event loop

  for (int ich = 0; ich < NCHAN; ++ich)
    hCount->SetBinContent(ich + 1, ncount[ich]);
  // summary
  printf("****** generated %i events.\nphoton count:\n", ngen);
  for (int ih = 1; ih < NCHAN; ++ih)
  {
    printf(" chan %i photons %i\n", ih, (int)hCount->GetBinContent(ih));
  }
  // fout->ls();
  hEventPass->Print("all");
  printf("********* end of btb with ngen %i triggers %i ********\n", ngen, nTrigger);
  // simRun->print();
}

// static TBRun *theTBRun;
int main(int argc, char *argv[])
{
  /* setup minuit fit*/
  gMinuit = new TMinuit(NPAR); // initialize TMinuit nphotons + event position vector
  gMinuit->SetFCN(fcn);
  gMinuit->SetPrintLevel(-1);

  /* define fit error */
  Int_t ierflg = 0;
  arglist[0] = 0.5; // UP for likelihood
  gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

  int ngen = 100000;

  std::cout << "  usage: btbSim <ngen> default 1000000  " << argv[0] << std::endl;
  printf("\n ");
  if (argc > 1)
  {
    ngen = atoi(argv[1]);
  }
  /* dont need this prints all the time
  double par[NPAR] = {numPhotons, 0, 0, 0};
  peakFitQsum[0] = 1.;
  peakFitQsum[1] = 2.;
  peakFitQsum[2] = 3.;
  triggerPeakFitShow = true;
  // set total photon yield
  peakFit(par);
  triggerPeakFitShow = false;
  */

  printf("***** START of btb ngen = %.0E *****\n", double(ngen));
  btb(ngen);
  printf("... %s ngen %i passed %lld total efficiency %.3f  file %s exit\n", argv[0], ngen, ntFit->GetEntries(), totalEventEffiency, fout->GetName());
  // fout->ls();
  fout->Write();
  fout->Close();
  exit(0);
}
