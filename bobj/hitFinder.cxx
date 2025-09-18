// revised Jan 27 2023 -- simplified findiang
// derivative paak finding
// class to make hits from vector data
// P. ugec et al. Pulse processing routines for neutron time-of-flight data. Nucl. Instrum. Meth., A812:134–144, 2016.
//////////////////////////////////////////////////////////
#include <sstream>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <complex> //includes std::pair, std::make_pair
#include <valarray>
//
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TROOT.h>
#include <TVirtualFFT.h>
#include <TChain.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <Rtypes.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFormula.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <algorithm> // std::sort
#include "TSpectrum.h"
#include "TRandom3.h"

#include "TBRun.hxx"
#include "hitFinder.hxx"

hitFinder::hitFinder(TFile *theFile, TBRun *brun, TString theTag, int nSamples, vector<int> vchan, vector<double> sigmaValue, vector<double> theGains)
{
  tbrun = brun;
  verboseB = false;
  verbose = false;
  isCAEN = false;
  doFFT = false;
  fFFT = NULL;
  fInverseFFT = NULL;
  // store gains
  for (unsigned ig = 0; ig < theGains.size(); ++ig)
    detGains.push_back(theGains[ig]);

  if (nSamples == CAENLENGTH)
    isCAEN = true;
  channelSigmaValue = sigmaValue;
  doPeakCorrection = true;
  TString templateDir = TString(getenv("BOBJ"));
  templateFileName = templateDir + TString("/templates-2023-05-01-15-06.root");

  // save vchan
  vChannel = vchan;
  if (verbose)
    cout << "INSTANCE OF HITFINDER using only PUP and NUP type crossings "
         << " vchan.size " << vchan.size() << endl;
  smoothing = true;
  fout = theFile;

  finderDir = (TDirectory *)fout->FindObject("finderDir");
  if (!finderDir)
  {
    printf("no finderDir\n");
  }
  splitDir = (TDirectory *)fout->FindObject("splitDir");
  if (!splitDir)
  {
    printf("no split dir\n");
  }
  sumWaveDir = (TDirectory *)fout->FindObject("sumWaverDir");
  if (!sumWaveDir)
  {
    printf("no sum wave  dir\n");
  }
  fitSingletDir = (TDirectory *)fout->FindObject("fitSingletDir");
  if (!fitSingletDir)
  {
    printf("no fit singlet dir\n");
  }

  tag = theTag;
  nsamples = nSamples;
  int nSize = nsamples + 100;
  // initialize fft

  microSec = 1.0E-3;
  timeUnit = 4.0; // ns per count
  maxPeakLength = 10000;
  thresholdStepSize = 1;

  fout->cd();
  htemplate = new TH1D("template", "template", nsamples, 0, nsamples);
  hPeakCount = new TH1D("PeakCount", " peaks by det ", vchan.size(), 0, vchan.size());
  hHitLength = new TH1I("HitLength", " hit length", 1000, 0, 1000);
  hPeakNWidth = new TH1I("PeakNWidth", "PeakNWidth", 1000, 0, 1000);
  hPeakValue = new TH1D("PeakValue", "Peak value (not trigger)", 1000, 0, 5000);
  hPeakCrossingBin = new TH1D("PeakCrossingBin", "peak Crossing Bin", 100, 0, 100);
  hPeakCrossingRatio = new TH1D("PeakCrossingRatio", "peak Crossing Ratio", 100, 0., 1.);
  hOverlap = new TH1D("Overlap", "peak overlap (samples ) ", 300, 0., 300.);
  hPeakCorrectionOffset = new TH1D("PeakCorrectionOffset", "PeakCorrectionOffset", 1100, -100., 10000.);
  if (doFFT)
  {
    if (verbose)
      cout << "line101 initialize  FFT  " << endl;
    fFFT = TVirtualFFT::FFT(1, &nSize, "R2C M K");
    fInverseFFT = TVirtualFFT::FFT(1, &nSize, "C2R M K");
    // make this one directory here
    fftDir = fout->mkdir("fftDir");
    fftDir->cd();
    for (unsigned index = 0; index < vchan.size(); ++index)
    {
      int id = vchan[index];
      TDet *deti = tbrun->getDet(id);
      chanMap.insert(std::pair<int, int>(id, index));
      hFFT.push_back(new TH1D(Form("FFTDET%i", id), Form("FFT Channel %i ", id), nsamples / 2, 0, nsamples / 2));
      hInvFFT.push_back(new TH1D(Form("InvFFTDET%i", id), Form("Inverse FFT Channel %i ", id), nsamples, 0, nsamples));
      hFFT[index]->SetDirectory(nullptr);
      hInvFFT[index]->SetDirectory(nullptr);
      hFFTFilt.push_back(new TH1D(Form("FFTFiltDET%i", id), Form("filtered FFT Channel %i ", id), nsamples / 2, 0, nsamples / 2));
      // hFFTFilt[index]->SetDirectory(nullptr);
      printf(" create  index %i vchan %i %s %s \n", index, id, hFFT[index]->GetName(), hFFT[index]->GetTitle());
    }

    htemplate = new TH1D("template", "template", nsamples, 0, nsamples);
    htemplateFFT = new TH1D("templateFFT", "templateFFT", nsamples / 2, 0, nsamples / 2);
    hWFilter = new TH1D("WFilter", "WFilter", nsamples / 2, 0, nsamples / 2);
  }

  splitDir->cd();
  for (unsigned index = 0; index < vchan.size(); ++index)
  {
    int id = vchan[index];
    TDet *deti = tbrun->getDet(id);
    chanMap.insert(std::pair<int, int>(id, index));
    hCrossingBinA.push_back(new TH1D(Form("CrossingBinA%i", id), Form("Crossing Bin chan  %i ", id), nsamples / 2, 0, nsamples / 2));
    hCrossingBinB.push_back(new TH1D(Form("CrossingBinB%i", id), Form("Crossing Bin chan  %i ", id), nsamples / 2, 0, nsamples / 2));
    hCrossingBinC.push_back(new TH1D(Form("CrossingBinC%i", id), Form("Crossing Bin chan  %i ", id), nsamples / 2, 0, nsamples / 2));
    hCrossingMaxBin.push_back(new TH1D(Form("CrossingMaxBin%i", id), Form("Crossing max bin chan  %i ", id), nsamples / 2, 0, nsamples / 2));
    hMaxBinVal.push_back(new TH1D(Form("MaxBinVal%i", id), Form(" max bin val/gain chan  %i ", id), 1000, 0., 10.));
  }

  fout->cd();
  hDeriv8 = new TH1D("Deriv8", "Deriv8", nsamples, 0, nsamples);
  finderDir->cd();
  for (unsigned index = 0; index < vchan.size(); ++index)
  {
    int id = vchan[index];
    TDet *deti = tbrun->getDet(id);
    hEvWave.push_back(new TH1D(Form("EvWave%s", deti->GetName()), Form("Wave%s", deti->GetName()), nsamples, 0, nsamples));
    hEvHitPeakWave.push_back(new TH1D(Form("EvHitPeakWave%s", deti->GetName()), Form("HitPeakWave%s", deti->GetName()), nsamples, 0, nsamples));
    hEvSmooth.push_back(new TH1D(Form("EvSmooth%s", deti->GetName()), Form("Smooth%s", deti->GetName()), nsamples, 0, nsamples));
    hEvCross.push_back(new TH1D(Form("EvCross%s", deti->GetName()), Form("Cross%s", deti->GetName()), nsamples, 0, nsamples));
    hEvPeakCross.push_back(new TH1D(Form("EvPeakCross%s", deti->GetName()), Form("PeaCross%s", deti->GetName()), nsamples, 0, nsamples));
    hEvDerWave.push_back(new TH1D(Form("EvDerWave%s", deti->GetName()), Form("DerWave%s", deti->GetName()), nsamples, 0, nsamples));
    hEvFiltWave.push_back(new TH1D(Form("EvFiltWave%s", deti->GetName()), Form("FiltWave%s", deti->GetName()), nsamples, 0, nsamples));
    hEvHitWave.push_back(new TH1D(Form("EvHitWave%s", deti->GetName()), Form("HitWave%s", deti->GetName()), nsamples, 0, nsamples));
    hDigiVal.push_back(new TH1D(Form("DigiVal%i", id), Form("digi value chan %id", id), 2000, -1000., 1000.));
    hDerivativeVal.push_back(new TH1D(Form("DerivativeVal%i", id), Form("derivative value chan %i", id), 2000, -1000., 1000.));
    hDerivativeValTime.push_back(new TH2D(Form("DerivativeValTime%i", id), Form("derivative value vs sample chan %i", id), nsamples, 0, nsamples, 2000, -1000., 1000.));
    hPeakCut.push_back(new TH1D(Form("PeakCut%i", id), Form("peak cut chan %i", id), 1000, 0., 1000.));
    hPeakCutAndTime.push_back(new TH2D(Form("PeakCutAndTime%i", id), Form("peak cut ADC vs time chan %i", id),
                                       30, 0, 7500, 100, 0, 1000));
    hEvWave[index]->SetDirectory(nullptr);
    hEvHitPeakWave[index]->SetDirectory(nullptr);
    hEvCross[index]->SetDirectory(nullptr);
    hEvSmooth[index]->SetDirectory(nullptr);
    hEvDerWave[index]->SetDirectory(nullptr);
    hEvHitWave[index]->SetDirectory(nullptr);
    hEvFiltWave[index]->SetDirectory(nullptr);
    hHitSum.push_back(new TH1D(Form("HitSum%s", deti->GetName()), Form("HitSum%s", deti->GetName()), nsamples, 0, nsamples));
    printf(" create  index %i vchan %i %s %s \n", index, id, hEvWave[index]->GetName(), hEvWave[index]->GetTitle());
  }

  hEvAllSumWave = new TH1D("EvAllSumWave", "EvAllSumWave", nsamples, 0, nsamples);
  hEvAllSumWave->SetDirectory(nullptr);

  fout->cd("sumDir");
  for (unsigned index = 0; index < vchan.size(); ++index)
  {
    int id = vchan[index];
    TDet *deti = tbrun->getDet(id);
    hUnFilteredSummedWave.push_back(new TH1D(Form("UnFilteredSummedWave%s", deti->GetName()), Form(" un filtered summed wave%s", deti->GetName()), nsamples, 0, nsamples));
    hFilteredSummedWave.push_back(new TH1D(Form("FilteredSummedWave%s", deti->GetName()), Form("filtered summed wave%s", deti->GetName()), nsamples, 0, nsamples));
  }
  fout->cd();
  ntFinder = new TNtuple("ntFinder", " hit finder ", "event:chan:nhit:startt:peakBin:lastBin:qpeak");
  ntSplit = new TNtuple("ntSplit", " split for finder ", "event:chan:cross:nsplit:bin:ratio:batr:width");
  ntPeakFix = new TNtuple("ntPeakFix", "peak fix for singlet", "detHits:idet:singlett:peakt:qpeak:qpeakFix");

  int templateChan = 8;
  gotTemplate = getTemplate(templateChan);

  cout << " created hitFinder with " << tbrun->GetName() << " nsamples =  " << nsamples << " ndet " << hEvWave.size() << " ";
  if (gotTemplate)
    cout << " totSumSPE Template " << htemplate->GetName() << endl;
  else
    cout << " SPE Template not found ! " << endl;

  // set wfilter size
  wfilter.resize(nsamples);
  for (int i = 0; i < nsamples; ++i)
    wfilter[i] = 1.;
  //
  if (gotTemplate && doFFT)
  {
    // make transorm
    templateTransform = forwardFFT(SPEdigi);
    // fill htemplateFFT start with first nonzero bin;
    printf(" ********   complex transform  size %lu ******** \n", templateTransform.size());
    // make filter
    fillWFilter(templateChan);
    for (int i = 0; i < nsamples / 2; ++i)
    {
      hWFilter->SetBinContent(i, wfilter[i]);
      // printf(" wfilter %i %f \n", i, wfilter[i]);
      htemplateFFT->SetBinContent(i, std::abs(templateTransform[i]));
    }
  }
  printf(" channel mapping \n");
  for (unsigned index = 0; index < vchan.size(); ++index)
  {
    int id = chanMap.at(vchan[index]);
    printf("index %i chan %i mapped to index  %i %s %s\n", index, vchan[index], id,
           hEvWave[id]->GetName(), hEvWave[id]->GetTitle());
  }
  printf("GAINS: \n");
  for (unsigned ichan = 0; ichan < detGains.size(); ++ichan)
    printf("chan %i gain %f ; ", ichan, detGains[ichan]);
  printf("\n");
  printf("\t HHHHHHHH INSTANCE of hitFinder verbose %i verboseB %i\n", verbose, verboseB);
}
//
void hitFinder::fillWFilter(int ichan)
{
  printf("hitFinder::fillWFilter called %i \n", gotTemplate);
  if (!gotTemplate)
    return;
  double noiseVal = channelSigmaValue[ichan];
  for (int i = 0; i < nsamples; ++i)
  {
    double val = std::abs(templateTransform[i]);
    wfilter[i] = val / (val + noiseVal);
  }
}

bool hitFinder::getTemplate(int ichan)
{
  printf(" hitFinder::getTemplate looking for  %s \n", templateFileName.Data());
  bool exists = false;
  FILE *aFile;
  aFile = fopen(templateFileName.Data(), "r");
  if (aFile)
  {
    fclose(aFile);
    exists = true;
  }
  if (!exists)
  {
    printf(" template file %s does not exist \n", templateFileName.Data());
    return false;
  }
  TH1D *hist = NULL;
  TFile *f1 = new TFile(templateFileName, "readonly");
  if (f1->IsZombie())
  {
    printf(" no  file for %s \n", templateFileName.Data());
    return false;
  }
  f1->GetObject(Form("QPEShapeChan%i", ichan), hist);
  if (!hist)
    return false;

  printf(" got template %s  from file %s \n", hist->GetName(), templateFileName.Data());

  // fill SPEdigi;
  SPEdigi.resize(nsamples);
  if (0)
  {
    int maxBin = hist->GetMaximumBin();
    for (int ibin = 0; ibin < hist->GetNbinsX(); ++ibin)
    {
      if (hist->GetBinContent(ibin) == 0)
      {
        continue;
      }
      if (ibin >= maxBin)
        SPEdigi[ibin - maxBin] = hist->GetBinContent(ibin);
      else
      {
        printf(" %i %i %i \n", ibin, -maxBin + ibin, int(SPEdigi.size()) - maxBin + ibin);
        SPEdigi[int(SPEdigi.size()) - maxBin + ibin] = hist->GetBinContent(ibin);
      }
    }
  }
  else
  {
    int fillBin = 0;
    for (int ibin = 0; ibin < hist->GetNbinsX(); ++ibin)
    {
      if (hist->GetBinContent(ibin) == 0)
        continue;
      SPEdigi[fillBin++] = hist->GetBinContent(ibin);
    }
  }

  // fill template
  for (int ibin = 0; ibin < SPEdigi.size(); ++ibin)
    htemplate->SetBinContent(ibin, SPEdigi[ibin]);

  return true;
}

void hitFinder::printPeakList(std::string message)
{
  cout << message << " peakList size " << peakList.size() << endl;
  if (peakList.size() < 1)
    return;
  for (unsigned ip = 0; ip < peakList.size(); ++ip)
  {
    unsigned peakStart = std::get<0>(peakList[ip]);
    unsigned peakEnd = std::get<1>(peakList[ip]);
    printf("\t peak %i (%i,%i) \n", ip, peakStart, peakEnd);
  }
}

void hitFinder::event(int ichan, Long64_t ievent, vector<double> inputDigi, double theDerivativeThreshold, double theHitThreshold, unsigned step)
{
  fSinglet = NULL;
  /*if (ichan == 9 && ievent == 14)
    verbose = true;
  else
    verbose = false;
    */
  /////  copy to internal class vector////////
  digi = inputDigi;
  if (digi.size() != CAENLENGTH)
  {
    printf("line343 hitFinder BAD DIGI SIZE %lu \n", digi.size());
    return;
  }

  // check validity of digi
  /*
  bool validDigi = true;
  unsigned badi = 0;
  for (unsigned idigi = 0; idigi < digi.size(); ++idigi)
  {
    if (digi[idigi] > 1.E9)
    {
      validDigi = false;
      badi = idigi;
      break;
    }
  }
  if (!validDigi)
  {
    printf("line361 hitFinder INVALID DIGI chan %i event %lld badi %u  \n", ichan, ievent, badi);
    return;
  }
  */

  bool trig = ichan == 9 || ichan == 10 || ichan == 11;
  theEvent = ievent;
  hitThreshold = theHitThreshold;
  derivativeThreshold = theDerivativeThreshold;
  diffStep = step;
  int idet = chanMap.at(ichan);
  splitCount.clear();
  for (int i = 0; i < vChannel.size(); ++i)
    splitCount.push_back(0);

  if (verbose)
    printf("line340 HHHHHH hitFinder START ievent %llu ichan %i idet %i derivative threshold %.1f digi size %lu \n", ievent, ichan, idet, derivativeThreshold, digi.size());

  double triggerTime = 0;
  double firstCharge = 0;

  for (int i = 0; i < nsamples; ++i)
  {

    hDigiVal[idet]->Fill(digi[i]);
  }
  // FFT and convolution
  if (doFFT)
  {
    if (verbose)
      printf("line357 HHHHHH hitFinder do FFT \n");
    std::vector<std::complex<double>> inputWaveTransform = forwardFFT(digi);
    for (int i = 0; i < nsamples / 2; ++i)
    {
      hFFT[idet]->SetBinContent(i, std::abs(inputWaveTransform[i]));
      hFFTFilt[idet]->SetBinContent(i, hFFTFilt[idet]->GetBinContent(i) + std::abs(inputWaveTransform[i]));
    }

    unsigned maxFrequency = inputWaveTransform.size();
    if (verbose)
      printf("line359 max frequency  %u  \n", maxFrequency);
    // apply FFT convolution here
    if (gotTemplate)
    {
      fillWFilter(ichan); // use channel noise
      for (unsigned iw = 1; iw < maxFrequency; ++iw)
      {
        // divide out the SPE shape
        inputWaveTransform[iw] = wfilter[iw] * inputWaveTransform[iw]; // templateTransform[iw];
      }
    }

    fdigi = backwardFFT(inputWaveTransform);
    for (unsigned isample = 0; isample < digi.size(); isample++)
    {
      hUnFilteredSummedWave[idet]->SetBinContent(isample + 1, digi[isample] + hUnFilteredSummedWave[idet]->GetBinContent(isample + 1));
      hFilteredSummedWave[idet]->SetBinContent(isample + 1, fdigi[isample] + hFilteredSummedWave[idet]->GetBinContent(isample + 1));
    }
  } // if doFFT
  else
    fdigi = digi;

  // use filtered waveforms
  // for (unsigned isample = 0; isample < 20; isample++)
  // printf(" wfilter ??? %i %f %f ?? %f \n", isample, wfilter[isample], digi[isample], fdigi[isample]);
  // if (gotTemplate) {
  //   digi = fdigi;
  //}
  hEvAllSumWave->Reset("ICESM");
  hEvWave[idet]->Reset("ICESM");
  hEvDerWave[idet]->Reset("ICESM");
  // fill wave for smoothing
  for (unsigned isample = 0; isample < digi.size(); isample++)
  {
    hEvSmooth[idet]->SetBinContent(isample + 1, digi[isample]);
    if (doFFT)
      hEvFiltWave[idet]->SetBinContent(isample + 1, fdigi[isample]);
    if (doFFT)
      hInvFFT[idet]->SetBinContent(isample + 1, fdigi[isample]);
    // sum all waves for this event
  }
  // smooth and fill vector
  sdigi.resize(digi.size());
  /* Based on algorithm 353QH twice presented by J. Friedman in Proc. of the 1974 CERN School of Computing, Norway, 11-24 August, 1974.
   See also Section 4.2 in J. Friedman, Data Analysis Techniques for High Energy Physics.
  hEvSmooth[idet]->Smooth(1); // one time
  for (unsigned ibin = 1; ibin < hEvSmooth[idet]->GetNbinsX(); ibin++)
    sdigi[ibin - 1] = hEvSmooth[idet]->GetBinContent(ibin);
   */

  // SGFitler smoothing
  //  smooth
  int nwindowSG = 10; // in samples
  int npoly = 3;
  sdigi = sgfilt->SavGolFilter(digi, nwindowSG, npoly);

  // use smooth wave if smoothing
  if (verbose)
    printf("line401  smoothing ? %i  digi size %lu \n", smoothing, digi.size());
  if (smoothing)
    digi = sdigi;

  // fill event wave after smoothing
  for (unsigned isample = 0; isample < digi.size(); isample++)
  {
    hEvWave[idet]->SetBinContent(isample + 1, digi[isample]);
    hEvAllSumWave->SetBinContent(isample + 1, hEvAllSumWave->GetBinContent(isample + 1) + digi[isample]);
  }

  ddigi.clear();
  // check validity of sdigi
  /*
  validDigi = true;
  badi = 0;
  for (unsigned idigi = 0; idigi < sdigi.size(); ++idigi)
  {
    if (sdigi[idigi] > 1.E9)
    {
      validDigi = false;
      badi = idigi;
      break;
    }
  }
  if (!validDigi)
  {
    printf("line466 hitFinder after smoothing INVALID DIGI chan %i event %lld badi %u  sdigi %E %E \n", ichan, ievent, badi,
           sdigi[badi - 1], sdigi[badi]);
    return;
  }
    */

  differentiate(ichan, ievent);
  for (unsigned isample = 0; isample < ddigi.size(); isample++)
  {
    hDerivativeVal[idet]->Fill(ddigi[isample]);
    hDerivativeValTime[idet]->Fill(double(isample), ddigi[isample]);
    hEvDerWave[idet]->SetBinContent(isample + 1, ddigi[isample]);
  }
  // find peaks
  // for derivativePeaks, window in time is timeUnit*windowSize (ns) . timeUnit = 2
  // min, max width in time bins for simple peaks
  Int_t windowSize = 10;
  unsigned maxWidth = 100000;
  unsigned minWidth = 10;
  findDerivativeCrossings(idet);
  // findThresholdCrossings(idet, hitThreshold);
  makePeaks(idet, digi);
  /*
     if (peakList.size() > 0)
       fitSinglet(idet, ievent);
  */
  // added back July 21 2025 removed aug 6
  // splitPeaks(idet);
  makeHits(idet, triggerTime, firstCharge);
  hPeakCount->Fill(idet, peakList.size());
  // fill hits
  if (verbose)
    cout << "line429 finished makePeaks  event " << ievent << " chan " << ichan << " det " << idet
         << "  ddigi size " << ddigi.size()
         << "  crossings size " << crossings.size()
         << "  peakList size " << peakList.size()
         << "  detHits size " << detHits.size()
         << endl;
  hdigi.clear();
  hdigi.resize(digi.size());

  // store threshold
  tbrun->detList[idet]->thresholds = derivativeThreshold;

  //  for (const auto &[key, value] : m)
  //    std::cout << '[' << key << "] = " << value << "; "
  // push hits to tbrun
  int icount = 0;
  bool triggerChannel = false;
  double startTimeCut = 800.0; // cut for singlet
  if (ichan == 9 || ichan == 10 || ichan == 11)
    triggerChannel = true;

  hEvHitPeakWave[idet]->Reset("ICESM");
  int hitNumber = 0;
  TString hitTitle;

  /*
    loop over detHits and add hits to TDet in tbrun
  */

  for (hitMapIter hitIter = detHits.begin(); hitIter != detHits.end(); ++hitIter)
  {
    TDetHit hiti = hitIter->second;
    tbrun->detList[idet]->hits.push_back(hiti);

    // fill hit digi
    for (unsigned iv = 0; iv < digi.size(); ++iv)
      if (iv >= hiti.firstBin && iv <= hiti.lastBin)
        hdigi[iv] = digi[iv];
    // if (hiti.qsum > 7000 && hiti.qsum < 10000) // FILL ONLY SINGLE PE

    // fill hit peak wave first bin is number one!
    hEvHitPeakWave[idet]->SetBinContent(hiti.peakBin + 1, hiti.qpeak);
    if (idet == 9 && verbose)
      printf("line461 size %lu hit%i idet %i time %u peakBin %i qpeak  %f \n", detHits.size(), hitNumber++, idet, hitIter->first, hiti.peakBin, hiti.qpeak);
    // make sums with cut
    if (hiti.qsum > hitThreshold)
    {
      tbrun->detList[idet]->qarea += hiti.qsum;
      tbrun->detList[idet]->qpeak += hiti.qpeak;
      if (hiti.startTime < startTimeCut)
      {
        // tbrun->detList[idet]->qPrompt += hiti.qsum;
        tbrun->detList[idet]->hitPrompt += hiti.qpeak;
      }
    }
    if (verbose)
      printf("line473 hitFinder event %lld chan %i thres %f qpeak sum %f\n", ievent, ichan, hitThreshold, tbrun->detList[idet]->qpeak);
    if (!triggerChannel)
      hPeakValue->Fill(hiti.qpeak);
    hitTitle.Form("TDetHit %i event %llu chan %i index %i ", icount++, ievent, ichan, idet);
    if (verbose)
      cout << hitTitle << endl;
    hiti.SetTitle(hitTitle);
  }
  // save the hits
  // cout << "btree entries " << tbrun->btree->GetEntries() << endl;
  // tbrun->fill();

  // save some split histograms
  for (unsigned idet = 0; idet < tbrun->detList.size(); ++idet)
    if (splitCount[idet] > 0 && splitDir->GetList()->GetEntries() < 500)
    {
      // printf("line486 plot SplitEvent %llu %i \n", theEvent, idet);
      plotEvent(splitDir, idet, theEvent);
    }

  for (unsigned isample = 0; isample < hdigi.size(); isample++)
  {
    hEvHitWave[idet]->SetBinContent(isample + 1, hdigi[isample]);
    hHitSum[idet]->SetBinContent(isample + 1, hdigi[isample] + hHitSum[idet]->GetBinContent(isample + 1));
  }

  //
  if (tbrun->detList[idet]->hits.size() > 1 && (verbose || verboseB))
  {
    TDet *tdet = tbrun->detList[idet];
    cout << "HHHH  END hitFinder::event " << theEvent << " idet= " << idet << " " << tdet->channel << " hits.size " << tdet->hits.size() << endl;

    for (unsigned ip = 0; ip < peakList.size(); ++ip)
    {
      unsigned peakStart = std::get<0>(peakList[ip]);
      unsigned peakEnd = std::get<1>(peakList[ip]);
      printf("line612  hitFinder::  event %lld det %i peak %i  peakStart %i peakEnd %i\n", theEvent, idet, ip, peakStart, peakEnd);
    }

    for (unsigned ihit = 0; ihit < tdet->hits.size(); ++ihit)
    {
      cout << " \t finder hit number  " << ihit << " peak bin " << tdet->hits[ihit].peakBin << " qpeak " << tdet->hits[ihit].qpeak << endl;
    }
  }
  if (verbose || verboseB)
    cout << "HHHH  END hitFinder::event ichan " << ichan << " event " << ievent << "   " << tbrun->detList[idet]->hits.size() << "  " << detHits.size() << endl;
}

// revised derivative Jan 27 2023 MG
void hitFinder::differentiate(int ichan, Long64_t ievent)
{
  if (verbose)
    printf("line516 hitFinder::differentiate nsamples %lu step %u\n", digi.size(), diffStep);

  // check validity of digi
  /*
  bool validDigi = true;
  unsigned badi = 0;
  for (unsigned idigi = 0; idigi < digi.size(); ++idigi)
    if (digi[idigi] > 1.E9)
    {
      validDigi = false;
      badi = idigi;
    }
  if (!validDigi)
  {
    printf("line601 hitFinder INVALID DIGI chan %i event %lld bad %u \n", ichan, ievent, badi);
  }
    */

  ddigi.clear();
  ddigi.resize(digi.size());
  Double_t sump = 0;
  Double_t summ = 0;
  unsigned nsamples = digi.size();
  ddigi[0] = 0; // first entry is zero
  for (unsigned i = 1; i < nsamples; ++i)
  {
    // sum limit
    int maxSum = diffStep;
    if (i < diffStep)
      maxSum = i;
    if (nsamples - 1 - i < diffStep)
      maxSum = nsamples - 1 - i;
    //
    sump = 0;
    for (unsigned j = 0; j < maxSum; ++j)
    {
      // f (i + 1 + j > 7480)
      //   printf("line596 hitFinder::differentiate event %lld ichan %i i %i j%i i+1+j %i digi %E  \n", ievent, ichan, i, j, i + 1 + j, digi[i + 1 + j]);
      sump += digi[i + 1 + j];
    }
    summ = 0;
    for (unsigned j = 0; j < maxSum; ++j)
    {
      // if (j > 7480)
      //   printf("line603 hitFinder::differentiate event %lld ichan %i i %i j%i i+1+j %i digi %E  \n", ievent, ichan, i, j, i - 1 - j, digi[i - 1 - j]);
      summ += digi[i - 1 - j];
    }
    ddigi[i] = sump - summ;
    // if (i > 7480)
    // printf("line608 hitFinder::differentiate event %lld ichan %i step %i bin %i maxSum %u sump %E summ %E ddigi %E \n", ievent, ichan, diffStep, i, maxSum, sump, summ, ddigi[i]);
  }
}

// revised derivative Jun22  2023 MG
vector<double> hitFinder::differentiate(int step, vector<double> pdigi)
{
  if (verbose)
    printf("line551 hitFinder::differentiate step %i size %lu \n", step, pdigi.size());
  vector<double> pddigi;
  pddigi.clear();
  pddigi.resize(pdigi.size());
  if (pdigi.size() == 0)
    return pddigi;
  Double_t sump = 0;
  Double_t summ = 0;
  unsigned nsamples = pdigi.size();
  pddigi[0] = 0; // first entry is zero
  for (unsigned i = 1; i < nsamples; ++i)
  {
    // sum limit
    int maxSum = step;
    if (i < step)
      maxSum = i;
    summ = 0;
    if (verbose)
      printf("line569 hitFinder::differentiate ind %i maxSum %i \n", i, maxSum);
    for (unsigned j = 0; j < maxSum; ++j)
      summ += pdigi[i - 1 - j];

    if (nsamples - 1 - i < step)
      maxSum = nsamples - 1 - i;
    //
    sump = 0;
    for (unsigned j = 0; j < maxSum; ++j)
      sump += pdigi[i + 1 + j];
    //

    pddigi[i] = sump - summ;
  }
  return pddigi;
}

// threshold crossings
void hitFinder::findThresholdCrossings(Int_t idet, double thresh)
{
  crossings.clear();
  crossingBin.clear();
  crossingTime.clear();
  unsigned vsize = digi.size();
  // Double_t cut = tbrun->detList[idet]->sigma * threshold;
  //  fixed cut value
  Double_t cut = hitThreshold;
  for (unsigned ibin = 0; ibin < digi.size(); ++ibin)
  {
    Double_t u = double(ibin) * timeUnit;
    if (digi[ibin] < cut && digi[ibin + 1] > cut)
    {
      crossings.push_back(PUP);
      crossingBin.push_back(ibin + 1);
      crossingTime.push_back(u);
      if (verbose)
        printf("line605 PUP det %i  bin %i %f %f  \n", idet, ibin, digi[ibin], digi[ibin + 1]);
    }
  }
  if (verbose)
    printf("line609 findTresholdCrossings det %i  crossings %lu \n", idet, crossings.size());
}
//
void hitFinder::findDerivativeCrossings(Int_t idet)
{
  unsigned step = 1;
  Double_t cut = derivativeThreshold;
  if (verbose)
    printf(" line617 findDerivativeCrossings  det = %i ddigi size %lu step %u cut %f \n", idet, ddigi.size(), step, cut);
  crossings.clear(); // crossing type
  crossingBin.clear();
  crossingTime.clear();
  unsigned vsize = ddigi.size();
  // find all crossings
  for (unsigned ibin = 0; ibin < vsize - step; ++ibin)
  {
    Double_t u = double(ibin) * timeUnit;
    Double_t vi = ddigi[ibin];
    Double_t vj = ddigi[ibin + step];
    unsigned ctype = 10;
    // crossing types
    if (vi < cut && vj > cut)
    {
      if (verbose)
        printf("line633 PUP det %i  bin %i vi %f vj %f  \n", idet, ibin, vi, vj);

      // if (ibin > 7488)
      //   printf("line698 PUP ddigi size %lu det %i  bin %i vi %f vj %f  \n", ddigi.size(), idet, ibin, vi, vj);
      //  if(idet==13&& ibin>1040&&ibin<1070)
      //    printf("line635  PUP det %i  bin %i %f %f  \n", idet, ibin, digi[ibin], digi[ibin + 1]);
      crossings.push_back(PUP);
      ctype = PUP;
      crossingBin.push_back(ibin + 1);
      crossingTime.push_back(u);
    }
    else if (vi > cut && vj < -cut)
    {
      if (verbose)
        printf("line644 UPDOWN det %i  bin %i vi %f vj %f  \n", idet, ibin, vi, vj);
      crossings.push_back(UPDOWN);
      ctype = UPDOWN;
      crossingBin.push_back(ibin + 1);
      crossingTime.push_back(u);
    }
    else if (vi > cut && vj < cut)
    {
      if (verbose)
        printf("line653  NUP det %i  bin %i vi %f vj %f  \n", idet, ibin, vi, vj);
      crossings.push_back(NUP);
      ctype = NUP;
      crossingBin.push_back(ibin + 1);
      crossingTime.push_back(u);
    }
    else if (vi < -cut && vj > cut)
    {
      if (verbose)
        printf("line662 DOWNUP det %i  bin %i vi %f vj %f  \n", idet, ibin, vi, vj);
      crossings.push_back(DOWNUP);
      ctype = DOWNUP;
      crossingBin.push_back(ibin + 1);
      crossingTime.push_back(u);
    }
    else if (vi < -cut && vj > -cut)
    {
      if (verbose)
        printf("line671 PDOWN det %i  bin %i %f %f  \n", idet, ibin, vi, vj);
      crossings.push_back(PDOWN);
      ctype = PDOWN;
      crossingBin.push_back(ibin + 1);
      crossingTime.push_back(u);
    }
    else if (vi > -cut && vj < -cut)
    {
      crossings.push_back(NDOWN);
      ctype = NDOWN;
      crossingBin.push_back(ibin + 1);
      crossingTime.push_back(u);
    }
    // if (idet==5&&ctype<10)  printf("....... %u vj %f vi %f cut %f cross type %u \n", ibin, vj, vi, cut, ctype );
    // if (idet==1&&ibin>2350&&ibin<2450)  printf("\t %u vj %f vi %f ctype %u  \n", ibin, vj, vi, ctype );
  }

  if (verbose)
    printf("line689  findDerivativeCrossings >> finished det = %i crossings found %lu  \n", idet, crossings.size());

  return;
}
// make peaks to zero of waveform from PUP crossing type
// peaks are kept in peakList
void hitFinder::makePeaks(int idet, std::vector<Double_t> v)
{
  if (verbose)
    printf("line697 hitFinder::makePeaks det %i crossings %lu \n", idet, crossings.size());
  double sigma = tbrun->detList[idet]->sigma;
  peakList.clear();
  peakKind.clear();
  double nominalGain = detGains[idet];
  hEvCross[idet]->Reset("ICESM");
  // loop over crossings using  PUP or NUP
  for (int icross = 0; icross < crossings.size(); ++icross)
  {
    // if (!(crossings[icross] == PUP || crossings[icross] == NUP))
    //   continue;
    //  find local max
    hCrossingBinA[idet]->Fill(crossingBin[icross]);
    unsigned imax = 0;

    // for (unsigned ibin = crossingBin[icross]; ibin < v.size(); ++ibin)
    double maxVal = -99999.;
    if (crossings[icross] == PUP) // case PUP
    {
      for (unsigned ibin = crossingBin[icross]; ibin < v.size(); ++ibin) //
      {
        if (v[ibin] < maxVal) // passed the peak
          break;
        imax = ibin;
        maxVal = v[ibin];
      }
    }
    /* sept 11 use only PUP */
    /*else if (crossings[icross] == NUP) // case NUP
    {                                  // NUP is other side of derivative going through zero
      maxVal = -99999.;
      for (unsigned ibin = crossingBin[icross]; ibin > 0; --ibin)
      {
        if (v[ibin] < maxVal)
          break;
        imax = ibin;
        maxVal = v[ibin];
      }
    }
    */
    hMaxBinVal[idet]->Fill(maxVal / nominalGain);
    hCrossingMaxBin[idet]->Fill(imax);

    // too small is garbage
    if (maxVal < hitThreshold)
      continue;

    // not to close to start of wave
    if (imax < 50)
      continue;

    if (verbose)
      printf("line740 hitFinder::makePeaks cross det %i icross %i cross bin %i  maxVal/nominal %f  \n", idet, icross, crossingBin[icross], maxVal / nominalGain);
    // find limits of peak
    /*
        just use a fixed window around the maximum
    */

    unsigned ilow = imax - 20;
    unsigned ihigh = imax + 50;
    if (ihigh > unsigned(v.size() - 1))
      ihigh = unsigned(v.size() - 1);
    if (ilow < 0)
      ilow = 0;
    // printf(" line776 \t\t ilow Diff! in makePeaks %lli add det %i icross %u  imax %u ilow %u \n",theEvent, idet, icross, imax, ilow );
    /*
    unsigned ilow = 0; //crossingBin[icross];
    //unsigned ilow = v.size();
    unsigned ihigh = v.size();
    // LLLLLLL low will be 10% of peak value need to correct for rise time
    // double nominalLowCut = 0.1 * maxVal;
    double nominalLowCut = 0.5 * maxVal;
    unsigned searchWindoe = 3.*8
    for (unsigned ibin = imax; ibin<v.size(); ++ibin)
    {
      if (v[ibin] < nominalLowCut) // fixed may 2 2024
        break;
      ihigh = ibin;
    }
    hCrossingBinB[idet]->Fill(ihigh);
    // look high
    // if ibin < trigEnd use fSinglet for baseline
    for (unsigned ibin = imax; ibin > 0; --ibin)
    {
      if (v[ibin] < nominalLowCut)  // fixed may 2 2024
        break;
      ilow = ibin;
    }:
    */

    hCrossingBinB[idet]->Fill(ilow);
    double nominalLowCut = 0; // no longer used

    if (verbose)
      printf("line767  hitFinder::makePeaks cross det %i lowCut %f imax %i val %f icross %u from (%u,%u) \n", idet, nominalLowCut, imax, maxVal, crossingBin[icross], ilow, ihigh);
    // if (idet == 13 && crossingBin[icross] > 1040 && crossingBin[icross] < 1070)
    //   printf("line 769 hitFinder::makePeaks cross det %i lowCut %f imax %i val %f icross %u from (%u,%u) \n", idet, nominalLowCut, imax, maxVal, crossingBin[icross], ilow, ihigh);
    //  if (ihigh - ilow > maxPeakLength)
    //    ihigh = ilow + maxPeakLength;

    // check that this hit has not already been found.
    bool found = false;
    for (unsigned ip = 0; ip < peakList.size(); ++ip)
    {
      unsigned peakStart = std::get<0>(peakList[ip]);
      unsigned peakEnd = std::get<1>(peakList[ip]);
      // if ((peakStart > ilow && peakStart < ihigh) || (peakEnd > ilow && peakEnd < ihigh))
      if (peakStart == ilow && peakEnd == ihigh)
      {
        if (verbose)
          printf("line852  hitFinder::makePeaks REMOVE OVERLAP event %lld det %i ilow,ihigh (%i %i) # %i ( %i %i )\n", theEvent, idet, ilow, ihigh, ip, peakStart, peakEnd);
        found = true;
      }
    }
    if (!found)
    {
      if (verbose)
        printf("line786  hitFinder::makePeaks add det %i  imax %i val %f icross %u from (%u,%u) \n", idet, imax, maxVal, crossingBin[icross], ilow, ihigh);
      peakList.push_back(std::make_pair(ilow, ihigh));
      peakKind.push_back(crossings[icross]);
      hCrossingBinC[idet]->Fill(ilow);
      if (imax - ilow > 30)
        printf("line830 ilow Diff! in makePeaks %lli add det %i  imax %i val %f icross %u from (%u,%u) %i\n", theEvent, idet, imax, maxVal, crossingBin[icross], ilow, ihigh, icross);
    }

    if (ilow >= 7500 && ihigh >= 7500)
      printf("line855 makePeaks ERROR!! LATE %lli add det %i  imax %i val %f icross %u from (%u,%u) %i\n", theEvent, idet, imax, maxVal, crossingBin[icross], ilow, ihigh, icross);

    /** fix peak list to remove overlaps look at next peak in list */
    for (unsigned ip = 0; ip < peakList.size() - 1; ++ip)
    {
      unsigned peakStart = std::get<0>(peakList[ip]);
      unsigned peakEnd = std::get<1>(peakList[ip]);
      unsigned peakStartNext = std::get<0>(peakList[ip + 1]);
      if (peakEnd > peakStartNext)
      {
        peakList.at(ip) = std::make_pair(peakStart, peakStartNext - 1);
        if (verboseB)
          printf("line971 peak %i start %u end %u new end %u \n", ip, std::get<0>(peakList[ip]), peakEnd, std::get<1>(peakList[ip]));
      }
    }
  }
}

void hitFinder::makeHits(int idet, Double_t &triggerTime, Double_t &firstCharge)
{
  double nominalGain = detGains[idet];
  double sigma = tbrun->detList[idet]->sigma;
  if (verbose)
    printf("line799 hitFinder::makeHits: AT event %lli det %i sigma %f peakList size %lu digi size %lu \n", theEvent, idet, sigma, peakList.size(), digi.size());
  triggerTime = 1E9;
  firstCharge = 0;
  detHits.clear();
  if (peakList.size() < 1)
    return;
  Double_t qmax = 0;
  if (isCAEN)
    qmax = 50; // about 5x CAEN noise
  // double hitThreshold = 5.0 * channelSigmaValue[idet];

  unsigned minLength = 3;
  if (peakList.size() < 1)
    return;
  for (unsigned ip = 0; ip < peakList.size(); ++ip)
  {
    unsigned klow = std::get<0>(peakList[ip]);
    unsigned khigh = std::get<1>(peakList[ip]);

    // protect against end of array MG Aug. 19 2025
    // if (idet == 12)
    //   printf("line881 hitFinder::makeHits event %lli det %i hit  %u (%u,%u) kind %i length %u \n", theEvent, idet, ip, klow, khigh, peakKind[ip], khigh - klow);
    if (verbose)
      printf("line895 hitFinder::makeHits event %lli det %i hit  %u (%u,%u) kind %i length %u \n", theEvent, idet, ip, klow, khigh, peakKind[ip], khigh - klow);
    hHitLength->Fill(khigh - klow + 1);
    if (khigh - klow + 1 < minLength)
    {
      continue;
    }
    Double_t qhit = 0;
    UInt_t peakt = 0;
    Double_t qpeak = 0;
    Double_t qsum = 0;
    for (unsigned k = klow; k < khigh; ++k)
    {
      double qdigik = digi[k];
      qsum += qdigik;
      if (qdigik > qpeak)
      {
        peakt = k;
        qpeak = qdigik;
      }
    }

    /***  second peak will be later in time then first.
           first search for local minimum
    *****/
    UInt_t minSample = 0;
    Double_t minAdc = qpeak; // must be less than this
    for (unsigned k = peakt; k < khigh; ++k)
    {
      double qdigik = digi[k];
      if (qdigik < minAdc)
      {
        minSample = k;
        minAdc = qdigik;
      }
      // break if next is greatger
      UInt_t lastk = CAENLENGTH - 1;
      double qnext = digi[min(k + 1, lastk)]; // dont go over the edge of array
      if (qnext > qdigik)
        break;
    }
    if (verboseB)
      printf("line1048 det %i  peak %u klow %u khigh %u minSample %u peakt %u qpeak %f\n ", idet, ip, klow, khigh, minSample, peakt, qpeak);
    /* if minimum is not at the end, then the peak we want is after the minimum */
    if (minSample > 0 && minSample < khigh - 1)
    {
      qsum = 0;
      qpeak = 0;
      for (unsigned k = minSample; k < khigh; ++k)
      {
        double qdigik = digi[k];
        qsum += qdigik;
        if (qdigik > qpeak)
        {
          peakt = k;
          qpeak = qdigik;
        }
      }
      if (verboseB)
        printf("line1057 NEW PEAK det %i  peak %u klow %u khigh %u minSample %u peakt %u qpeak %f\n ", idet, ip, klow, khigh, minSample, peakt, qpeak);
    }
    /*
    else
    {
      if (verboseB)
        printf("line1048 KEEP PEAK det %i  peak %u klow %u khigh %u minSample %u peakt %u qpeak %f\n ", idet, ip, klow, khigh, minSample, peakt, qpeak);
    }
    */

    // if (idet == 12)
    //   printf("line905 HitFinderMakeHits ihit %i qpeak %f time %f \n ", int(detHits.size()), qpeak, double(peakt));

    // redefine low, ihgh relative to this peak
    // this is a bug because of possible negative unisgned!
    int diff = peakt - peakWidth;
    unsigned kstart = 0;
    if (diff > 0)
      kstart = unsigned(diff);
    // BUG FIX MG August 19, 2025 cannot past end of digi!
    unsigned kend = TMath::Min(unsigned(digi.size() - 1), peakt + peakWidth);

    if (kstart >= 7500 || kend >= 7500)
      printf("line941 hitFinder::makeHit LATE %u width %u (start %u,end %u) ip %u \n", peakt, peakWidth, kstart, kend, ip);
    // cut small peaks below hitThreshold
    hPeakCut[idet]->Fill(qpeak);
    hPeakCutAndTime[idet]->Fill(peakt, qpeak);

    // if (qpeak < hitThreshold && idet == 12)
    //   printf("line911 HitFinderMakeHits ihit %i qpeak %f thresh %f \n ", int(detHits.size()), qpeak, hitThreshold);

    if (qpeak < hitThreshold)
      continue;

    TDetHit dhit;
    if (vChannel[idet] < 9)
      for (unsigned k = kstart; k < kend; ++k)
        dhit.digi.push_back(digi[k]);

    // redo the qsum here as kstart to kend;

    if (verbose)
      printf("line854 hitFinder::makeHits hit chan %i (%i,%i) size %lu \n ", vChannel[idet], klow, khigh, dhit.digi.size());

    /* bug fix */
    if (dhit.qpeak > 1.E5)
      printf("line964 hitFinder BUG very large qpeak ch %i klow %i khigh %i firstbin %i val %E \n", vChannel[idet], klow, khigh, dhit.firstBin, dhit.qpeak);

    dhit.peakBin = Int_t(peakt);
    dhit.qsum = qsum;
    dhit.qpeak = qpeak;
    dhit.firstBin = kstart;
    dhit.lastBin = kend;
    dhit.peakMaxTime = peakt;
    dhit.peakt = peakt;
    dhit.startTime = kstart;
    dhit.peakWidth = kstart - kend + 1;
    // this is N= q/qnorm and delta q = root(n)*qnorm;
    dhit.qerr = sqrt(pow(sigma * Double_t(dhit.peakWidth), 2) + qnorm * qsum);
    dhit.kind = peakKind[ip];
    // printf("line854 hitFinder::makeHits hit chan %i (%i,%i) size %lu kind %i \n ", vChannel[idet], klow, khigh, dhit.digi.size(), dhit.kind);

    // just use the biggest pulse
    if (qsum > qmax)
    {
      qmax = qsum;
      triggerTime = dhit.startTime * timeUnit * microSec;
      firstCharge = qsum;
    }
    // make  map key just startTime
    // Double_t hitTime = dhit.startTime * timeUnit * microSec;

    if (dhit.startTime >= CAENLENGTH || dhit.lastBin >= CAENLENGTH)
    {
      printf("line993 hitFinder::makeHits !!!LATE HIT TIME!!! %llu insert hit idet %i  time %f (%u,%u) peak bin %i kind %i length %u qpeak %f detHit size %lu  \n", theEvent, idet, dhit.startTime, dhit.firstBin, dhit.lastBin, dhit.peakBin, peakKind[ip], khigh - klow + 1, qpeak, detHits.size());
    }

    // if (idet == 12)
    //   printf("line951HitFinderMakeHits ihit %i time %f qpeak %f \n ", int(detHits.size()), double(dhit.startTime), dhit.qpeak);
    /*****************************************************************************
     *          duplicate peak cut
     *          ensure new hit it does not have peak bin too close to another hit
     ******************************************************************************/
    bool used = false;
    for (hitMapIter hitIter = detHits.begin(); hitIter != detHits.end(); ++hitIter)
    {
      TDetHit hiti = hitIter->second;
      // if (hiti.peakBin == dhit.peakBin)
      //  check in 5 sample range
      // if (hiti.peakBin > dhit.peakBin - 5 && hiti.peakBin < dhit.peakBin + 5)
      if (hiti.peakBin > dhit.peakBin - 2 && hiti.peakBin < dhit.peakBin + 2)
      {
        used = true;
        if (hiti.peakBin != dhit.peakBin)
          printf("line963 hitFinder::makeHit REMOVE CLOSE HIT det %i this (%i,%i,%i)  last peak (%i,%i,%i) dethit size %lu \n", idet, hiti.firstBin, hiti.peakBin, hiti.lastBin, dhit.firstBin, dhit.peakBin, dhit.lastBin, detHits.size());
      }
    }

    // if (idet == 12)
    //  printf("line968  hitFinder::makeHit  det %i last peak (%i,%i,%i) dethit size %lu used %i \n", idet, dhit.firstBin, dhit.peakBin, dhit.lastBin, detHits.size(), int(used));
    if (used)
      continue;

    /*fix peak if after singlet */
    if (fSinglet != NULL && peakt > singletPeakTime && peakt < trigEnd)
    {
      double xbin = hEvWave[idet]->GetBinLowEdge(peakt);
      double offset = fSinglet->Eval(xbin);
      if (ntPeakFix->GetEntries() < 1.E6)
        ntPeakFix->Fill(float(detHits.size()), float(idet), float(singletPeakTime), float(peakt), qpeak, qpeak - offset);
      if (verbose)
        printf("line919 list size %lu idet %i singlett %u peakt %u qpeak %f fixed %f\n", detHits.size(), idet, singletPeakTime, peakt, dhit.qpeak, dhit.qpeak - offset);
      // fix here
      dhit.qpeak = dhit.qpeak - offset;
    }
    else if (verbose && peakt > singletPeakTime && peakt < trigEnd)
    {
      printf("line925 detHits %lu  %i singlett %u peakt %u  \n", detHits.size(), idet, singletPeakTime, peakt);
    }

    // fill tFinder
    // ntFinder = new TNtuple("ntFinder", " hit finder ", "event:chan:nhit:startt:peakBin:lastBin:qpeak");
    // if (idet == 12)
    //  printf("line990  detHits fill ntFinder %lu  %i  peakt %f qpeak %f  \n", detHits.size(), idet, dhit.startTime, dhit.qpeak);
    if (ntFinder->GetEntries() < 1.E6)
      ntFinder->Fill(float(theEvent), float(idet), float(detHits.size()), float(dhit.startTime), float(dhit.peakBin), float(dhit.lastBin), dhit.qpeak);

    if (dhit.qpeak < hitThreshold && idet == 12)
      printf("line975HitFinderMakeHits ihit %i qpeak %f thresh %f \n ", int(detHits.size()), dhit.qpeak, hitThreshold);

    // cheak after peak fix
    if (dhit.qpeak < hitThreshold)
      continue;

    detHits.insert(std::pair<Double_t, TDetHit>(dhit.peakt, dhit));
    hPeakNWidth->Fill(dhit.lastBin - dhit.firstBin + 1);
    if (verbose)
    {
      printf("line941 hitFinder::makeHits %llu insert hit idet %i  time %u (%u,%u) peak bin %i kind %i length %u qpeak %f detHit size %lu  \n", theEvent, idet, dhit.peakt, dhit.firstBin, dhit.lastBin, dhit.peakBin, peakKind[ip], khigh - klow + 1, qpeak, detHits.size());
    }

    /* debugging but should never print */
    if (dhit.startTime > CAENLENGTH)
    {
      printf("line1035 hitFinder::makeHits !!!LATE HIT TIME!!! %llu insert hit idet %i  time %i (%u,%u) peak bin %i kind %i length %u qpeak %f detHit size %lu  \n", theEvent, idet, int(dhit.startTime), dhit.firstBin, dhit.lastBin, dhit.peakBin, peakKind[ip], khigh - klow + 1, qpeak, detHits.size());
    }
    /* debugging but should never print
    if (dhit.startTime == 0)
    {
      printf("line1040 hitFinder::makeHits !!!ZERO time hit!!! %llu insert hit idet %i  klow %u time %i (%u,%u) peak bin %i kind %i length %u qpeak %f detHit size %lu  \n", theEvent, idet, klow, int(dhit.startTime), dhit.firstBin, dhit.lastBin, dhit.peakBin, peakKind[ip], khigh - klow + 1, qpeak, detHits.size());
    }
    */
    /* debugging
    if (idet == 12)
    {
      printf("line1030hitFinder::makeHits event %llu insert hit idet %i  first,last (%u,%u) peak bin %i  ADC kind %i length %u qpeak %f ADC %f %f %f  detHit size %lu  \n", theEvent, idet, dhit.firstBin, dhit.lastBin, dhit.peakBin, peakKind[ip], khigh - klow + 1, qpeak, digi[dhit.peakBin - 1], digi[dhit.peakBin], digi[dhit.peakBin + 1], detHits.size());
    }
    */
  }

  /* debugging
  unsigned jjhit = 0;
  for (hitMapIter hitIter1 = detHits.begin(); hitIter1 != detHits.end(); ++hitIter1)
  {
    printf("line1078 ALL HITS hitFinder  channel %i hit %u bin %i last bin %i  qpeak  %E \n", idet, jjhit++, hitIter1->second.firstBin, hitIter1->second.lastBin, hitIter1->second.qpeak);
  }
  */

  int nhit = 0;

  /*
  unsigned jhit = 0;
  for (hitMapIter hitIter1 = detHits.begin(); hitIter1 != detHits.end(); ++hitIter1)
  {
    if (hitIter1->second.qpeak > 1.E5)
      printf("line1078 hitFinder  channel %i hit %u bin %i qpeak  %E \n", idet, jhit++, hitIter1->second.firstBin, hitIter1->second.qpeak);
  }
      */

  // this messes ip yaxis on chan13 EvWave??
  // do this differently with very short hits
  /* do subraction for overlapping hits  only correct immediate preceeding hit*/
  if (idet != 13 && doPeakCorrection)
  {

    // make a list of pointers for this detector
    std::vector<TDetHit> detHitList;
    for (hitMapIter hitIter1 = detHits.begin(); hitIter1 != detHits.end(); ++hitIter1)
    {
      detHitList.push_back(hitIter1->second);
      // printf(" line 1076 event %llu det %i uncorrected peak value peak bin %i qpeak %f \n", theEvent, idet, hitIter1->second.peakBin, hitIter1->second.qpeak);
    }

    /* correct qoeak for hit overlap */
    if (detHitList.size() > 1)
    {
      vector<unsigned> peakTimeList;
      vector<unsigned> indexList;
      // loop over this list of hits
      for (unsigned j = 0; j < detHitList.size() - 1; ++j)
      {
        TDetHit hitj = detHitList[j]; // earlier hit
        // map is ordered by hit time, so only need next hit
        int nextIndex = j + 1;
        TDetHit hitNext = detHitList[nextIndex]; // later hit
        int overLap = hitNext.peakBin - hitj.peakBin;
        hOverlap->Fill(overLap);
        // insure this peak is after previous, and separation is greater than minOverlap
        if (overLap < minOverlap)
        {
          // printf("line1094  hitFinder::makeHit event hit OVERLAP  %llu det %i  this hit (%i,%i,%i) last peak (%i,%i,%i) overlap %i  this qpeak  %f last qpeak %f \n", theEvent, idet, hitNext.firstBin, hitNext.peakBin, hitNext.lastBin, hitj.firstBin, hitj.peakBin, hitj.lastBin, overLap, hitNext.qpeak, hitj.qpeak);
          peakTimeList.push_back(hitNext.peakt);
          indexList.push_back(nextIndex);
          splitCount[idet] += 1;
          // fit is in axis value fit to previous peak
          double fitStart = hEvWave[idet]->GetBinCenter(hitj.peakBin);
          double fitEnd = hEvWave[idet]->GetBinCenter(hitj.lastBin);
          hEvWave[idet]->Fit("landau", "Q", "", fitStart, fitEnd); // Q for quiet
          //  switch to landau ?? offset too small!
          // TFitResultPtr fitptr = hEvWave[idet]->Fit("landau", "QS0", "", hitj.peakBin, hitj.lastBin); // was 20
          TF1 *expFit = (TF1 *)hEvWave[idet]->GetListOfFunctions()->FindObject("landau");
          // printf("line1092 expFit !!!!!!!!!!!!!!");
          // expFit->Print("all");
          if (!expFit->IsValid())
          {
            printf("line1092 fit to expo fails");
            continue;
          }
          // double slope = expFit->GetParameter(1);
          double slope = expFit->GetParameter(1);
          double offSet = expFit->Eval(hitNext.peakBin); // evaluate at peak of hit to correct

          if (offSet > 0. && offSet < nominalGain)
          { // this is hack for bad fit
            double qpeakBefore = hitNext.qpeak;
            detHitList[nextIndex].qpeak -= offSet;
            hPeakCorrectionOffset->Fill(offSet);
            if (verbose)
              printf("line1121  hitFinder::makeHit event %llu det %i hit %i found overlap this hit (%i,%i,%i) last peak (%i,%i,%i) fit range (%fi,%f) slope  %f offset %f  peak was %f corrected %f \n", theEvent, idet, ++nhit, hitNext.firstBin, hitNext.peakBin, hitNext.lastBin, hitj.firstBin, hitj.peakBin, hitj.lastBin, fitStart, fitEnd, slope, offSet, qpeakBefore, detHitList[nextIndex].qpeak);
            //   having corrected this peak,
          }
          else
          {
            if (verbose)
              printf("line1126 hitFinder::makeHit event BAD OFFSET  %llu det %i  found overlap this hit (%i,%i,%i) last peak (%i,%i,%i) sigma %f offset %f nominalGain %f \n", theEvent, idet, hitNext.firstBin, hitNext.peakBin, hitNext.lastBin, hitj.firstBin, hitNext.peakBin, hitNext.lastBin, slope, offSet, nominalGain);
          }
          // correct
          // overlap fix hitj is the first
          // fill histogram from hit
        }
      }

      // correct structure
      for (unsigned j = 0; j < peakTimeList.size(); ++j)
      {
        double oldPeak = detHits.at(peakTimeList[j]).qpeak;
        // new peak
        // this line does not wpork
        // detHits.at(peakTimeList[j]).qpeak = detHits.at(peakTimeList[j]).qpeak;
        detHits.erase(peakTimeList[j]);
        detHits.insert(std::pair<Double_t, TDetHit>(detHitList[indexList[j]].peakt, detHitList[indexList[j]]));
        if (verbose)
          printf(" line1144 event %llu det %i corected peak value peak bin %i time %u qpeak %f to %f  \n", theEvent, idet, detHits.at(peakTimeList[j]).peakBin, peakTimeList[j], oldPeak, detHits.at(peakTimeList[j]).qpeak);
      }

      /*
      int iterNumber = 0;
      for (hitMapIter hitIter1 = detHits.begin(); hitIter1 != detHits.end(); ++hitIter1)
      {
        printf(" line1151  event %llu det %i corected peak %i value %f  \n", theEvent, idet, ++iterNumber, hitIter1->second.qpeak);
        // hitIter1->second.qpeak = detHitList[hitNumber++].qpeak;
      }
        */
    }
  }
  // first time, charge from map
  /*
  hitMapIter hitIter;
  hitIter=detHits.begin();
  TDetHit dhit0 = hitIter->second;
  triggerTime = dhit0.startTime*microSec;
  firstCharge = dhit0.qsum;
  */
  /* bug fix */
  // printf("line1184  hitFinder::makeHits return event %lld det %i with %lu hits \n", theEvent, idet, detHits.size());
  unsigned jhit = 0;
  for (hitMapIter hitIter1 = detHits.begin(); hitIter1 != detHits.end(); ++hitIter1)
  {
    if (hitIter1->second.qpeak > 1.E5)
      printf("line1188 BUG very large qpeak hitFinder  channel %i hit %u bin %i last bin %i  qpeak  %E \n", idet, jhit++, hitIter1->second.firstBin, hitIter1->second.lastBin, hitIter1->second.qpeak);
  }

  if (verbose)
    printf(" hitFinder::makeHits return event %lld det %i with %lu made \n", theEvent, idet, detHits.size());
  return;
}

//
void hitFinder::findPeakCrossings(Int_t idet, unsigned peakStart, unsigned peakEnd)
{
  peakCrossings.clear();
  peakCrossingBin.clear();
  peakCrossingTime.clear();
  unsigned vsize = ddigi.size();
  peakThreshold = 7.0;
  if (vsize < peakEnd)
    return;
  Double_t cut = tbrun->detList[idet]->sigma * peakThreshold;
  unsigned step = 1;
  // find all crossings
  for (unsigned ibin = peakStart; ibin < peakEnd; ++ibin)
  {
    Double_t u = double(ibin) * timeUnit;
    Double_t vi = ddigi[ibin];
    Double_t vj = ddigi[ibin + step];
    unsigned ctype = 10;
    // crossing types
    if (vi < cut && vj > cut)
    {
      if (verbose)
        printf("line963 PUP det %i  bin %i %f %f  \n", idet, ibin, vi, vj);
      peakCrossings.push_back(PUP);
      ctype = PUP;
      peakCrossingBin.push_back(ibin + 1);
      peakCrossingTime.push_back(u);
    }
    else if (vi > cut && vj < -cut)
    {
      if (verbose)
        printf("line972 UPDOWN det %i  bin %i %f %f  \n", idet, ibin, vi, vj);
      peakCrossings.push_back(UPDOWN);
      ctype = UPDOWN;
      peakCrossingBin.push_back(ibin + 1);
      peakCrossingTime.push_back(u);
    }
    else if (vi > cut && vj < cut)
    {
      if (verbose)
        printf("line981 NUP det %i  bin %i %f %f  \n", idet, ibin, vi, vj);
      peakCrossings.push_back(NUP);
      ctype = NUP;
      peakCrossingBin.push_back(ibin + 1);
      peakCrossingTime.push_back(u);
    }
    else if (vi < -cut && vj > cut)
    {
      if (verbose)
        printf("line990  DOWNUP det %i  bin %i %f %f  \n", idet, ibin, vi, vj);
      peakCrossings.push_back(DOWNUP);
      ctype = DOWNUP;
      peakCrossingBin.push_back(ibin + 1);
      peakCrossingTime.push_back(u);
    }
    else if (vi < -cut && vj > -cut)
    {
      if (verbose)
        printf("line999  PDOWN det %i  bin %i %f %f  \n", idet, ibin, vi, vj);
      peakCrossings.push_back(PDOWN);
      ctype = PDOWN;
      peakCrossingBin.push_back(ibin + 1);
      peakCrossingTime.push_back(u);
    }
    else if (vi > -cut && vj < -cut)
    {
      peakCrossings.push_back(NDOWN);
      ctype = NDOWN;
      peakCrossingBin.push_back(ibin + 1);
      peakCrossingTime.push_back(u);
    }
    // if (idet==5&&ctype<10)  printf("....... %u vj %f vi %f cut %f cross type %u \n", ibin, vj, vi, cut, ctype );
    // if (idet==1&&ibin>2350&&ibin<2450)  printf("\t %u vj %f vi %f ctype %u  \n", ibin, vj, vi, ctype );
  }

  return;
}

void hitFinder::fitSinglet(int idet, Long64_t ievent)
{
  fSinglet = NULL;
  if (idet == 12)
    return;
  // first find peak max for fit range
  double ymax = 0;
  int maxBin = nominalTrigger - 30;
  for (int ibin = nominalTrigger - 30; ibin < trigEnd; ++ibin)
  {
    if (hEvWave[idet]->GetBinContent(ibin) > ymax)
    {
      ymax = hEvWave[idet]->GetBinContent(ibin);
      maxBin = ibin;
    }
  }
  // need to save this so we dont subtract from this peak
  singletPeakTime = unsigned(maxBin);
  // printf("line1079 fitSinglet  idet %i event %lld %s \n", idet, ievent, hEvWave[idet]->GetName());
  //  this prevents crash!!!

  /* do not do this: memory leak first clone*/
  // TH1D* hEvClone = (TH1D*) hEvWave[idet]->Clone("EvClone");
  hEvWave[idet]->GetListOfFunctions()->Clear();
  // do not make TCanvas // MINUIT error matrix not postive def. switch to Likleihood
  double fitStart = hEvWave[idet]->GetBinLowEdge(maxBin - 10);
  double fitEnd = hEvWave[idet]->GetBinLowEdge(maxBin + 20);
  TFitResultPtr fitptr = hEvWave[idet]->Fit("landau", "QS0", "", fitStart, fitEnd); // was 20
  int fitStatus = fitptr;
  // check its(int) value which is 0 if ok, -1 if not .
  //  status = 0 : the fit has been performed successfully(i.e no error occurred).
  // hEvWave[idet]->GetListOfFunctions()->ls();
  if (fitStatus == 0)
    fSinglet = (TF1 *)hEvWave[idet]->GetListOfFunctions()->FindObject("landau");
  // status = migradStatus + 10*minosStatus + 100*hesseStatus + 1000*improveStatus.
  if (!fSinglet)
  {
    printf("line1095  hitFinder::fitSinglet fSinglet NULL for det %i event %lld  range (%0.f,%0.f) fitStatus %i \n", idet, ievent, fitStart, fitEnd, fitStatus);
    if (fitSingletDir->GetList()->GetEntries() < 100)
      plotEvent(fitSingletDir, idet, ievent);
  }
}

// split peak based on derivaive
// requires digi, ddigi vectors
void hitFinder::splitPeaks(int idet)
{
  peakType addPeak;
  std::vector<Int_t> addKind;
  std::vector<unsigned> erasePeak;
  addPeak.clear();
  addKind.clear();
  erasePeak.clear();
  splitVerbose = true;
  vector<int> splitAt; // list of peaks to be split
  if (peakList.size() < 1 || digi.size() < 1 || ddigi.size() < 1)
    return;
  hEvPeakCross[idet]->Reset("ICESM");
  //
  // hPeakCrossingBin->Fill(0);
  Double_t cut = tbrun->detList[idet]->sigma * peakThreshold;
  vector<unsigned> indexSplit;
  vector<unsigned> isplit;

  // loop over peaks
  for (unsigned ip = 0; ip < peakList.size(); ++ip)
  {
    // access value in the memory to which the pointer
    // is referencing
    unsigned peakStart = std::get<0>(peakList[ip]);
    unsigned peakEnd = std::get<1>(peakList[ip]);
    // max for this peak
    double peakMax = 0;
    unsigned nsplits = 0;
    indexSplit.clear();
    isplit.clear();
    for (unsigned k = peakStart; k < peakEnd; ++k)
    {
      if (digi[k] > peakMax)
        peakMax = digi[k];
    }

    // use peak crossings from derivative
    findPeakCrossings(idet, peakStart, peakEnd);
    if (splitVerbose)
      printf(" \t peak  = %u   max %f crossings %lu cut  %f  \n", ip, peakMax, peakCrossings.size(), cut);

    for (unsigned ipc = 0; ipc < peakCrossings.size(); ++ipc)
    {
      bool pickit = theEvent == 0 && vChannel[idet] == 6;

      if (peakCrossings[ipc] == PDOWN || peakCrossings[ipc] == NDOWN || pickit)
      {
        double ratio = digi[peakCrossingBin[ipc]] / peakMax;
        double binAtr = 20. - (20. / .3) * ratio;
        double subBin = double(peakCrossingBin[ipc] - peakStart);
        // study splitting
        if (splitVerbose && (ratio < 0.5 || (theEvent == 0 && vChannel[idet] == 6)))
        {
          printf("line1082 event %llu idet %i  peak %u npeaks = %lu peakThreshold %.2f ratio %f binAtr %f \n", theEvent, vChannel[idet], ip, peakList.size(), peakThreshold, ratio, binAtr);
          printf("line1083\t\t crossing  %i type %i bin %i peakStart %i  peakEnd %i ddigi %f digi %f max %f ratio to max %f  subBin %.0f \n", ipc, peakCrossings[ipc], peakCrossingBin[ipc], peakStart, peakEnd, ddigi[peakCrossingBin[ipc]], digi[peakCrossingBin[ipc]], peakMax, ratio, subBin);
        }

        // split peak at largest subBin
        if (subBin > binAtr && binAtr > 0 && ratio < 0.5)
        {
          isplit.push_back(peakCrossingBin[ipc]);
          indexSplit.push_back(ip);
          if (splitVerbose)
            printf("line1092 sssssssssss splitting peak number %u nsplits %lu \n", ip, isplit.size());
        }
        /*
        hPeakCrossingBin->Fill(peakCrossingBin[ipc] - peakStart);
        hEvPeakCross[idet]->SetBinContent(peakCrossingBin[ipc], digi[peakCrossingBin[ipc]]);
        hPeakCrossingRatio->Fill(digi[peakCrossingBin[ipc]] / peakMax);
        */
        if (ntSplit->GetEntries() < 1.E6)
          ntSplit->Fill(theEvent, float(vChannel[idet]), float(ipc), float(indexSplit.size()), float(peakCrossingBin[ipc] - peakStart), float(ratio), float(binAtr), float(peakEnd - peakStart));
      }
    } // peak crossing loop

    if (indexSplit.size() > 0)
    {
      printPeakList("before");
      for (unsigned index = 0; index < indexSplit.size(); ++index)
      {
        erasePeak.push_back(indexSplit[index]);
        addPeak.push_back(std::make_pair(peakStart, isplit[index]));
        addKind.push_back(0);
        addPeak.push_back(std::make_pair(isplit[index], peakEnd));
        addKind.push_back(0);
      }
      printPeakList("after");
    }

    splitCount[idet] += indexSplit.size();
  } // peakList loop

  // remove old
  if (addPeak.size() > 0)
  {
    printf(" event %llu idet %i  npeaks = %lu peakThreshold %.2f \n", theEvent, vChannel[idet], peakList.size(), peakThreshold);
    printf(" before %lu erase %lu  starting peak list \n", peakList.size(), erasePeak.size());
    printPeakList("before");
    //
    for (unsigned jp = 0; jp < min(erasePeak.size(), peakList.size()); ++jp)
    {
      if (peakList.begin() + erasePeak[jp] < peakList.end())
      {
        printf("PEAK ERASE peak index %u  at %u  \n", jp, erasePeak[jp]);
        peakList.erase(peakList.begin() + erasePeak[jp]);
      }
      else
        printf("line1389 BAD PEAK ERASE POSITION  sizze %lu index %u at %u !!!! \n", peakList.size(), jp, erasePeak[jp]);
    }

    if (splitVerbose && addPeak.size() > 0)
    {
      printf(" ADD PEAKS event %llu chan %i add %lu \n ", theEvent, vChannel[idet], addPeak.size());
      for (unsigned jp = 0; jp < addPeak.size(); ++jp)
      {
        printf("line 1398 add peak %i (%i,%i) \n", jp, std::get<0>(addPeak[jp]), std::get<1>(addPeak[jp]));
      }
    }

    // I dont think the order matters so put them at the end
    for (unsigned jp = 0; jp < addPeak.size(); ++jp)
    {
      peakList.push_back(addPeak[jp]);
      peakKind.push_back(addKind[jp]);
    }
    printPeakList("after");
  }
}

void hitFinder::trimPeaks(int idet, std::vector<Double_t> v)
{

  if (peakList.size() < 1)
    return;
  for (unsigned ip = 0; ip < peakList.size(); ++ip)
  {
    // trim first peak
    unsigned peakStart = std::get<0>(peakList[ip]);
    unsigned peakEnd = std::get<1>(peakList[ip]);
    for (unsigned kp = peakEnd; kp > peakStart; --kp)
    {
      double vp = v[kp];
      if (vp > 0)
        break;
      std::get<1>(peakList[ip]) = kp;
    }

    for (unsigned kp = peakStart; kp < peakEnd; ++kp)
    {
      double vp = v[kp];
      if (vp > 0)
        break;
      std::get<0>(peakList[ip]) = kp;
    }
  }
}

std::vector<std::complex<double>> hitFinder::forwardFFT(std::vector<double> rdigi)
{
  unsigned nsamples = rdigi.size();
  std::vector<std::complex<double>> VectorComplex;
  for (unsigned is = 0; is < nsamples; ++is)
    fFFT->SetPoint(is, rdigi[is]);
  fFFT->Transform();

  std::vector<Double_t> realVec, imVec;
  for (unsigned i = 0; i < nsamples; ++i)
  {
    double rl, im;
    fFFT->GetPointComplex(i, rl, im);
    std::complex<double> c(rl, im); //.real or .imag accessors
    VectorComplex.push_back(c);
  }
  return VectorComplex;
}

std::vector<Double_t> hitFinder::backwardFFT(std::vector<std::complex<double>> VectorComplex)
{
  unsigned nsamples = VectorComplex.size();
  std::vector<Double_t> Signal;
  for (int is = 0; is < nsamples; ++is)
  {
    fInverseFFT->SetPoint(is, VectorComplex[is].real(), VectorComplex[is].imag());
  }
  fInverseFFT->Transform();

  for (unsigned i = 0; i < nsamples; ++i)
  {
    double rl = fInverseFFT->GetPointReal(i);
    Signal.push_back(rl);
  }

  // normalize
  for (unsigned i = 0; i < Signal.size(); ++i)
    Signal[i] /= double(nsamples);

  return Signal;
}

void hitFinder::plot1Wave(TDirectory *dir, int idet, Long64_t jentry)
{
  dir->cd();
  TString histName;
  TString detName = tbrun->detList[idet]->GetName();
  histName.Form("EvWave%lli%s", jentry, detName.Data());
  TH1D *hwave = (TH1D *)hEvWave[idet]->Clone(histName);
  hwave->SetTitle(histName);
}

void hitFinder::plotWave(int idet, Long64_t jentry)
{

  printf(" \t plotWave idet %i event %lld  %lu %lu %lu %lu \n", idet, jentry, digi.size(), ddigi.size(), hdigi.size(), fdigi.size());

  TString hname;
  hname.Form("raw-det-%i-event-%lli", idet, jentry);
  TH1S *hraw = new TH1S(hname, hname, nsamples, 0, nsamples);

  hname.Form("der-det-%i-event-%lli", idet, jentry);
  TH1S *hder = new TH1S(hname, hname, nsamples, 0, nsamples);

  hname.Form("hit-det-%i-event-%lli", idet, jentry);
  TH1S *hhit = new TH1S(hname, hname, nsamples, 0, nsamples);

  hname.Form("filter-det-%i-event-%lli", idet, jentry);
  TH1S *hfilt = new TH1S(hname, hname, nsamples, 0, nsamples);

  for (int i = 0; i < rdigi.size(); ++i)
    hraw->SetBinContent(i + 1, rdigi[i]);
  for (int i = 0; i < ddigi.size(); ++i)
    hder->SetBinContent(i + 1, ddigi[i]);
  for (int i = 0; i < hdigi.size(); ++i)
    hhit->SetBinContent(i + 1, hdigi[i]);
  for (int i = 0; i < fdigi.size(); ++i)
    hfilt->SetBinContent(i + 1, fdigi[i]);

  TString cname;
  cname.Form("det-%i-event-%lli-nhits-%ld", idet, jentry, detHits.size());
  TCanvas *can = new TCanvas(cname, cname);
  can->Divide(1, 4);
  can->cd(1);
  hraw->Draw();
  can->cd(2);
  hfilt->Draw();
  can->cd(3);
  hder->Draw();
  can->cd(4);
  hhit->Draw();
  can->Print(".gif");
}

void hitFinder::plotEvent(TDirectory *dir, unsigned ichan, Long64_t ievent)
{
  int idet = chanMap.at(ichan);
  int nhits = tbrun->detList[idet]->hits.size();
  // evDir->cd();
  dir->cd();
  // printf("hitFinder::plotEvent %s ichan %i event %lld \n",dir->GetName(),ichan, ievent);
  TString histName;
  TString histTitle;
  TString detName = tbrun->detList[idet]->GetName();

  histName.Form("EvWave%lli%s", ievent, detName.Data());
  TH1D *hwave = (TH1D *)hEvWave[idet]->Clone(histName);
  hwave->SetTitle(histName);

  /*
  // cout << " det " << idet << " "  << hwave->GetName() << " ," << hwave->GetTitle() << endl;
  histName.Form("EvSmooth%lli_%s", ievent, detName.Data());
  TH1D *hsmooth = (TH1D *)hEvSmooth[idet]->Clone(histName);
  hsmooth->SetTitle(histName);
  */

  histName.Form("EvDerWave%lli%s", ievent, detName.Data());
  TH1D *hdwave = (TH1D *)hEvDerWave[idet]->Clone(histName);
  hdwave->SetTitle(histName);

  histTitle.Form("EvHitPeakWave%lli%s hits %i", ievent, detName.Data(), nhits);
  histName.Form("EvHitPeakWave%lli%s", ievent, detName.Data());
  TH1D *hhitPeakWave = (TH1D *)hEvHitPeakWave[idet]->Clone(histName);
  hhitPeakWave->SetTitle(histTitle);

  /*
  histName.Form("EvCross%lli%s", ievent, detName.Data());
  TH1D *hcross = (TH1D *)hEvCross[idet]->Clone(histName);
  hcross->SetTitle(histName);

  histName.Form("EvPeakCross%lli%s", ievent, detName.Data());
  TH1D *hpeakCross = (TH1D *)hEvPeakCross[idet]->Clone(histName);
  hpeakCross->SetTitle(histName);

  histName.Form("EvFiltWave%lli_%s", ievent, detName.Data());
  TH1D *hfiltwave = (TH1D *)hEvFiltWave[idet]->Clone(histName);

   histName.Form("EvFFT%lli_%s", ievent, detName.Data());
   TH1D *hfft = (TH1D *)hFFT[idet]->Clone(histName);

   histName.Form("EvInvFFT%lli_%s", ievent, detName.Data());
   TH1D *hinvfft = (TH1D *)hInvFFT[idet]->Clone(histName);

  histName.Form("EvHitWave%lli_%s", ievent, detName.Data());
  TH1D *hhitWave = (TH1D *)hEvHitWave[idet]->Clone(histName);

  histName.Form("EvFFTFilt%lli_DET%1i_%s", ievent, idet, detName.Data());
  TH1D *hfftfilt = (TH1D *)hFFTFilt[idet]->Clone(histName);

  histName.Form("EvBase%lli_DET%1i_%s", ievent,idet,detName.Data());
  TH1D* hbase = (TH1D*)hBaselineWMA[idet]->Clone(histName);

  histName.Form("EvDWave%lli_DET%1i_%s", ievent,idet,detName.Data());
  TH1D* hdwave = (TH1D*)hEvDWave[idet]->Clone(histName);

  histName.Form("EvSignal%lli_DET%1i_%s", ievent,idet,detName.Data());
  TH1D* hsignal = (TH1D*)hEvSignal[idet]->Clone(histName);

  histName.Form("EvPeaks%lli_DET%1i_%s", ievent,idet,detName.Data());
  TH1D* hpeaks = (TH1D*)hEvPeaks[idet]->Clone(histName);

  histName.Form("EvDPeaks%lli_DET%1i_%s", ievent,idet,detName.Data());
  TH1D* hdpeaks = (TH1D*)hEvDPeaks[idet]->Clone(histName);

  histName.Form("EvWeight%lli_DET%1i_%s", ievent,idet,detName.Data());
  TH1D* hweight = (TH1D*)hEvWeight[idet]->Clone(histName);
  */

  fout->cd();
}
/*
if (ibin>nominalTrigger && ibin < trigEnd && fSinglet)
      {
        double xbin = hEvWave[idet]->GetBinLowEdge(ibin);
        lowCut = fSinglet->Eval(xbin);
        if (lowCut < nominalLowCut)
          lowCut = nominalLowCut;
        // printf("line751 hitFinder::makePeaks imax %i ibin %i xbin %f lowCut %f v %f \n",imax,ibin,xbin,lowCut,v[ibin]);
      }
*/