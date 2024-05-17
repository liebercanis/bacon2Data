//  M.Gold June 2022
// revised Jan 27 2023 -- simplified finding
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

hitFinder::hitFinder(TFile *theFile, TBRun *brun, TString theTag, int nSamples, vector<int> vchan, vector<double> sigmaValue)
{
  isCAEN = false;
  if (nSamples == CAENLENGTH)
    isCAEN = true;
  channelSigmaValue = sigmaValue;
  verbose = false;
  if (isCAEN)
    QPEPeak = 200;
  else
    QPEPeak = 50;
  for (unsigned i = 0; i < vchan.size(); ++i)
    QPEnominal.push_back(QPEPeak);
  templateFileName = TString("$HOME/bacon2Data/bobj/templates-2023-05-01-15-06.root");
  // CAEN casN
  if (nSamples == CAENLENGTH)
  {
    templateFileName = TString("$HOME/bacon2Data/bobj/templatesCaen-2023-05-17-12-00.root");
    for (unsigned i = 0; i < 9; ++i)
      QPEnominal[i] = 10.0;
    QPEnominal[9] = 50.0;
    QPEnominal[10] = 50.0;
    QPEnominal[11] = 50.0;
    QPEnominal[12] = 3.5;
  }
  // save vchan
  vChannel = vchan;
  if (verbose)
    cout << "INSTANCE OF HITFINDER "
         << " vchan.size " << vchan.size() << endl;
  smoothing = false;
  fout = theFile;
  fftDir = fout->mkdir("fftDir");
  finderDir = fout->mkdir("finderDir");
  splitDir = fout->mkdir("splitDir");
  sumWaveDir = fout->mkdir("sumWaveDir");
  tag = theTag;
  tbrun = brun;
  nsamples = nSamples;
  int nSize = nsamples + 100;
  // initialize fft
  if (verbose)
    cout << " initialize  FFT  " << endl;
  fFFT = TVirtualFFT::FFT(1, &nSize, "R2C M K");
  fInverseFFT = TVirtualFFT::FFT(1, &nSize, "C2R M K");

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
    hMaxBinVal.push_back(new TH1D(Form("MaxBinVal%i", id), Form(" max bin val/gain chan  %i ", id),210, -1, 20));
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
    hDerivativeValTime.push_back(new TH2D(Form("DerivativeValTime%i", id), Form("derivative value vs sample chan %i", id), nsamples, 0, nsamples,  2000, -1000., 1000.));
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
  //ntFinder->Fill(float(theEvent), float(idet), float(detHits.size()), float(dhit.firstBin), float(dhit.startTime), float(dhit.peakBin), float(dhit.lastBin));
  ntFinder = new TNtuple("ntFinder", " hit finder ", "event:chan:nhit:first:start:peak:last:qpeak");
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
  if (gotTemplate)
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
  printf("QPE: \n");
  for (unsigned ichan = 0; ichan < QPEnominal.size(); ++ichan)
    printf("chan %i QPEnominal %f ; ", ichan, QPEnominal[ichan]);
  printf("\n");
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
    printf(" couldnt open template file %s\n", templateFileName.Data());
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

void hitFinder::printPeakList()
{
  cout << "peakList size " << peakList.size() << endl;
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
  // PMT bad
  if (ichan == 12)
    return;
  /*if(ichan == 7 && ievent==0 )
    verbose = true;
    else
      verbose = false;
      */
  /////  copy to internal class vector////////
  digi = inputDigi;
  QPEPeak = QPEnominal[ichan];
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
  // use filtered waveforms
  // for (unsigned isample = 0; isample < 20; isample++)
  // printf(" wfilter ??? %i %f %f ?? %f \n", isample, wfilter[isample], digi[isample], fdigi[isample]);
  // if (gotTemplate) {
  //   digi = fdigi;
  //}
  hEvAllSumWave->Reset("ICESM");
  for (unsigned isample = 0; isample < digi.size(); isample++)
  {
    hEvWave[idet]->SetBinContent(isample + 1, digi[isample]);
    hEvSmooth[idet]->SetBinContent(isample + 1, digi[isample]);
    hEvFiltWave[idet]->SetBinContent(isample + 1, fdigi[isample]);
    hInvFFT[idet]->SetBinContent(isample + 1, fdigi[isample]);
    // sum all waves for this event
    hEvAllSumWave->SetBinContent(isample + 1, hEvAllSumWave->GetBinContent(isample + 1) + digi[isample]);
  }
  // smooth and fill vector
  hEvSmooth[idet]->Smooth(1); // one time
  sdigi.resize(digi.size());
  for (unsigned ibin = 1; ibin < hEvSmooth[idet]->GetNbinsX(); ibin++)
    sdigi[ibin - 1] = hEvSmooth[idet]->GetBinContent(ibin);

  // use smooth wave if smoothing
  if (verbose)
    printf("line401  smoothing ? %i  digi size %lu \n", smoothing, digi.size());
  if (!smoothing)
    digi = sdigi;

  ddigi.clear();
  differentiate();
  for (unsigned isample = 0; isample < ddigi.size(); isample++)
  {
    hDerivativeVal[idet]->Fill(ddigi[isample]);
    hDerivativeValTime[idet]->Fill(double(isample),ddigi[isample]);
    hEvDerWave[idet]->SetBinContent(isample + 1, ddigi[isample]);
  }
  // find peaks
  // for derivativePeaks, window in time is timeUnit*windowSize (ns) . timeUnit = 2
  // min, max width in time bins for simple peaks
  Int_t windowSize = 10;
  unsigned maxWidth = 100000;
  unsigned minWidth = 10;
  findDerivativeCrossings(idet);
  // findThresholdCrossings(idet, threshold);
  fitSinglet(idet,ievent);
  makePeaks(idet, digi);
  //splitPeaks(idet);
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
  for (hitMapIter hitIter = detHits.begin(); hitIter != detHits.end(); ++hitIter)
  {
    TDetHit hiti = hitIter->second;
    tbrun->detList[idet]->hits.push_back(hiti);

    // fill hit digi
    for (unsigned iv = 0; iv < digi.size(); ++iv)
      if (iv >= hiti.firstBin && iv <= hiti.lastBin)
        hdigi[iv] = digi[iv];
    // if (hiti.qsum > 7000 && hiti.qsum < 10000) // FILL ONLY SINGLE PE

    // fill hit peak wave
    hEvHitPeakWave[idet]->SetBinContent(hiti.peakBin, hiti.qpeak);
    if(verbose) printf("line461 size %lu hit%i idet %i time %f peakBin %i qpeak  %f \n", detHits.size(),hitNumber++,idet, hitIter->first, hiti.peakBin, hiti.qpeak);
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
    if(verbose) printf("line473 hitFinder event %lld chan %i thres %f qpeak sum %f\n", ievent, ichan, hitThreshold, tbrun->detList[idet]->qpeak);
    if (!triggerChannel)
      hPeakValue->Fill(hiti.qpeak);
    hitTitle.Form("TDetHit %i event %llu chan %i index %i ", icount++, ievent, ichan, idet);
    if(verbose)
      cout << hitTitle << endl;
    hiti.SetTitle(hitTitle);
  }
  // save the hits
  //tbrun->fill();

  // save some split histograms
  for (unsigned idet = 0; idet < tbrun->detList.size(); ++idet)
    if (splitCount[idet] > 0 && splitDir->GetList()->GetEntries() < 500)
    {
      printf("line486 plot SplitEvent %llu %i \n", theEvent, idet);
      plotEvent(splitDir, idet, theEvent);
    }

  for (unsigned isample = 0; isample < hdigi.size(); isample++)
  {
    hEvHitWave[idet]->SetBinContent(isample + 1, hdigi[isample]);
    hHitSum[idet]->SetBinContent(isample + 1, hdigi[isample] + hHitSum[idet]->GetBinContent(isample + 1));
  }

  //
  if (1)
  {
    TDet *tdet = tbrun->detList[idet];
    if (tdet->hits.size() > 1 && verbose)
    {
      cout << "HHHH  END hitFinder::event " << theEvent << " idet= " << idet << " " << tdet->channel << " hits.size " << tdet->hits.size() << endl;
      for (unsigned ihit = 0; ihit < tdet->hits.size(); ++ihit)
      {
        cout << " \t finder hit number  " << ihit << " peak bin " << tdet->hits[ihit].peakBin << endl;
      }
    }
  }
  if(1) cout << "HHHH  END hitFinder::event ichan " << ichan << " event " << ievent << "  " << detHits.size() << endl;
}

// revised derivative Jan 27 2023 MG
void hitFinder::differentiate()
{
  if (verbose)
    printf("line516 hitFinder::differentiate nsamples %lu step %u\n", digi.size(), diffStep);
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
      sump += digi[i + 1 + j];
    }
    summ = 0;
    for (unsigned j = 0; j < maxSum; ++j)
    {
      summ += digi[i - 1 - j];
    }
    // if(verbose) printf(" hitFinder::differentiate bin %i maxSum %u sump %E summ %E \n",i,maxSum,sump,summ);
    ddigi[i] = sump - summ;
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
  crossings.clear();
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
      //if(idet==13&& ibin>1040&&ibin<1070) 
      //  printf("line635  PUP det %i  bin %i %f %f  \n", idet, ibin, digi[ibin], digi[ibin + 1]);
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
void hitFinder::makePeaks(int idet, std::vector<Double_t> v)
{
  if (verbose)
    printf("line697 hitFinder::makePeaks det %i crossings %lu \n", idet, crossings.size());
  double sigma = tbrun->detList[idet]->sigma;
  peakList.clear();
  peakKind.clear();
  hEvCross[idet]->Reset("ICESM");
  // loop over crossings using  PUP or NUP
  for (int icross = 0; icross < crossings.size(); ++icross)
  {
    if (!(crossings[icross] == PUP || crossings[icross] == NUP))
      continue;
    // find local max
    hCrossingBinA[idet]->Fill(crossingBin[icross]);
    unsigned imax = 0;

    //for (unsigned ibin = crossingBin[icross]; ibin < v.size(); ++ibin)
    double maxVal = -99999.;
    if (crossings[icross] == PUP)  // case PUP
    {
      for (unsigned ibin = crossingBin[icross]; ibin < v.size(); ++ibin)
      {
        if (v[ibin] < maxVal)
          break;
        imax = ibin;
        maxVal = v[ibin];
      }
    } else {  // NUP is other side of derivative going through zero
      maxVal = -99999.;
      for (unsigned ibin = crossingBin[icross]; ibin > 0; --ibin)
      {
        if (v[ibin] < maxVal)
          break;
        imax = ibin;
        maxVal = v[ibin];
      }
    }
    // too small is garbage
    hMaxBinVal[idet]->Fill(maxVal/nominalGain);
    if(maxVal<hitThreshold)
      continue;

    // not to close to start of wave
    if (imax<50)
      continue;

    hCrossingMaxBin[idet]->Fill(imax);

    if (verbose)
      printf("line740 hitFinder::makePeaks cross det %i icross %i maxVal %f  \n", idet, icross, maxVal);
    // find limits of peak
    /* 
    just use a fixed window around the maximum so look for peak  in peak -40  to peak +50*/
    
    
    unsigned ilow =  imax-20;
    unsigned ihigh = imax+50;
    if(ihigh>unsigned(v.size()-1))
      ihigh = unsigned(v.size()-1);
    if(ilow<0)
      ilow = 0;

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
    }
    */
    
    
    hCrossingBinB[idet]->Fill(ilow);
    double nominalLowCut = 0; // no longer used

    if (verbose)
      printf("line767  hitFinder::makePeaks cross det %i lowCut %f imax %i val %f icross %u from (%u,%u) \n", idet, nominalLowCut, imax, maxVal, crossingBin[icross], ilow, ihigh);
    //if (idet == 13 && crossingBin[icross] > 1040 && crossingBin[icross] < 1070)
    //  printf("line 769 hitFinder::makePeaks cross det %i lowCut %f imax %i val %f icross %u from (%u,%u) \n", idet, nominalLowCut, imax, maxVal, crossingBin[icross], ilow, ihigh);
    // if (ihigh - ilow > maxPeakLength)
    //   ihigh = ilow + maxPeakLength;

    // check that this hit has not already been found.
    bool found = false;
    for (unsigned ip = 0; ip < peakList.size(); ++ip)
    {
      // trim first peak
      unsigned peakStart = std::get<0>(peakList[ip]);
      unsigned peakEnd = std::get<1>(peakList[ip]);
      if (peakStart == ilow && peakEnd == ihigh)
        found = true;
    }
    if (!found)
    {
      if (verbose)
        printf("line786  hitFinder::makePeaks add det %i  imax %i val %f icross %u from (%u,%u) \n", idet, imax, maxVal, crossingBin[icross], ilow, ihigh);
      peakList.push_back(std::make_pair(ilow, ihigh));
      peakKind.push_back(0);
      hCrossingBinC[idet]->Fill(ilow);
      // if(vChannel[idet]==6) printf(" %lli add det %i  imax %i val %f icross %u from (%u,%u) %i\n",theEvent, idet, imax, maxVal, crossingBin[icross], ilow, ihigh,icross);
    }
  }
}

void hitFinder::makeHits(int idet, Double_t &triggerTime, Double_t &firstCharge)
{
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
    if (verbose)
      printf("line819 hitFinder::makeHits event %lli det %i hit  %u (%u,%u) kind %i length %u \n", theEvent, idet, ip, klow, khigh, peakKind[ip], khigh - klow);
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
    // cut small peaks below hitThreshold
    hPeakCut[idet]->Fill(qpeak);
    hPeakCutAndTime[idet]->Fill(peakt,qpeak);
    if (qpeak < hitThreshold)
      continue;

    unsigned kstart = TMath::Max(unsigned(0), klow - 20);
    unsigned kend = TMath::Min(unsigned(digi.size()), khigh + 20);

    TDetHit dhit;
    if (vChannel[idet] < 9)
      for (unsigned k = kstart; k < kend; ++k)
        dhit.digi.push_back(digi[k]);

    if (verbose)
      printf("line854 hitFinder::makeHits hit chan %i (%i,%i) size %lu \n ", vChannel[idet], klow, khigh, dhit.digi.size());

  
    dhit.peakBin = Int_t(peakt);
    dhit.qsum = qsum;
    dhit.qpeak = qpeak;
    dhit.firstBin = klow;
    dhit.lastBin = khigh;
    dhit.peakMaxTime = peakt;
    dhit.peakt = peakt;
    dhit.startTime = klow;
    dhit.peakWidth = khigh - klow + 1;
    // this is N= q/qnorm and delta q = root(n)*qnorm;
    dhit.qerr = sqrt(pow(sigma * Double_t(dhit.peakWidth), 2) + qnorm * qsum);
    dhit.kind = peakKind[ip];
    
    // just use the biggest pulse
    if (qsum > qmax)
    {
      qmax = qsum;
      triggerTime = dhit.startTime * timeUnit * microSec;
      firstCharge = qsum;
    }
    Double_t hitTime = dhit.startTime * timeUnit * microSec;

    // ensure new hit it does not have same peak bin and check
    bool used = false;
    for (hitMapIter hitIter = detHits.begin(); hitIter != detHits.end(); ++hitIter)
    {
      TDetHit hiti = hitIter->second;
      if (hiti.peakBin == dhit.peakBin)
      {
        used = true;
        if (verbose)
          printf("line887 hitFinder::makeHit found multiple hits det %i this (%i,%i,%i)  last peak (%i,%i,%i) dethit size %lu \n", idet, hiti.firstBin, hiti.peakBin, hiti.lastBin,dhit.firstBin, dhit.peakBin, dhit.lastBin, detHits.size());
      }
    }
    ntFinder->Fill(float(theEvent), float(idet), float(detHits.size()), float(dhit.firstBin), float(dhit.startTime), float(dhit.peakBin), float(dhit.lastBin), dhit.qpeak);

    if (used)
      continue;

    /*fix peak if after singlet */
    if (fSinglet != NULL && peakt > singletPeakTime && peakt < trigEnd)
    {
      double xbin = hEvWave[idet]->GetBinLowEdge(peakt);
      double offset = fSinglet->Eval(xbin);
      ntPeakFix->Fill(float(detHits.size()), float(idet), float(singletPeakTime), float(peakt), qpeak, qpeak - offset);
      if(verbose )printf("line919 list size %lu idet %i singlett %u peakt %u qpeak %f fixed %f\n", detHits.size(), idet, singletPeakTime, peakt, dhit.qpeak, dhit.qpeak - offset);
      // fix here 
      dhit.qpeak = dhit.qpeak - offset;
    }
    else if (verbose&&peakt > singletPeakTime && peakt < trigEnd)
    {
      printf("line925 detHits %lu  %i singlett %u peakt %u  \n", detHits.size(), idet, singletPeakTime, peakt);
    }

    // cheak after peak fix
    if (dhit.qpeak < hitThreshold)
      continue;

    detHits.insert(std::pair<Double_t, TDetHit>(hitTime, dhit));
    hPeakNWidth->Fill(dhit.lastBin - dhit.firstBin + 1);
    if (verbose) {
      printf("line941 hitFinder::makeHits %llu insert hit idet %i  time %f (%u,%u) peak bin %i kind %i length %u qpeak %f detHit size %lu  \n", theEvent, idet, hitTime, dhit.firstBin, dhit.lastBin, dhit.peakBin, peakKind[ip], khigh - klow + 1, qpeak, detHits.size());
    }
  }

    int nhit = 0;
    // this messes ip yaxis on chan13 EvWave??
    /* do this differently with very short hits
    if(idet!=13) for (hitMapIter hitIter1 = detHits.begin(); hitIter1 != detHits.end(); ++hitIter1)
    {
      TDetHit hitj = hitIter1->second;
      for (hitMapIter hitIter2 = detHits.begin(); hitIter2 != detHits.end(); ++hitIter2)
      {
        if (hitIter2 == hitIter1)
          continue;
        TDetHit hiti = hitIter2->second;
        if (hiti.peakBin > hitj.firstBin && hiti.peakBin < hitj.firstBin + 150 && hiti.peakBin != hitj.peakBin)
        {
          splitCount[idet] += 1;
          hEvWave[idet]->Fit("expo", "", "", hitj.peakBin, hitj.lastBin);
          TF1 *expFit = (TF1 *)hEvWave[idet]->GetListOfFunctions()->FindObject("expo");
          if(!expFit->IsValid())
            continue;
          double slope = expFit->GetParameter(1);
          double offSet = expFit->Eval(hiti.peakBin);
          if (offSet > 0. && offSet < nominalGain  ) { // this is hack for bad fit
            double qpeakBefore = hiti.qpeak;
          hitIter2->second.qpeak -= offSet;
          printf("line919 hitFinder::makeHit event %llu det %i hit %i found overlap this hit (%i,%i,%i) last peak (%i,%i,%i) slope %f offset %f  peak was %f corrected %f \n", theEvent, idet, ++nhit, hiti.firstBin, hiti.peakBin, hiti.lastBin, hitj.firstBin, hitj.peakBin, hitj.lastBin, slope, offSet, qpeakBefore, hitIter2->second.qpeak);
          }
          // correct
          // overlap fix hitj is the first
          // fill histogram from hit

          //
        }
      }
    }
    */
  // first time, charge from map
  /*
  hitMapIter hitIter;
  hitIter=detHits.begin();
  TDetHit dhit0 = hitIter->second;
  triggerTime = dhit0.startTime*microSec;
  firstCharge = dhit0.qsum;
  */
  if (verbose)
      printf(" hitFinder::makeHits return event %lld det %i with %lu made \n", theEvent, idet, detHits.size());
    return;
}

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
  //printf("line1079 fitSinglet  idet %i event %lld %s \n", idet, ievent, hEvWave[idet]->GetName());
  // this prevents crash!!!
  hEvWave[idet]->GetListOfFunctions()->Clear();
  hEvWave[idet]->Fit("landau", "RQS", "", maxBin, maxBin + 20);
  //check its(int) value which is 0 if ok, -1 if not .
  // status = 0 : the fit has been performed successfully(i.e no error occurred).
  //hEvWave[idet]->GetListOfFunctions()->ls();
  fSinglet = (TF1 *)hEvWave[idet]->GetListOfFunctions()->FindObject("landau");
  if (!fSinglet)
    printf("line1095  hitFinder::fitSinglet fSinglet NULL so returning det %i event %lld \n", idet, ievent);
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
  splitVerbose = false;
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
        hPeakCrossingBin->Fill(peakCrossingBin[ipc] - peakStart);
        hEvPeakCross[idet]->SetBinContent(peakCrossingBin[ipc], digi[peakCrossingBin[ipc]]);
        hPeakCrossingRatio->Fill(digi[peakCrossingBin[ipc]] / peakMax);
        ntSplit->Fill(theEvent, float(vChannel[idet]), float(ipc), float(indexSplit.size()), float(peakCrossingBin[ipc] - peakStart), float(ratio), float(binAtr), float(peakEnd - peakStart));
      }
    } // peak crossing loop

    for (unsigned index = 0; index < indexSplit.size(); ++index)
    {
      erasePeak.push_back(indexSplit[index]);
      addPeak.push_back(std::make_pair(peakStart, isplit[index]));
      addKind.push_back(0);
      addPeak.push_back(std::make_pair(isplit[index], peakEnd));
      addKind.push_back(0);
    }

    splitCount[idet] += indexSplit.size();
  } // peakList loop

  // remove old
  if (splitVerbose && addPeak.size() > 0)
  {
    printf(" event %llu idet %i  npeaks = %lu peakThreshold %.2f \n", theEvent, vChannel[idet], peakList.size(), peakThreshold);
    printf(" before %lu erase %lu \n", peakList.size(), erasePeak.size());
  }
  for (unsigned jp = 0; jp < min(erasePeak.size(), peakList.size()); ++jp)
  {
    if (peakList.begin() + erasePeak[jp] < peakList.end())
    {
      peakList.erase(peakList.begin() + erasePeak[jp]);
      peakKind.erase(peakKind.begin() + erasePeak[jp]);
    }
    else
      printf("BAD PEAK ERASE POSITION  %u !!!! \n", erasePeak[jp]);
  }

  if (splitVerbose && addPeak.size() > 0)
  {
    printPeakList();
    printf(" ADD PEAKS event %llu chan %i %lu \n ", theEvent, vChannel[idet], addPeak.size());
    for (unsigned jp = 0; jp < addPeak.size(); ++jp)
    {
      printf("\t peak %i (%i,%i) \n", jp, std::get<0>(addPeak[jp]), std::get<1>(addPeak[jp]));
    }
  }

  // I dont think the order matters so put them at the end
  for (unsigned jp = 0; jp < addPeak.size(); ++jp)
  {
    peakList.push_back(addPeak[jp]);
    peakKind.push_back(addKind[jp]);
  }

  if (splitVerbose && addPeak.size() > 0)
  {
    printf("------------------------\n");
    printf(" after %lu  \n", peakList.size());
    printPeakList();
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

void hitFinder::plot1Wave(TDirectory* dir,int  idet, Long64_t jentry)
{
  dir->cd();
  TString histName;
  TString detName = tbrun->detList[idet]->GetName();
  histName.Form("EvWave%lli%s",jentry, detName.Data());
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

  void hitFinder::plotEvent(TDirectory * dir, unsigned ichan, Long64_t ievent)
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