////////////////////////////////////////////////////////
//  M.Gold June 2022
// revised Jan 27 2023 -- simplified finding
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
  channelSigmaValue = sigmaValue;
  verbose = false;
  QPEPeak = 50;
  for (unsigned i = 0; i < vchan.size(); ++i)
    QPEnominal.push_back(QPEPeak);
  templateFileName = TString("../bobj/templates-2023-05-01-15-06.root");
  // CAEN casN
  if (nSamples == CAENLENGTH)
  {
    templateFileName = TString("../bobj/templatesCaen-2023-05-17-12-00.root");
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
  tag = theTag;
  tbrun = brun;
  nsamples = nSamples;
  int nSize = nsamples + 100;
  // initialize fft
  fFFT = TVirtualFFT::FFT(1, &nSize, "R2C M K");
  fInverseFFT = TVirtualFFT::FFT(1, &nSize, "C2R M K");

  microSec = 1.0E-3;
  timeUnit = 8.0; // ns per count
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
  hWFilter = new TH1D("WFilter","WFilter", nsamples / 2, 0, nsamples / 2);

  fout->cd();
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

  fout->cd("sumDir");
  for (unsigned index = 0; index < vchan.size(); ++index)
  {
    int id = vchan[index];
    TDet *deti = tbrun->getDet(id);
    hUnFilteredSummedWave.push_back(new TH1D(Form("UnFilteredSummedWave%s", deti->GetName()), Form(" un filtered summed wave%s", deti->GetName()), nsamples, 0, nsamples));
    hFilteredSummedWave.push_back(new TH1D(Form("FilteredSummedWave%s", deti->GetName()), Form("filtered summed wave%s", deti->GetName()), nsamples, 0, nsamples));
  }
  fout->cd();

  ntFinder = new TNtuple("finder", " hit finder ", "chan:ncross:npeak");
  ntSplit = new TNtuple("split", " split for finder ", "event:chan:cross:nsplit:bin:ratio:batr:width");

  int templateChan = 8;
  gotTemplate = getTemplate(templateChan);

  cout << " created hitFinder with " << tbrun->GetName() << " nsamples =  " << nsamples << " ndet " << hEvWave.size() << " ";
  if (gotTemplate)
    cout << " SPE Template " << htemplate->GetName() << endl;
  else
    cout << " SPE Template not found ! " << endl;

  // set wfilter size
  wfilter.resize(nsamples);
  for (int i = 0; i < nsamples; ++i) wfilter[i]=1.;
  //
  if (gotTemplate)
  {
    // make transorm
    templateTransform = forwardFFT(SPEdigi);
    // fill htemplateFFT start with first nonzero bin;
    printf(" ********   complex transform  size %lu ******** \n", templateTransform.size());
    
    // make filter
    fillWFilter(templateChan);
  }
  for (int i = 0; i < nsamples / 2; ++i)
  {
    hWFilter->SetBinContent(i, wfilter[i]);
    //printf(" wfilter %i %f \n", i, wfilter[i]);
    htemplateFFT->SetBinContent(i, std::abs(templateTransform[i]));
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
void hitFinder::fillWFilter(int ichan){
  double noiseVal = channelSigmaValue[ichan];
  for (int i = 0; i < nsamples; ++i)
  {
    double val = std::abs(templateTransform[i]);
    wfilter[i] = val / (val + noiseVal);
  }
}

bool hitFinder::getTemplate(int ichan)
{
  TH1D *hist = NULL;
  TFile *f1 = new TFile(templateFileName, "readonly");
  if (f1->IsZombie())
  {
    printf(" no  file for %s \n",templateFileName.Data());
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

void hitFinder::event(int ichan, Long64_t ievent, vector<double> eventDigi, double thresh, unsigned step)
{
  /*if(ichan == 7 && ievent==0 )
    verbose = true;
    else
      verbose = false;
      */
  QPEPeak = QPEnominal[ichan];
  digi = eventDigi;
  bool trig = ichan == 9 || ichan == 10 || ichan == 11;

  theEvent = ievent;
  threshold = thresh;
  diffStep = step;
  int idet = chanMap.at(ichan);
  splitCount.clear();
  for (int i = 0; i < vChannel.size(); ++i)
    splitCount.push_back(0);

  if (verbose)
    printf(" HHHHHH hitFinder ievent %llu ichan %i idet %i threshold %.1f \n", ievent, ichan, idet, tbrun->detList[idet]->sigma * threshold);

  double triggerTime = 0;
  double firstCharge = 0;

  for (int i = 0; i < nsamples; ++i)
  {
    hDigiVal[idet]->Fill(digi[i]);
  }
  // FFT and convolution
  std::vector<std::complex<double>> inputWaveTransform = forwardFFT(digi);
  for (int i = 0; i < nsamples / 2; ++i) {
    hFFT[idet]->SetBinContent(i, std::abs(inputWaveTransform[i]));
    hFFTFilt[idet]->SetBinContent(i, hFFTFilt[idet]->GetBinContent(i)+std::abs(inputWaveTransform[i]));
  }

  unsigned maxFrequency = inputWaveTransform.size(); 
  // apply FFT convolution here
  fillWFilter(ichan); // use channel noise
  for (unsigned iw = 1; iw < maxFrequency; ++iw)
  {
    
    // divide out the SPE shape
    inputWaveTransform[iw] = wfilter[iw] * inputWaveTransform[iw]; // templateTransform[iw];
  }
  fdigi = backwardFFT(inputWaveTransform);
  for (unsigned isample = 0; isample < digi.size(); isample++)
  {
    hUnFilteredSummedWave[idet]->SetBinContent(isample + 1, digi[isample] + hUnFilteredSummedWave[idet]->GetBinContent(isample + 1));
    hFilteredSummedWave[idet]->SetBinContent(isample + 1, fdigi[isample] + hFilteredSummedWave[idet]->GetBinContent(isample+1));
  }
  // use filtered waveforms
  //for (unsigned isample = 0; isample < 20; isample++)
    //printf(" wfilter ??? %i %f %f ?? %f \n", isample, wfilter[isample], digi[isample], fdigi[isample]);
  //if (gotTemplate) {
 //   digi = fdigi;
  //}

  for(unsigned isample = 0; isample < digi.size(); isample++)
  {
    hEvWave[idet]->SetBinContent(isample + 1, digi[isample]);
    hEvSmooth[idet]->SetBinContent(isample + 1, digi[isample]);
    hEvFiltWave[idet]->SetBinContent(isample + 1, fdigi[isample]);
    hInvFFT[idet]->SetBinContent(isample + 1, fdigi[isample]);
  }
  // smooth and fill vector
  hEvSmooth[idet]->Smooth(1); // one time
  sdigi.resize(digi.size());
  for (unsigned ibin = 1; ibin < hEvSmooth[idet]->GetNbinsX(); ibin++)
    sdigi[ibin - 1] = hEvSmooth[idet]->GetBinContent(ibin);

  // use smooth wave if smoothing
  if (!smoothing)
    digi = sdigi;

  differentiate();
  for (unsigned isample = 0; isample < ddigi.size(); isample++)
  {
    hEvDerWave[idet]->SetBinContent(isample + 1, ddigi[isample]);
  }
  // find peaks
  // for derivativePeaks, window in time is timeUnit*windowSize (ns) . timeUnit = 2
  // min, max width in time bins for simple peaks
  Int_t windowSize = 10;
  unsigned maxWidth = 100000;
  unsigned minWidth = 10;
  findThresholdCrossings(idet,threshold);
  makePeaks(idet, digi);
  // splitPeaks(idet);
  detHits = makeHits(idet, triggerTime, firstCharge);
  hPeakCount->Fill(idet, peakList.size());
  // fill hits
  if (verbose)
    cout << " finished make hits event " << ievent << " chan " << ichan << " det " << idet
         << "  ddigi size " << ddigi.size()
         << "  crossings size " << crossings.size()
         << "  peakList size " << peakList.size()
         << "  detHits size " << detHits.size() << endl
         << endl;
  hdigi.clear();
  hdigi.resize(digi.size());

  //  for (const auto &[key, value] : m)
  //    std::cout << '[' << key << "] = " << value << "; "
  // push hits to tbrun
  int icount = 0;
  bool triggerChannel = false;
  double startTimeCut = 70.0;
  if (ichan == 9 || ichan == 10 || ichan == 11)
    triggerChannel = true;
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
    // make sums with cut
    if (hiti.qsum > hitQThreshold)
    {
      tbrun->detList[idet]->qSum += hiti.qsum;
      tbrun->detList[idet]->hitSum += hiti.qpeak;
      if (hiti.startTime < startTimeCut)
      {
        tbrun->detList[idet]->qPrompt += hiti.qsum;
        tbrun->detList[idet]->hitPrompt += hiti.qpeak;
      }
    }
    if (!triggerChannel)
      hPeakValue->Fill(hiti.qpeak);
    hiti.SetTitle(Form("TDetHit %i event %llu chan %i index %i ", icount++, ievent, ichan, idet));
  }
  // save the hits
  tbrun->fill();

  // save some split histograms
  for (unsigned idet = 0; idet < tbrun->detList.size(); ++idet)
    if (splitCount[idet] > 0 && splitDir->GetList()->GetEntries() < 500)
      plotSplitEvent(idet, theEvent);

  for (unsigned isample = 0; isample < hdigi.size(); isample++)
  {
    hEvHitWave[idet]->SetBinContent(isample + 1, hdigi[isample]);
    hHitSum[idet]->SetBinContent(isample + 1, hdigi[isample] + hHitSum[idet]->GetBinContent(isample + 1));
  }

  //
  ntFinder->Fill(float(idet), float(crossings.size()), float(peakList.size()));
  if (verbose)
  {
    for (unsigned idet = 0; idet < tbrun->detList.size(); ++idet)
    {
      TDet *tdet = tbrun->detList[idet];
      cout << tdet->channel << " hits.size " << tdet->hits.size() << endl;
      for (unsigned ihit = 0; ihit < tdet->hits.size(); ++ihit)
      {
          cout << theEvent << " finder hit number  " << ihit << " bin " << tdet->hits[ihit].peakBin << endl;
      }
    }
  }
}

// revised derivative Jan 27 2023 MG
void hitFinder::differentiate()
{
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
      sump += digi[i + 1 + j];
    summ = 0;
    for (unsigned j = 0; j < maxSum; ++j)
      summ += digi[i - 1 - j];
    ddigi[i] = sump - summ;
  }
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
  Double_t cut = threshold;
  for (unsigned ibin = 0; ibin < digi.size(); ++ibin)
  {
    Double_t u = double(ibin) * timeUnit;
    if (digi[ibin] < cut && digi[ibin + 1] > cut)
    {
      crossings.push_back(PUP);
      crossingBin.push_back(ibin + 1);
      crossingTime.push_back(u);
      if (verbose)
        printf(" PUP det %i  bin %i %f %f  \n", idet, ibin, digi[ibin], digi[ibin + 1]);
    }
  }
  if (verbose)
    printf(" findTresholdCrossings det %i  crossings %lu \n", idet, crossings.size());
}

void hitFinder::findDerivativeCrossings(Int_t idet)
{

  crossings.clear();
  crossingBin.clear();
  crossingTime.clear();
  unsigned vsize = ddigi.size();
  Double_t cut = tbrun->detList[idet]->sigma * threshold;
  unsigned step = 1;
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
        printf(" PUP det %i  bin %i %f %f  \n", idet, ibin, vi, vj);
      crossings.push_back(PUP);
      ctype = PUP;
      crossingBin.push_back(ibin + 1);
      crossingTime.push_back(u);
    }
    else if (vi > cut && vj < -cut)
    {
      if (verbose)
        printf(" UPDOWN det %i  bin %i %f %f  \n", idet, ibin, vi, vj);
      crossings.push_back(UPDOWN);
      ctype = UPDOWN;
      crossingBin.push_back(ibin + 1);
      crossingTime.push_back(u);
    }
    else if (vi > cut && vj < cut)
    {
      if (verbose)
        printf(" NUP det %i  bin %i %f %f  \n", idet, ibin, vi, vj);
      crossings.push_back(NUP);
      ctype = NUP;
      crossingBin.push_back(ibin + 1);
      crossingTime.push_back(u);
    }
    else if (vi < -cut && vj > cut)
    {
      if (verbose)
        printf(" DOWNUP det %i  bin %i %f %f  \n", idet, ibin, vi, vj);
      crossings.push_back(DOWNUP);
      ctype = DOWNUP;
      crossingBin.push_back(ibin + 1);
      crossingTime.push_back(u);
    }
    else if (vi < -cut && vj > -cut)
    {
      if (verbose)
        printf(" PDOWN det %i  bin %i %f %f  \n", idet, ibin, vi, vj);
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

  return;
}
// make peaks to zero of waveform from PUP crossing type
void hitFinder::makePeaks(int idet, std::vector<Double_t> v)
{
  peakList.clear();
  peakKind.clear();
  hEvCross[idet]->Reset("ICESM");
  // loop over crossings using only PUP
  for (int icross = 0; icross < crossings.size(); ++icross)
  {
    if (crossings[icross] != PUP)
      continue;
    // find local max
    unsigned imax = 0;
    double maxVal = -99999.;
    for (unsigned ibin = crossingBin[icross]; ibin < v.size(); ++ibin)
    {
      if (v[ibin] < maxVal)
        break;
      imax = ibin;
      maxVal = v[ibin];
    }

    // no local max found
    if (verbose)
      printf(" cross det %i  imax %u \n", idet, imax);
    if (imax >= v.size() - 1 || imax == 0)
      continue;

    // histogram local max
    hEvCross[idet]->SetBinContent(imax, v[imax]);

    // extend to baseline
    unsigned ilow = v.size();
    unsigned ihigh = 0;
    // look low
    for (unsigned ibin = imax; ibin > 0; --ibin)
    {
      if (v[ibin] < 0)
        break;
      ilow = ibin;
    }
    // look high
    for (unsigned ibin = imax; ibin < v.size(); ++ibin)
    {
      if (v[ibin] < 0)
        break;
      ihigh = ibin;
    }
    if (verbose)
      printf(" cross det %i  imax %i val %f icross %u from (%u,%u) \n", idet, imax, maxVal, crossingBin[icross], ilow, ihigh);
    //if (ihigh - ilow > maxPeakLength)
    //  ihigh = ilow + maxPeakLength;

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
        printf(" add det %i  imax %i val %f icross %u from (%u,%u) \n", idet, imax, maxVal, crossingBin[icross], ilow, ihigh);
      peakList.push_back(std::make_pair(ilow, ihigh));
      peakKind.push_back(0);
      // if(vChannel[idet]==6) printf(" %lli add det %i  imax %i val %f icross %u from (%u,%u) %i\n",theEvent, idet, imax, maxVal, crossingBin[icross], ilow, ihigh,icross);
    }
  }
}

hitMap hitFinder::makeHits(int idet, Double_t &triggerTime, Double_t &firstCharge)
{
  double sigma = tbrun->detList[idet]->sigma;
  if (verbose)
    printf(" makeHits: event %lli det %i sigma %f peakList size %lu \n", theEvent, idet, sigma, peakList.size());
  triggerTime = 1E9;
  firstCharge = 0;
  hitMap detHits;
  if (peakList.size() < 1)
    return detHits;
  Double_t qmax = 0;

  unsigned minLength = 3;
  if (peakList.size() < 1)
    return detHits;
  for (unsigned ip = 0; ip < peakList.size(); ++ip)
  {
    unsigned klow = std::get<0>(peakList[ip]);
    unsigned khigh = std::get<1>(peakList[ip]);
    if (verbose)
      printf(" event %lli det %i hit  %u (%u,%u) kind %i length %u \n", theEvent, idet, ip, klow, khigh, peakKind[ip], khigh - klow);
    hHitLength->Fill(khigh - klow + 1);
    if (khigh - klow + 1 < minLength)
    {
      continue;
    }
    Double_t qhit = 0;
    UInt_t peakt = 0;
    Double_t qpeak = 0;
    Double_t qsum = 0;
    TDetHit dhit;
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

    unsigned kstart=TMath::Max(unsigned(0), klow-20);
    unsigned kend=TMath::Min(unsigned(digi.size()), khigh+20);

    if(vChannel[idet]<9) for(unsigned k = kstart; k < kend;  ++k)
      dhit.digi.push_back(digi[k]);

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

    if (verbose)
      printf(" makeHits %llu insert hit idet %i  number %lu time %f (%u,%u) kind %i length %u qpeak %f digi size %lu  \n", theEvent, idet, detHits.size(), hitTime, dhit.firstBin, dhit.lastBin, peakKind[ip], khigh - klow + 1, qpeak, dhit.digi.size());
    // for (unsigned k=klow; k<khigh;  ++k) printf(" \t %u %f ; ", k, ddigi[k]);
    // cout << endl;
    detHits.insert(std::pair<Double_t, TDetHit>(hitTime, dhit));
    hPeakNWidth->Fill(dhit.lastBin - dhit.firstBin + 1);
  }

  // first time, charge from map
  /*
  hitMapIter hitIter;
  hitIter=detHits.begin();
  TDetHit dhit0 = hitIter->second;
  triggerTime = dhit0.startTime*microSec;
  firstCharge = dhit0.qsum;
  */
  // printf(" return from makeHits with %lu made \n",detHits.size());
  return detHits;
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
        printf(" PUP det %i  bin %i %f %f  \n", idet, ibin, vi, vj);
      peakCrossings.push_back(PUP);
      ctype = PUP;
      peakCrossingBin.push_back(ibin + 1);
      peakCrossingTime.push_back(u);
    }
    else if (vi > cut && vj < -cut)
    {
      if (verbose)
        printf(" UPDOWN det %i  bin %i %f %f  \n", idet, ibin, vi, vj);
      peakCrossings.push_back(UPDOWN);
      ctype = UPDOWN;
      peakCrossingBin.push_back(ibin + 1);
      peakCrossingTime.push_back(u);
    }
    else if (vi > cut && vj < cut)
    {
      if (verbose)
        printf(" NUP det %i  bin %i %f %f  \n", idet, ibin, vi, vj);
      peakCrossings.push_back(NUP);
      ctype = NUP;
      peakCrossingBin.push_back(ibin + 1);
      peakCrossingTime.push_back(u);
    }
    else if (vi < -cut && vj > cut)
    {
      if (verbose)
        printf(" DOWNUP det %i  bin %i %f %f  \n", idet, ibin, vi, vj);
      peakCrossings.push_back(DOWNUP);
      ctype = DOWNUP;
      peakCrossingBin.push_back(ibin + 1);
      peakCrossingTime.push_back(u);
    }
    else if (vi < -cut && vj > -cut)
    {
      if (verbose)
        printf(" PDOWN det %i  bin %i %f %f  \n", idet, ibin, vi, vj);
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
    // if (splitVerbose)
    // printf(" \t peak  = %u   max %f crossings %lu cut  %f  \n", ip, peakMax, peakCrossings.size(), cut);

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
          printf(" event %llu idet %i  peak %u npeaks = %lu peakThreshold %.2f ratio %f binAtr %f \n", theEvent, vChannel[idet], ip, peakList.size(), peakThreshold, ratio, binAtr);
          printf(" \t\t crossing  %i type %i bin %i peakStart %i  peakEnd %i ddigi %f digi %f max %f ratio to max %f  subBin %.0f \n", ipc, peakCrossings[ipc], peakCrossingBin[ipc], peakStart, peakEnd, ddigi[peakCrossingBin[ipc]], digi[peakCrossingBin[ipc]], peakMax, ratio, subBin);
        }

        // split peak at largest subBin
        if (subBin > binAtr && binAtr > 0 && ratio < 0.5)
        {
          isplit.push_back(peakCrossingBin[ipc]);
          indexSplit.push_back(ip);
          if (splitVerbose)
            printf("sssssssssss splitting peak number %u nsplits %lu \n", ip, isplit.size());
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
  for (unsigned  is = 0; is < nsamples; ++is)
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

  for (unsigned  i = 0; i < nsamples; ++i)
  {
    double rl = fInverseFFT->GetPointReal(i);
    Signal.push_back(rl);
  }

  // normalize
  for (unsigned i = 0; i < Signal.size(); ++i)
    Signal[i] /= double(nsamples);

  return Signal;
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

void hitFinder::plotSplitEvent(unsigned idet, Long64_t ievent)
{
  // evDir->cd();
  splitDir->cd();
  TString histName;
  TString detName = tbrun->detList[idet]->GetName();

  histName.Form("EvSplitWave%lli_%s", ievent, detName.Data());
  TH1D *hwave = (TH1D *)hEvWave[idet]->Clone(histName);
  hwave->SetTitle(histName);

  histName.Form("EvSplitPeakWave%lli_%s", ievent, detName.Data());
  TH1D *hhitPeakWave = (TH1D *)hEvHitPeakWave[idet]->Clone(histName);
  hhitPeakWave->SetTitle(histName);

  fout->cd();
}

void hitFinder::plotEvent(unsigned ichan, Long64_t ievent)
{
  int idet = chanMap.at(ichan);
  // evDir->cd();
  fftDir->cd();
  TString histName;
  TString detName = tbrun->detList[idet]->GetName();

  histName.Form("EvWave%lli_%s", ievent, detName.Data());
  TH1D *hwave = (TH1D *)hEvWave[idet]->Clone(histName);
  hwave->SetTitle(histName);

  /*
  // cout << " det " << idet << " "  << hwave->GetName() << " ," << hwave->GetTitle() << endl;
  histName.Form("EvSmooth%lli_%s", ievent, detName.Data());
  TH1D *hsmooth = (TH1D *)hEvSmooth[idet]->Clone(histName);
  hsmooth->SetTitle(histName);

  histName.Form("EvDerWave%lli_%s", ievent, detName.Data());
  TH1D *hdwave = (TH1D *)hEvDerWave[idet]->Clone(histName);
  hdwave->SetTitle(histName);

  histName.Form("EvCross%lli_%s", ievent, detName.Data());
  TH1D *hcross = (TH1D *)hEvCross[idet]->Clone(histName);

  histName.Form("EvPeakCross%lli_%s", ievent, detName.Data());
  TH1D *hpeakCross = (TH1D *)hEvPeakCross[idet]->Clone(histName);
  */

  histName.Form("EvFiltWave%lli_%s", ievent, detName.Data());
  TH1D *hfiltwave = (TH1D *)hEvFiltWave[idet]->Clone(histName);

 /*
  histName.Form("EvFFT%lli_%s", ievent, detName.Data());
  TH1D *hfft = (TH1D *)hFFT[idet]->Clone(histName);

  histName.Form("EvInvFFT%lli_%s", ievent, detName.Data());
  TH1D *hinvfft = (TH1D *)hInvFFT[idet]->Clone(histName);
  */

  histName.Form("EvHitWave%lli_%s", ievent, detName.Data());
  TH1D *hhitWave = (TH1D *)hEvHitWave[idet]->Clone(histName);

  /*
      histName.Form("EvHitPeakWave%lli_%s", ievent, detName.Data());
      TH1D *hhitPeakWave = (TH1D *)hEvHitPeakWave[idet]->Clone(histName);


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