////////////////////////////////////////////////////////
//  M.Gold June 2022
// revised Jan 27 2023 -- simplified finding
// class to make hits from vector data
//P. ugec et al. Pulse processing routines for neutron time-of-flight data. Nucl. Instrum. Meth., A812:134–144, 2016.

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

hitFinder::hitFinder(TFile *theFile, TBRun *brun, TString theTag, int nSamples, vector<int> vchan)
{
  verbose = false;
  if(verbose)
    cout << "INSTANCE OF HITFINDER " <<  " vchan.size " << vchan.size() <<endl;
  smoothing = false;
  fout = theFile;
  fftDir = fout->mkdir("fftDir");
  finderDir = fout->mkdir("finderDir");
  tag = theTag;
  tbrun = brun;
  nsamples = nSamples;
  // initialize fft
  fFFT = TVirtualFFT::FFT(1, &nsamples, "R2C M K");
  fInverseFFT = TVirtualFFT::FFT(1, &nsamples, "C2R M K");

  microSec = 1.0E-3;
  timeUnit = 8.0; // ns per count
  maxPeakLength = 100;
  thresholdStepSize=1;
  fout->cd();
  hPeakCount = new TH1D("PeakCount", " peaks by det ", vchan.size(),0,vchan.size());
  hHitLength = new TH1I("HitLength", " hit length", 1000, 0, 1000);
  hPeakNWidth = new TH1I("PeakNWidth", "PeakNWidth", 1000, 0, 1000);
  hPeakValue = new TH1D("PeakValue", "Peak value (not trigger)", 1000, 0, 5000);

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
    printf(" create  index %i vchan %i %s %s \n",index, id, hFFT[index]->GetName(),hFFT[index]->GetTitle());
  }

  fout->cd();
  finderDir->cd();
  for (unsigned index = 0; index < vchan.size(); ++index)
  {
    int id = vchan[index];
    TDet *deti = tbrun->getDet(id);
    hEvWave.push_back(new TH1D(Form("EvWave%s", deti->GetName()), Form("Wave%s", deti->GetName()), nsamples, 0, nsamples));
    hEvSmooth.push_back(new TH1D(Form("EvSmooth%s", deti->GetName()), Form("Smooth%s", deti->GetName()), nsamples, 0, nsamples));
    hEvCross.push_back(new TH1D(Form("EvCross%s", deti->GetName()), Form("Cross%s", deti->GetName()), nsamples, 0, nsamples));
    hEvDerWave.push_back(new TH1D(Form("EvDerWave%s", deti->GetName()), Form("DerWave%s", deti->GetName()), nsamples, 0, nsamples));
    hEvFiltWave.push_back(new TH1D(Form("EvFiltWave%s", deti->GetName()), Form("FiltWave%s", deti->GetName()), nsamples, 0, nsamples));
    hEvHitWave.push_back(new TH1D(Form("EvHitWave%s", deti->GetName()), Form("HitWave%s", deti->GetName()), nsamples, 0, nsamples));
    hDigiVal.push_back(new TH1D(Form("DigiVal%i", id), Form("digi value chan %id", id), 2000, -1000., 1000.));
    hEvWave[index]->SetDirectory(nullptr);
    hEvCross[index]->SetDirectory(nullptr);
    hEvSmooth[index]->SetDirectory(nullptr);
    hEvDerWave[index]->SetDirectory(nullptr);
    hEvHitWave[index]->SetDirectory(nullptr);
    hEvFiltWave[index]->SetDirectory(nullptr);
    hHitSum.push_back(new TH1D(Form("HitSum%s", deti->GetName()), Form("HitSum%s", deti->GetName()), nsamples, 0, nsamples));
    printf(" create  index %i vchan %i %s %s \n", index, id, hEvWave[index]->GetName(), hEvWave[index]->GetTitle());
  }
  fout->cd();
  ntFinder = new TNtuple("finder"," hit finder ","chan:ncross:npeak");

  gotTransforms = getTransforms();

  cout << " created hitFinder with " << tbrun->GetName() << " nsamples =  " << nsamples << " ndet " << hEvWave.size() << " gotTransforms= " << gotTransforms << endl;

  printf(" channel mapping \n");
  for (unsigned index  = 0; index < vchan.size(); ++index){
    int id = chanMap.at(vchan[index]);
    printf("index %i chan %i mapped to index  %i %s %s\n", index, vchan[index],id,
      hEvWave[id]->GetName(), hEvWave[id]->GetTitle());
  }
}

void hitFinder::event(int ichan, Long64_t ievent, vector<double> eventDigi,double thresh, unsigned step)
{
  /*if(ichan == 7 && ievent==0 )
    verbose = true;
    else
      verbose = false;
      */
    threshold = thresh;
    diffStep = step;
    int idet = chanMap.at(ichan);
    digi = eventDigi;
    if (verbose)
    printf(" HHHHHH hitFinder ievent %llu ichan %i idet %i threshold %.1f \n", ievent, ichan, idet, tbrun->detList[idet]->sigma * threshold);

    double triggerTime = 0;
    double firstCharge = 0;

    for (int i = 0; i < nsamples; ++i)
    {
    hDigiVal[idet]->Fill(digi[i]);
  }
  // FFT and filter
    fdigi.resize(digi.size());
    std::vector<std::complex<double>> complexTransform;
    complexTransform = FFT(idet, digi);
    unsigned maxFrequency = complexTransform.size();
    //unsigned maxN = gTransform[0]->GetN();
    double xpoint;
    double ypoint=1.0;
    for (unsigned iw = 1; iw < maxFrequency; ++iw)
    {
      int ipoint = iw - 1;
      if (iw > maxFrequency / 2)
        ipoint = maxFrequency - iw - 1;
      if (gotTransforms)
        gTransform[0]->GetPoint(ipoint, xpoint, ypoint);

      complexTransform[iw] *= ypoint;
    }
    fdigi = inverseFFT(idet, complexTransform, digi);

  for (unsigned isample = 0; isample < digi.size(); isample++)
  {
    hEvWave[idet]->SetBinContent(isample + 1, digi[isample]);
    hEvSmooth[idet]->SetBinContent(isample + 1, digi[isample]);
    if (gotTransforms)
      hEvFiltWave[idet]->SetBinContent(isample + 1, fdigi[isample]);
  }
  // smooth and fill vector
  hEvSmooth[idet]->Smooth(1); // one time
  sdigi.resize(digi.size());
  for (unsigned ibin = 1; ibin < hEvSmooth[idet]->GetNbinsX() ; ibin++)
    sdigi[ibin - 1] = hEvSmooth[idet]->GetBinContent(ibin);

  // use filtered waveforms
  if (gotTransforms)
    digi = fdigi;

  if(!smoothing)
    sdigi = digi;

  // 7;
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
  findThresholdCrossings(idet);
  makePeaks(idet, digi);
  //splitPeaks(idet, ddigi);
  detHits = makeHits(idet, triggerTime, firstCharge);
  hPeakCount->Fill(idet, peakList.size());
  // fill hits
  if (verbose)
    cout << " finished make hits event " << ievent << " chan " << ichan << " det " << idet
         << "  ddigi size " << ddigi.size()
         << "  crossings size " << crossings.size()
         << "  peakList size " << peakList.size()
         << "  detHits size " << detHits.size() << endl << endl;
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
    tbrun->detList[idet]->qSum += hiti.qsum;
    tbrun->detList[idet]->hitSum += hiti.qpeak;
    if(hiti.startTime<startTimeCut) {
      tbrun->detList[idet]->qPrompt += hiti.qsum;
      tbrun->detList[idet]->hitPrompt += hiti.qpeak;
    }
    if(!triggerChannel)
      hPeakValue->Fill(hiti.qpeak);
    hiti.SetTitle(Form("TDetHit %i event %llu chan %i index %i ", icount++, ievent, ichan, idet));
    tbrun->detList[idet]->hits.push_back(hiti);
    for (unsigned iv = 0; iv < digi.size(); ++iv)
      if (iv >= hiti.firstBin && iv <= hiti.lastBin){
        hdigi[iv] = digi[iv];
      }
  }
  tbrun->fill();

  for (unsigned isample = 0; isample < hdigi.size(); isample++){
    hEvHitWave[idet]->SetBinContent(isample + 1, hdigi[isample]);
    hHitSum[idet]->SetBinContent(isample + 1, hdigi[isample] + hHitSum[idet]->GetBinContent(isample + 1));
  }

  //
  ntFinder->Fill(float(idet),float(crossings.size()), float(peakList.size()));
  if(verbose) {
  for (unsigned idet = 0; idet < tbrun->detList.size(); ++idet)
  {
    TDet *tdet = tbrun->detList[idet];
    cout << tdet->channel << " hits.size " << tdet->hits.size() << endl;
    for (unsigned ihit = 0; ihit < tdet->hits.size(); ++ihit)
    {
      TDetHit thit = tdet->hits[ihit];
      cout << "finder hit number  " << ihit << " bin " << thit.peakBin << endl;
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
  ddigi[0]=0; // first entry is zero
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
    for (unsigned j = 1; j < maxSum; ++j)
      sump += digi[i + j];
    summ = 0;
    for (unsigned j = 1; j < maxSum; ++j)
      sump -= digi[i - j];
    ddigi[i]=sump - summ;
  }
}

// threshold crossings
void hitFinder::findThresholdCrossings(Int_t idet)
{
  crossings.clear();
  crossingBin.clear();
  crossingTime.clear();
  unsigned vsize = digi.size();
  Double_t cut = tbrun->detList[idet]->sigma * threshold;
  unsigned step = thresholdStepSize;
  unsigned ibin = 0;
  while (ibin < digi.size() - step)
  {
    Double_t vi = digi[ibin];
    Double_t vj = digi[ibin + step];
    Double_t u = double(ibin) * timeUnit;
    if (vi < cut && vj > cut)
    { 
      crossings.push_back(PUP);
      crossingBin.push_back(ibin+step);
      crossingTime.push_back(u);
      if (verbose)
        printf(" PUP det %i  bin %i %f %f  \n", idet, ibin, vi, vj);
    }
    ibin += step;
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
  for (unsigned ibin = 0; ibin < vsize-step;  ++ibin)
  {
    Double_t u = double(ibin) * timeUnit;
    Double_t vi = ddigi[ibin];
    Double_t vj = ddigi[ibin + step];
    unsigned ctype = 10;
    // crossing types
    if (vi < cut && vj > cut)
    {
      if(verbose)  printf(" PUP det %i  bin %i %f %f  \n", idet, ibin, vi, vj);
      crossings.push_back(PUP);
      ctype = PUP;
      crossingBin.push_back(ibin+1);
      crossingTime.push_back(u);
    }
    else if (vi > cut  && vj < -cut)
    {
      if(verbose)  printf(" UPDOWN det %i  bin %i %f %f  \n", idet, ibin, vi, vj);
      crossings.push_back(UPDOWN);
      ctype = UPDOWN;
      crossingBin.push_back(ibin+1);
      crossingTime.push_back(u);
    }
    else if (vi > cut && vj < cut)
    {
      if(verbose)  printf(" NUP det %i  bin %i %f %f  \n", idet, ibin, vi, vj);
      crossings.push_back(NUP);
      ctype = NUP;
      crossingBin.push_back(ibin+1);
      crossingTime.push_back(u);
    }
    else if (vi < -cut && vj > cut)
    {
      if(verbose)  printf(" DOWNUP det %i  bin %i %f %f  \n", idet, ibin, vi, vj);
      crossings.push_back(DOWNUP);
      ctype = DOWNUP;
      crossingBin.push_back(ibin+1);
      crossingTime.push_back(u);
    }
    else if (vi < -cut && vj > -cut)
    {
      if(verbose)  printf(" PDOWN det %i  bin %i %f %f  \n", idet, ibin, vi, vj);
      crossings.push_back(PDOWN);
      ctype = PDOWN;
      crossingBin.push_back(ibin+1);
      crossingTime.push_back(u);
    }
    else if (vi > -cut && vj < -cut)
    {
      crossings.push_back(NDOWN);
      ctype = NDOWN;
      crossingBin.push_back(ibin+1);
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
    if(imax >= v.size() -1 || imax ==0 )
      continue;
    
    // histogram local max
    hEvCross[idet]->SetBinContent(imax, v[imax]);

    // extend to baseline
    unsigned  ilow = v.size();
    unsigned  ihigh = 0;
    // look low
    for(unsigned ibin = imax; ibin > 0; --ibin) {
        if (v[ibin] < 0)
          break;
        ilow = ibin;
    }
    // look high
    for (unsigned ibin = imax; ibin< v.size() ; ++ ibin) {
        if (v[ibin] < 0)
        break;
        ihigh  = ibin;
    }
    if (verbose)
        printf(" cross det %i  imax %i val %f icross %u from (%u,%u) \n", idet, imax, maxVal, crossingBin[icross],ilow, ihigh );
    if (ihigh - ilow > maxPeakLength)
        ihigh = ilow + maxPeakLength;

    // check that this hit has not already been found.
    bool found = false;
    for (unsigned ip = 0; ip < peakList.size(); ++ip)
    {
        // trim first peak
        unsigned peakStart = std::get<0>(peakList[ip]);
        unsigned peakEnd = std::get<1>(peakList[ip]);
        if(peakStart==ilow && peakEnd == ihigh )
        found = true;
    }
    if(!found) {
        if (verbose)
        printf(" add det %i  imax %i val %f icross %u from (%u,%u) \n", idet, imax, maxVal, crossingBin[icross], ilow, ihigh);
        peakList.push_back(std::make_pair(ilow, ihigh));
        peakKind.push_back(0);
    }
  }
}

hitMap hitFinder::makeHits(int idet, Double_t &triggerTime, Double_t &firstCharge)
{
    double sigma = tbrun->detList[idet]->sigma;
    if (verbose)
    printf(" makeHits: det %i sigma %f peakList size %lu \n", idet, sigma, peakList.size());
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
        printf(" idet %i hit  %u (%u,%u) kind %i length %u \n", idet, ip, klow, khigh, peakKind[ip], khigh - klow + 1);
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

    TDetHit dhit;
    dhit.peakBin = Int_t(peakt);
    dhit.qsum = qsum;
    dhit.qpeak = qpeak;
    dhit.firstBin = klow;
    dhit.lastBin = khigh;
    dhit.peakMaxTime = peakt;
    dhit.peakt = peakt;
    dhit.startTime = klow;
    dhit.peakWidth = khigh - klow;
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
      printf(" makeHits insert hit idet %i  number %lu time %f (%u,%u) kind %i length %u qpeak %f  \n", idet, detHits.size(), hitTime, dhit.firstBin, dhit.lastBin, peakKind[ip], khigh - klow + 1, qpeak);
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

// split peak based on derivaive
void hitFinder::splitPeaks(int idet, std::vector<Double_t> v)
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
      double vp =  v[kp];
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

std::vector<std::complex<double>> hitFinder::FFT(Int_t id, std::vector<double> rdigi, bool first)
{
  std::vector<std::complex<double>> VectorComplex;
  for (int is = 0; is < nsamples; ++is)
    fFFT->SetPoint(is, rdigi[is]);
  fFFT->Transform(); //
  // cout << " FFT histogram " << hFFT[id]->GetName() << endl;

  std::vector<Double_t> realVec, imVec;
  for (int i = 0; i < nsamples; ++i)
  {
    double rl, im;
    fFFT->GetPointComplex(i, rl, im);
    std::complex<double> c(rl, im); //.real or .imag accessors
    VectorComplex.push_back(c);
    // skip first bin which is pedestal
    if (i < nsamples / 2)
    {
      if (first)
        hFFT[id]->SetBinContent(i, hFFT[id]->GetBinContent(i) + std::abs(c));
      hFFTFilt[id]->SetBinContent(i, hFFTFilt[id]->GetBinContent(i) + std::abs(c));
    }

    realVec.push_back(VectorComplex[i].real());
    imVec.push_back(VectorComplex[i].imag());
  }
  return VectorComplex;
}

std::vector<Double_t> hitFinder::inverseFFT(Int_t id, std::vector<std::complex<double>> VectorComplex, std::vector<double> rdigi)
{
  std::vector<Double_t> Signal;
  double sum = 0;
  for (int is = 0; is < nsamples; ++is)
  {
    fInverseFFT->SetPoint(is, VectorComplex[is].real(), VectorComplex[is].imag());
    sum += rdigi[is];
  }
  fInverseFFT->Transform();
  // cout << " FFT histogram " << hInvFFT[id]->GetName() << " rdigi sum " << sum << endl;

  Double_t norm = 0;
  for (int i = 0; i < nsamples; ++i)
  {
    double rl = fInverseFFT->GetPointReal(i);
    norm += rl;
    Signal.push_back(rl);
  }
  for (int i = 0; i < nsamples; i++)
  {
    Signal[i] = Signal[i] * sum / norm;
    hInvFFT[id]->SetBinContent(i + 1, hInvFFT[id]->GetBinContent(i + 1) + Signal[i]);
  }
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

void hitFinder::plotEvent(unsigned ichan, Long64_t ievent)
{
  int idet = chanMap.at(ichan);
  // evDir->cd();
  fftDir->cd();
  TString histName;
  TString detName = tbrun->detList[idet]->GetName();

  histName.Form("EvWave%lli_%s", ievent, detName.Data());
  TH1D *hwave = (TH1D *)hEvWave[idet]->Clone(histName);
  //cout << " det " << idet << " "  << hwave->GetName() << " ," << hwave->GetTitle() << endl;
  histName.Form("EvSmooth%lli_%s", ievent, detName.Data());
  TH1D *hsmooth = (TH1D *)hEvSmooth[idet]->Clone(histName);

  histName.Form("EvDerWave%lli_%s",ievent, detName.Data());
  TH1D *hdwave = (TH1D *)hEvDerWave[idet]->Clone(histName);
  hdwave->SetTitle(histName);

  histName.Form("EvCross%lli_%s", ievent, detName.Data());
  TH1D *hcross = (TH1D *)hEvCross[idet]->Clone(histName);

  histName.Form("EvFiltWave%lli_%s", ievent,detName.Data());
  TH1D *hfiltwave = (TH1D *)hEvFiltWave[idet]->Clone(histName);

  histName.Form("EvFFT%lli_%s", ievent,detName.Data());
  TH1D *hfft = (TH1D *)hFFT[idet]->Clone(histName);

  histName.Form("EvInvFFT%lli_%s", ievent, detName.Data());
  TH1D *hinvfft = (TH1D *)hInvFFT[idet]->Clone(histName);

  histName.Form("EvHitWave%lli_%s", ievent, detName.Data());
  TH1D *hhitWave = (TH1D *)hEvHitWave[idet]->Clone(histName);

  /*
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
bool hitFinder::getTransforms()
{
  TGraph *gw = NULL;
  TString fileName = TString("../bobj/WienerTransforms.root");
  TFile *f1 = new TFile(fileName, "readonly");
  if (f1->IsZombie())
  {
    printf(" no  file for %s \n", tag.Data());
    return false;
  }

  bool val = true;
  int idet = 0;
  while (val)
  {
    gTransform.clear();
    f1->GetObject(Form("WienerTransDet%i", idet), gw);
    if (!gw)
    {
      val = false;
      break;
    }
    gTransform.push_back((TGraph *)gw->Clone(Form("WeinerForDet%i-%s", idet, tag.Data())));
    unsigned itran = gTransform.size()-1;
    printf(" got transform %s points %i for run %s \n", gTransform[itran]->GetName(), gTransform[itran]->GetN(), tag.Data());
    fftDir->Append(gTransform[itran]);
  }
  if(val) printf(" got transforms for file %s %i \n", fileName.Data(), val);
  return val;
}