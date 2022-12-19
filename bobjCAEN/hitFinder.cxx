////////////////////////////////////////////////////////
//  M.Gold June 2022
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
  fout = theFile;
  diffStep = 3;
  fftDir = fout->mkdir("fftDir");
  tag = theTag;
  tbrun = brun;
  nsamples = nSamples;
  // initialize fft
  fFFT = TVirtualFFT::FFT(1, &nsamples, "R2C M K");
  fInverseFFT = TVirtualFFT::FFT(1, &nsamples, "C2R M K");

  microSec = 1.0E-3;
  timeUnit = 1.0; // ns per count

  fout->cd();
  hPeakCount = new TH1D("PeakCount", " peaks by det ", vchan.size(),0,vchan.size());
  hHitLength = new TH1I("HitLength", " hit length", 100, 0, 100);
  hPeakNWidth = new TH1I("PeakNWidth", "PeakNWidth", 100, 0, 100);
  // tbrun->print();

  for (unsigned index = 0; index < vchan.size(); ++index)
  {
    int id = vchan[index];
    TDet *deti = tbrun->getDet(id);
    chanMap.insert(std::pair<int, int>(id, index));
    hFFT.push_back(new TH1D(Form("FFTDET%i", id), Form("FFT Channel %i ", id), nsamples / 2, 0, nsamples / 2));
    hInvFFT.push_back(new TH1D(Form("InvFFTDET%i", id), Form("Inverse FFT Channel %i ", id), nsamples, 0, nsamples));
    hEvWave.push_back(new TH1D(Form("EvWave%s", deti->GetName()), Form("Wave%s", deti->GetName()), nsamples, 0, nsamples));
    hEvSmooth.push_back(new TH1D(Form("EvSmooth%s", deti->GetName()), Form("Smooth%s", deti->GetName()), nsamples, 0, nsamples));
    hEvCross.push_back(new TH1D(Form("EvCross%s", deti->GetName()), Form("Cross%s", deti->GetName()), nsamples, 0, nsamples));
    hEvDerWave.push_back(new TH1D(Form("EvDerWave%s", deti->GetName()), Form("DerWave%s", deti->GetName()), nsamples, 0, nsamples));
    hEvFiltWave.push_back(new TH1D(Form("EvFiltWave%s", deti->GetName()), Form("FiltWave%s", deti->GetName()), nsamples, 0, nsamples));
    hEvHitWave.push_back(new TH1D(Form("EvHitWave%s", deti->GetName()), Form("HitWave%s", deti->GetName()), nsamples, 0, nsamples));
    hDigiVal.push_back(new TH1D(Form("DigiVal%i",id),Form("digi value chan %id",id),2000,-1000.,1000.));
    hFFT[index]->SetDirectory(nullptr);
    hInvFFT[index]->SetDirectory(nullptr);
    hEvWave[index]->SetDirectory(nullptr);
    hEvCross[index]->SetDirectory(nullptr);
    hEvSmooth[index]->SetDirectory(nullptr);
    hEvDerWave[index]->SetDirectory(nullptr);
    hEvHitWave[index]->SetDirectory(nullptr);
    hEvFiltWave[index]->SetDirectory(nullptr);
    hHitSum.push_back(new TH1D(Form("HitSum%s", deti->GetName()), Form("HitSum%s", deti->GetName()), nsamples, 0, nsamples));
    fftDir->cd();
    hFFTFilt.push_back(new TH1D(Form("FFTFiltDET%i", id), Form("filtered FFT Channel %i ", id), nsamples / 2, 0, nsamples / 2));
    // hFFTFilt[index]->SetDirectory(nullptr);
    fout->cd();
    printf(" create  index %i vchan %i %s %s \n",index, id, hEvWave[index]->GetName(),hEvWave[index]->GetTitle());
  }

  gotTransforms = getTransforms();

  cout << " created hitFinder with " << tbrun->GetName() << " nsamples =  " << nsamples << " ndet " << hEvWave.size() << " gotTransforms= " << gotTransforms << endl;

  printf(" channel mapping \n");
  for (unsigned index  = 0; index < vchan.size(); ++index){
    int id = chanMap.at(vchan[index]);
    printf("index %i chan %i mapped to index  %i %s %s\n", index, vchan[index],id,
      hEvWave[id]->GetName(), hEvWave[id]->GetTitle());
  }
}

void hitFinder::event(int ichan, Long64_t ievent, vector<double> eventDigi)
{

  int idet = chanMap.at(ichan);
  digi = eventDigi;
  printf( " HHHHHH hitFinder ievent %llu ichan %i idet %i ndigi %lu \n",ievent, ichan,idet , digi.size() );

  double triggerTime = 0;
  double firstCharge = 0;

  for (int i = 0; i < nsamples; ++i) {
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

  // 7;
  differentiate(diffStep);
  for (unsigned isample = 0; isample < ddigi.size(); isample++)
  {
    hEvDerWave[idet]->SetBinContent(isample + 1, ddigi[isample]);
  }
  findHits(idet, ievent);

  detHits = makeHits(idet, triggerTime, firstCharge);
  // fill hits
  cout << " made hits " << ievent << " det " << idet << " << detHits size " << detHits.size() << endl;
  hdigi.clear();
  hdigi.resize(digi.size());

//  for (const auto &[key, value] : m)
//    std::cout << '[' << key << "] = " << value << "; "
  int icount = 0;
  for (hitMapIter hitIter = detHits.begin(); hitIter != detHits.end(); ++hitIter)
  {
    TDetHit hiti = hitIter->second;
    tbrun->detList[idet]->qSum += hiti.qsum;
    //if (hiti.startTime > 720 && hiti.startTime > 770)
    //  tbrun->detList[idet]->hitPrompt += hiti.qsum;
    //
    hiti.SetTitle(Form("TDetHit %i event %llu chan %i index %i ",icount++, ievent, ichan, idet));
    //printf(" new hit title = %s \n", hiti.GetTitle());
    tbrun->detList[idet]->hits.push_back(hiti);
    // Double_t phitQErr = phiti.qerr*timeUnit*1E9;
    // fill hdigi vector
    //cout << hiti.firstBin << " , " << hiti.lastBin << endl;
    for (unsigned iv = 0; iv < digi.size(); ++iv)
      if (iv >= hiti.firstBin && iv <= hiti.lastBin){
        hdigi[iv] = digi[iv];
      }
  }

  for (unsigned isample = 0; isample < hdigi.size(); isample++){
    hEvHitWave[idet]->SetBinContent(isample + 1, hdigi[isample]);
    hHitSum[idet]->SetBinContent(isample + 1, hdigi[isample] + hHitSum[idet]->GetBinContent(isample + 1));
  }
  //
  tbrun->fill();
  printf(" HHHHHH hitFinder finished \n");
}

void hitFinder::findHits(int idet, Long64_t jentry)
{

  // find peaks
  // for derivativePeaks, window in time is timeUnit*windowSize (ns) . timeUnit = 2
  // min, max width in time bins for simple peaks
  Int_t windowSize = 10;
  unsigned maxWidth = 100000;
  unsigned minWidth = 10;
  printf(" call derivative peaks idet % i \n ",idet);
  derivativePeaks(idet, 10.0);
  hPeakCount->Fill(idet, peakList.size());
}

/*
   void hitFinder::getAverage(std::vector<Double_t> digi, Double_t& ave, Double_t& sigma, unsigned max)
   {

   if(max==0) max = digi.size();
   double sum=0;
  double sum2=0;
  for (unsigned is=0; is< max; ++is) {
    sum += digi[is];
    sum2 += pow(digi[is],2.);
  }
  ave = sum/double(max);
  sigma = sqrt( sum2/double(max) - pow(ave,2.));
}
*/

// revised derivative Dec 8 2022 MG
void hitFinder::differentiate(unsigned diffStep)
{
  ddigi.clear();
  Double_t sump = 0;
  Double_t summ = 0;
  ddigi.push_back(0); // first entry is zero
  for (unsigned i = 1; i < nsamples; ++i)
  {
    // sum limit 
    int maxSum = diffStep;
    if(i<diffStep)
      maxSum = i;
    if(nsamples -1 -i < diffStep)
      maxSum = nsamples - 1 - i;
    //
    sump = 0;
    for (unsigned j = 1; j < maxSum; ++j)
      sump += sdigi[i + j];
    summ = 0;
    for (unsigned j = 1; j < maxSum; ++j)
      sump -= sdigi[i - j];
    ddigi.push_back(sump - summ);
  }
}

void hitFinder::derivativePeaks(Int_t idet, Double_t rms)
{
  peakList.clear();
  peakKind.clear();
  std::vector<unsigned> crossings;
  std::vector<unsigned> crossingBin;
  std::vector<double> crossingTime;
  unsigned vsize = ddigi.size();
  Double_t cut = tbrun->detList[idet]->sigma * rms;
  unsigned step = 1;
  // cout << " for det " << idet << " in derivative peaks >>>> rms " << rms << " cut " << cut << endl;
  Double_t ncut = -cut;
  // find all crossings
  for (unsigned ibin = step; ibin < vsize;  ++ibin)
  {
    Double_t u = double(ibin) * timeUnit * microSec;
    Double_t vi = ddigi[ibin];
    Double_t vj = ddigi[ibin - step];
    unsigned ctype = 10;
    //if (idet == 5)
     //printf(" det %i  bin %i %f %f  \n", idet, ibin, vj, vi);
    if (vj > cut && vi < ncut)
    {
      crossings.push_back(DOUBLEUPCROSS);
      ctype = DOUBLEUPCROSS;
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    }
    else if (vj < ncut && vi > cut)
    {
      crossings.push_back(DOUBLEDOWNCROSS);
      ctype = DOUBLEDOWNCROSS;
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    }
    else if (vi > cut && vj < cut)
    {
      crossings.push_back(UPCROSS);
      ctype = UPCROSS;
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    }
    else if (vi < cut && vj > cut)
    {
      crossings.push_back(UPCROSS);
      ctype = UPCROSS;
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    }
    else if (vi < ncut && vj > ncut)
    {
      crossings.push_back(DOWNCROSS);
      ctype = DOWNCROSS;
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    }
    else if (vi > ncut && vj < ncut)
    {
      crossings.push_back(DOWNCROSS);
      ctype = DOWNCROSS;
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    }
    // if (idet==5&&ctype<10)  printf("....... %u vj %f vi %f cut %f cross type %u \n", ibin, vj, vi, cut, ctype );
    // if (idet==1&&ibin>2350&&ibin<2450)  printf("\t %u vj %f vi %f ctype %u  \n", ibin, vj, vi, ctype );
  }

  if (crossings.size() < 4)
    return;
  ;

  // label found crossings, intially all false
  std::vector<bool> crossingFound;
  crossingFound.resize(crossings.size());
  for (unsigned jc = 0; jc < crossings.size(); ++jc)
  {
    //printf(" det %i crossing %i bin %i time %f type %i \n",idet,jc,crossingBin[jc],crossingTime[jc],crossings[jc]);
    crossingFound[jc] = false;
  }

  // histogram crossings
  for (int icross = 0; icross < crossings.size(); ++icross) {
    double crossVal = 0;
    if (crossings[icross] == UPCROSS)
      crossVal = 0.4 * sdigi[crossingBin[icross]-step];
    else if(crossings[icross] == DOWNCROSS)
      crossVal = -0.4 * sdigi[crossingBin[icross]-step];
    else if(crossings[icross] == DOUBLEUPCROSS)
      crossVal = 0.9*sdigi[crossingBin[icross-step]];
    else if(crossings[icross] == DOUBLEDOWNCROSS)
      crossVal = -0.9*sdigi[crossingBin[icross-step]];
    //
    hEvCross[idet]->SetBinContent(crossingBin[icross], crossVal);
  }
  // parse crossings to make pairs;
  /* first find sequence of 0,2,x x>0 crossings */
  unsigned ip = 0;
  while (ip <= crossings.size() - 3)
  {
    int ibin = crossingBin[ip];
    if (crossings[ip] == UPCROSS && crossings[ip + 1] == DOUBLEUPCROSS && crossings[ip + 2] > UPCROSS)
    {
    if (idet==5)  printf("\t peak %lu ibin %i type (%i,%i,%i) \n", peakList.size() ,ibin, crossings[ip],crossings[ip+1],crossings[ip+2] );
      peakList.push_back(std::make_pair(crossingBin[ip], crossingBin[ip + 2]));
      peakKind.push_back(0);
      /*
      std::pair<unsigned,unsigned>  pp = std::make_pair(crossingBin[ip], crossingBin[ip + 2]);
      unsigned ip = peakList.size() - 1;
      printf(" derivativePeaks ip = %lu \n", peakList.size());
      unsigned klow = pp.first;
      unsigned khigh = pp.second;
      printf(" derivativePeaks   (%u,%u) \n ", klow, khigh);
      printf(" det %i make peak  (%i,%i) kind %i  \n", idet, crossingBin[ip], crossingBin[ip + 2], peakKind[peakKind.size() - 1]);
      */
      // ntDer->Fill(rms, v[crossingBin[ip]], double(crossingBin[ip + 2] - crossingBin[ip]), double(0)); // sigma:d0:step:dstep
      crossingFound[ip] = true;
      crossingFound[ip + 1] = true;
      crossingFound[ip + 2] = true;
      ip = ip + 3;
    }
    else if (crossings[ip] == UPCROSS && crossings[ip + 1] == UPCROSS )
    {
      crossingFound[ip] = true;
      crossingFound[ip + 1] = true;
      ip = ip + 2;
      peakList.push_back(std::make_pair(crossingBin[ip], crossingBin[ip + 1]));
      peakKind.push_back(1);
    } else 
        ++ip;
    }
  return;
}

hitMap hitFinder::makeHits(int idet, Double_t &triggerTime, Double_t &firstCharge)
{
  double sigma = tbrun->detList[idet]->sigma;
  printf(" makeHits: det %i sigma %f \n",idet, sigma );
  triggerTime = 1E9;
  firstCharge = 0;
  hitMap detHits;
  if (peakList.size() < 1)
    return detHits;
  Double_t qmax = 0;

  unsigned minLength = 3;
  cout << " ..... makeHits peak size " << peakList.size() << " peakKind  " << peakKind.size() << endl;
  if (peakList.size() < 1)
    return detHits;
  for (unsigned ip = 0; ip < peakList.size(); ++ip)
  {
    unsigned klow = std::get<0>(peakList[ip]);
    unsigned khigh = std::get<1>(peakList[ip]);
    printf(" hit  %u (%u,%u) kind %i length %u \n", ip, klow, khigh, peakKind[ip], khigh - klow + 1);
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
      double qdigik =  digi[k];
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
    printf("  insert hit  %lu time %f (%u,%u) kind %i length %u qsum %f  \n", detHits.size(), hitTime, dhit.firstBin, dhit.lastBin, peakKind[ip], khigh - klow + 1,qsum);
    //for (unsigned k=klow; k<khigh; ++k) printf(" \t %u %f ; ", k, ddigi[k]);
    //cout << endl;
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

// extend peaks to zero of waveform
void hitFinder::extendPeaks(int idet, std::vector<Double_t> v)
{
  if (peakList.size() < 1)
    return;
  for (unsigned ip = 0; ip < peakList.size(); ++ip)
  {
    // high direction
    unsigned high = std::get<1>(peakList[ip]);
    unsigned next = v.size();
    if (ip < peakList.size() - 1)
      next = std::get<1>(peakList[ip + 1]);
    for (unsigned kp = high; kp < next; ++kp)
    {
      double vp =  v[kp];
      if (vp < 0)
        break;
      std::get<1>(peakList[ip]) = kp;
    }

    // low direction
    unsigned low = std::get<0>(peakList[ip]);
    unsigned prev = 1;
    if (ip > 0)
      prev = TMath::Max(prev, std::get<0>(peakList[ip - 1]));
    for (unsigned kp = low; kp > prev; --kp)
    {
      double vp = v[kp];
      if (vp < 0)
        break;
      std::get<0>(peakList[ip]) = kp;
    }
    // printf("\t  extend peak %u from (%u,%u) to  (%u,%u)  \n",ip,low,high,std::get<0>(peakList[ip]),std::get<1>(peakList[ip])    );
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

  histName.Form("EvDer%uWave%lli_%s", diffStep, ievent, detName.Data());
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