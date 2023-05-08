//////////////////////////////////////////////////////////
//  M.Gold May 2022
//////////////////////////////////////////////////////////
#include <sstream>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <complex> //includes std::pair, std::make_pair
#include <valarray>
#include <numeric>
#include <algorithm> // std::sort
// root
#include <TROOT.h>
#include <TKey.h>
#include <TBranch.h>
#include <TBranchElement.h>
#include <TVirtualFFT.h>
#include <TChain.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TBranchElement.h>
#include <TFile.h>
#include <Rtypes.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFormula.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
// bobj classes
#include "TBWave.hxx"
#include "TBRun.hxx"
#include "TBRawEvent.hxx"
#include "hitFinder.hxx"
#include "TBFile.hxx"

class anaRun
{
public:
  enum
  {
    UPCROSS,
    DOWNCROSS,
    DOUBLEUPCROSS,
    DOUBLEDOWNCROSS
  };
  TFile *fout;
  TFile *fin;
  TDirectory *fftDir;
  vector<TTree *> treeList;
  vector<int> chanList;
  std::map<int, int> chanMap;
  vector<TBRawEvent *> rawBr;
  TBRun *tbrun;
  TNtuple *ntChan;
  TNtuple *ntChanSum;
  vector<TH1D *> noiseHist;
  vector<TH1D *> skewHist;
  vector<TH1D *> sumWave;
  vector<TH1D *> sumHitWave;
  vector<TH1D *> valHist;
  vector<TH1D *> sumWaveB;
  vector<TH1D *> valHistB;
  vector<TH1D *> threshHist;
  vector<TH1D *> crossHist;
  vector<TH1D *> hQSum;
  vector<TH1D *> hQPeak;
  vector<TH1D *> hQPEShape;
  TH1D *hEvBaseWave;
  vector<TH1D *> hEvGaus;
  TH1D *evCount;
  TH1D *histQSum;
  TH1D *histHitCount;
  TH1D *hEventPass;
  // TH1D *histQPE;
  TH1D *histQPrompt;
  vector<double> digi;
  vector<double> ddigi;
  vector<double> hdigi;
  std::vector<unsigned> thresholds;
  std::vector<unsigned> crossings;
  std::vector<unsigned> crossingBin;
  std::vector<double> crossingTime;
  vector<double> slope;
  vector<double> eslope;
  vector<double> chan;
  vector<double> echan;
  ofstream dumpFile;
  // vector<TBWave *> waveList;
  hitFinder *finder;
  TString tag;
  int currentBuffer;
  Long64_t currentBufferCount;
  anaRun(TString theTag = TString("dirName"));
  ~anaRun(){}
  Long64_t anaRunFile(TString theFile, Long64_t maxEntries);
  void clear();
  bool openFile(TString fileName);
  unsigned getListOfFiles(TString dir);
  unsigned getTrees();
  void getSummedHists();
  void getBranches();
  void getEvent(Long64_t entry);
  bool anaEvent(Long64_t entry);
  void differentiate(unsigned diffStep);          // not used
  void derivativeCount(TDet *idet, Double_t rms); // not used
  void negativeCrossingCount(double thresh);
  void thresholdCrossingCount(double thresh);

  TDirectory *rawSumDir;
  TDirectory *badDir;
  TDirectory *badTrigDir;
  TDirectory *evDir;
  TDirectory *sumDir;
  TDirectory *anaDir;
  Long64_t nentries;
  double QPEPeak;
};

void anaRun::clear()
{
  hQSum.clear();
  hQPEShape.clear();
  hQPeak.clear();
  treeList.clear();
  chanList.clear();
  chanMap.clear();
  rawBr.clear();
  noiseHist.clear();
  skewHist.clear();
  sumWave.clear();
  sumHitWave.clear();
  valHist.clear();
  sumWaveB.clear();
  valHistB.clear();
  threshHist.clear();
  crossHist.clear();
  hEvGaus.clear();
  digi.clear();
  ddigi.clear();
  hdigi.clear();
  thresholds.clear();
  crossings.clear();
  crossingBin.clear();
  crossingTime.clear();
  slope.clear();
  eslope.clear();
  chan.clear();
  echan.clear();
}

bool anaRun::openFile(TString theFile)
{
  // open input file and make some histograms
  TString fileName;
  fileName.Form("rootData/%s", theFile.Data());
  printf(" looking for file %s\n", fileName.Data());

  bool exists = false;
  FILE *aFile;
  aFile = fopen(fileName.Data(), "r");
  if (aFile)
  {
    fclose(aFile);
    exists = true;
  }
  if (!exists)
  {
    printf(" couldnt open file %s\n", fileName.Data());
    return false;
  }
  fin = new TFile(fileName, "readonly");
  printf(" opened file %s\n", fileName.Data());
  return true;
}

unsigned anaRun::getTrees()
{
  treeList.clear();
  chanList.clear();
  chanMap.clear();
  //
  TIter next(fin->GetListOfKeys());
  TKey *key;
  while (TKey *key = (TKey *)next())
  {
    TClass *cl = gROOT->GetClass(key->GetClassName());

    if (!cl->InheritsFrom("TTree"))
      continue;
    TTree *tree = (TTree *)key->ReadObj();
    // cout << tree->GetName() << " " << key->GetCycle() << endl;
    if (tree->GetEntries() < 1)
      continue;
    // pick up only last in cycle
    int ichan = TString(TString(tree->GetName())(5, 2)).Atoi();
    if (find(chanList.begin(), chanList.end(), ichan) != chanList.end())
      continue;
    cout << "tree name " << tree->GetName() << " chan " << ichan << endl;
    treeList.push_back(tree);
    chanList.push_back(ichan);
    chanMap.insert(std::pair<int, int>(ichan, chanList.size() - 1));
  }
  return treeList.size();
}
// get summed histos
void anaRun::getSummedHists()
{
  rawSumDir = fout->mkdir("rawSumDir");
  rawSumDir->cd();
  TIter next(fin->GetListOfKeys());
  TKey *key;
  while (TKey *key = (TKey *)next())
  {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1D"))
      continue;
    TH1D *h = (TH1D *)key->ReadObj();
    TString name;
    name.Form("SumWave-%s-%s", h->GetName(), tag.Data());
    TH1D *hsave = (TH1D *)h->Clone(name);
    // rawSumDir->Add(hsave);
  }
  fout->cd();
  return;
}

/* get rawBr */
void anaRun::getBranches()
{
  for (unsigned it = 0; it < treeList.size(); ++it)
  {
    rawBr[it] = NULL;
    TObjArray *brList = treeList[it]->GetListOfBranches();
    TIter next(brList);
    TBranchElement *aBranch = NULL;
    while ((aBranch = (TBranchElement *)next()))
    {
      // cout << "tree " << it << " " << aBranch->GetName() << endl;
      TBRawEvent *chan = (TBRawEvent *)aBranch->GetObject();
      if (!chan)
      {
        cout << " \t\t !!!! NULL " << it << endl;
        continue;
      }
      rawBr[it] = chan;
      // cout << "rawBr name  ==>  " << rawBr[it]->GetName() << endl;
    }
  }
}

void anaRun::getEvent(Long64_t entry)
{
  for (unsigned it = 0; it < treeList.size(); ++it)
    treeList[it]->GetEntry(entry);
  getBranches();
}

/* analyze rawBr */
bool anaRun::anaEvent(Long64_t entry)
{
  QPEPeak = 100;
  tbrun->clear();
  // loop over channels
  for (int ib = 0; ib < rawBr.size(); ++ib)
  {
    unsigned ichan = chanList[ib];
    bool trig = ichan == 9 || ichan == 10 || ichan == 11;
    int nbins = rawBr[ib]->rdigi.size();
    //cout << ib << " nbins " << nbins << " max hist " << hEvGaus.size() << " rawBr.size() " << rawBr.size() << endl;
    if (nbins != 1024)
      continue;

    // simple baseline
    std::vector<unsigned short> orderDigi = rawBr[ib]->rdigi;
    std::sort(orderDigi.begin(), orderDigi.end());
    unsigned baseLength = orderDigi.size() / 4;
    double base = 0;
    for (unsigned j = 0; j < baseLength; ++j)
    {
      base += orderDigi[j];
    }
    base /= double(baseLength);

    // baseline correction from fitted Gaussian
    hEvGaus[ib]->Reset("ICES");
    for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
      hEvGaus[ib]->Fill((double)rawBr[ib]->rdigi[j] - base);
    hEvGaus[ib]->Fit("gaus", "QO", "", hEvGaus[ib]->GetMean() - 100, hEvGaus[ib]->GetMean() + 100);
    TF1 *gfit = (TF1 *)hEvGaus[ib]->GetListOfFunctions()->FindObject("gaus");
    double ave = hEvGaus[ib]->GetMean();
    double sigma = hEvGaus[ib]->GetRMS();
    double skew = 0;
    double gausFitMean = 0;
    if (!isnan(hEvGaus[ib]->GetSkewness()))
      skew = hEvGaus[ib]->GetSkewness();
    if (gfit != nullptr)
    {
      ave = gfit->GetParameter(1);
      gausFitMean = ave;
      sigma = gfit->GetParameter(2);
    }

    evDir->cd();
    if (evDir->GetList()->GetEntries() < 100)
    {
      TH1D *EvGaussEvent = (TH1D *)hEvGaus[ib]->Clone(Form("EvGaussEvent%lld-Ch%i", entry, ichan));
    }
    fout->cd();

    // fitptr->Print();
    noiseHist[ib]->Fill(sigma);
    skewHist[ib]->Fill(skew);
    double sign = TMath::Sign(1., skew);
    TDet *idet = tbrun->getDet(ichan);
    if (idet == NULL)
    {
      printf("!!!!!NULL idet br %u ichan %i\n", ib, ichan);
      continue;
    }

    idet->ave = ave;
    idet->sigma = sigma;
    idet->skew = skew;
    idet->event = entry;
    idet->trigger = rawBr[ib]->trigger;
    idet->base = base+gausFitMean;

    /* fill baseline subtracted vector */
    digi.clear();
    for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
    {
      double val = double(rawBr[ib]->rdigi[j]) - idet->base;
      digi.push_back(val);
    }

    // find first peak maximum
    double peakMax = 0;
    for (unsigned j = 0; j < digi.size(); ++j)
    {
      if (digi[j] > peakMax)
        peakMax = digi[j];
    }

    // some sums
    double sum2 = 0;
    double sum = 0;
    unsigned peakWidth = 100;
    for (unsigned j = 0; j < digi.size(); ++j)
    {
      if (j > peakWidth)
        sum2 += digi[j];
      sum += digi[j];
    }
    sum2 /= double(digi.size() - peakWidth);
    sum /= double(digi.size());

   
    idet->peak = peakMax;
    idet->sum2 = sum2;
    idet->sum = sum;

    //cout << "tree " << ib << " ch " << ichan << " "
    //     << " ave " << ave << " sigma " << sigma << " skew " << skew << " base " << base << endl;


    // these 2 functions use didi and crossings,threshold vectors
    // i admit ugly
    negativeCrossingCount(10.);
    idet->crossings = int(crossings.size());
    crossHist[ib]->Fill(idet->crossings);
    if (trig)
      thresholdCrossingCount(15.);
    else
      thresholdCrossingCount(1.);
    idet->thresholds = int(thresholds.size());
    threshHist[ib]->Fill(idet->thresholds);

    ntChan->Fill(float(rawBr[ib]->trigger), float(ichan), float(ave), float(sigma), float(skew), float(base), float(peakMax), float(sum2), float(sum), float(crossings.size()), float(thresholds.size()));
    // set cut
    unsigned maxCrossings = 1;
    unsigned maxThresholds = 1;
    idet->pass = true;
    if (idet->crossings > maxCrossings)
      idet->pass = false;
    if (trig && idet->thresholds > maxThresholds)
      idet->pass = false;

      // plot some events
    TString hname;
    hEvBaseWave = NULL;
    if (!trig && !(idet->pass) && badDir->GetList()->GetEntries() < 100)
    {
      badDir->cd();
      hname.Form("EvBaseWaveEv%ich%ithresh%icross%i", int(entry), ichan, idet->thresholds, idet->crossings);
      // cout << " failed " << hname << endl;
      hEvBaseWave = new TH1D(hname, hname, nbins, 0, nbins);
      for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
        hEvBaseWave->SetBinContent(j + 1, digi[j]);
      }
    fout->cd();

    // plot some events
    if (trig && !(idet->pass) && badTrigDir->GetList()->GetEntries() < 100)
    {
      badTrigDir->cd();
      hname.Form("EvBaseWaveEv%ich%ithresh%icross%i", int(entry), ichan, idet->thresholds, idet->crossings);
      // cout << " failed " << hname << endl;
      hEvBaseWave = new TH1D(hname, hname, nbins, 0, nbins);
      for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
        hEvBaseWave->SetBinContent(j + 1, digi[j]);
    }
    fout->cd();

    if (idet->pass && evDir->GetList()->GetEntries() < 100)
    {
      evDir->cd();
      hname.Form("EvBaseWaveEv%ich%ithresh%icross%i", int(entry), ichan, idet->thresholds, idet->crossings);
      hEvBaseWave = new TH1D(hname, hname, nbins, 0, nbins);
      for (unsigned j = 0; j < digi.size(); ++j)
        hEvBaseWave->SetBinContent(j + 1, digi[j]);
    }
    fout->cd();

  } // channel loop
  // fill output Tree for each event
  // cout << "finished " << entry << endl;
  tbrun->fill();
  vector<float> fsum;
  fsum.resize(14);
  int passBit = 0;
  bool eventPass = true;
  for (int ib = 0; ib < rawBr.size(); ++ib)
  {
    unsigned ichan = chanList[ib];
    TDet *idet = tbrun->getDet(ichan);
    //cout << ichan << " " << idet->sum << endl;
    fsum[ichan] = idet->sum;
    if (!idet->pass && ichan != 12)
    {
      // printf(" %llu bad chan %u thresh %u crossing %u \n ", entry,ichan, idet->thresholds, idet->crossings);
      eventPass = false;
      passBit |= 0x1;
    }
    if(ichan == 12 && idet->sum>20){  // cosmic cut 
      passBit |= 0x2;
      eventPass = false;
    }
  }
  fsum[13] = float(eventPass);
  ntChanSum->Fill(&fsum[0]);

  evCount->Fill(-1); // underflow bin
  hEventPass->Fill(passBit);
  if (!eventPass)
    return eventPass;

  for (int ib = 0; ib < rawBr.size(); ++ib)
  {
    unsigned ichan = chanList[ib];
    bool trig = ichan == 9 || ichan == 10 || ichan == 11;
    int nbins = rawBr[ib]->rdigi.size();
    // if (nbins != 1024)
    //   continue;
    TDet *idet = tbrun->getDet(ichan);
    if (idet == NULL)
    {
      printf("!!!!!NULL idet br %u ichan %i\n", ib, ichan);
      continue;
    }
    // construct digi for channel
    digi.clear();
    for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
    {
      double val = (double)rawBr[ib]->rdigi[j] - idet->base;
      digi.push_back(val);
    }

    /* fill the summed histogram*/
    // histogram bad events
    if (!idet->pass)
    {
      for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
      {
        double val = digi[j];
        sumWaveB[ib]->SetBinContent(j + 1, sumWaveB[ib]->GetBinContent(j + 1) + val);
        valHistB[ib]->Fill(val);
      }
    }

    for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
    {
      double val = digi[j];
      sumWave[ib]->SetBinContent(j + 1, sumWave[ib]->GetBinContent(j + 1) + val);
      valHist[ib]->Fill(val);
    }

    // if (idet->thresholds < 1)
    //      continue;

    evCount->Fill(ichan); // chan 0 from GetBinContent(0)

    /* pulse finding
      hitFinder::event(int ichan, Long64_t ievent, vector<double> eventDigi,double thresh, unsigned step)
    */
    double hitThreshold = 10.0;
    if (ichan == 12)
      hitThreshold = 100;
    finder->event(ichan, entry, digi, hitThreshold, 1); // DEG suggests 10
    fftDir = (TDirectory *)fout->FindObject("fftDir");
    if (!fftDir)
    {
      cout << " Error no fftDir" << endl;
      fout->ls();
      return false;
    }
   
  } // second channel loop

  // fill sumHitWave and Q sums
  // loop over detector channels
  bool plotThisEvent=false;
  for (unsigned idet = 0; idet < tbrun->detList.size(); ++idet)
  {
    plotThisEvent = false;
    TDet *tdet = tbrun->detList[idet];
    // hit threshold cut done in hitFinder
    histQSum->Fill(tdet->channel, tdet->qSum);
    histQPrompt->Fill(tdet->channel, tdet->qPrompt);
    int ichan = tdet->channel;
    // add some event plots
    bool trig = tdet->channel == 9 || tdet->channel == 10 || tdet->channel == 11;
    TDirectory *fftDir = (TDirectory *)fout->FindObject("fftDir");
    if (!trig && tdet->hits.size() > 0 && fftDir->GetList()->GetEntries() < 2000){
      finder->plotEvent(tdet->channel, entry);
      printf("ev %llu chan %i nhits %lu\n", entry, tdet->channel, tdet->hits.size());
    }

    //if (ichan == 6)
      //printf("ev %llu chan %i nhits %lu\n", entry, tdet->channel, tdet->hits.size());
    histHitCount->SetBinContent(ichan, histHitCount->GetBinContent(ichan) + tdet->hits.size());
    // loop over hits
    for (unsigned ihit = 0; ihit < tdet->hits.size(); ++ihit)
    {
      TDetHit thit = tdet->hits[ihit];
      // fill QPE shape
      if(thit.qsum>7000&&thit.qsum<10000&&thit.startTime>65) {
        plotThisEvent = true;
        hQSum[idet]->Fill(thit.qsum);
        hQPeak[idet]->Fill(thit.qpeak);
        TH1D *hitPlot = thit.plot();
        if(hitPlot) fftDir->Append(hitPlot);
        printf(" \t %lli 1QPE hit %i start,peak,end (%i,%i,%i)  peak val %f \n", entry, ihit, thit.firstBin, thit.peakBin, thit.lastBin, thit.qpeak);
        //  loop over ADC values
        for (unsigned jbin = thit.firstBin; jbin < thit.lastBin; ++jbin)
        {
          int fillBin = hQPEShape[idet]->GetNbinsX() / 2 + jbin - thit.peakBin;
          double val = double(rawBr[idet]->rdigi[jbin]) - tdet->base;
          hQPEShape[idet]->SetBinContent(fillBin, hQPEShape[idet]->GetBinContent(fillBin) + val);
          //if (ichan == 6)
            //printf("\t\t  fill val at %i val %f binval %f \n", jbin - thit.peakBin, val, hQPEShape[idet]->GetBinContent(fillBin));
          }
      /*} else {
        if (ichan == 6) printf(" \t other hit %i start %i end %i peak bin %i peak val %f \n", ihit, thit.firstBin, thit.lastBin, thit.peakBin, thit.qpeak);
        */
      } 
      // do threshold for summed waveform
      if (thit.qsum > hitQThreshold)
        sumHitWave[idet]->SetBinContent(thit.firstBin + 1, sumHitWave[idet]->GetBinContent(thit.firstBin + 1) + thit.qsum);
    }
    //if (fftDir->GetList()->GetEntries() < 2000)
    //    finder->plotEvent(ichan, entry);
  }

  //printf(" event %llu  pass %i fail 1 %i cosmic only %i fail both %i \n",entry, int(hEventPass->GetBinContent(1)), int(hEventPass->GetBinContent(2)), int(hEventPass->GetBinContent(3)), int(hEventPass->GetBinContent(4)));

  return eventPass;
} // anaEvent

// revised derivative Dec 8 2022 MG
void anaRun::differentiate(unsigned diffStep)
{
  ddigi.clear();
  Double_t sump = 0;
  Double_t summ = 0;
  unsigned nsamples = digi.size();
  ddigi.push_back(0); // first entry is zero
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
      summ = digi[i - 1 - j];
    ddigi.push_back(sump - summ);
  }
}

void anaRun::negativeCrossingCount(double thresh)
{
  crossings.clear();
  Double_t cut = QPEPeak * thresh;
  for (unsigned ibin = 0; ibin < digi.size(); ++ibin)
  {
    Double_t vi = digi[ibin];
    if (vi < -cut)
      crossings.push_back(ibin);
  }
}
// count threshold crossings
void anaRun::thresholdCrossingCount(double thresh)
{
  thresholds.clear();
  Double_t cut = QPEPeak * thresh;
  for (unsigned ibin = 0; ibin < digi.size(); ++ibin)
  {
    if (digi[ibin] < cut && digi[ibin + 1] > cut)
      thresholds.push_back(ibin);
  }
}

void anaRun::derivativeCount(TDet *idet, Double_t rms)
{
  crossings.clear();
  crossingBin.clear();
  crossingTime.clear();
  unsigned vsize = ddigi.size();
  double microSec = 1.0E-3;
  double timeUnit = 8.0;
  Double_t cut = idet->sigma * rms;
  unsigned step = 1;
  // cout << " for det " << idet->channel  << " in derivative peaks >>>> rms " << rms << " cut " << cut << endl;
  Double_t ncut = -cut;
  // find all crossings
  for (unsigned ibin = step; ibin < vsize; ++ibin)
  {
    Double_t u = double(ibin) * timeUnit;
    Double_t vi = ddigi[ibin];
    Double_t vj = ddigi[ibin - step];
    unsigned ctype = 10;
    // if (idet == 5)
    // printf(" det %i  bin %i %f %f  \n", idet, ibin, vj, vi);
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
  }
  return;
}

Long64_t anaRun::anaRunFile(TString theFile, Long64_t maxEntries)
{
  clear();
  currentBuffer = -1;
  currentBufferCount = 0;
  cout << " starting anaRun entries = " << maxEntries << " file =  " << theFile << endl;
  if (!openFile(theFile))
  {
    printf("cannot open file!\n");
    return 0;
  }
  getTrees();
  cout << "treeList size " << treeList.size() << endl;

  if (treeList.size() < 1)
  {
    printf("EXIT no trees\n");
    return 0;
  }
  rawBr.clear();
  rawBr.resize(treeList.size());
  getEvent(0);
  string sfilename(theFile.Data());
  string shortName = sfilename.substr(0, sfilename.find_last_of("."));
  cout << " instance of anaRun with shortName= " << shortName << endl;
  fout = new TFile(Form("myData/anaRun-%s-%llu.root", shortName.c_str(), maxEntries), "recreate");
  evDir = fout->mkdir("evDir");
  badDir = fout->mkdir("badDir");
  badTrigDir = fout->mkdir("badTrigDir");
  fout->cd();
  cout << " opened output file " << fout->GetName() << endl;
  getSummedHists();
  // fout->ls();

  // fin->ls();
  TBFile *bf;
  fin->GetObject("tbfile", bf);
  if (bf)
  {
    bf->print();
    fout->Append(bf);
  }
  else
    printf("anaRun ERROR!! no tbfile in input file \n");

  tbrun = new TBRun(tag);
  for (unsigned it = 0; it < chanList.size(); ++it)
    tbrun->addDet(chanList[it]);

  // chanlist check
  cout << "**** chanList check **** size " << chanList.size() << endl;
  for (unsigned ic = 0; ic < chanList.size(); ++ic)
    printf(" %u list/key  %u mapped to %u treeList %s det in tbrun %s \n ",
           ic, chanList[ic], chanMap.at(chanList[ic]), treeList[ic]->GetName(), tbrun->detList[ic]->GetName());

  // fout->append(tbrun->btree);
  int totalchannels = 13;
  ntChan = new TNtuple("ntchan", "channel ntuple", "trig:chan:ave:sigma:skew:base:peakmax:sum2:sum:negcrossings:thresholds");
  ntChanSum = new TNtuple("ntchansum", "channel ntuple", "sum0:sum1:sum2:sum3:sum4:sum5:sum6:sum7:sum8:sum9:sum10:sum11:sum12:pass");
  hEventPass = new TH1D("EventPass", " event failures",4,0,4);
  evCount = new TH1D("eventcount", "event count", totalchannels, 0, totalchannels);
  histQSum = new TH1D("histqsum", "qsum by channel", totalchannels, 0, totalchannels);
  histHitCount = new TH1D("hitCount", "hit count by channel", totalchannels, 0, totalchannels);
  // nn/histqpe = new th1d("histqpe", "qpe by channel", totalchannels, 0, totalchannels);
  histQPrompt = new TH1D("histqprompt", "qprompt by channel", totalchannels, 0, totalchannels);
  histQSum->Sumw2();
  histQPrompt->Sumw2();

  //
  anaDir = fout->mkdir("anadir");
  anaDir->cd();
  for (unsigned i = 0; i < rawBr.size(); ++i)
  {
    unsigned ichan = chanList[i];
    noiseHist.push_back(new TH1D(Form("noiseChan%i", ichan), Form("noiseChan%i", ichan), 1000, 0, 1000));
    skewHist.push_back(new TH1D(Form("skewChan%i", ichan), Form("skewChan%i", ichan), 200, -3, 7));
    threshHist.push_back(new TH1D(Form("threshChan%i", ichan), Form("threshChan%i", ichan), 20, 0, 20));
    crossHist.push_back(new TH1D(Form("crossChan%i", ichan), Form("crossChan%i", ichan), 100, 0, 100));
    if (ichan > 8 && ichan < 12)
    {
      valHist.push_back(new TH1D(Form("valChan%i", ichan), Form("valChan%i", ichan), 1500, -500, 1000));
      valHistB.push_back(new TH1D(Form("valBadChan%i", ichan), Form("valBadChan%i", ichan), 1500, -500, 1000));
      hEvGaus.push_back(new TH1D(Form("evGaus%i", ichan), Form("evGaus%i", ichan), 200, -500, 500));
    }
    else
    {
      valHist.push_back(new TH1D(Form("valChan%i", ichan), Form("valChan%i", ichan), 1000, -200, 200));
      valHistB.push_back(new TH1D(Form("valBadChan%i", ichan), Form("valBadChan%i", ichan), 1000, -200, 200));
      hEvGaus.push_back(new TH1D(Form("evGaus%i", ichan), Form("evGaus%i", ichan), 200, -100, 100));
    }
  }
  sumDir = fout->mkdir("sumDir");
  sumDir->cd();
  double limit;
  double plimit;
  for (unsigned i = 0; i < rawBr.size(); ++i)
  {
    unsigned ichan = chanList[i];

    bool trigger = ichan == 9 || ichan == 10 || ichan == 11;
    if (trigger)
    {
      limit = 200000;
      plimit = 10000;
    }
    else
    {
      plimit = 2000;
    }

    // else if(ichan==12) limit = 10000;
    if (ichan == 12)
    {
      limit = 10000;
    }

    hQSum.push_back(new TH1D(Form("QSumChan%i", ichan), Form("QSumChan%i", ichan), 1000, 0, limit));
    hQPEShape.push_back(new TH1D(Form("QPEShapeChan%i", ichan), Form("QPEShapeChan%i", ichan), 200, -100, 100));
    hQPEShape[hQPEShape.size() - 1]->SetMarkerStyle(20);
    hQPeak.push_back(new TH1D(Form("QPeakChan%i", ichan), Form("QPeakChan%i", ichan), 1000, 0, plimit));
    sumWave.push_back(new TH1D(Form("sumWave%i", ichan), Form("sumWave%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));
    sumWaveB.push_back(new TH1D(Form("sumWaveBad%i", ichan), Form("sumWaveBad%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));
    sumHitWave.push_back(new TH1D(Form("sumHitWave%i", ichan), Form("sumHitWave%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));
  }
  fout->cd();

  // fout->ls();
  cout << " make hitFinder dets = " << rawBr.size() << endl;
  finder = new hitFinder(fout, tbrun, tag, rawBr[0]->rdigi.size(), chanList);
  Long64_t nentries = treeList[0]->GetEntries();
  if (maxEntries > 0)
    nentries = TMath::Min(maxEntries, nentries);
  printf("... total entries  %llu looping over %llu \n ", treeList[0]->GetEntries(), nentries);
  unsigned npass = 0;
  unsigned nfail = 0;
  for (Long64_t entry = 0; entry < nentries; ++entry)
  {
    if (entry / 1000 * 1000 == entry)
      printf("... %llu pass %u fail %u \n", entry, npass, nfail);
    getEvent(entry);
    /* test
    for (unsigned it = 0; it < treeList.size(); ++it)
      printf("tree %u chan %i event %i time %lld \n", it, rawBr[it]->channel, rawBr[it]->trigger, rawBr[it]->time);
    */
    if (anaEvent(entry))
      ++npass;
    else
      ++nfail;
  }
  // normailize to number of nentries.
  // loop over detector channels
  double scaleFactor = 1. / double(nentries);
  printf(" \n \n At END OF FILE scale by %E\n", scaleFactor); 
  histQSum->Scale(scaleFactor);
  histQPrompt->Scale(scaleFactor);
  for (unsigned i = 0; i < sumHitWave.size(); ++i)
  {
    sumWave[i]->Scale(scaleFactor); // DEF added 
    sumHitWave[i]->Scale(scaleFactor);
    hQSum[i]->Scale(scaleFactor);
    hQPeak[i]->Scale(scaleFactor);
  }

  // finder->hPeakCount->Print("all");
  // tbrun->btree->GetListOfBranches()->ls();
  // fout->ls();
  /* do fits at end of run */

  TString graphName = TString("slope-graph");
  TString graphTitle = TString(Form("slope-graph-%s", shortName.c_str()));
  printf(" making slope graph %s \n", graphName.Data());

  for (unsigned i = 0; i < sumHitWave.size(); ++i)
  {
    sumHitWave[i]->Fit("expo", "Q", "", 100, 300); // DEG suggests
    TF1 *g = (TF1 *)sumHitWave[i]->GetListOfFunctions()->FindObject("expo");
    chan.push_back(chanList[i]);
    echan.push_back(0);
    if (g)
    {
      printf("%s %E %E \n", sumHitWave[i]->GetName(), g->GetParameter(1), g->GetParError(1));
      slope.push_back(g->GetParameter(1));
      eslope.push_back(g->GetParError(1));
    }
    else
    {
      slope.push_back(0);
      eslope.push_back(0);
    }
  }
  TGraphErrors *grslope = new TGraphErrors(chan.size() - 4, &chan[3], &slope[3], &eslope[3], &echan[3]);
  grslope->SetName(graphName);
  grslope->SetTitle(graphTitle);
  fout->Append(grslope);

  printf(" FINISHED pass %i fail 1 %i cosmic only %i fail both %i \n", int(hEventPass->GetBinContent(1)), int(hEventPass->GetBinContent(2)), int(hEventPass->GetBinContent(3)), int(hEventPass->GetBinContent(4)));

  printf(" FINISHED npass %u nfail %u output file  %s \n", npass, nfail, fout->GetName());
  printf(" hit count \n");
  for (int ibin = 0; ibin <= histHitCount->GetNbinsX(); ++ ibin)
    printf(" bin %i count %i ;", ibin, int(histHitCount->GetBinContent(ibin)));
  printf("  \n");

  fout->Write();
  fout->Close();
  return nentries;
}

anaRun::anaRun(TString theTag)
{
  tag = theTag;
  cout << " instance of anaRun with tag= " << tag << endl;
}
