/////////////////////////////////////////////////////////
//  M.Gold April 2023 read CAEN files
/////////////////////////////////////////////////////////
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
#include "TBEventData.hxx"
#include "TBRawEvent.hxx"
#include "hitFinder.hxx"
#include "TBFile.hxx"

class anaCRun
{
public:
  enum
  {
    UPCROSS,
    DOWNCROSS,
    DOUBLEUPCROSS,
    DOUBLEDOWNCROSS
  };
  enum
  {
    CHANNELS = 13
  };
  enum
  {
    WAVELENGTH = 7500
  };
  TFile *fout;
  TFile *fin;
  TTree *rawTree;
  std::map<int, int> chanMap;
  vector<TBRawEvent *> rawBr;
  TBEventData *eventData;
  TBRun *tbrun;
  TNtuple *ntChan;
  TNtuple *ntChanSum;
  vector<TH1D *> baseHist;
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
  TH1D *hEvBaseWave;
  vector<TH1D *> hEvGaus;
  TH1D *evCount;
  TH1D *histQSum;
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
  anaCRun(TString theTag = TString("dirName"));
  ~anaCRun() {}
  Long64_t anaCRunFile(TString theFile, Long64_t maxEntries);
  void clear();
  bool openFile(TString fileName);
  unsigned getListOfFiles(TString dir);
  void getSummedHists();
  unsigned getBranches();
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

void anaCRun::clear()
{
  hQSum.clear();
  hQPeak.clear();
  chanMap.clear();
  baseHist.clear();
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

bool anaCRun::openFile(TString theFile)
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
  rawTree = NULL;
  fin->GetObject("RawTree", rawTree);
  if (!rawTree)
  {
    printf(" no RawTree in file %s\n", fileName.Data());
    return false;
  }
  cout << "  RawTree has " << rawTree->GetEntries() << " entries " << endl;
  eventData = new TBEventData();
  rawTree->SetBranchAddress("eventData", &eventData);
  if (!eventData)
  {
    printf(" eventData not found in file  %s\n", fileName.Data());
    return false;
  }
  printf(" rawTree has %u channels stored in rawBr \n", getBranches());
  for (unsigned i = 0; i < rawBr.size(); ++i)
    printf(" branch %s chan %i \n", rawBr[i]->GetName(), i);

  return true;
}

/* get rawBr */
unsigned anaCRun::getBranches()
{
  TObjArray *brList = rawTree->GetListOfBranches();
  TString cname;
  TIter next(brList);
  TBranch *aBranch = NULL;
  while ((aBranch = (TBranch *)next()))
    {
      TString s(aBranch->GetName());
      if (s != TString("eventData"))
      {
        int ichan = TString(s(s.Last('n') + 1, s.Length())).Atoi();
        // rawTree->GetBranch(aBranch->GetName())->SetAutoDelete(kTRUE);
        cout << s << "  " << aBranch->GetName() << " return val =  " << rawTree->SetBranchAddress(aBranch->GetName(), &rawBr[ichan]) << endl;
      }
    }
    return rawBr.size();
}

// get summed histos
void anaCRun::getSummedHists()
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
  cout << " found " << rawSumDir->GetList()->GetEntries() << " summed histos " << endl;
  fout->cd();
  return;
}

/* analyze rawBr */
bool anaCRun::anaEvent(Long64_t entry)
{
  double vsign = 1.0;
  QPEPeak = 100;
  tbrun->clear();
  // loop over channels
  for (unsigned ib = 0; ib < rawBr.size(); ++ib)
  {
    unsigned ichan = ib;
    // define trigger sipms
    bool trig = ichan == 9 || ichan == 10 || ichan == 11;
    // set sign of waveform
    if (ib > 8)
      vsign = -1.0;
    int nbins = rawBr[ib]->rdigi.size();
    //cout << ib << " nbins " << nbins << " max hist " << hEvGaus.size() << " rawBr.size() " << rawBr.size() << endl;

    // sanity check
    if (rawBr[ib]->rdigi.size() != WAVELENGTH)
      continue;
    // simple baseline
    std::vector<unsigned short> orderDigi = rawBr[ib]->rdigi;
    std::sort(orderDigi.begin(), orderDigi.end());
    unsigned baseLength = orderDigi.size() /2;
    double base = 0;
    for (unsigned j = 0; j < baseLength; ++j)
    {
      base += orderDigi[j];
    }
    base /= double(baseLength);

    // baseline correction from fitted Gaussian
    hEvGaus[ib]->Reset("ICES");
    for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j){
      double val = vsign * (double(rawBr[ib]->rdigi[j])-base);
      hEvGaus[ib]->Fill(val);
    }
    
    hEvGaus[ib]->Fit("gaus", "QO", "", hEvGaus[ib]->GetMean()-100, hEvGaus[ib]->GetMean()+100);
    TF1 *gfit = (TF1 *)hEvGaus[ib]->GetListOfFunctions()->FindObject("gaus");
    double ave = hEvGaus[ib]->GetMean();
    double sigma = hEvGaus[ib]->GetRMS();
    double skew = 0;
    double fitMean = 0;
    if (!isnan(hEvGaus[ib]->GetSkewness()))
      skew = hEvGaus[ib]->GetSkewness();
    if (gfit != nullptr)
    {
      ave = gfit->GetParameter(1);
      fitMean = ave;  //fit mean
      sigma = gfit->GetParameter(2);
    }

    evDir->cd();
    if (evDir->GetList()->GetEntries() < 100)
    {
      TH1D *EvGaussEvent = (TH1D *)hEvGaus[ib]->Clone(Form("EvGaussEvent%lld-Ch%i", entry, ib));
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
    idet->base = base+fitMean; // add in fit mean if fit succeeded
    
    /*********
     *  fill baseline subtracted vector digi
     *********/
    digi.clear();
    for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
    {
      double val = vsign * (double(rawBr[ib]->rdigi[j]) - idet->base);
      baseHist[ichan]->Fill(val);
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
    // add some other variables
    idet->peak = peakMax;
    idet->sum2 = sum2;
    idet->sum = sum;

    //cout << "tree " << ib << " ch " << ichan << " "
    //     << " ave " << ave << " sigma " << sigma << " skew " << skew << " base " << base << endl;

    /*** 
     * these 2 functions use didi and crossings,threshold vectors  i admit ugly 
    ****/
    negativeCrossingCount(10.);
    idet->crossings = int(crossings.size());
    crossHist[ib]->Fill(idet->crossings);
    if (trig)
      thresholdCrossingCount(15.);
    else
      thresholdCrossingCount(1.);
    idet->thresholds = int(thresholds.size());
    threshHist[ib]->Fill(idet->thresholds);

    // set cut
    unsigned maxCrossings = 1;
    unsigned maxThresholds = 1;
    idet->pass = true;
    if (idet->crossings > maxCrossings)
      idet->pass = false;
    if (trig && idet->thresholds > maxThresholds)
      idet->pass = false;

    ntChan->Fill(float(rawBr[ib]->trigger), float(ichan), float(ave), float(sigma), float(skew), float(base), float(peakMax), float(sum2), float(sum), float(crossings.size()), float(thresholds.size()), float(idet->pass));

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

    if (idet->pass && idet->thresholds>0  && evDir->GetList()->GetEntries() < 100)
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
  for (unsigned ib = 0; ib < rawBr.size(); ++ib)
  {
    unsigned ichan = ib;
    
    TDet *idet = tbrun->getDet(ichan);
    // cout << ichan << " " << idet->sum << endl;
    fsum[ichan] = idet->sum;
    //if (!idet->pass && ichan != 12)
    if (!idet->pass && ichan<9)
    {
      // printf(" %llu bad chan %u thresh %u crossing %u \n ", entry,ichan, idet->thresholds, idet->crossings);
      eventPass = false;
      passBit |= 0x1;
    }
    if (ichan == 12 && idet->sum > 20)
    { // cosmic cut
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
  vsign = 1.0;
  for (unsigned ib = 0; ib < rawBr.size(); ++ib)
  {
    if (ib > 8)
      vsign = -1.0;
    unsigned ichan = ib;
    int nbins = rawBr[ib]->rdigi.size();
    // if (nbins != 1024)
    //   continue;
    TDet *idet = tbrun->getDet(ichan);
    if (idet == NULL)
    {
      printf("!!!!!NULL idet br %u ichan %i\n", ib, ichan);
      continue;
    }

    /* fill the summed histogram*/
    digi.clear();

    // histogram bad events
    if (!idet->pass)
    {
      for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
      {
        double val = vsign * (double(rawBr[ib]->rdigi[j]) - idet->base);
        sumWaveB[ib]->SetBinContent(j + 1, sumWaveB[ib]->GetBinContent(j + 1) + val);
        valHistB[ib]->Fill(val);
      }
    }

    for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
    {
      double val = vsign * (double(rawBr[ib]->rdigi[j]) - idet->base);
      digi.push_back(val);
      sumWave[ib]->SetBinContent(j + 1, sumWave[ib]->GetBinContent(j + 1) + val);
      valHist[ib]->Fill(val);
    }

    // if (idet->thresholds < 1)
    //      continue;

    evCount->Fill(ichan); // chan 0 from GetBinContent(0)

    //pulse finding
    // hitFinder::event(int ichan, Long64_t ievent, vector<double> eventDigi,double thresh, unsigned step)
    
    //return true;
    finder->event(ichan, entry, digi, 10., 1); // DEG suggests 10
    TDirectory *fftDir = (TDirectory *)fout->FindObject("fftDir");
    if (!fftDir)
    {
      cout << " Error no fftDir" << endl;
      fout->ls();
      return false;
    }
    // add some event plots
    if (fftDir->GetList()->GetEntries() < 1000)
      finder->plotEvent(ichan, entry);
  } // second channel loop

  // fill sumHitWave and Q sums
  // loop over detector channels
  for (unsigned idet = 0; idet < tbrun->detList.size(); ++idet)
  {
    TDet *tdet = tbrun->detList[idet];
    // hit threshold cut done in hitFinder
    histQSum->Fill(tdet->channel, tdet->qSum);
    histQPrompt->Fill(tdet->channel, tdet->qPrompt);

    for (unsigned ihit = 0; ihit < tdet->hits.size(); ++ihit)
    {
      TDetHit thit = tdet->hits[ihit];
      hQSum[idet]->Fill(thit.qsum);
      hQPeak[idet]->Fill(thit.qpeak);
      // do threshold for summed waveform
      if (thit.qsum > hitQThreshold)
        sumHitWave[idet]->SetBinContent(thit.firstBin + 1, sumHitWave[idet]->GetBinContent(thit.firstBin + 1) + thit.qsum);
    }
  }

  // printf(" event %llu  pass %i fail 1 %i cosmic only %i fail both %i \n",entry, int(hEventPass->GetBinContent(1)), int(hEventPass->GetBinContent(2)), int(hEventPass->GetBinContent(3)), int(hEventPass->GetBinContent(4)));

  return eventPass;
} // anaEvent

// revised derivative Dec 8 2022 MG
void anaCRun::differentiate(unsigned diffStep)
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

void anaCRun::negativeCrossingCount(double thresh)
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
void anaCRun::thresholdCrossingCount(double thresh)
{
  thresholds.clear();
  Double_t cut = QPEPeak * thresh;
  for (unsigned ibin = 0; ibin < digi.size(); ++ibin)
  {
    if (digi[ibin] < cut && digi[ibin + 1] > cut)
      thresholds.push_back(ibin);
  }
}

void anaCRun::derivativeCount(TDet *idet, Double_t rms)
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

Long64_t anaCRun::anaCRunFile(TString theFile, Long64_t maxEntries)
{
  clear();
  currentBuffer = -1;
  currentBufferCount = 0;
  cout << " starting anaCRun entries = " << maxEntries << " file =  " << theFile << endl;
  if (!openFile(theFile)) // and get branches
  {
    printf("cannot open file!\n");
    return 0;
  }

  string sfilename(theFile.Data());
  string shortName = sfilename.substr(0, sfilename.find_last_of("."));
  cout << " instance of anaCRun with shortName= " << shortName << endl;

  // need to fill rawBr[0]->rdigi.size()
  printf("Read zeroth entry from tree \n");
  rawTree->GetEntry(0);
  printf("\t\t start of file %i %i %i : %i\n", eventData->day, eventData->mon, eventData->year, eventData->hour);
  printf("\t\t SIZE OF WAVEFORM = %lu \n", rawBr[0]->rdigi.size());
  if (rawBr[0]->rdigi.size() != WAVELENGTH)
  {
    printf(" \n\n\n\n ERROR rdigi size %lu !!! \n", rawBr[0]->rdigi.size());
    //return 0;
  }
  Long64_t nentries = rawTree->GetEntries();
  if (maxEntries > 0)
    nentries = TMath::Min(maxEntries, nentries);
  printf("... total entries  %llu looping over %llu \n ", rawTree->GetEntries(), nentries);

/* test 
  for (Long64_t entry = 0; entry < nentries; ++entry)
  {
    cout << " entry " << entry << " ";
    rawTree->GetEntry(entry);
    for (unsigned i = 0; i < rawBr.size(); ++i)
      cout << i << " " << rawBr[i]->rdigi[0] << " ; ";
    cout << endl;
  }

  return nentries;
  */

  // open outout file and make histograms

  fout = new TFile(Form("myData/anaCRun-%s-%llu.root", shortName.c_str(), maxEntries), "recreate");
  evDir = fout->mkdir("evDir");
  badDir = fout->mkdir("badDir");
  badTrigDir = fout->mkdir("badTrigDir");
  fout->cd();
  cout << " opened output file " << fout->GetName() << endl;
  getSummedHists();
  // fout->ls();

  // make output tree
  tbrun = new TBRun(tag);
  for (unsigned it = 0; it < rawBr.size(); ++it)
  {
    tbrun->addDet(it);
  }

  // fout->append(tbrun->btree);
  int totalchannels = rawBr.size() + 1;
  ntChan = new TNtuple("ntchan", "channel ntuple", "trig:chan:ave:sigma:skew:base:peakmax:sum2:sum:negcrossings:thresholds:pass");
  ntChanSum = new TNtuple("ntchansum", "channel ntuple", "sum0:sum1:sum2:sum3:sum4:sum5:sum6:sum7:sum8:sum9:sum10:sum11:sum12:pass");
  hEventPass = new TH1D("EventPass", " event failures", 4, 0, 4);
  evCount = new TH1D("eventcount", "event count", totalchannels, 0, totalchannels);
  histQSum = new TH1D("histqsum", "qsum by channel", totalchannels, 0, totalchannels);
  // nn/histqpe = new th1d("histqpe", "qpe by channel", totalchannels, 0, totalchannels);
  histQPrompt = new TH1D("histqprompt", "qprompt by channel", totalchannels, 0, totalchannels);
  histQSum->Sumw2();
  histQPrompt->Sumw2();

  //
  anaDir = fout->mkdir("anadir");
  anaDir->cd();
  for (unsigned i = 0; i < rawBr.size(); ++i)
  {
    unsigned ichan = i;
    baseHist.push_back(new TH1D(Form("noiseChan%i", ichan), Form("noiseChan%i", ichan), 2000,-1000, 1000));
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
  double limit = 10000;
  double plimit = 2000;
  for (unsigned i = 0; i < rawBr.size(); ++i)
  {
    unsigned ichan = i;

    bool trigger = ichan == 9 || ichan == 10 || ichan == 11;
    if (trigger)
    {
      limit = 200000;
      plimit = 10000;
    }

    // else if(ichan==12) limit = 10000;

    hQSum.push_back(new TH1D(Form("QSumChan%i", ichan), Form("QSumChan%i", ichan), 1000, 0, limit));
    hQPeak.push_back(new TH1D(Form("QPeakChan%i", ichan), Form("QPeakChan%i", ichan), 1000, 0, plimit));
    sumWave.push_back(new TH1D(Form("sumWave%i", ichan), Form("sumWave%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));
    sumWaveB.push_back(new TH1D(Form("sumWaveBad%i", ichan), Form("sumWaveBad%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));
    sumHitWave.push_back(new TH1D(Form("sumHitWave%i", ichan), Form("sumHitWave%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));
  }
  fout->cd();

  // fout->ls();
  cout << " make hitFinder dets = " << rawBr.size() << "  size " << rawBr[0]->rdigi.size() << endl;
  //finder = NULL;
  vector<int> chanList;
  for (int ichan = 0; ichan<rawBr.size(); ++ichan)
    chanList.push_back(ichan);
  finder = new hitFinder(fout, tbrun, tag, rawBr[0]->rdigi.size(), chanList);

  unsigned npass = 0;
  unsigned nfail = 0;
  printf("... total entries  %llu looping over %llu \n ", rawTree->GetEntries(), nentries);
  for (Long64_t entry = 0; entry < nentries; ++entry)
  {
    if (entry / 1000 * 1000 == entry)
      printf("... %llu pass %u fail %u \n", entry, npass, nfail);
    rawTree->GetEntry(entry);

    /*
    for (unsigned it = 0; it < rawBr.size(); ++it)
      if (rawBr[it])
        printf("branch %u chan %i trigger %i time %lld digi0 %u \n",
               it, rawBr[it]->channel, rawBr[it]->trigger, rawBr[it]->time, rawBr[it]->rdigi[0]);
      else
        printf("NULL branch %i \n", it);
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
    chan.push_back(i);
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

  fout->Write();
  fout->Close();
  printf(" FINISHED npass %u nfail %u output file  %s \n", npass, nfail, fout->GetName());

  return nentries;
}

anaCRun::anaCRun(TString theTag)
{
  tag = theTag;
  cout << " instance of anaCRun with tag= " << tag << "CHANNELS" << CHANNELS << endl;

  rawBr.clear();
  for (int ichan = 0; ichan < CHANNELS; ++ichan)
  {
    TBRawEvent *rawEv = new TBRawEvent(ichan);
    rawEv->rdigi.resize(7500);
    rawEv->SetName(Form("rawChan%i", ichan));
    rawBr.push_back(rawEv);
  }
}