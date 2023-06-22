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
  int passBit;
  std::map<int, int> chanMap;
  vector<TBRawEvent *> rawBr;
  TBEventData *eventData;
  TBEventData *rawEventData;
  TBRun *tbrun;
  TNtuple *ntChan;
  TNtuple *ntChanSum;
  vector<TH1D *> baseHist;
  vector<TH1D *> noiseHist;
  vector<TH1D *> skewHist;
  vector<TH1D *> sumWave;
  vector<TH1D *> sumHitWave;
  vector<TH1D *> sumPeakWave;
  vector<TH1D *> valHist;
  vector<TH1D *> sumWaveB;
  vector<TH1D *> valHistB;
  
  vector<TH1D *> hQSum;
  vector<TH1D *> hQPeak;
  vector<TH1D *> hQPEShape;
  TH1D *hEvBaseWave;
  vector<TH1D *> hEvGaus;
  vector<TH1D *> hChannelGaus;
  TH1D *evCount;
  TH1D *histQSum;
  TH1D *hEventPass;
  TH1D *histHitCount;
  TH1D *hSumPMT;
  TH1D *threshHist;
  TH1D *crossHist;
  TH1D *cosmicCut1;
  TH1D *cosmicCut2;
  // TH1D *histQPE;
  TH1D *histQPrompt;
  vector<double> channelSigmaValue;
  vector<double> channelSigma;
  vector<double> channelSigmaErr;
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
  void negativeCrossingCount(int ichan);
  void thresholdCrossingCount(double thresh);


  TDirectory *rawSumDir;
  TDirectory *badDir;
  TDirectory *badTrigDir;
  TDirectory *evDir;
  TDirectory *sumDir;
  TDirectory *anaDir;
  TDirectory *pmtDir;
  Long64_t nentries;
  double QPEPeak;
};


void anaCRun::clear()
{
  hQSum.clear();
  hQPEShape.clear();
  hQPeak.clear();
  chanMap.clear();
  baseHist.clear();
  noiseHist.clear();
  skewHist.clear();
  sumWave.clear();
  sumHitWave.clear();
  sumPeakWave.clear();
  valHist.clear();
  sumWaveB.clear();
  valHistB.clear();
  hEvGaus.clear();
  hChannelGaus.clear();
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
  // fill channel sigma in order of branches
  channelSigmaValue.resize(13);
  channelSigmaValue[0] = 12.6634;
  channelSigmaValue[1] = 11.7833;
  channelSigmaValue[2] = 11.9779;
  channelSigmaValue[3] = 0.247759;
  channelSigmaValue[4] = 14.1827;
  channelSigmaValue[5] = 12.4958;
  channelSigmaValue[6] = 12.346;
  channelSigmaValue[7] = 11.9752;
  channelSigmaValue[8] = 12.4379;
  channelSigmaValue[9] = 474.388;
  channelSigmaValue[10] = 450.223;
  channelSigmaValue[11] = 426.53;
  channelSigmaValue[12] = 3.50614;
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
  rawEventData = new TBEventData();
  rawTree->SetBranchAddress("eventData", &rawEventData);
  if (!rawEventData)
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
  // copy event data
  eventData->evtime = rawEventData->evtime;
  eventData->sec  = rawEventData->sec;
  eventData->min  = rawEventData->min;
  eventData->hour = rawEventData->hour;
  eventData->day  = rawEventData->day;
  eventData->mon  = rawEventData->mon;
  eventData->year = rawEventData->year;
  eventData->isdst = rawEventData->isdst;
  QPEPeak = 100;
  tbrun->clear();
  // loop over channels
  double vsign = 1.0;
  for (unsigned ib = 0; ib < rawBr.size(); ++ib)
  {
    unsigned ichan = ib;
    // define trigger sipms
    bool trig = ichan == 9 || ichan == 10 || ichan == 11;
    // set sign of waveform
    vsign = 1.0;
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
      double val = vsign * (double(rawBr[ib]->rdigi[j])-base); // base is > digi value!
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
    idet->base = base+vsign*fitMean; // add in fit mean if fit succeeded

    /*********
     *  fill baseline subtracted vector digi
     *********/
    digi.clear();
    for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
    {
      double val = vsign * (double(rawBr[ib]->rdigi[j])-idet->base);
      baseHist[ichan]->Fill(val);
      digi.push_back(val);
      if (hChannelGaus.size()>0) hChannelGaus[ib]->Fill(val);
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

    if (ichan == 12)
      hSumPMT->Fill(sum);

      //printf("chan %i %f sum %f \n", ichan, fitMean, sum);

      // cout << "tree " << ib << " ch " << ichan << " "
      //      << " ave " << ave << " sigma " << sigma << " skew " << skew << " base " << base << endl;

      /***
       * these 2 functions use didi and crossings,threshold vectors  i admit ugly
       ****/
    negativeCrossingCount(ichan); // multiples of QPEPeak
    idet->crossings = int(crossings.size());
    thresholdCrossingCount(1000.);
    idet->thresholds = int(thresholds.size());

    if (!trig)  crossHist->Fill(idet->crossings);
    if (trig)   threshHist->Fill(idet->thresholds);
    // set cut
    unsigned maxCrossings = 1;
    unsigned maxThresholds = 1;
    idet->pass = true;
    if (!trig && idet->crossings > maxCrossings)
      idet->pass = false;
    if (trig && idet->thresholds > maxThresholds)
      idet->pass = false;
      /*
    if (!trig && idet->crossings > maxCrossings)
      printf("fail crossings %i %i %i\n",ichan,idet->crossings, idet->pass);
    if (trig && idet->thresholds > maxThresholds)
      printf("fail threshold  %i %i %i\n",ichan,idet->thresholds,idet->pass);
      */

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
    // plot PMT
    if (ichan == 12 && pmtDir->GetList()->GetEntries() < 100)
    {
      pmtDir->cd();
      hname.Form("EvBaseWaveEv%ich%ithresh%icross%i", int(entry), ichan, idet->thresholds, idet->crossings);
      // cout << " failed " << hname << endl;
      hEvBaseWave = new TH1D(hname, hname, nbins, 0, nbins);
      for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
        hEvBaseWave->SetBinContent(j + 1, digi[j]);
    }

    fout->cd();

    // events with hits only data
    if (idet->pass && idet->thresholds>0  && evDir->GetList()->GetEntries() < 500)
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

  
  // defined in class as int
  passBit = 0;
  bool eventPass = true;
  for (unsigned ib = 0; ib < rawBr.size(); ++ib)
  {
    unsigned ichan = ib;
    bool trig = ichan == 9 || ichan == 10 || ichan == 11;

    TDet *idet = tbrun->getDet(ichan);
    if (trig && !idet->pass)
    {
      //printf(" %llu bad chan %u thresh %u crossing %u \n ", entry,ichan, idet->thresholds, idet->crossings);
      eventPass = false;
      passBit |= 0x1;
      //printf(".... set bit theshold  %i %i %i\n",ichan,idet->thresholds,passBit);
    }
    if(!trig && !idet->pass)
    {
      // printf(" %llu bad chan %u thresh %u crossing %u \n ", entry,ichan, idet->thresholds, idet->crossings);
      eventPass = false;
      passBit |= 0x2;
      //printf(".... set bit crossings  %i %i %i\n",ichan,idet->crossings,passBit);
    }

    if (ichan == 12 )
      cosmicCut1->Fill(idet->sum);

    if (ichan == 12 && idet->sum > 20)
    { // cosmic cut
      passBit |= 0x4;
      eventPass = false;
      //printf(".... set bit cosmic1  %i %f %i\n",ichan,idet->sum,passBit);
    }
  }

  evCount->Fill(-1); // underflow bin
  if (!eventPass)
    return eventPass;
  // continue if event passes 
  for (unsigned ib = 0; ib < rawBr.size(); ++ib)
  {
    vsign = 1.0;
    if (ib > 8)
      vsign = -1.0;
    unsigned ichan = ib;
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
    double hitThreshold = 5.0*channelSigmaValue[ib];
    finder->event(ichan, entry, digi,hitThreshold, 1); // DEG suggests 10
    TDirectory *fftDir = (TDirectory *)fout->FindObject("fftDir");
    if (!fftDir)
    {
      cout << " Error no fftDir" << endl;
      fout->ls();
      return false;
    } 
  } // second channel loop

  /**** second cosmic cut.****/
  TDet *tdet = tbrun->detList[12];
  bool cosmicTwo = true;
  if (tdet->channel != 12)
    printf("ERROR! detlist miss match to channel 12!!!\n");
  else {
    // loop over hits
    for (unsigned ihit = 0; ihit < tdet->hits.size(); ++ihit){
      TDetHit thit = tdet->hits[ihit];
      if(thit.startTime > 1000) cosmicCut2->Fill(thit.qpeak);
      if (thit.startTime > 1000 && thit.qpeak > 1000)
        eventPass = false;
    }
  }
  if (!eventPass){
    passBit |= 0x8;
    return eventPass;
  }
  // fill total light
  vector<float> fsum;
  fsum.resize(tbrun->detList.size());

  // loop over detector channels
  for (unsigned idet = 0; idet < tbrun->detList.size(); ++idet)
    {
      TDet *tdet = tbrun->detList[idet];
      fsum[tdet->channel] = tdet->sum;
      // hit threshold cut done in hitFinder
      if(tdet->qSum>0){
         histQSum->Fill(tdet->channel, tdet->qSum);
         histQPrompt->Fill(tdet->channel, tdet->qPrompt);
      }

      // add some event plots
      bool trig = tdet->channel == 9 || tdet->channel == 10 || tdet->channel == 11;
      TDirectory *fftDir = (TDirectory *)fout->FindObject("fftDir");
      if (!trig && tdet->hits.size() > 0 && fftDir->GetList()->GetEntries() < 2000)
        finder->plotEvent(tdet->channel, entry);

      // loop over hits
      //double hitThreshold = 5.0 * channelSigmaValue[idet];
      double hitThreshold = 2000; // value for hit.qsum threshold
      for (unsigned ihit = 0; ihit < tdet->hits.size(); ++ihit)
      {
        TDetHit thit = tdet->hits[ihit];
        hQSum[idet]->Fill(thit.qsum);
        hQPeak[idet]->Fill(thit.qpeak);
        // do threshold for summed waveform
        if (thit.qsum > hitThreshold)
        {
          sumHitWave[idet]->SetBinContent(thit.firstBin + 1, sumHitWave[idet]->GetBinContent(thit.firstBin + 1) + thit.qsum);
          sumPeakWave[idet]->SetBinContent(thit.firstBin + 1, sumPeakWave[idet]->GetBinContent(thit.firstBin + 1) + thit.qpeak);
          // only count hits passing cut 
          histHitCount->SetBinContent(tdet->channel, histHitCount->GetBinContent(tdet->channel) + 1 );
        }

        if (thit.qpeak > 100 && thit.qpeak < 300 && thit.startTime > 800)
        {
          for (unsigned jbin = thit.firstBin; jbin < thit.lastBin; ++jbin)
          {
            int fillBin = hQPEShape[idet]->GetNbinsX() / 2 + jbin - thit.peakBin;
            double val = double(rawBr[idet]->rdigi[jbin]) - tdet->base;
            hQPEShape[idet]->SetBinContent(fillBin, hQPEShape[idet]->GetBinContent(fillBin) + val);
          }
        }
      }
    }
    // printf(" event %llu  pass %i fail 1 %i cosmic only %i fail both %i \n",entry, int(hEventPass->GetBinContent(1)), int(hEventPass->GetBinContent(2)), int(hEventPass->GetBinContent(3)), int(hEventPass->GetBinContent(4)));
    ntChanSum->Fill(&fsum[0]); // fill sumHitWave and Q sums
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

void anaCRun::negativeCrossingCount(int ichan)
{
  crossings.clear();
  Double_t cut = 10. * channelSigmaValue[ichan];
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
  Double_t cut = thresh;
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
  printf("\t\t start of file %i %i %i : %i\n", rawEventData->day, rawEventData->mon, rawEventData->year, rawEventData->hour);
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

  fout = new TFile(Form("caenData/anaCRun-%s-%llu.root", shortName.c_str(), maxEntries), "recreate");
  evDir = fout->mkdir("evDir");
  pmtDir = fout->mkdir("pmtDir");
  badDir = fout->mkdir("badDir");
  badTrigDir = fout->mkdir("badTrigDir");
  fout->cd();
  cout << " opened output file " << fout->GetName() << endl;
  getSummedHists();
  // fout->ls();

  // make output tree
  tbrun = new TBRun(tag);
  // and event time
  eventData = new TBEventData();
  tbrun->btree->Branch("eventData",&eventData);

  for (unsigned it = 0; it < rawBr.size(); ++it)
  {
    tbrun->addDet(it);
  }

  // fout->append(tbrun->btree);
  int totalchannels = rawBr.size() + 1;
  ntChan = new TNtuple("ntchan", "channel ntuple", "trig:chan:ave:sigma:skew:base:peakmax:sum2:sum:negcrossings:thresholds:pass");
  ntChanSum = new TNtuple("ntchansum", "channel ntuple", "sum0:sum1:sum2:sum3:sum4:sum5:sum6:sum7:sum8:sum9:sum10:sum11:sum12:pass");
  hEventPass = new TH1D("EventPass", " event failures", 16, 0, 16);
  evCount = new TH1D("eventcount", "event count", totalchannels, 0, totalchannels);
  histHitCount = new TH1D("hitCount", "hit count by channel", totalchannels, 0, totalchannels);
  histQSum = new TH1D("histqsum", "qsum by channel", totalchannels, 0, totalchannels);
  // nn/histqpe = new th1d("histqpe", "qpe by channel", totalchannels, 0, totalchannels);
  histQPrompt = new TH1D("histqprompt", "qprompt by channel", totalchannels, 0, totalchannels);
  histQSum->Sumw2();
  histQPrompt->Sumw2();
  hSumPMT = new TH1D("SumPMT","SumPMT",1100,-100,2000);

  //
  anaDir = fout->mkdir("anadir");
  anaDir->cd();
  for (unsigned i = 0; i < rawBr.size(); ++i)
  {
    unsigned ichan = i;
    hChannelGaus.push_back(new TH1D(Form("channelGaus%i", ichan), Form("channelGaus%i", ichan ), 600, -100, 500));
    noiseHist.push_back(new TH1D(Form("noiseChan%i", ichan), Form("noiseChan%i", ichan), 1000, 0, 1000));
    skewHist.push_back(new TH1D(Form("skewChan%i", ichan), Form("skewChan%i", ichan), 200, -3, 7));
    if (ichan > 8 && ichan<12)
    {
      valHist.push_back(new TH1D(Form("valChan%i", ichan), Form("valChan%i", ichan), 1500, -500, 1000));
      valHistB.push_back(new TH1D(Form("valBadChan%i", ichan), Form("valBadChan%i", ichan), 1500, -500, 1000));
      hEvGaus.push_back(new TH1D(Form("evGaus%i", ichan), Form("evGaus%i", ichan), 200, -500, 500));
      baseHist.push_back(new TH1D(Form("baseChan%i", ichan), Form("baseChan%i", ichan), 200,-10000, 1000));
    }
    else
    {
      baseHist.push_back(new TH1D(Form("baseChan%i", ichan), Form("baseChan%i", ichan), 200,-100, 100));
      valHist.push_back(new TH1D(Form("valChan%i", ichan), Form("valChan%i", ichan), 1000, -200, 200));
      valHistB.push_back(new TH1D(Form("valBadChan%i", ichan), Form("valBadChan%i", ichan), 1000, -200, 200));
      hEvGaus.push_back(new TH1D(Form("evGaus%i", ichan), Form("evGaus%i", ichan), 200, -100, 100));
    }
  }
  cosmicCut1 = new TH1D("cosmicCut1", " cosmic total sum chan 12 ", 100, 0, 100);
  cosmicCut2 = new TH1D("cosmicCut2", " cosmic late large hit chan 12 ", 200, 0, 2000);
  threshHist = new TH1D("threshHist"," threshold crossings trig channels ", 20, 0, 20);
  crossHist = new TH1D("crossHist", "  negative crossings non trigger channels", 100, 0, 100);
  sumDir = fout->mkdir("sumDir");
  sumDir->cd();
  double limit = 100000;
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
    hQPEShape.push_back(new TH1D(Form("QPEShapeChan%i", ichan), Form("QPEShapeChan%i", ichan), 400, -100, 300));
    hQPEShape[hQPEShape.size() - 1]->SetMarkerStyle(20);

    hQSum.push_back(new TH1D(Form("QSumChan%i", ichan), Form("QSumChan%i", ichan), 1000, 0, limit));
    hQPeak.push_back(new TH1D(Form("QPeakChan%i", ichan), Form("QPeakChan%i", ichan), 1000, 0, plimit));
    sumWave.push_back(new TH1D(Form("sumWave%i", ichan), Form("sumWave%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));
    sumWaveB.push_back(new TH1D(Form("sumWaveBad%i", ichan), Form("sumWaveBad%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));
    sumHitWave.push_back(new TH1D(Form("sumHitWave%i", ichan), Form("sumHitWave%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));
    sumPeakWave.push_back(new TH1D(Form("sumPeakWave%i", ichan), Form("sumPeakWave%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));
  }
  fout->cd();

  // fout->ls();
  cout << " make hitFinder dets = " << rawBr.size() << "  size " << rawBr[0]->rdigi.size() << endl;
  //finder = NULL;
  vector<int> chanList;
  for (int ichan = 0; ichan<rawBr.size(); ++ichan)
    chanList.push_back(ichan);
  finder = new hitFinder(fout, tbrun, tag, rawBr[0]->rdigi.size(), chanList, channelSigmaValue);

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

    hEventPass->Fill(-1);
    if (anaEvent(entry))
      ++npass;
    else
      ++nfail;
    //if(passBit!=0)
    //  printf("event fails with passBit = %x npass %i nfail %i \n", passBit,npass,nfail);
    hEventPass->Fill(passBit);
  }
  // normailize to number of nentries.
  // loop over detector channels
  //double scaleFactor = 1. / double(npass);
  printf(" \n \n At END OF FILE total pass  = %i  \n",npass);
  // do not scale here
  /*
  histQSum->Scale(scaleFactor);
  histQPrompt->Scale(scaleFactor);
  for (unsigned i = 0; i < sumHitWave.size(); ++i)
  {
    sumWave[i]->Scale(scaleFactor); // DEF added
    sumHitWave[i]->Scale(scaleFactor);
    sumPeakWave[i]->Scale(scaleFactor);
    hQSum[i]->Scale(scaleFactor);
    hQPeak[i]->Scale(scaleFactor);
  }
  */

  // finder->hPeakCount->Print("all");
  // tbrun->btree->GetListOfBranches()->ls();
  // fout->ls();
  /* do fits at end of run */

  TString graphName = TString("slopeGraph");
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

  // get channel sigma

  channelSigma.resize(chanList.size());
  channelSigmaErr.resize(chanList.size());
  for (int index = 0; index < hChannelGaus.size(); ++index)
  {
    printf(" fit to %s  %i ", hChannelGaus[index]->GetName(), int(hChannelGaus[index]->GetEntries()));
    hChannelGaus[index]->Fit("gaus", "QO", "", hChannelGaus[index]->GetMean() - 100, hChannelGaus[index]->GetMean() + 100);
    TF1 *gfit = (TF1 *)hChannelGaus[index]->GetListOfFunctions()->FindObject("gaus");
    double sigma = hChannelGaus[index]->GetRMS();
    double sigmaErr = 0;
    if (gfit != nullptr)
    {
      sigma = gfit->GetParameter(2);
      sigmaErr = gfit->GetParError(2);
    }
    printf(" chan point ind %i channe %i %f %f \n", index, chanList[index], sigma, sigmaErr);
    channelSigma[index] = sigma;
    channelSigmaErr[index] = sigmaErr;
  }
  for (int index = 0; index < hChannelGaus.size(); ++index)
  {
    printf(" chan %f %i channelSigma %f %f \n", chan[index], chanList[index], channelSigma[index], channelSigmaErr[index]);
  }

  if (hChannelGaus.size() > 0)
  {
    graphName = TString("channelSigmaGraph");
    graphTitle = TString(Form("channel sigma graph-%s", shortName.c_str()));
    printf(" making channel sigma graph %s %lu %lu \n", graphName.Data(), chan.size(), channelSigma.size());
    TGraphErrors *grChannelSigma = new TGraphErrors(channelSigma.size(), &chan[0], &channelSigma[0], &echan[0], &channelSigmaErr[0]);
    grChannelSigma->SetName(graphName);
    grChannelSigma->SetTitle(graphTitle);
    grChannelSigma->SetMarkerStyle(21);
    grChannelSigma->SetLineStyle(0);
    fout->Append(grChannelSigma);
  }
  printf(" \t\t *** hit count by channel ***  \n");
  for (int ibin = 0; ibin < histHitCount->GetNbinsX(); ++ibin)
    printf(" bin %i count %i ;", ibin, int(histHitCount->GetBinContent(ibin)));
  printf("  \n");

  printf(" ******* FINISHED***** \n \t hits by channel  \n");
  for (int i = 0; i < 13; ++i)
  {
    printf("channel %i %f  \n", i, hQSum[i]->GetEntries());
  }
    // printf(" FINISHED npass %u nfail %u output file  %s \n", npass, nfail, fout->GetName());
    printf(" FINISHED pass %i  fail %i thresh only %i cross only %i 1st cosmic only %i 2nd cosmic only %i output file %s  \n",
           int(hEventPass->GetBinContent(1)),
           int(nfail),
           int(hEventPass->GetBinContent(1+1)),
           int(hEventPass->GetBinContent(2+1)),
           int(hEventPass->GetBinContent(4+1)),
           int(hEventPass->GetBinContent(8+1)),
           fout->GetName());

    //hEventPass->Print("all");

    fout->Write();
    fout->Close();
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
