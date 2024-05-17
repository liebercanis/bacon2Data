////////////////////////////////////////////////////g
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
  // add two for summed waveform 0-12 channels, 13 summed trigger 14 summed nontrigger
  // NONSUMCANNELS are nonsummed
  enum
  {
    CHANNELS = 14,
    NONSUMCHANNELS = CHANNELS - 1
  };
  enum
  {
    WAVELENGTH = 7500
  };

  TFile *fout;
  TFile *fin;
  TTree *rawTree;
  TNtuple *ntHit;
  // vectors for gains
  std::vector<double> sipmGain;
  std::vector<double> sipmGainError;
  //
  std::map<int, int> chanMap;
  vector<TBRawEvent *> rawBr;
  TBEventData *eventData;
  TBEventData *rawEventData;
  TBRun *tbrun;
  TNtuple *ntChan;
  TNtuple *ntChanSum;
  TNtuple *ntTrigTime;
  TNtuple *ntSpeYield;
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
  vector<TH1D *> hSingletShape;
  TH1D *hEvBaseWave;
  vector<TH1D *> hEvGaus;
  vector<TH1D *> hChannelGaus;

  // for sums needed for gains
  std::vector<TH1D *> hTotSum;
  std::vector<TH1D *> hPreSum;
  std::vector<TH1D *> hTrigSum;
  std::vector<TH1D *> hLateSum;

  TH1D *hPreQpeak;
  TH1D *hLateQpeak;
  TH1D *hCountPre;
  TH1D *hCountLate;
  TH1D *hCountLateTime;
  TH2D *hCountLateTimeQpeak;
  TH1D *evCount;
  TH1D *histQSum;
  TH1D *hEventPass;
  TH1D *histHitCount;
  TH1D *hNoPeak;
  TH1D *hSumPMT;
  TH1D *threshHist;
  TH2D *threshValueHist;
  TH1D *crossHist;
  TH1D *cosmicCut1;
  TH1D *cosmicCut2;
  // TH1D *histQPE;
  TH1D *histQPrompt;
  TH1D *hTriggerTimeAll;
  TH1D *hTriggerHitTimeAll;
  TH1D *hTriggerTime;
  TH1D *hTriggerTimeTrig;
  TH1D *hTriggerShift;
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
  Long64_t anaCRunFile(TString theFile, Long64_t maxEntries, Long64_t firstEntry = 0);
  void clear();
  bool openFile(TString fileName);
  unsigned getListOfFiles(TString dir);
  bool readGains(TString fileName);
  void getSummedHists();
  unsigned getBranches();
  int anaEvent(Long64_t entry);                   // return passBit
  void differentiate(unsigned diffStep);          // not used
  void derivativeCount(TDet *idet, Double_t rms); // not used
  void negativeCrossingCount(int ichan);
  void thresholdCrossingCount(double thresh);
  std::vector<double> sumDigi();
  unsigned getTriggerTime(int ichan, double &adc);
  void getTriggerTimeStats(unsigned *timeArray, double &ave, double &sigma);
  unsigned fixedTriggerTime(int ichan, double &adc);
  void doTimeShiftAndNorm();
  std::vector<std::vector<double>> fixedDigi; // all the fixed waveforms
  std::vector<unsigned> trigTimes;
  std::vector<unsigned> sTrigTimes; // after correction
  std::vector<double> adcBin;
  std::vector<double> speCount;

  TDirectory *earlyPeakDir;
  TDirectory *rawSumDir;
  TDirectory *badDir;
  TDirectory *badTrigDir;
  TDirectory *evDir;
  TDirectory *sumDir;
  TDirectory *anaDir;
  TDirectory *pmtDir;
  Long64_t nentries;
  double QPEPeak;
  //
  unsigned trigStart = 600;
  unsigned trigEnd = 800;
  int nominalTrigger = 729;
  double nominalGain = 227.4; // average
  unsigned firstTime;         // corrected trigger time for event
  unsigned timeOffset = 13; // changed from 17 may 13, 2024
  ULong_t taveEarlyCut = 710;
  ULong_t taveLateCut = 740;
  ULong_t timeEarlyCut = 690;
  ULong_t timeLateCut = 890;
  double prePeakCut = 0.5;
  double latePeakCut = 3.5; // march 18 2024 2.5;
  double diffStep = 1.;     // back to one from 4 May 4 2024
  };

/** get trigger time **/
unsigned anaCRun::getTriggerTime(int ic, double &adc)
{
  TDet *idet = tbrun->getDet(ic);
  unsigned time = 0;
  for (unsigned j = 0; j < 801; ++j)
  {
    double val = double(rawBr[ic]->rdigi[j]) - idet->base;
    val *= nominalGain / sipmGain[ic];
    if (val > 0.5 * nominalGain)
    {
      adc = val;
      time = j;
      break;
    }
  }
  return time;
}

void anaCRun::getTriggerTimeStats(unsigned *timeArray, double& ave, double& sigma){
  unsigned nave = 0;
  ave = 0;
  sigma = 0;
  // calculate ave
  for (unsigned ic = 0; ic < 3; ++ic)
  {
    // 690 < time < 890
    if (timeArray[ic] < timeLateCut && timeArray[ic]>timeEarlyCut)
    {
      ave += double(timeArray[ic]);
      ++nave;
    }
  }
  // will cast as unsigned 
  if(nave>0&&ave>0) 
    ave /= double(nave);
  else 
    ave = nominalTrigger;

  if(nave==0)
    return;
  // calculate sigma
  for (unsigned ic = 0; ic < 3; ++ic)
  {
      sigma  += pow(double(timeArray[ic]-ave),2.);
  }
  sigma = sqrt(sigma) / double(nave);
}

/** get trigger time **/
unsigned anaCRun::fixedTriggerTime(int ic, double &adc)
{
  unsigned time = 0;
  for (unsigned j = 0; j < 801; ++j)
  {
    double val = fixedDigi[ic][j];
    if (val > 0.5 * nominalGain)
    {
      adc = val;
      time = j;
      break;
    }
  }
  return time;
}

// shift the times and apply gain norm void anaCRun::doTimeShiftAndNorm()
// bug fixed
void anaCRun::doTimeShiftAndNorm()
{
  fixedDigi.clear();
  // timeShift = 0; // for debugging!!
  // loop over all branches
  std::vector<double> fDigi;
  fDigi.clear();
  fDigi.resize(rawBr[0]->rdigi.size());
  std::fill(fDigi.begin(), fDigi.end(), ULong_t(0));
  for (unsigned ib = 0; ib < NONSUMCHANNELS; ++ib)
  {
    // printf(" first %u  shift %i off %u chan %u \n",firstTime, timeShift,timeOffset,ib);
    int timeShift = nominalTrigger - firstTime;
    if (ib < 9)
      timeShift += int(timeOffset); // amplifier delay
    // after doing time shift set jstsart,jstop,absShift
    ULong_t jstart = TMath::Max(-timeShift, 0);
    ULong_t jstop = TMath::Min(int(rawBr[0]->rdigi.size()), int(rawBr[0]->rdigi.size()) - timeShift);
    int absShift = TMath::Abs(timeShift);
    hTriggerShift->Fill(timeShift);
    TDet *idet = tbrun->getDet(ib);
    // printf(" chan %u nominal %i first %i shift %i\n",ib,nominalTrigger,firstTime,timeShift);
    /* take care here for summed ib=CHANNELS-2 and set appropriate hitThreshold */
    // which nominal gain, hit threshold id the samea
    for (ULong_t j = jstart; j < jstop; ++j)
    {
      double val = double(rawBr[ib]->rdigi[j]) - idet->base;
      // scale all channels by nominal gain
      val *= nominalGain / sipmGain[ib];
      if (timeShift > 0)
        fDigi[j + ULong_t(absShift)] = val;
      else
        fDigi[j - ULong_t(absShift)] = val;
    }
    fixedDigi.push_back(fDigi);
  }
}

// we must also do the time shift here
std::vector<double> anaCRun::sumDigi()
{
  std::vector<double> digiSum;
  // fix all the waveforms
  // loop over summed times
  for (unsigned j = 0; j < rawBr[0]->rdigi.size(); ++j)
  {
    double summedDigi = 0;
    // loop over branches < 12 to sum digi at this sample time
    for (unsigned ic = 0; ic < NONSUMCHANNELS - 1; ++ic)
    {
      summedDigi += fixedDigi[ic][j];
    }
    digiSum.push_back(summedDigi);
  }
  return digiSum;
}

bool anaCRun::readGains(TString fileName)
{
  TFile *fin = new TFile(fileName, "readonly");
  if (fin->IsZombie())
  {
    std::cout << "Error opening file" << fileName << std::endl;
    return false;
  }
  cout << " opened sipm gain file " << fileName << endl;
  TGraphErrors *gGain = NULL;
  fin->GetObject("gGain", gGain);
  if (gGain == NULL)
  {
    cout << "no gGain in file " << endl;
    return false;
  }
  cout << "found graph named " << gGain->GetName() << " in file " << fileName << endl;
  sipmGain.clear();
  sipmGainError.clear();
  sipmGain.resize(NONSUMCHANNELS);
  sipmGainError.resize(NONSUMCHANNELS);
  for (unsigned long j = 0; j < sipmGain.size(); ++j)
  {
    sipmGain[j] = nominalGain;
    sipmGainError[j] = sqrt(nominalGain);
  }
  for (int i = 0; i < gGain->GetN(); ++i)
  {
    int index = int(gGain->GetPointX(i));
    sipmGain[index] = gGain->GetPointY(i);
    sipmGainError[index] = gGain->GetErrorY(i);
  }

  printf("stored gains %lu \n", sipmGain.size());
  for (unsigned long j = 0; j < sipmGain.size(); ++j)
  {
    printf(" %lu  gain %.4f error %.4f   \n", j, sipmGain[j], sipmGainError[j]);
  }
  return true;
}

void anaCRun::clear()
{
  hQSum.clear();
  hQPEShape.clear();
  hSingletShape.clear();
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
  sipmGain.clear();
  sipmGainError.clear();
  // fill channel sigma in order of branches
  channelSigmaValue.resize(CHANNELS);
  for (unsigned long j = 0; j < channelSigmaValue.size(); ++j)
    channelSigmaValue[j] = 10;
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
int anaCRun::anaEvent(Long64_t entry)
{
  //printf("start  %lld \n",entry);
  // clear
  speCount.clear();
  speCount.resize(NONSUMCHANNELS);
  std::fill(speCount.begin(), speCount.end(), 0);
  int passBit = 0;
  // previously 40 but channel 9 was missing peaks
  double derivativeThreshold;
  double hitThreshold;
  
  eventData->evtime = rawEventData->evtime;
  eventData->sec = rawEventData->sec;
  eventData->min = rawEventData->min;
  eventData->hour = rawEventData->hour;
  eventData->day = rawEventData->day;
  eventData->mon = rawEventData->mon;
  eventData->year = rawEventData->year;
  eventData->isdst = rawEventData->isdst;
  QPEPeak = 100;
  // loop over channels but not including summed channel
  for (unsigned ib = 0; ib < NONSUMCHANNELS; ++ib)
  {
    unsigned ichan = ib;
    // define trigger sipms
    bool trig = ichan == 9 || ichan == 10 || ichan == 11;
    // deal with trigger channel sign by overwriting rdigi
    if (trig)
    {
      for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
      {
        rawBr[ib]->rdigi[j] = -1. * (rawBr[ib]->rdigi[j] - pow(2, 14)); // base is > digi value!
      }
    }

    int nbins = rawBr[ib]->rdigi.size();
    // cout << ib << " nbins " << nbins << " max hist " << hEvGaus.size() << " rawBr.size() " << rawBr.size() << endl;

    // sanity check
    if (rawBr[ib]->rdigi.size() != WAVELENGTH)
      continue;
    // simple baseline
    std::vector<unsigned short> orderDigi = rawBr[ib]->rdigi;
    std::sort(orderDigi.begin(), orderDigi.end());
    unsigned baseLength = orderDigi.size() / 2;
    double base = 0;
    for (unsigned j = 0; j < baseLength; ++j)
    {
      base += orderDigi[j];
    }
    base /= double(baseLength);

    // baseline correction from fitted Gaussian
    hEvGaus[ib]->Reset("ICES");
    for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
    {
      double val = double(rawBr[ib]->rdigi[j]) - base; // base is > digi value!
      hEvGaus[ib]->Fill(val);
    }
    // get the distribution mode
    double mode = hEvGaus[ib]->GetBinLowEdge(hEvGaus[ib]->GetMaximumBin()) + 0.5 * hEvGaus[ib]->GetBinWidth(hEvGaus[ib]->GetMaximumBin());

    hEvGaus[ib]->Fit("gaus", "QO", "", hEvGaus[ib]->GetMean() - 100, hEvGaus[ib]->GetMean() + 100);
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
      fitMean = ave; // fit mean
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
    idet->base = base + fitMean; // add in fit mean if fit succeeded
    idet->mode = mode;

    /*********
     *  fill baseline subtracted vector digi
     *********/
    digi.clear();
    for (unsigned j = 0; j < rawBr[ib]->rdigi.size() - 1; ++j)
    {
      double val = double(rawBr[ib]->rdigi[j]) - idet->base;
      baseHist[ichan]->Fill(val);
      digi.push_back(val);
      if (hChannelGaus.size() > 0)
        hChannelGaus[ib]->Fill(val);
    }
    idet->pass = true;
    // waveform sums on baseline subtracted waveform
    idet->totSum = 0;
    idet->preSum = 0;
    idet->trigSum = 0;
    idet->lateSum = 0;
    idet->totPeakSum = 0;
    idet->prePeakSum = 0;
    idet->trigPeakSum = 0;
    idet->latePeakSum = 0;
    double peakMax = 0;
    for (unsigned j = 0; j < digi.size(); ++j)
    {
      idet->totSum += digi[j];
      if (j < trigStart)
        idet->preSum += digi[j];
      else if (j < trigEnd)
      {
        idet->trigSum += digi[j];
        if (digi[j] > peakMax)
          peakMax = digi[j];
      }
      else
        idet->lateSum += digi[j];
    }
    // add some other variables
    idet->peak = peakMax;

    ntChan->Fill(float(rawBr[ib]->trigger), float(ichan), float(ave), float(sigma), float(skew), float(base), float(peakMax), float(idet->trigSum), float(idet->totSum), float(crossings.size()), float(thresholds.size()), float(idet->pass));

  } // channel loop

  /* find trigger time from trigger sipms */
  trigTimes.resize(12);
  sTrigTimes.resize(12);
  adcBin.resize(12);
  for (unsigned ic = 0; ic < 12; ++ic)
  {
    double val = 0;
    unsigned time = getTriggerTime(ic, val);
    trigTimes[ic] = time;
    adcBin[ic] = val;
    // if (time < nominalTrigger - 5 * triggerSigma)
    if (ic < 9 && time > 0) // add in amplifier delay zero not found
      time += timeOffset;
    if (time > 0 && time < unsigned(taveEarlyCut))
    {
      //printf(" failed timeEarlyCut event %llu chan %i cut %lu time %u \n", entry, ic, taveEarlyCut, time);
      passBit |= 0x1;
    }
    if (time > 0)
      hTriggerTimeAll->Fill(double(time));
    if (time > 0 && ic>8 )
      hTriggerTime->Fill(double(time));
  }

  /* ave trigger sipm times before shift */
  double trigTimeAve = 0;
  double trigTimeSigma = 0;
  getTriggerTimeStats(&trigTimes[9], trigTimeAve, trigTimeSigma);
  firstTime = unsigned(trigTimeAve);

  /* ave non trig  sipm times before shift */
  double nonTimeAve = 0;
  double nonTimeSigma = 0;
  getTriggerTimeStats(&trigTimes[6], nonTimeAve, nonTimeSigma);

  if (firstTime > taveLateCut)
  {
    printf(" failed timeLateCut event %llu cut %lu time %u \n", entry, timeLateCut, firstTime);
    passBit |= 0x2;
  }

  
  hTriggerTime->Fill(double(firstTime));
  int timeShift = nominalTrigger - firstTime;
  hTriggerShift->Fill(timeShift);

  /*******
   * now that we have the firstTime
        align to nominalTrigger
        defined as timeShift>0 shift right
        normalize to nominal gain
   ********/
  //printf("doTimeSiftAndNorm %lld \n",entry);
  doTimeShiftAndNorm();
  /* make ntuple of before and after shift */
  for (unsigned ic = 6; ic < 12; ++ic)
  {
    double val;
    sTrigTimes[ic]=fixedTriggerTime(ic, val);
  }
  /* ave trigger sipm times after shift */
  double trigTimeAve2 = 0;
  double trigTimeSigma2 = 0;
  getTriggerTimeStats(&sTrigTimes[9], trigTimeAve2, trigTimeSigma2);

  /* ave non trig  sipm times after shift*/
  double nonTimeAve2 = 0;
  double nonTimeSigma2 = 0;
  getTriggerTimeStats(&sTrigTimes[6], nonTimeAve2, nonTimeSigma2);

  ntTrigTime->Fill(entry, 
                   double(firstTime), trigTimeSigma, nonTimeAve,nonTimeSigma,
                   trigTimeAve2, trigTimeSigma2, nonTimeAve2, nonTimeSigma2 );

  /* *******
        do pulse finding on summed line for event cuts
  ****   */
  digi.clear();
  // hitThreshold = 0.74 * nominalGain;
  digi = sumDigi();
  TDet *tdet = tbrun->getDet(NONSUMCHANNELS);
  tdet->hits.clear();
  derivativeThreshold = 30; // for summed waveform
  //printf("call to finder for chan 13 %lld \n",entry);
  hitThreshold = 0.25 * nominalGain; // for summed waveform
  finder->event(NONSUMCHANNELS, entry, digi, derivativeThreshold, hitThreshold, diffStep); // DEG suggests 10
  // which nominal gain, hit threshold id the same

  // event cuts
  int nPreHits = 0;
  int nLateHits = 0;
  for (unsigned ihit = 0; ihit < tdet->hits.size(); ++ihit)
  {
    TDetHit hiti = tdet->hits[ihit];
    ULong_t hitStartTime = ULong_t(tdet->hits[ihit].startTime);
    hCountLateTimeQpeak->Fill(hitStartTime, tdet->hits[ihit].qpeak / nominalGain);
    if (hitStartTime < timeEarlyCut)
      hPreQpeak->Fill(tdet->hits[ihit].qpeak / nominalGain);
    if (hitStartTime < timeEarlyCut && tdet->hits[ihit].qpeak / nominalGain > prePeakCut)
    {
      ++nPreHits;
      // printf("event preHits %llu cut %lu hitStartTime %lu  qpeak %.2f nPreHits %i \n", entry, timeEarlyCut, hitStartTime, hiti.qpeak, nPreHits);
    }
    if (hitStartTime > timeLateCut)
      hLateQpeak->Fill(tdet->hits[ihit].qpeak / nominalGain);
    if (hitStartTime > timeLateCut && tdet->hits[ihit].qpeak / nominalGain > latePeakCut)
    {
      ++nLateHits;
      hCountLateTime->Fill(hitStartTime);
      // printf("event lateHits %llu cut %lu hitStartTime %lu  qpeak %.2f nLateHits %i \n", entry, timeLateCut, hitStartTime, hiti.qpeak, nLateHits);
      //  hCountLateTime->Fill(tdet->hits[ihit].startTime);
    }
    // if (hitStartTime > 600 && hitStartTime < 800 && hitStartTime < firstTime)
    //   firstTime = hitStartTime;
  }
  hCountPre->Fill(nPreHits);
  hCountLate->Fill(nLateHits);

  if (nPreHits > 0)
    passBit |= 0x4;
  if (nLateHits > 0)
  {
    // printf(" failed nLate event %llu nLate %i \n", entry,nLateHits);
    passBit |= 0x8;
  }

  evCount->Fill(-1); // underflow bin
  if (passBit != 0)
  {
    printf("det %i nhits %u \n", NONSUMCHANNELS, tbrun->detList[NONSUMCHANNELS]->nhits());
    //return passBit;
  }

  // continue if event passes

  /*******
        start of second channel loop doing pulse finding
   ********/

  for (unsigned ib = 0; ib < NONSUMCHANNELS; ++ib)
  {
    unsigned ichan = ib;
    TDet *tdet = tbrun->getDet(ib);
    tdet->hits.clear();
    bool trig = ichan == 9 || ichan == 10 || ichan == 11;
    int nbins = rawBr[ib]->rdigi.size();

    /* take care here for summed ib=CHANNELS-2 and set appropriate hitThreshold */
    digi.clear();
    digi = fixedDigi[ib];
    
    // used time aligned and gain normaized waveforms
    // fill histograms
    for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
    {
      if (tdet->pass)
      {
        sumWave[ib]->SetBinContent(j + 1, sumWave[ib]->GetBinContent(j + 1) + digi[j]);
        valHist[ib]->Fill(digi[j]);
        // histogram//  bad events
      }
      else
      {
        sumWaveB[ib]->SetBinContent(j + 1, sumWaveB[ib]->GetBinContent(j + 1) + digi[j]);
        valHistB[ib]->Fill(digi[j]);
      }
    }

    evCount->Fill(ib); // chan 0 from GetBinContent(0)

    
    // tried to fix gap but didnt work
    //if (ichan >8 ) //lower for trigger sipms
    //  derivativeThreshold = 10.;
    derivativeThreshold = 20;              // for non summed
    hitThreshold = 0.25 * nominalGain;     // for non summed
    finder->event(ichan, entry, digi, derivativeThreshold, hitThreshold,diffStep); // DEG suggests 10

    TDirectory *fftDir = (TDirectory *)fout->FindObject("fftDir");
    if (!fftDir)
    {
      cout << " Error no fftDir" << endl;
      fout->ls();
      return false;
    }
    TDirectory *sumWaveDir = (TDirectory *)fout->FindObject("sumWaveDir");
    if (!sumWaveDir)
    {
      cout << " Error no sumWaveDir" << endl;
      return false;
    }
  } // second channel loop after pulse finding
  if (passBit != 0) {
    printf("event %lld det %i nhits %u \n", entry,NONSUMCHANNELS, tbrun->detList[NONSUMCHANNELS]->nhits());
    return passBit;
  }
  // fill total light
  vector<float> fsum;
  fsum.resize(tbrun->detList.size());
  // loop over detector channels
    for (unsigned idet = 0; idet < tbrun->detList.size(); ++idet)
    {
      TDet *tdet = tbrun->detList[idet];
      // printf(" anaCRuna::event at event %llu idet %i chan %i hits %lu \n", entry, idet, tdet->channel, tdet->hits.size());
      fsum[tdet->channel] = tdet->totSum;

      // add some event plots
      bool trig = tdet->channel == 9 || tdet->channel == 10 || tdet->channel == 11;
      TDirectory *fftDir = (TDirectory *)fout->FindObject("fftDir");
      TDirectory *finderDir = (TDirectory *)fout->FindObject("finderDir");
      if (finderDir->GetList()->GetEntries() < 2000)
      {
        // count late hits
        int lateHits = 0;
        int earlyHits = 0;
        int thitStartTime = 0;
        TDetHit tlateHit;
        for (unsigned ihit = 0; ihit < tdet->hits.size(); ++ihit)
        {
          TDetHit thit = tdet->hits[ihit];
          if (thit.startTime > 4000)
          {
            ++lateHits;
            tlateHit = thit;
          }
          if (thit.startTime > 700 && thit.startTime < 720)
            ++earlyHits;
          thitStartTime = thit.startTime;
        }
        if (earlyHits > 0)
        {
          // printf(" found %i earlyHits det %i time %i event %lld \n", earlyHits, tdet->channel,thitStartTime, entry);
          if (earlyPeakDir->GetList()->GetEntries() < 100)
          {
            finder->plot1Wave(earlyPeakDir, tdet->channel, entry);
            finder->plot1Wave(earlyPeakDir, 9, entry);
          }
        }

        /*
        if (lateHits > 0 && tdet->channel < 9)
        {
          finder->plotEvent(finderDir, tdet->channel, entry);
          // printf("\t plotEvent %lld late %i start time %f peak %f \n", entry,lateHits,tlateHit.startTime, tlateHit.qpeak);
        }
        if (tdet->channel ==9)
        {
          finder->plotEvent(finderDir, tdet->channel, entry);
        }
        */
      }
      // finder->plotEvent(fftDir, 8, entry);

      TDirectory *sumWaveDir = (TDirectory *)fout->FindObject("sumWaveDir");
      if (idet == 13 && tdet->hits.size() > 1 && sumWaveDir->GetList()->GetEntries() < 5000)
      {
        // printf("xxxxxxx anaCRun::event event %llu chan %i hits %lu der thresh %f hit thresh %f \n", entry, tdet->channel, tdet->hits.size(), derivativeThreshold, hitThreshold);
        finder->plotEvent(sumWaveDir, tdet->channel, entry);
      }

      if (trig && tdet->hits.size() == 0 && fftDir->GetList()->GetEntries() < 2000)
      {
        // printf("!!!!!! anaCRuna::event plot event %llu idet %i chan %i hits %lu \n", entry, idet, tdet->channel, tdet->hits.size());
        finder->plotEvent(fftDir, tdet->channel, entry);
      }

      // loop over hits
      // printf(" event %llu  det %u nhits %lu \n", entry, idet, tdet->hits.size());
      // add peak sums
      if (tdet->hits.size() == 0)
        hNoPeak->SetBinContent(tdet->channel + 1, hNoPeak->GetBinContent(tdet->channel + 1) + 1);
      int firstHitTime = rawBr[NONSUMCHANNELS]->rdigi.size();
      histQSum->SetBinContent(tdet->channel + 1, histQSum->GetBinContent(tdet->channel + 1) + tdet->qpeak);
      histQPrompt->SetBinContent(tdet->channel + 1, histQPrompt->GetBinContent(tdet->channel + 1) + tdet->hitPrompt);
      // printf(" event %lld det %i sum qpeak %f sum qprompt %f\n", entry, idet, tdet->qpeak, tdet->hitPrompt);
      for (unsigned ihit = 0; ihit < tdet->hits.size(); ++ihit)
      {
        TDetHit thit = tdet->hits[ihit];
        if (thit.qpeak < 1)
          printf("line822 chan %i ihit %i startTime %i  peak %f\n", tdet->channel, ihit, int(thit.startTime), thit.qpeak);
        hQSum[idet]->Fill(thit.qsum);
        hQPeak[idet]->Fill(thit.qpeak);
        unsigned hitTime = unsigned(thit.startTime);
        // do peak sums
        tdet->totPeakSum += thit.qpeak;

        if (hitTime > timeEarlyCut && hitTime < timeLateCut && hitTime < firstHitTime)
          firstHitTime = hitTime;
        //
        if (hitTime < trigStart)
          tdet->prePeakSum += thit.qpeak;
        else if (hitTime < trigEnd)
        {
          tdet->trigPeakSum += thit.qpeak;
        }
        else
          tdet->latePeakSum += thit.qpeak;
        // fill here for gains
        if (hitTime > trigEnd)
          hLateSum[idet]->Fill(thit.qpeak);

        // do threshold for summed waveform
        // if (thit.qsum > hitThreshold)
        sumHitWave[idet]->SetBinContent(thit.firstBin + 1, sumHitWave[idet]->GetBinContent(thit.firstBin + 1) + thit.qsum);
        sumPeakWave[idet]->SetBinContent(thit.firstBin + 1, sumPeakWave[idet]->GetBinContent(thit.firstBin + 1) + thit.qpeak);
        histHitCount->SetBinContent(tdet->channel + 1, histHitCount->GetBinContent(tdet->channel + 1) + 1);
        ntHit->Fill(double(entry), double(passBit), double(idet), thit.startTime, thit.peakt, thit.qpeak / nominalGain);
        // sum of photons in SPE for this channel
        speCount[idet] += thit.qpeak / nominalGain;

        // singlet shapes
        if (thit.qpeak < 0.5 * nominalGain)
        {
          if (thit.startTime > trigEnd)
          {
            for (unsigned jbin = thit.firstBin; jbin < thit.lastBin; ++jbin)
            {
              int fillBin = hQPEShape[idet]->GetNbinsX() / 2 + jbin - thit.peakBin;
              double val = double(rawBr[idet]->rdigi[jbin]) - tdet->base;
              hQPEShape[idet]->SetBinContent(fillBin, hQPEShape[idet]->GetBinContent(fillBin) + val);
            }
          }

          if (thit.startTime > trigStart && thit.startTime < trigEnd)
          {
            for (unsigned jbin = thit.firstBin; jbin < thit.lastBin + 100; ++jbin)
            {
              int fillBin = hSingletShape[idet]->GetNbinsX() / 2 + jbin - thit.peakBin;
              double val = double(rawBr[idet]->rdigi[jbin]) - tdet->base;
              hSingletShape[idet]->SetBinContent(fillBin, hSingletShape[idet]->GetBinContent(fillBin) + val);
            }
          }
        }

        /*
        if (trig)
          printf(" anaCRun::event %llu det %i nhits %lu , tot %f pre %f trig %f late %f\n", entry, tdet->channel, tdet->hits.size(),
                 tdet->totPeakSum, tdet->prePeakSum, tdet->trigPeakSum, tdet->latePeakSum);
         */
      } // hit loop
      hTriggerHitTimeAll->Fill(firstHitTime);
      // fill sums do not fill for zero
      if (tdet->totPeakSum > 0)
        hTotSum[idet]->Fill(tdet->totPeakSum);
      if (tdet->prePeakSum > 0)
        hPreSum[idet]->Fill(tdet->prePeakSum);
      if (tdet->trigPeakSum > 0)
        hTrigSum[idet]->Fill(tdet->trigPeakSum);
      // for each peak not summed
      // if (tdet->latePeakSum>0)
      //  hLateSum[idet]->Fill(tdet->latePeakSum);
      // printf(" anaCRun::event %llu det %i nhits %lu , tot %f pre %f trig %f late %f\n", entry, tdet->channel, tdet->hits.size(),
      //       tdet->totPeakSum, tdet->prePeakSum, tdet->trigPeakSum, tdet->latePeakSum);

    } // det loop

    // printf(" event %llu  pass %i fail 1 %i cosmic only %i fail both %i \n",entry, int(hEventPass->GetBinContent(1)), int(hEventPass->GetBinContent(2)), int(hEventPass->GetBinContent(3)), int(hEventPass->GetBinContent(4)));
    ntChanSum->Fill(&fsum[0]); // fill sumHitWave and Q sums
    ntSpeYield->Fill(entry,
                     speCount[0], speCount[1], speCount[2],
                     speCount[3], speCount[4], speCount[5],
                     speCount[6], speCount[7], speCount[8],
                     speCount[9], speCount[10], speCount[11]);
    // printf(" !!!!! end of event %llu event returns with pass bit  %x \n",entry, eventPass);
    //  was debugging
    /*
    for (unsigned idet = 0; idet < tbrun->detList.size(); ++idet)
    {
      TDet *tdet = tbrun->detList[idet];
      if(idet>8&&idet<12)
        printf(" anaCRun::event %llu det %i nhits %lu , tot %f pre %f trig %f late %f\n",
        entry, tdet->channel, tdet->hits.size(),
        tdet->totPeakSum, tdet->prePeakSum, tdet->trigPeakSum, tdet->latePeakSum);
    }
   */
    // printf("finished  %lld pass %i \n",entry, passBit);
    return passBit;
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
  for (unsigned ibin = 1500; ibin < digi.size(); ++ibin)
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

Long64_t anaCRun::anaCRunFile(TString theFile, Long64_t maxEntries, Long64_t firstEntry)
{
  clear();
  currentBuffer = -1;
  currentBufferCount = 0;
  printf(" anaCRun::anaCRunFile starting anaCRun file %s maxEntries %llu firstEntry %llu \n",
         theFile.Data(), maxEntries, firstEntry);
  if (!openFile(theFile)) // and get branches
  {
    printf("cannot open file!\n");
    return 0;
  }

  string sfilename(theFile.Data());
  string shortName = sfilename.substr(0, sfilename.find_last_of("."));
  cout << " anaCRunFile  with shortName= " << shortName << endl;
  // new gain file
  TString gainFileName = TString("gains-2024-02-15-17-26-save.root");
  cout << "read gains from file " << gainFileName << endl;
  readGains(gainFileName);

  // need to fill rawBr[0]->rdigi.size()
  printf("Read zeroth entry from tree \n");
  if (!rawTree)
  {
    printf("EEEEEE rawTree is null!!!!!\n");
    return 0;
  }
  cout << " RawTree still has has " << rawTree->GetEntries() << " entries " << endl;
  cout << " rawTree return " << rawTree->GetEntry(0) << endl;
  printf("got rawTree entry 0 \n");
  printf("\n\n\t\t >>>>>>>>> start of file %i %i %i : %i <<<<<<<<<<<< \n", rawEventData->day, rawEventData->mon, rawEventData->year, rawEventData->hour);
  printf("\t\t SIZE OF WAVEFORM = %lu \n", rawBr[0]->rdigi.size());
  if (rawBr[0]->rdigi.size() != WAVELENGTH)
  {
    printf(" \n\n\n\n ERROR rdigi size %lu !!! \n", rawBr[0]->rdigi.size());
    // return 0;
  }
  Long64_t nentries = rawTree->GetEntries();
  if (maxEntries > 0)
    nentries = TMath::Min(maxEntries, nentries);
  printf("... total entries  %llu looping over %llu starting from %llu \n ", rawTree->GetEntries(), nentries, firstEntry);

  // open outout file and make histograms

  fout = new TFile(Form("caenData/anaCRun-%s-%llu.root", shortName.c_str(), maxEntries), "recreate");
  evDir = fout->mkdir("evDir");
  pmtDir = fout->mkdir("pmtDir");
  badDir = fout->mkdir("badDir");
  badTrigDir = fout->mkdir("badTrigDir");
  earlyPeakDir = fout->mkdir("earlyPeakDir");
  fout->cd();
  cout << " opened output file " << fout->GetName() << endl;
  getSummedHists();
  // fout->ls();

  // make output tree
  tbrun = new TBRun(tag);
  // and event time
  eventData = new TBEventData();
  tbrun->btree->Branch("eventData", &eventData);

  for (unsigned it = 0; it < rawBr.size(); ++it)
  {
    tbrun->addDet(it);
  }

  // fout->append(tbrun->btree);e("ntHit", " hits
  ntHit = new TNtuple("ntHit", "hit ntuple", "event:flag:chan:time:peakTime:qpeak");
  ntChan = new TNtuple("ntchan", "channel ntuple", "trig:chan:ave:sigma:skew:base:peakmax:sum2:sum:negcrossings:thresholds:pass");
  ntSpeYield = new TNtuple("ntSpeYield", "spe per sipm",
                           "event:spe0:spe1:spe2:spe3:spe4:spe5:spe6:spe7:spe8:spe9:spe10:spe11");
  ntTrigTime = new TNtuple("ntTrigTime", "trigger time ntuple",
    "entry:trigTime:trigTimeSigma:nonTimeAve:nonTimeSigma:trigTimeAve2:trigTimeSigma2:nonTimeAve2:nonTimeSigma2");
  ntChanSum = new TNtuple("ntchansum", "channel ntuple", "sum0:sum1:sum2:sum3:sum4:sum5:sum6:sum7:sum8:sum9:sum10:sum11:sum12:pass");
  hEventPass = new TH1D("EventPass", " event failures", 16, 0, 16);
  evCount = new TH1D("eventcount", "event count", CHANNELS, 0, CHANNELS);
  hNoPeak = new TH1D("noPeak", "no peak events count by channel", CHANNELS, 0, CHANNELS);
  histHitCount = new TH1D("hitCount", "hit count by channel", CHANNELS, 0, CHANNELS);
  histQSum = new TH1D("histqsum", "qsum by channel", CHANNELS, 0, CHANNELS);
  // nn/histqpe = new th1d("histqpe", "qpe by channel", CHANNELS, 0, CHANNELS);
  histQPrompt = new TH1D("histqprompt", "qprompt by channel", CHANNELS, 0, CHANNELS);
  histQSum->Sumw2();
  histQPrompt->Sumw2();

  //
  anaDir = fout->mkdir("anadir");
  anaDir->cd();
  hTriggerTimeTrig = new TH1D("TriggerTimeTrig", " first hit time in trigger Sipm", 800, 0, 800);
  hTriggerTime = new TH1D("TriggerTime", " ave of trigger Sipm times ", 800, 0, 800);
  hTriggerShift = new TH1D("TriggerShift", " ave trigger time shift ", 40, -20, 20);
  hTriggerTimeAll = new TH1D("TriggerTimeAll", " first time all channels ", 800, 00, 800);
  hTriggerHitTimeAll = new TH1D("TriggerHitTimeAll", " first hit time all channels ", 80, 0, 800);
  TString htitle;
  htitle.Form(" pre time < %lu normalized qpeak", timeEarlyCut);
  hPreQpeak = new TH1D("PreQpeak", htitle, 100, 0, 10);
  htitle.Form(" pre time > %lu normalized qpeak", timeLateCut);
  hLateQpeak = new TH1D("LateQpeak", htitle, 100, 0, 10);
  hCountPre = new TH1D("CountPre", " hits sample<600 in sum", 20, 0, 20);
  htitle.Form("hits qpeak>%.2f SPE sample>%luin sum", latePeakCut, timeLateCut);
  hCountLate = new TH1D("CountLate", htitle, 20, 0, 20);
  htitle.Form("number of late time hits with qpeak>%.2f", latePeakCut);
  hCountLate->GetXaxis()->SetTitle(htitle);
  htitle.Form("hits qpeak>%.2f SPE sample>%lu in sum", latePeakCut, timeLateCut);
  hCountLateTime = new TH1D("CountLateTime ", htitle, 30, 0, 7500);
  hCountLateTime->GetXaxis()->SetTitle("sample time");
  hCountLateTime->Sumw2();
  hCountLateTimeQpeak = new TH2D("CountLateTimeQpeak", " sample>890 in sum qpeak vs time ", 30, 0, 7500, 20, 0, 20);
  hCountLateTimeQpeak->GetXaxis()->SetTitle("sample time");
  hCountLateTimeQpeak->GetYaxis()->SetTitle("qpeak [SPE]");

  for (unsigned i = 0; i < rawBr.size(); ++i)
  {
    unsigned ichan = i;
    hChannelGaus.push_back(new TH1D(Form("channelGaus%i", ichan), Form("channelGaus%i", ichan), 600, -100, 500));
    noiseHist.push_back(new TH1D(Form("noiseChan%i", ichan), Form("noiseChan%i", ichan), 1000, 0, 1000));
    skewHist.push_back(new TH1D(Form("skewChan%i", ichan), Form("skewChan%i", ichan), 200, -3, 7));
    if (ichan > 8 && ichan < 12)
    {
      valHist.push_back(new TH1D(Form("valChan%i", ichan), Form("valChan%i", ichan), 1500, -500, 1000));
      valHistB.push_back(new TH1D(Form("valBadChan%i", ichan), Form("valBadChan%i", ichan), 1500, -500, 1000));
      hEvGaus.push_back(new TH1D(Form("evGaus%i", ichan), Form("evGaus%i", ichan), 200, -500, 500));
      baseHist.push_back(new TH1D(Form("baseChan%i", ichan), Form("baseChan%i", ichan), 200, -10000, 1000));
    }
    else
    {
      baseHist.push_back(new TH1D(Form("baseChan%i", ichan), Form("baseChan%i", ichan), 200, -100, 100));
      valHist.push_back(new TH1D(Form("valChan%i", ichan), Form("valChan%i", ichan), 1000, -200, 200));
      valHistB.push_back(new TH1D(Form("valBadChan%i", ichan), Form("valBadChan%i", ichan), 1000, -200, 200));
      hEvGaus.push_back(new TH1D(Form("evGaus%i", ichan), Form("evGaus%i", ichan), 200, -100, 100));
    }

    // for summary //
    double limit = 5.*nominalGain;
    int nbins = 400.;
    hTotSum.push_back(new TH1D(Form("TotPeakSumChan%i", i), Form("tot peak sum chan %i", i), nbins, 0, limit));
    hPreSum.push_back(new TH1D(Form("PrePeakSumChan%i", i), Form("pre peak sum chan %i", i), nbins, 0, limit));
    hTrigSum.push_back(new TH1D(Form("TrigPeakSumChan%i", i), Form("trig peak sum chan %i", i), nbins, 0, limit));
    hLateSum.push_back(new TH1D(Form("LatePeakSumChan%i", i), Form("late peak sum chan %i", i), nbins, 0, limit));
  }
  cosmicCut1 = new TH1D("cosmicCut1", " cosmic total sum chan 12 ", 100, 0, 100);
  cosmicCut2 = new TH1D("cosmicCut2", " cosmic late large hit chan 12 ", 200, 0, 2000);
  threshValueHist = new TH2D("threshValueHist", " threshold crossings value channels by time  ", 7500, 0, 7500, 1000, 0, 100000);
  threshHist = new TH1D("threshHist", " threshold crossings trig channels ", 20, 0, 20);
  crossHist = new TH1D("crossHist", "  negative crossings non trigger channels", 100, 0, 100);
  sumDir = fout->mkdir("sumDir");
  sumDir->cd();
  double limit;
  double plimit;
  for (unsigned i = 0; i < rawBr.size(); ++i)
  {
    unsigned ichan = i;
    limit = 100000;
    plimit = 2000;

    bool trigger = ichan == 9 || ichan == 10 || ichan == 11;
    if (trigger)
    {
      limit = 200000;
      plimit = 10000;
    }
    if (ichan == 12)
    {
      limit = 5.E3;
      plimit = 1.E3;
    }

    // else if(ichan==12) limit = 10000;
    hQPEShape.push_back(new TH1D(Form("QPEShapeChan%i", ichan), Form("QPEShapeChan%i", ichan), 400, -100, 300));
    hQPEShape[hQPEShape.size() - 1]->SetMarkerStyle(20);

    hSingletShape.push_back(new TH1D(Form("SingletShapeChan%i", ichan), Form("SingletShapeChan%i", ichan), 400, -100, 300));
    hSingletShape[hSingletShape.size() - 1]->SetMarkerStyle(20);

    hQSum.push_back(new TH1D(Form("QSumChan%i", ichan), Form("QSumChan%i", ichan), 1000, 0, limit));
    hQPeak.push_back(new TH1D(Form("QPeakChan%i", ichan), Form("QPeakChan%i", ichan), 1000, 0, plimit));
    sumWave.push_back(new TH1D(Form("sumWave%i", ichan), Form("sumWave%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));
    sumWaveB.push_back(new TH1D(Form("sumWaveBad%i", ichan), Form("sumWaveBad%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));
    sumHitWave.push_back(new TH1D(Form("sumHitWave%i", ichan), Form("sumHitWave%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));
    sumPeakWave.push_back(new TH1D(Form("sumPeakWave%i", ichan), Form("sumPeakWave%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));
  }
  fout->cd();

  // fout->ls();
  cout << " make hitFinder dets = " << CHANNELS << "  size " << rawBr[0]->rdigi.size() << endl;
  vector<int> chanList;
  for (int ichan = 0; ichan < CHANNELS; ++ichan)
    chanList.push_back(ichan);

  finder = NULL;
  finder = new hitFinder(fout, tbrun, tag, rawBr[0]->rdigi.size(), chanList, channelSigmaValue);

  int npass = 0;
  int nfail = 0;
  printf("... total entries  %llu looping over %llu firstEntry %llu \n ", rawTree->GetEntries(), nentries, firstEntry);
  for (Long64_t entry = firstEntry; entry < nentries; ++entry)
  {
    tbrun->clear();
    if (entry / 1000 * 1000 == entry)
    {
      printf("... entry %llu pass %u fail %u \n", entry, npass, nfail);
      //hEventPass->Print("all"); bacondaq seg violated here
      if (npass > 0)
      {
        printf(" \t hits by channel  \n");
        for (int ibin = 0; ibin < histHitCount->GetNbinsX() - 1; ++ibin)
          printf(" chan %i count %i frac %f ; zero %i \n", ibin,
                 int(histHitCount->GetBinContent(ibin + 1)), double(histHitCount->GetBinContent(ibin + 1)) / double(npass), int(hNoPeak->GetBinContent(ibin + 1)));
        printf("  \n");
      }
    }
    rawTree->GetEntry(entry);

    int passBit = anaEvent(entry);
    if (passBit == 0)
    {
      ++npass;
    }
    else
    {
      ++nfail;
      printf(" event %llu fails with pass bit  %x pass %i fail %i \n", entry, passBit, npass, nfail);
    }
    // hEventPass->Fill(-1);
    //  use total entries for all and bin 0 for passing
    hEventPass->SetBinContent(passBit, hEventPass->GetBinContent(passBit) + 1);
    printf(" line1353 event %lld passbit %x bin 0 %f \n",entry, passBit,hEventPass->GetBinContent(2));
    //  if(eventPass!=0)
    //    printf("event fails with eventPass = %x npass %i nfail %i \n", eventPass,npass,nfail);
    //tbrun->print();
    //printf(" %s %lu \n", tbrun->detList[13]->GetName(), tbrun->detList[13]->hits.size());
    /*
    for (unsigned it = 0; it < tbrun->detList[13]->hits.size();++it)
     tbrun->detList[13]->hits[it].print();

    tbrun->detList[13]->clear();
  */
    tbrun->fill();
  }
  printf(" \n \n At END OF FILE total pass  = %i fail %i  \n", npass, nfail);

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

  printf(" ******* FINISHED***** \n \t hits by channels %i   \n", histHitCount->GetNbinsX());
  for (int ibin = 0; ibin < histHitCount->GetNbinsX() - 1; ++ibin)
    printf(" chan %i count %i frac %f ; zero %i \n", ibin,
           int(histHitCount->GetBinContent(ibin + 1)), double(histHitCount->GetBinContent(ibin + 1)) / double(npass), int(hNoPeak->GetBinContent(ibin + 1)));
  printf("  \n");

  printf(" \n \t ***** sums by channel with entries %.0f \n", hTotSum[0]->GetEntries());
  for (int idet = 0; idet < hTotSum.size() ; ++idet){
    printf(" \t chan %i means: tot %.2f pre %.2f trig %.2f late %.2f\n",idet,
           hTotSum[idet]->GetMean(),
           hPreSum[idet]->GetMean(),
           hTrigSum[idet]->GetMean(),
           hLateSum[idet]->GetMean());
  }

  // printf(" FINISHED npass %u nfail %u output file  %s \n", npass, nfail, fout->GetName());
  printf(" FINISHED %i ( %i ) pass %i (%i) fail %i ( frac %0.3f ) output file %s  \n",
         npass + nfail,
         int(hEventPass->GetEntries()),
         npass, int(hEventPass->GetBinContent(0)),
         nfail,
         double(nfail) / double(npass + nfail),
         // int(hEventPass->GetBinContent(8 + 1)),
         fout->GetName());

  // hEventPass->Print("all");
  printf("pass fractions total = %.0f \n", hEventPass->GetBinContent(0));
  for (int ibin = 1; ibin < hEventPass->GetNbinsX(); ++ibin)
    printf(" bin %i fail %.f frac %.3f \n", ibin, hEventPass->GetBinContent(ibin), hEventPass->GetBinContent(ibin) / hEventPass->GetBinContent(0));

  hEventPass->Print("all");

  fout->Write();
  fout->Close();
  return nentries;
}

anaCRun::anaCRun(TString theTag)
{
  tag = theTag;
  cout << " anaCRun::anaCRun instance of anaCRun with tag= " << tag << " CHANNELS = " << CHANNELS - 1 << endl;

  rawBr.clear();

  for (int ichan = 0; ichan < CHANNELS; ++ichan)
  {
    TBRawEvent *rawEv = new TBRawEvent(ichan);
    rawEv->rdigi.resize(7500);
    rawEv->SetName(Form("rawChan%i", ichan));
    rawBr.push_back(rawEv);
  }
}