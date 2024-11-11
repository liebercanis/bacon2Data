// ***This is GAMMA version Sept 25 2024 * **
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
// root/chan
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

  // pass bit failures
  enum FAILURECODES
  {
    PASS = 0,
    GAUSFAIL = 0x1,
    EARLYCUT = 0x2,
    FIRSTTIME = 0x4,
    COSMIC = 0x8,
    GAMMA = 0x10,
    TOTALCODES = 2 * GAMMA
  };

  std::vector<TString> codeNames;

  int badEvent = 5671;

  bool doNotOverWrite = true;
  bool theFirstFile = true;
  int badEventDirMax = 1000;
  bool reportFailures = false;
  TBRun *tbrun;
  TFile *fout;
  TFile *fin;
  TTree *rawTree;
  TNtuple *ntHit;
  unsigned orderFraction = 10;
  // vectors for gains
  std::vector<double> sipmGain;
  std::vector<double> sipmGainError;
  //
  std::map<int, int> chanMap;
  vector<int> nSpeSum;
  vector<TBRawEvent *> rawBr;
  TBEventData *eventData;
  TBEventData *rawEventData;
  TNtuple *ntThresholdAll;
  TNtuple *ntThresholdAdc;
  TNtuple *ntThreshold;
  TNtuple *ntChan;
  TNtuple *ntChanSum;
  TNtuple *ntTrigTime;
  TNtuple *ntSetTrigTime;
  TNtuple *ntSpeYield;
  TNtuple *ntAdc;
  vector<TH1D *> baseHist;
  vector<TH1D *> noiseHist;
  vector<TH1D *> skewHist;
  vector<TH1D *> sumWave;
  vector<TH1D *> sumHitWave;
  vector<TH1D *> sumPeakWave;
  vector<TH1D *> valHist;
  vector<TH1D *> sumWaveB;
  vector<TH1D *> valHistB;

  vector<TH1D *> hMult;
  vector<TH1D *> hQSum;
  vector<TH1D *> hQPeak;
  vector<TH1D *> hQSpe;
  TH1D *hEvBaseWave;
  vector<TH1D *> hEvGaus;
  vector<TH1D *> hEvRawWave;
  vector<TH1D *> hChannelGaus;
  std::vector<std::vector<TH1D *>> hSPEShape; // 4 shapes per channel
  std::vector<TH1D *> hSPEShapeLate;
  // for sums needed for gains
  std::vector<TH1D *> hTotSum;
  std::vector<TH1D *> hPreSum;
  std::vector<TH1D *> hTrigSum;
  std::vector<TH1D *> hLateSum;
  std::vector<TH1D *> hWave;

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
  TH1D *hGammaCut;
  TH1D *hCosmicCut;
  TH1D *hCosmicMult;
  // TH1D *histQPE;
  TH1D *histQPrompt;
  TH1D *hTriggerTimeDiff;
  TH1D *hTriggerTimeAllVal;
  TH1D *hTriggerTimeAllValPmt;
  TH1D *hTriggerHitTimeAll;
  TH1D *hTriggerTime;
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
  vector<double> chanThreshold;
  // sums
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
  bool outFileCheck(TString outFileName);
  unsigned getListOfFiles(TString dir);
  void printGains();
  bool readGains(TString fileName);
  void getSummedHists();
  unsigned getBranches();
  int anaEvent(Long64_t entry);                   // return passBit
  void differentiate(double step);                //
  void derivativeCount(TDet *idet, Double_t rms); // not used
  void negativeCrossingCount(int ichan);
  void thresholdCrossingCount(double thresh);
  std::vector<double> sumDigi();
  unsigned getTriggerTime(int ichan, double &adc);
  void getTriggerTimeStats(unsigned *timeArray, double &ave, double &sigma, unsigned &ichan, double &dmax);
  unsigned fixedTriggerTime(int ichan, double &adc);
  void doTimeShiftAndNorm();
  void getMaxRawAdc(int ichan, double base, double &maxAdc, int &maxSample);
  void setTBRun(TBRun *theTBRun)
  {
    tbrun = theTBRun;
  }
  std::vector<std::vector<double>> fixedDigi; // all the fixed waveforms
  std::vector<unsigned> trigTimes;
  std::vector<unsigned> sTrigTimes; // after correction
  std::vector<double> adcBin;
  std::vector<double> speCount;
  TDirectory *threshDir;
  TDirectory *earlyPeakDir;
  TDirectory *rawSumDir;
  TDirectory *badDir;
  TDirectory *badTrigDir;
  TDirectory *baseDir;
  TDirectory *sumDir;
  TDirectory *anaDir;
  TDirectory *badEventDir;
  TDirectory *pmtDir;
  Long64_t nentries;
  double QPEPeak;
  //
  int MaxSPEShape = 4;
  unsigned trigStart = 600;
  int nominalTrigger = 753; // was 729;
  double nominalGain = 160.0;
  // 227.4; // average
  //  double nominalGain = 160.0; // average
  unsigned firstTime;       // corrected trigger time for event
  unsigned timeOffset = 13; // changed from 17 may 13, 2024
  double passValEarlyCut = 100.0;
  double passValEarlyPmtCut = 225.0;
  ULong_t triggerEnd = 800; // 740;
  ULong_t lateTimeStart = 900;
  ULong_t triggerStart = 730; // 740;
  ULong_t timeVeryLateCut = 3500;
  double prePeakCut = 0.5;
  double latePeakCut = 3.5; // march 18 2024 2.5;
  double diffStepSipm = 3.; // 6 ns steps for SIPM
  double diffStepPmt = 1.;  // back to one on Oct 15 2024
  double cosmicCut = 3.E3;  // set Nov 3 2024
  double gammaCut = 1.E5;   // set Nov 3 2024
  double qpeakCosmicCut = 200.;
  double hitThresholdPmt = 30.;
};

void anaCRun::getMaxRawAdc(int ichan, double base, double &maxAdc, int &maxSample)
{
  maxAdc = -1.E0;
  for (unsigned j = 0; j < rawBr[ichan]->rdigi.size(); ++j)
  {
    double adc = double(rawBr[ichan]->rdigi[j]) - base;
    if (adc > maxAdc)
    {
      maxAdc = adc;
      maxSample = int(j);
    }
  }
}

/** get trigger time
 * do not make the cut here but after this call
 * **/
unsigned anaCRun::getTriggerTime(int ic, double &adc)
{
  TDet *idet = tbrun->getDet(ic);
  unsigned utime = 0;
  unsigned off = 0;
  if (ic < 9)
    off = timeOffset;
  adc = 0;
  for (unsigned j = 0; j < triggerEnd; ++j)
  {
    double val = double(rawBr[ic]->rdigi[j]) - idet->base;
    // here I want the pure ADC count
    // val *= nominalGain / sipmGain[ic];
    if (val > adc)
    {
      adc = val;
      utime = j + off;
    }
  }
  return utime;
}

void anaCRun::getTriggerTimeStats(unsigned *timeArray, double &ave, double &sigma, unsigned &ichan, double &dmax)
{
  unsigned nave = 0;
  ave = 0;
  sigma = 0;
  // calculate ave
  for (unsigned ic = 0; ic < 3; ++ic)
  {
    if (timeArray[ic] < triggerEnd && timeArray[ic] > triggerStart)
    {
      ave += double(timeArray[ic]);
      ++nave;
    }
  }
  // will cast as unsigned
  if (nave > 0 && ave > 0)
    ave /= double(nave);
  else
    ave = nominalTrigger;

  if (nave == 0)
    return;
  // calculate sigma
  for (unsigned ic = 0; ic < 3; ++ic)
  {
    sigma += pow(double(timeArray[ic] - ave), 2.);
  }
  sigma = sqrt(sigma) / double(nave);
  // find channel that is outlier

  ichan = 0;
  dmax = 0;
  for (unsigned ic = 0; ic < 3; ++ic)
  {
    double delta = abs(double(timeArray[ic]) - ave);
    if (delta > dmax)
    {
      dmax = delta;
      ichan = ic;
    }
  }
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
    int timeShift = nominalTrigger - firstTime;
    // printf(" first %u  shift %i off %u chan %u \n", firstTime, timeShift, timeOffset, ib);
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

void anaCRun::printGains()
{
  printf("got %lu gains \n", sipmGain.size());
  for (unsigned long j = 0; j < sipmGain.size(); ++j)
  {
    printf(" %lu  gain %.4f error %.4f   \n", j, sipmGain[j], sipmGainError[j]);
  }
}

bool anaCRun::readGains(TString fileName)
{
  /* define nominal */
  sipmGain.clear();
  sipmGainError.clear();
  sipmGain.resize(NONSUMCHANNELS);
  sipmGainError.resize(NONSUMCHANNELS);
  /*    preliminaru gains
        no trig 158
        trig 760
        PMT 380
  */
  for (unsigned long j = 0; j < sipmGain.size(); ++j)
  {
    if (j < 9)
    {
      sipmGain[j] = 158.;
      sipmGainError[j] = sqrt(158.);
    }
    else if (j < 12)
    {
      sipmGain[j] = 755.;
      sipmGainError[j] = sqrt(755.);
    }
    else
    {
      sipmGain[j] = nominalGain;
      sipmGainError[j] = sqrt(nominalGain);
    }
  }
  /* look for gain file */
  bool exists = false;
  FILE *aFile;
  aFile = fopen(fileName.Data(), "r");
  if (aFile)
  {
    fclose(aFile);
    exists = false;
  }

  if (!exists)
  {
    printf(" couldnt open template file %s\n", fileName.Data());
    return false;
  }

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

  for (int i = 0; i < gGain->GetN(); ++i)
  {
    int index = int(gGain->GetPointX(i));
    sipmGain[index] = gGain->GetPointY(i);
    sipmGainError[index] = gGain->GetErrorY(i);
  }
  return true;
}

void anaCRun::clear()
{
  nSpeSum.clear();
  hQSum.clear();
  hSPEShape.clear();
  hSPEShapeLate.clear();
  hQPeak.clear();
  hQSpe.clear();
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
  hEvRawWave.clear();
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
  chanThreshold.resize(CHANNELS);
  channelSigmaValue.resize(CHANNELS);
  // based on step 3 for SIPMs and 1 for PMT Oct 16,2024
  for (unsigned long j = 0; j < channelSigmaValue.size(); ++j)
  {
    chanThreshold[j] = 2. * 36.4;
  }
  // special threshold values
  chanThreshold[9] = 2. * 21.6;
  chanThreshold[10] = 2. * 21.6;
  chanThreshold[11] = 2. * 21.6;
  chanThreshold[12] = 2. * 2.1;  // based on histogram sigma, had been 5.3;
  chanThreshold[13] = 2. * 36.4; // should be same as trigger sipm
  nSpeSum.resize(CHANNELS);
}

bool anaCRun::outFileCheck(TString outFileName)
{
  // does file exist?
  printf(" check for existing output file %s\n", outFileName.Data());
  bool exists = false;
  FILE *aFile;
  aFile = fopen(outFileName.Data(), "r");
  if (aFile)
  {
    fclose(aFile);
    exists = true;
  }
  if (!exists)
    return false;
  // check that file was closed properly
  TFile *fcheck = new TFile(outFileName, "readonly");
  TTree *tree = nullptr;
  fcheck->GetObject("RunTree", tree);
  if (tree == nullptr)
  {
    printf(" file not closed properly %s so run again \n", outFileName.Data());
    return false;
  }
  printf(" outFileCheck of %s returns true  \n", outFileName.Data());
  return true;
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
  fin->ls();
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
  rawSumDir->cd();
  TIter next(fin->GetListOfKeys());
  printf(" getSummedHists ........ list of fin \n");
  fin->GetListOfKeys()->ls();
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
  //  clear
  TTree *tree = NULL;
  fout->GetObject("RunTree", tree);
  if (!tree)
  {
    printf("line 531 ERROR!! anaEvent no tree event %lld \n", entry);
    fout->ls();
  }
  tbrun->clear(); // clear detList
  speCount.clear();
  speCount.resize(CHANNELS);
  std::fill(speCount.begin(), speCount.end(), 0);
  int passBit = PASS;
  // previously 40 but channel 9 was missing peaks
  // double hitThreshold = nominalGain - 3. * 16.; // this is 3 sigma of SPE peak June 3 2014
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
    TDet *idet = tbrun->getDet(ichan);
    // define trigger sipms
    bool trig = ichan == 9 || ichan == 10 || ichan == 11;
    // deal with trigger channel sign by overwriting rdigi
    // also invert pulse on PMT
    if (trig || ib == 12)
    {
      for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
      {
        rawBr[ib]->rdigi[j] = -1. * (rawBr[ib]->rdigi[j] - pow(2, 14)); // base is > digi value!
      }
    }

    int nbins = rawBr[ib]->rdigi.size();
    // cout << "@line611 " << ib << " nbins " << nbins << " max hist " << hEvGaus.size() << " rawBr.size() " << rawBr.size() << endl;

    // sanity check
    if (rawBr[ib]->rdigi.size() != WAVELENGTH)
      continue;
    // simple baseline
    /* this base is a smaller value and in the case of almost all noise, is biased
    std::vector<unsigned short> orderDigi = rawBr[ib]->rdigi;
    std::sort(orderDigi.begin(), orderDigi.end(), std::less<int>());
    unsigned baseLength = orderDigi.size() / orderFraction;
    double base = 0;
    for (unsigned j = 0; j < baseLength; ++j)
    {
      base += orderDigi[j];
    }
    base /= double(baseLength);
    */

    // baseline from pre trigger data
    double base = 0;
    for (unsigned j = 0; j < trigStart; ++j)
    {
      base += rawBr[ib]->rdigi[j];
    }
    base /= double(trigStart);

    // baseline correction from fitted Gaussian
    hEvGaus[ib]->Reset("ICES");
    hEvRawWave[ib]->Reset("ICES");
    for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
    {
      double val = double(rawBr[ib]->rdigi[j]) - base; // base is > digi value!
      hEvGaus[ib]->Fill(val);
      hEvRawWave[ib]->SetBinContent(j + 1, val);
    }
    // get the distribution mode
    double mode = hEvGaus[ib]->GetBinLowEdge(hEvGaus[ib]->GetMaximumBin()) + 0.5 * hEvGaus[ib]->GetBinWidth(hEvGaus[ib]->GetMaximumBin());

    /* dont do this memory leak first clone*/
    hEvGaus[ib]->GetListOfFunctions()->Clear();
    // TH1D* hEvClone = (TH1D*) hEvGaus[ib]->Clone("EvClone");
    TFitResultPtr fitptr = hEvGaus[ib]->Fit("gaus", "LQ0", "", -100, 100);
    int fitStatus = fitptr;
    TF1 *gfit = (TF1 *)hEvGaus[ib]->GetListOfFunctions()->FindObject("gaus");
    double ave = hEvGaus[ib]->GetMean();
    double sigma = hEvGaus[ib]->GetRMS();
    double skew = 0;
    double fitMean = 0;
    if (!isnan(hEvGaus[ib]->GetSkewness()))
      skew = hEvGaus[ib]->GetSkewness();
    if (gfit != nullptr && fitStatus == 0)
    {
      ave = gfit->GetParameter(1);
      fitMean = ave; // fit mean
      sigma = gfit->GetParameter(2);
    }
    // fit status = migradStatus + 10*minosStatus + 100*hesseStatus + 1000*improveStatus \n", fullResult);
    else
    {
      if (reportFailures)
        printf("@line750 failed Baseline Cut event %llu cut %lu time %u \n", entry, triggerEnd, firstTime);
      if (badEventDir->GetList()->GetEntries() < badEventDirMax)
      {
        badEventDir->cd();
        TH1D *EvRawWave = (TH1D *)hEvRawWave[ib]->Clone(Form("EvRawBaselineEvent%lld-Ch%i", entry, ib));
        EvRawWave->SetTitle(Form("EvRawBaselineEvent%lld-Ch%i", entry, ib));
      }
      passBit |= GAUSFAIL;
    }

    // printf("@line652 baseline %lld chan %u base %f ave %f  \n", entry, ib, base, ave);

    baseDir->cd();
    if (baseDir->GetList()->GetEntries() < 100 || fitStatus != 0)
    {
      TH1D *EvGaussEvent = (TH1D *)hEvGaus[ib]->Clone(Form("EvGaussEvent%lld-Ch%i", entry, ib));
    }
    fout->cd();

    // fitptr->Print();
    noiseHist[ib]->Fill(sigma);
    skewHist[ib]->Fill(skew);
    double sign = TMath::Sign(1., skew);

    if (idet == NULL)
    {
      printf("@line711!!!!!NULL idet br %u ichan %i\n", ib, ichan);
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
    /* add maxAdc */
    double maxAdc;
    int maxSample;
    getMaxRawAdc(ib, idet->base, maxAdc, maxSample);
    idet->maxAdc = maxAdc;
    idet->maxSample = maxSample;
    /* I have added this to the TDet as  maxSample maxAdc */
    if (ntAdc->GetEntries() < 1E9)
    {
      for (unsigned j = 0; j < rawBr[ichan]->rdigi.size(); ++j)
      {
        double adc = double(rawBr[ichan]->rdigi[j]) - idet->base;
        if (adc > 2. * idet->sigma)
          ntAdc->Fill(double(entry), double(ib), double(j), adc);
      }
    }

  } // channel loop
  /* find trigger time from trigger sipms */
  trigTimes.resize(NONSUMCHANNELS);
  sTrigTimes.resize(NONSUMCHANNELS);
  adcBin.resize(NONSUMCHANNELS);
  for (unsigned ib = 0; ib < NONSUMCHANNELS; ++ib)
  {
    double val = 0;
    // get time for maximim val before triggerEnd
    unsigned time = getTriggerTime(ib, val); // include timeOffset in routine
    ntSetTrigTime->Fill(double(entry), double(ib), double(time), double(val));
    trigTimes[ib] = time;
    adcBin[ib] = val;

    /* set passBit 1 val cut is different for PMT*/
    if (time < unsigned(triggerStart) && val > passValEarlyCut && ib < 12)
    {
      if (reportFailures)
        printf("@line757 failed passBit 1 triggerStart event %llu chan %i time %u val %f \n", entry, ib, time, val);
      passBit |= EARLYCUT;
      if (badEventDir->GetList()->GetEntries() < badEventDirMax)
      {
        badEventDir->cd();
        // printf("@line862 failed RawEarlyEvent event %llu chan%u cut %lu time %u val %f \n", entry, ib, triggerEnd, time, val);
        TH1D *EvRawWave = (TH1D *)hEvRawWave[ib]->Clone(Form("EvRawEarlyEvent%lld-Ch%i", entry, ib));
        EvRawWave->SetTitle(Form("EvRawEarlyEvent%lld-Ch%i", entry, ib));
      }
    }
    if (time < unsigned(triggerStart) && val > passValEarlyPmtCut && ib == 12)
    {
      if (reportFailures)
        printf("@line757 failed triggerStart event %llu chan %i time %u val %f \n", entry, ib, time, val);
      passBit |= EARLYCUT;
      if (badEventDir->GetList()->GetEntries() < badEventDirMax)
      {
        badEventDir->cd();
        // printf("@line862 failed RawEarlyEvent event %llu chan %u cut %lu time %u val %f \n", entry, ib, triggerEnd, time, val);
        TH1D *EvRawWave = (TH1D *)hEvRawWave[ib]->Clone(Form("EvRawEarlyEvent%lld-Ch%i", entry, ib));
        EvRawWave->SetTitle(Form("EvRawEarlyEvent%lld-Ch%i", entry, ib));
      }
    }

    if (ib != 12)
      hTriggerTimeAllVal->Fill(double(val));
    else
      hTriggerTimeAllValPmt->Fill(double(val));
    // fill an ntuple to monitor this cut
  }

  /* ave trigger sipm times before shift */
  double trigTimeAve = 0;
  double trigTimeSigma = 0;
  unsigned chanBad;
  double dmax;
  getTriggerTimeStats(&trigTimes[9], trigTimeAve, trigTimeSigma, chanBad, dmax);
  hTriggerTimeDiff->Fill(dmax);
  firstTime = unsigned(trigTimeAve); // av::e of trigger sipm times

  /* ave non trig  sipm times before shift */
  double nonTimeAve = 0;
  double nonTimeSigma = 0;
  unsigned chanBad2;
  double dmax2;
  getTriggerTimeStats(&trigTimes[6], nonTimeAve, nonTimeSigma, chanBad2, dmax2);

  hTriggerTime->Fill(double(firstTime));
  if (firstTime > triggerEnd && firstTime < triggerStart)
  {
    if (reportFailures)
      printf("@line893 failed triggerEnd event %llu cut %lu time %u \n", entry, triggerEnd, firstTime);
    passBit |= FIRSTTIME;
    if (badEventDir->GetList()->GetEntries() < badEventDirMax)
    {
      badEventDir->cd();
      TH1D *EvRawWave = (TH1D *)hEvRawWave[chanBad]->Clone(Form("EvRawFirstTimeEvent%lld-Ch%i", entry, chanBad));
      EvRawWave->SetTitle(Form("EvRawFirstTimeEvent%lld-Ch%i", entry, chanBad));
    }
  }
  /********************************************************
   * now that we have the firstTime
        align to nominalTrigger
        defined as timeShift>0 shift right
        normalize to nominal gain
   ********************************************************/
  // printf("doTimeSiftAndNorm %lld \n",entry);
  doTimeShiftAndNorm();

  /*  fill ntuple for threshold setting loop over channels */
  for (unsigned long ib = 0; ib < NONSUMCHANNELS; ++ib)
  {
    digi.clear();
    digi = fixedDigi[ib];
    TDet *idet = tbrun->getDet(ib);

    double peakMax = 0;
    // do digi sums on fixedDigi
    for (unsigned j = 0; j < digi.size(); ++j)
    {
      idet->totSum += digi[j];
      if (j < triggerStart)
        idet->preSum += digi[j];

      if (j > triggerStart && j < triggerEnd)
      {
        idet->trigSum += digi[j];
        if (digi[j] > peakMax)
          peakMax = digi[j];
      }

      if (j > lateTimeStart)
        idet->lateSum += digi[j];
    }
    // add some other variables
    idet->peakMax = peakMax;

    ntChan->Fill(float(rawBr[ib]->trigger), float(ib), float(idet->ave), float(idet->sigma), float(idet->skew), float(idet->base), float(peakMax), float(idet->totSum), float(idet->lateSum), float(crossings.size()), float(passBit));

    if (ib == 12)
      differentiate(diffStepPmt);
    else
      differentiate(diffStepSipm);
    unsigned long sampleLow = 0;
    double valLow = 1.E-9;
    unsigned long sampleHigh = 0;
    double valHigh = -1.E-9;
    unsigned long maxBin = 0;
    double adcMax = -1.E-9;
    // find high and low
    for (unsigned long idd = 0; idd < ddigi.size(); ++idd)
    {
      // limit size because otherwise this is huge
      if (ddigi[idd] < valLow)
      {
        valLow = ddigi[idd];
        sampleLow = idd;
      }
      if (ddigi[idd] > valHigh)
      {
        valHigh = ddigi[idd];
        sampleHigh = idd;
      }
      if (ntThresholdAll->GetEntries() < 1.0E7)
        ntThresholdAll->Fill(float(entry), float(ib), float(idd), float(ddigi[idd]));

      //  find max anywhere
      if (digi[idd] > adcMax)
      {
        maxBin = idd;
        adcMax = digi[idd];
      }
    } // ddigi loop
    // if (!((sampleLow - sampleHigh) > 0 && (sampleLow - sampleHigh) < 50))
    //   continue;

    ntThresholdAdc->Fill(entry, ib, sampleLow, sampleHigh, maxBin, adcMax);

    // plot to see a few of these
    if (abs(valLow) > chanThreshold[ib] && valHigh > chanThreshold[ib] && threshDir->GetList()->GetEntries() < 100)
    {
      threshDir->cd();
      for (int ibin = 0; ibin < hWave[ib]->GetNbinsX(); ++ibin)
        hWave[ib]->SetBinContent(ibin, digi[ibin]);
      TString histName;
      TString detName = tbrun->detList[ib]->GetName();
      histName.Form("Wave%lli%sMax%.0fDDigi%.0f", entry, detName.Data(), adcMax, valHigh);
      TH1D *hEventWave = (TH1D *)hWave[ib]->Clone(histName);
      hEventWave->SetTitle(histName);
    }

    // printf("@line808 %lld ichan %lu low %lu high %lu maxBin %lu adcMax %f  \n", entry, ib, sampleLow, sampleHigh, maxBin, adcMax);
    if (abs(valLow) > chanThreshold[ib] && valHigh > chanThreshold[ib])
      ntThreshold->Fill(entry, ib, sampleLow, ddigi[sampleLow], sampleHigh, ddigi[sampleHigh], maxBin, adcMax);
  }

  /* **** */
  /* make ntuple of before and after shift */
  for (unsigned ic = 0; ic < NONSUMCHANNELS; ++ic)
  {
    double val;
    sTrigTimes[ic] = fixedTriggerTime(ic, val);
    ntTrigTime->Fill(double(entry), double(ic), double(firstTime), trigTimes[ic], adcBin[ic], sTrigTimes[ic], val);
  }

  /*******
        start of second channel loop doing pulse finding
   ********/

  // also fill chan 13
  TDet *tdet13 = tbrun->getDet(NONSUMCHANNELS); // get channel 13 det
  tdet13->clear();
  hEvRawWave[NONSUMCHANNELS]->Reset("ICES");
  for (unsigned ib = 0; ib < NONSUMCHANNELS; ++ib)
  {
    unsigned ichan = ib;
    TDet *tdet = tbrun->getDet(ib);

    // fill summed wave
    for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
    {
      double val = double(rawBr[ib]->rdigi[j]) - tdet->base; // base is > digi value!
      hEvRawWave[NONSUMCHANNELS]->SetBinContent(j + 1, hEvRawWave[NONSUMCHANNELS]->GetBinContent(j + 1) + val);
    }

    // fill tdet13 sums
    tdet13->event = entry;
    tdet13->totSum += tdet->totSum;
    tdet13->preSum += tdet->preSum;
    tdet13->trigSum += tdet->trigSum;
    tdet13->lateSum += tdet->lateSum;

    /************************************/
    // make hit on channel ib
    tdet->hits.clear();
    bool trig = ichan == 9 || ichan == 10 || ichan == 11;
    int nbins = rawBr[ib]->rdigi.size();
    digi.clear();
    digi = fixedDigi[ib];

    evCount->Fill(ib);                       // chan 0 from GetBinContent(0)
    double hitThreshold = 0.2 * nominalGain; // 500.0;
    if (ib == 12)
      hitThreshold = hitThresholdPmt; // this is 5*(6 sigma noise)
    double theStep = diffStepSipm;
    if (ib == 12)
    {
      theStep = diffStepPmt;
    }
    finder->event(ichan, entry, digi, chanThreshold[ib], hitThreshold, theStep); // DEG suggests 10
    // add hits to channel 13
    for (unsigned ihit = 0; ihit < tdet->hits.size(); ++ihit)
    {
      tdet13->hits.push_back(tdet->hits[ihit]);
    }

    /* this was one cosmic
    if (badEvent == entry)
      finder->plot1Wave(badEventDir, tdet->channel, entry);
      */

    // if (entry / 1000 * 1000 == entry)
    //   printf("line974 badEventDir %lld det %i size %d \n ", entry, ib, badEventDir->GetList()->GetEntries());

    // look at PMT
    // TDirectory *badEventDir = (TDirectory *)fout->FindObject("badEventDir");

    /* do not fill here any more
    for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
    {
      double adc = double(rawBr[ib]->rdigi[j]) - tdet->base;
      if (ntAdc->GetEntries() < 1E9)
      {
        ntAdc->Fill(double(entry), double(ib), double(j), adc);
      }
    }
    */
  }

  // for cosmic cut, count large photons
  int nCosmicHits = 0;
  TDet *tdetPmt = tbrun->getDet(12);
  for (unsigned ihit = 0; ihit < tdetPmt->hits.size(); ++ihit)
  {
    if (tdetPmt->hits[ihit].qpeak > qpeakCosmicCut)
      ++nCosmicHits;
  }
  // do cosmic cut based on large pulse counting
  hCosmicMult->Fill(double(nCosmicHits));
  hCosmicCut->Fill(tdetPmt->totSum);
  if (nCosmicHits > 0)
  {
    if (reportFailures)
      printf("@line1077 failed cosmic event %llu cut %E totSum %E nCosmicHits %i \n", entry, cosmicCut, tdetPmt->totSum, nCosmicHits);
    passBit |= COSMIC;
    if (badEventDir->GetList()->GetEntries() < badEventDirMax)
    {
      badEventDir->cd();
      TH1D *EvRawWave = (TH1D *)hEvRawWave[12]->Clone(Form("EvRawCosmicEvent%lldTotSum%.0E-Ch%i", entry, tdetPmt->totSum, 12));
      EvRawWave->SetTitle(Form("EvRawCosmicEvent%lldTotSum%.3E-Ch%i", entry, tbrun->getDet(12)->totSum, 12));
    }
  }

  // do gamma cut
  hGammaCut->Fill(tbrun->getDet(13)->lateSum);
  if (tbrun->getDet(13)->lateSum > gammaCut)
  {
    if (reportFailures)
      printf("@line1090 failed gamma event %llu cut %E lateSum %E \n", entry, gammaCut, tbrun->getDet(13)->lateSum);
    passBit |= GAMMA;
    if (badEventDir->GetList()->GetEntries() < badEventDirMax)
    {
      badEventDir->cd();
      TH1D *EvRawWave = (TH1D *)hEvRawWave[NONSUMCHANNELS]->Clone(Form("EvRawGammaEvent%lldVal%.0E-Ch%i", entry, tbrun->getDet(13)->lateSum, 13));
      EvRawWave->SetTitle(Form("EvRawGammaEvent%lldVal%.3E-Ch%i", entry, tbrun->getDet(13)->lateSum, 13));
    }
  }

  // look at good PMT events
  if (passBit == 0 && tdetPmt->peakMax > hitThresholdPmt && pmtDir->GetList()->GetEntries() < 1000)
  {
    pmtDir->cd();
    printf("@line1120 print event %llu peakMax %E \n", entry, tdetPmt->peakMax);
    TH1D *EvRawWave = (TH1D *)hEvRawWave[12]->Clone(Form("EvRawPMTEvent%lldVal%.0E-Ch%i", entry, tdetPmt->peakMax, 12));
    EvRawWave->SetTitle(Form("EvRawPMTEvent%lldVal%.3E-Ch%i", entry, tdetPmt->peakMax, 12));
    if (tdetPmt->hits.size() > 0)
      finder->plotEvent(pmtDir, tdetPmt->channel, entry);
  }

  // printf("line975 chan 13 has %lu hits \n", tdet13->hits.size());

  // event cuts on summed line
  int nPreHits = 0;
  int nLateHits = 0;
  for (unsigned ihit = 0; ihit < tdet13->hits.size(); ++ihit)
  {
    TDetHit hiti = tdet13->hits[ihit];
    ULong_t hitStartTime = ULong_t(tdet13->hits[ihit].startTime);
    hCountLateTimeQpeak->Fill(hitStartTime, tdet13->hits[ihit].qpeak / nominalGain);
    if (hitStartTime < triggerStart)
      hPreQpeak->Fill(tdet13->hits[ihit].qpeak / nominalGain);
    if (hitStartTime < triggerStart && tdet13->hits[ihit].qpeak / nominalGain > prePeakCut)
    {
      ++nPreHits;
      // printf("event preHits %llu cut %lu hitStartTime %lu  qpeak %.2f nPreHits %i \n", entry, triggerStart, hitStartTime, hiti.qpeak, nPreHits);
    }
    if (hitStartTime > triggerEnd)
      hLateQpeak->Fill(tdet13->hits[ihit].qpeak / nominalGain);
    if (hitStartTime > triggerEnd && tdet13->hits[ihit].qpeak / nominalGain > latePeakCut)
    {
      ++nLateHits;
      // printf("event lateHits %llu cut %lu hitStartTime %lu  qpeak %.2f nLateHits %i \n", entry, firstTimeCut, hitStartTime, hiti.qpeak, nLateHits);
      //  hCountLateTime->Fill(tdet13->hits[ihit].startTime);
    }
    // if (hitStartTime > 600 && hitStartTime < 800 && hitStartTime < firstTime)
    //   firstTime = hitStartTime;
  }
  hCountPre->Fill(nPreHits);
  hCountLate->Fill(nLateHits);

  /*
  if (nPreHits > 0)
  {
    printf("@line901 failed nPreHits event %llu nPre %i \n", entry, nPreHits);
    passBit |= 0x4;
  }
  if (nLateHits > 0)
  {
    printf("@line904 failed nLate event %llu nLate %i \n", entry, nLateHits);
    passBit |= 0x8;
  }
  */

  evCount->Fill(-1); // underflow bin

  // fill histograms for good and bad events
  for (unsigned ib = 0; ib < NONSUMCHANNELS; ++ib)
  {
    /* take care here for summed ib=CHANNELS-2 and set appropriate hitThreshold */
    digi.clear();
    digi = fixedDigi[ib];
    for (unsigned j = 0; j < digi.size(); ++j)
    {
      if (passBit != 0)
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
  }

  // if (passBit != 0) return passBit;
  // printf("line818  event %lld passbit %i \n",entry,passBit);
  if (passBit != 0)
  {
    // printf("@line913 event %lld passBit %i det %i nhits %u \n",
    //        entry, int(passBit), NONSUMCHANNELS, tbrun->detList[NONSUMCHANNELS]->nhits());
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
    hMult[idet]->Fill(double(tdet->hits.size())); // hit multiplicity
    // add some event plots
    bool trig = tdet->channel == 9 || tdet->channel == 10 || tdet->channel == 11;
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
    if (tdet->hits.size() > 1 && tdet->channel == 12 && sumWaveDir->GetList()->GetEntries() < 5000)
    {
      // printf("line1083xxxxxxx anaCRun::event event %llu chan %i hits %lu der thresh %f hit thresh %f \n", entry, tdet->channel, tdet->hits.size(), chanThreshold[tdet->channel], hitThreshold);
      finder->plotEvent(sumWaveDir, tdet->channel, entry);
    }

    /*
    TDirectory *badEventDir = (TDirectory *)fout->FindObject("badEventDir");
    if (tdet->hits.size() > 1 && tdet->channel == 12 && badEventDir->GetList()->GetEntries() < 5000)
    {
      // printf("xxxxxxx anaCRun::event event %llu chan %i hits %lu der thresh %f hit thresh %f \n", entry, tdet->channel, tdet->hits.size(), derivativeThreshold, hitThreshold);
      finder->plot1Wave(badEventDir, tdet->channel, entry);
    }
    */

    TDirectory *fftDir = (TDirectory *)fout->FindObject("fftDir");
    if (fftDir)
    {
      if (trig && tdet->hits.size() == 0 && fftDir->GetList()->GetEntries() < 2000)
      {
        // printf("!!!!!! anaCRuna::event plot event %llu idet %i chan %i hits %lu \n", entry, idet, tdet->channel, tdet->hits.size());
        finder->plotEvent(fftDir, tdet->channel, entry);
      }
    }

    // loop over hits
    // if (tdet->hits.size() > 0 && idet == 12) // PMT
    //  printf("@line978 event %llu  det %u nhits %lu \n", entry, idet, tdet->hits.size());
    // add peak sums
    if (tdet->hits.size() == 0)
      hNoPeak->SetBinContent(tdet->channel + 1, hNoPeak->GetBinContent(tdet->channel + 1) + 1);
    int firstHitTime = rawBr[NONSUMCHANNELS]->rdigi.size();
    histQSum->SetBinContent(tdet->channel + 1, histQSum->GetBinContent(tdet->channel + 1) + tdet->qpeak);
    histQPrompt->SetBinContent(tdet->channel + 1, histQPrompt->GetBinContent(tdet->channel + 1) + tdet->hitPrompt);
    // printf(" event %lld det %i sum qpeak %f sum qprompt %f\n", entry, idet, tdet->qpeak, tdet->hitPrompt);
    if (tdet->hits.size() == 0)
    {
      hQPeak[idet]->Fill(-1);
      hQSpe[idet]->Fill(0);
    }

    // printf("@line1065 event %llu  det %u nhits %lu \n", entry, idet, tdet->hits.size());
    for (unsigned ihit = 0; ihit < tdet->hits.size(); ++ihit)
    {
      TDetHit thit = tdet->hits[ihit];
      if (thit.qpeak < 1)
        printf("line822 chan %i ihit %i startTime %i  peak %f\n", tdet->channel, ihit, int(thit.startTime), thit.qpeak);
      hQSum[idet]->Fill(thit.qsum / nominalGain);
      hQPeak[idet]->Fill(thit.qpeak / nominalGain);
      unsigned hitTime = unsigned(thit.startTime);
      // do peak sums
      tdet->totPeakSum += thit.qpeak;

      if (hitTime > triggerStart && hitTime < triggerEnd && hitTime < firstHitTime)
        firstHitTime = hitTime;
      //
      if (hitTime < trigStart)
        tdet->prePeakSum += thit.qpeak;
      else if (hitTime < triggerEnd)
      {
        tdet->trigPeakSum += thit.qpeak;
      }
      else if (hitTime > timeVeryLateCut)
        tdet->latePeakSum += thit.qpeak;
      // fill here for gains
      hTotSum[idet]->Fill(thit.qpeak);

      if (hitTime < trigStart)
        hPreSum[idet]->Fill(thit.qpeak);

      if (hitTime > trigStart && hitTime < triggerEnd)
        hTrigSum[idet]->Fill(thit.qpeak);

      if (hitTime > triggerEnd)
        hLateSum[idet]->Fill(thit.qpeak);

      // do threshold for summed waveform
      // if (thit.qsum > hitThreshold)
      // if(int(thit.peakt)-thit.firstBin > 30)
      //  printf("line980 in anaCRun event %lli  det %i  peak %u start %i \n",entry, idet, thit.peakt, thit.firstBin);

      // if(thit.qpeak > 7.5* nominalGain)  printf("line1008 idet %i qpeak %f \n",idet,thit.qpeak/nominalGain );
      sumHitWave[idet]->SetBinContent(thit.firstBin + 1, sumHitWave[idet]->GetBinContent(thit.firstBin + 1) + thit.qsum);
      sumPeakWave[idet]->SetBinContent(thit.firstBin + 1, sumPeakWave[idet]->GetBinContent(thit.firstBin + 1) + thit.qpeak);
      histHitCount->SetBinContent(tdet->channel + 1, histHitCount->GetBinContent(tdet->channel + 1) + 1);

      ntHit->Fill(double(entry), double(passBit), double(idet), thit.startTime, thit.peakt, thit.qpeak);
      // sum of photons in SPE for this channel
      speCount[idet] += thit.qpeak / nominalGain;

      /* fill the SPE histograms */
      // count SPE for this hit
      if (idet < NONSUMCHANNELS)
      { // exclude summed channel
        int nSPE = 0;
        if (thit.qpeak > 0.5 * nominalGain && thit.qpeak < 1.5 * nominalGain)
          nSPE = 1.;
        else if (thit.qpeak > 1.5 * nominalGain && thit.qpeak < 2.5 * nominalGain)
          nSPE = 2.;
        else if (thit.qpeak > 2.5 * nominalGain && thit.qpeak < 3.5 * nominalGain)
          nSPE = 3.;
        else if (thit.qpeak > 3.5 * nominalGain && thit.qpeak < 4.5 * nominalGain)
          nSPE = 4.;
        else if (thit.qpeak > 4.5 * nominalGain && thit.qpeak < 5.5 * nominalGain)
          nSPE = 5.;
        else if (thit.qpeak > 5.5 * nominalGain && thit.qpeak < 6.5 * nominalGain)
          nSPE = 6.;
        else if (thit.qpeak > 6.5 * nominalGain && thit.qpeak < 7.5 * nominalGain)
          nSPE = 7.;
        else
          nSPE = 8;
        // printf("line 995 SPE check %f %i \n",thit.qpeak,nSPE);
        hQSpe[idet]->Fill(nSPE);
        nSpeSum[idet] += nSPE;

        int thePeakBin = 100.; // set to bin 100 in the histogram
        // hSPE shape 1000 bins
        int sumStartBin = max(thit.peakBin - thePeakBin, 1);
        int sumEndBin = min(thit.peakBin - thePeakBin + hSPEShape[0][0]->GetNbinsX(), int(rawBr[NONSUMCHANNELS]->rdigi.size()));
        for (unsigned jbin = sumStartBin; jbin < sumEndBin; ++jbin)
        {
          int fillBin = thePeakBin - thit.peakBin + jbin;
          double val = fixedDigi[idet][jbin];
          // fill 1 SPE from late
          if (thit.startTime > triggerEnd && nSPE == 1)
            hSPEShapeLate[idet]->SetBinContent(fillBin, hSPEShapeLate[idet]->GetBinContent(fillBin) + val);
          // fill the right histogram for 1 SPE take from after trigger
          if (nSPE > 0 && nSPE < MaxSPEShape)
            hSPEShape[nSPE - 1][idet]->SetBinContent(fillBin, hSPEShape[nSPE - 1][idet]->GetBinContent(fillBin) + val);
        }
      }
    } // hit loop

    /* cross check on SPE */
    // double sumAfter = sumPeakWave[idet]->Integral();
    // printf(" line1051 chan %i spe sum  %i sumAfter %.3E  \n",
    //     idet, nSpeSum[idet], sumAfter/nominalGain);

    // hTriggerHitTimeAll->Fill(firstHitTime);
    //  fill sums do not fill for zero
    /* redefine these for gains
    if (tdet->totPeakSum > 0)
      hTotSum[idet]->Fill(tdet->totPeakSum);
    if (tdet->prePeakSum > 0)
      hPreSum[idet]->Fill(tdet->prePeakSum);
    if (tdet->trigPeakSum > 0)
      hTrigSum[idet]->Fill(tdet->trigPeakSum);
    */
    // for each peak not summed
    // if (tdet->latePeakSum>0)
    //  hLateSum[idet]->Fill(tdet->latePeakSum);
    // printf(" anaCRun::event %llu det %i nhits %lu , tot %f pre %f trig %f late %f\n", entry, tdet->channel, tdet->hits.size(),
    //       tdet->totPeakSum, tdet->prePeakSum, tdet->trigPeakSum, tdet->latePeakSum);

  } // det loop

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
// copied from hitFinder spt 18 2924
// revised derivative Jan 27 2023 MG
void anaCRun::differentiate(double diffStep)
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

/* not useed warning cut is based on idet-sigma which is now running sigma*/
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
  //
  string sfilename(theFile.Data());
  string shortName = sfilename.substr(0, sfilename.find_last_of("."));
  cout << " anaCRunFile  with shortName= " << shortName << endl;
  // open outout file
  TString outFileName;
  outFileName.Form("caenData/anaCRun-%s-%llu.root", shortName.c_str(), maxEntries);
  if (doNotOverWrite)
    if (outFileCheck(outFileName))
    {
      printf(" do not recreate %s file \n", outFileName.Data());
      return 0;
    }

  fout = new TFile(outFileName, "recreate");
  cout << " opened output file " << fout->GetName() << endl;

  rawSumDir = fout->mkdir("rawSumDir");
  baseDir = fout->mkdir("baseDir");
  badEventDir = fout->mkdir("badEventDir");
  pmtDir = fout->mkdir("pmtDir");
  badDir = fout->mkdir("badDir");
  badTrigDir = fout->mkdir("badTrigDir");
  threshDir = fout->mkdir("threshDir");
  earlyPeakDir = fout->mkdir("earlyPeakDir");
  anaDir = fout->mkdir("anadir");
  sumDir = fout->mkdir("sumDir");
  TDirectory *finderDir = fout->mkdir("finderDir");
  TDirectory *splitDir = fout->mkdir("splitDir");
  TDirectory *fitSingletDir = fout->mkdir("fitSingletDir");
  TDirectory *sumWaveDir = fout->mkdir("sumWaveDir");
  fout->ls();

  currentBuffer = -1;
  currentBufferCount = 0;
  printf(" anaCRun::anaCRunFile starting anaCRun file %s maxEntries %llu firstEntry %llu \n",
         theFile.Data(), maxEntries, firstEntry);
  if (!openFile(theFile)) // and get branches
  {
    printf("cannot open file!\n");
    return 0;
  }

  // new gain file
  TString gainFileName = TString(getenv("BOBJ")) + TString("/gains-2024-02-15-17-26-save.root");
  cout << "read gains from file " << gainFileName << endl;
  readGains(gainFileName);

  if (theFirstFile)
  {
    theFirstFile = false;
    printf("chanThreshold values \n");
    for (unsigned j = 0; j < chanThreshold.size(); ++j)
      printf("chan %u chanThreshold %.3f \n", j, chanThreshold[j]);
    printGains();
  }

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

  getSummedHists();
  // fout->ls();

  // make output tree
  tbrun = new TBRun(tag);
  fout->Append(tbrun->btree);
  // and event time
  eventData = new TBEventData();
  tbrun->btree->Branch("eventData", &eventData);

  for (unsigned it = 0; it < rawBr.size(); ++it)
  {
    tbrun->addDet(it);
  }

  // fout->append(tbrun->btree);e("ntHit", " hits
  ntAdc = new TNtuple("ntAdc", " ADC ntuple ", "event:chan:sample:digi");
  ntThresholdAll = new TNtuple("ntThresholdAll", "ntThreshold", "event:chan:sample:ddigi");
  ntThresholdAdc = new TNtuple("ntThresholdAdc", "ntThresholdAdc", "event:chan:sampleLow:sampleHigh:maxBin:adcMax");
  ntThreshold = new TNtuple("ntThreshold", "ntThreshold", "event:chan:sampleLow:ddigiLow:sampleHigh:ddigiHigh:maxBin:adcMax");
  ntHit = new TNtuple("ntHit", "hit ntuple", "event:flag:chan:time:peakTime:qpeak");
  ntChan = new TNtuple("ntchan", "channel ntuple", "trig:chan:ave:sigma:skew:base:peakmax:totSum:lateSum:negcrossings:pass");
  ntSpeYield = new TNtuple("ntSpeYield", "spe per sipm",
                           "event:spe0:spe1:spe2:spe3:spe4:spe5:spe6:spe7:spe8:spe9:spe10:spe11");
  ntSetTrigTime = new TNtuple("ntSetTrigTime", " trig time and val", "event:chan:time:val");
  ntTrigTime = new TNtuple("ntTrigTime", "trigger time ntuple", "entry:chan:firstTime:time:adc:ftime:fadc");
  ntChanSum = new TNtuple("ntchansum", "channel ntuple", "sum0:sum1:sum2:sum3:sum4:sum5:sum6:sum7:sum8:sum9:sum10:sum11:sum12:pass");
  hEventPass = new TH1D("EventPass", " event failures", TOTALCODES, 0, TOTALCODES);
  evCount = new TH1D("eventcount", "event count", CHANNELS, 0, CHANNELS);
  hNoPeak = new TH1D("noPeak", "no peak events count by channel", CHANNELS, 0, CHANNELS);
  histHitCount = new TH1D("hitCount", "hit count by channel", CHANNELS, 0, CHANNELS);
  histQSum = new TH1D("histqsum", "qsum by channel", CHANNELS, 0, CHANNELS);
  // nn/histqpe = new th1d("histqpe", "qpe by channel", CHANNELS, 0, CHANNELS);
  histQPrompt = new TH1D("histqprompt", "qprompt by channel", CHANNELS, 0, CHANNELS);
  histQSum->Sumw2();
  histQPrompt->Sumw2();
  hCosmicMult = new TH1D("CosmicMult", "CosmicMult", 10, 0, 10);
  hCosmicCut = new TH1D("CosmicCut", "cosmic total sum chan 12 ", 1000, 0, 10. * cosmicCut);
  hGammaCut = new TH1D("GammaCut", "gamma late sum chan 13 ", 1000, 0, 10. * gammaCut);

  //
  anaDir->cd();
  hTriggerTime = new TH1D("TriggerTime", " ave of trigger Sipm times ", 1000, 0, 1000);
  hTriggerTimeDiff = new TH1D("TriggerTimeDiff", " max trigger time diff ", 1000, 0, 1000);
  hTriggerShift = new TH1D("TriggerShift", " ave trigger time shift ", 200, -100, 100);
  hTriggerTimeAllVal = new TH1D("TriggerTimeAllVal", " first time val all channels ", 1000, 0, 1000);
  hTriggerTimeAllValPmt = new TH1D("TriggerTimeAllValPmt", " first time val Pmt ", 1000, 0, 1000);
  TString htitle;
  htitle.Form(" pre time < %lu normalized qpeak", triggerStart);
  hPreQpeak = new TH1D("PreQpeak", htitle, 100, 0, 10);
  htitle.Form(" pre time > %lu normalized qpeak", triggerEnd);
  hLateQpeak = new TH1D("LateQpeak", htitle, 100, 0, 10);
  hCountPre = new TH1D("CountPre", " hits sample<600 in sum", 20, 0, 20);
  htitle.Form("hits qpeak>%.2f SPE sample>%luin sum", latePeakCut, triggerEnd);
  hCountLate = new TH1D("CountLate", htitle, 20, 0, 20);
  htitle.Form("number of late time hits with qpeak>%.2f", latePeakCut);
  hCountLate->GetXaxis()->SetTitle(htitle);
  htitle.Form("hits qpeak>%.2f SPE sample>%lu in sum", latePeakCut, triggerEnd);
  hCountLateTime = new TH1D("CountLateTime ", htitle, 30, 0, 7500);
  hCountLateTime->GetXaxis()->SetTitle("sample time");
  hCountLateTime->Sumw2();
  hCountLateTimeQpeak = new TH2D("CountLateTimeQpeak", " sum qpeak vs time ", 750, 0, 7500, 80, 0, 20);
  hCountLateTimeQpeak->GetXaxis()->SetTitle("sample time");
  hCountLateTimeQpeak->GetYaxis()->SetTitle("qpeak [SPE]");

  printf("line1724 rawBr.size %lu \n", rawBr.size());
  for (unsigned i = 0; i < rawBr.size(); ++i)
  {
    unsigned ichan = i;
    hMult.push_back(new TH1D(Form("HitMultChan%i", ichan), Form("HitMultChan%i", ichan), 10, 0, 10));
    hWave.push_back(new TH1D(Form("waveChan%i", ichan), Form("WaveChan%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));
    hWave[hWave.size() - 1]->SetDirectory(nullptr);
    hChannelGaus.push_back(new TH1D(Form("channelGaus%i", ichan), Form("channelGaus%i", ichan), 600, -100, 500));
    noiseHist.push_back(new TH1D(Form("noiseChan%i", ichan), Form("noiseChan%i", ichan), 1000, 0, 1000));
    skewHist.push_back(new TH1D(Form("skewChan%i", ichan), Form("skewChan%i", ichan), 200, -3, 7));
    hEvRawWave.push_back(new TH1D(Form("evRawWave%i", ichan), Form("evRawWave%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));
    if (ichan > 8 && ichan < 12)
    {
      valHist.push_back(new TH1D(Form("valChan%i", ichan), Form("valChan%i", ichan), 1500, -500, 1000));
      valHistB.push_back(new TH1D(Form("valBadChan%i", ichan), Form("valBadChan%i", ichan), 1500, -500, 1000));
      hEvGaus.push_back(new TH1D(Form("evGaus%i", ichan), Form("evGaus%i", ichan), 200, -100, 100));
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
    double limit = 2000.;
    int nbins = 1000.;
    hTotSum.push_back(new TH1D(Form("TotPeakSumChan%i", i), Form("tot peak sum chan %i", i), nbins, 0, limit));
    hPreSum.push_back(new TH1D(Form("PrePeakSumChan%i", i), Form("pre peak sum chan %i", i), nbins, 0, limit));
    hTrigSum.push_back(new TH1D(Form("TrigPeakSumChan%i", i), Form("trig peak sum chan %i", i), nbins, 0, limit));
    hLateSum.push_back(new TH1D(Form("LatePeakSumChan%i", i), Form("late peak sum chan %i", i), nbins, 0, limit));
  }
  // one more for summed channel 13
  // hEvRawWave.push_back(new TH1D(Form("evRawWave%i", 13), Form("evRawWave%i", 13), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));

  threshValueHist = new TH2D("threshValueHist", " threshold crossings value channels by time  ", 7500, 0, 7500, 1000, 0, 100000);
  threshHist = new TH1D("threshHist", " threshold crossings trig channels ", 20, 0, 20);
  crossHist = new TH1D("crossHist", "  negative crossings non trigger channels", 100, 0, 100);
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

    hQSum.push_back(new TH1D(Form("QSumChan%i", ichan), Form("QSumChan%i", ichan), 1000, 0, 500.));
    hQPeak.push_back(new TH1D(Form("QPeakChan%i", ichan), Form("QPeakChan%i", ichan), 700, 0, 7.));
    hQSpe.push_back(new TH1D(Form("QSpeChan%i", ichan), Form("QSpeChan%i", ichan), 9, 0, 9.));
    sumWave.push_back(new TH1D(Form("sumWave%i", ichan), Form("sumWave%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));
    sumWaveB.push_back(new TH1D(Form("sumWaveBad%i", ichan), Form("sumWaveBad%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));
    sumHitWave.push_back(new TH1D(Form("sumHitWave%i", ichan), Form("sumHitWave%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));
    sumPeakWave.push_back(new TH1D(Form("sumPeakWave%i", ichan), Form("sumPeakWave%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));
  }
  // SPE Shapes
  for (unsigned ichan = 0; ichan < rawBr.size(); ++ichan)
  {
    hSPEShapeLate.push_back(new TH1D(Form("SPEShapeLateChan%i", ichan), Form("SPEShapeLateChan%i", ichan), 1000, 0, 1000));
    hSPEShapeLate[hSPEShapeLate.size() - 1]->SetMarkerStyle(20);
  }
  hSPEShape.resize(MaxSPEShape);
  for (int jspe = 0; jspe < MaxSPEShape; ++jspe)
  {
    for (unsigned ichan = 0; ichan < rawBr.size(); ++ichan)
    {
      hSPEShape[jspe].push_back(new TH1D(Form("SPE%iShapeChan%i", jspe + 1, ichan), Form("SPE%iShapeChan%i", jspe + 1, ichan), 1000, 0, 1000));
      hSPEShape[jspe][hSPEShape[jspe].size() - 1]->SetMarkerStyle(20);
    }
  }

  fout->cd();
  /// fout->ls();

  cout << " make hitFinder dets = " << CHANNELS << "  size " << rawBr[0]->rdigi.size() << endl;
  vector<int> chanList;
  for (int ichan = 0; ichan < CHANNELS; ++ichan)
    chanList.push_back(ichan);

  finder = NULL;
  finder = new hitFinder(fout, tbrun, tag, rawBr[0]->rdigi.size(), chanList, channelSigmaValue, nominalGain);
  if (!finder)
  {
    printf(" failed to make finder ");
    fout->Close();
    return 0;
  }

  int npass = 0;
  int nfail = 0;
  Long64_t lastEntry = firstEntry + nentries;
  printf("... total entries  %llu looping over %llu firstEntry %llu last %lld \n ", rawTree->GetEntries(), nentries, firstEntry, lastEntry);
  for (Long64_t entry = firstEntry; entry < lastEntry; ++entry)
  {
    tbrun->clear();
    rawTree->GetEntry(entry);

    // main ana routine
    int passBit = anaEvent(entry);
    if (passBit == 0)
    {
      ++npass;
    }
    else
    {
      ++nfail;
      // set pass bit and fill tbrun
      for (int idet = 0; idet < tbrun->detList.size(); ++idet)
      {
        tbrun->detList[idet]->pass = passBit;
      }
    }
    hEventPass->SetBinContent(passBit, hEventPass->GetBinContent(passBit) + 1);
    tbrun->fill();
    if (entry / 1000 * 1000 == entry)
    {
      printf("... entry %llu pass %u fail %u \n", entry, npass, nfail);
      hEventPass->Print("all");

      if (npass > 0)
      {
        printf(" \t hits by channel  \n");
        for (int ibin = 0; ibin < histHitCount->GetNbinsX() - 1; ++ibin)
          printf(" chan %i count %i frac %f ; zero %i \n", ibin,
                 int(histHitCount->GetBinContent(ibin + 1)), double(histHitCount->GetBinContent(ibin + 1)) / double(npass), int(hNoPeak->GetBinContent(ibin + 1)));
        printf("  \n");
      }
    }
  }
  printf(" \n \n At END OF FILE total pass  = %i fail %i  \n", npass, nfail);

  /*
  TString graphName = TString("slopeGraph");
  TString graphTitle = TString(Form("slope-graph-%s", shortName.c_str()));
  printf(" making slope graph %s \n", graphName.Data());

  for (unsigned i = 0; i < sumHitWave.size(); ++i)
  {
    sumHitWave[i]->Fit("expo", "Q0", "", 100, 300); // DEG suggests
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
    hChannelGaus[index]->Fit("gaus", "Q0", "", hChannelGaus[index]->GetMean() - 100, hChannelGaus[index]->GetMean() + 100);
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
  }*/

  printf(" ******* hit count summary ***** \n \t hits by channels %i   \n", histHitCount->GetNbinsX());
  for (int ibin = 0; ibin < histHitCount->GetNbinsX() - 1; ++ibin)
    printf(" chan %i count %i frac %f ; zero %i \n", ibin,
           int(histHitCount->GetBinContent(ibin + 1)), double(histHitCount->GetBinContent(ibin + 1)) / double(npass), int(hNoPeak->GetBinContent(ibin + 1)));
  printf("  \n");

  printf(" \n \t sums by channel with entries %.0f \n", hTotSum[0]->GetEntries());

  // calculate mean hits from waveforms
  std::vector<double> hitMean;
  std::vector<double> hitIntegral;
  for (int idet = 0; idet < sumHitWave.size(); ++idet)
  {
    double inte = sumPeakWave[idet]->Integral() / nominalGain;
    double mean = inte / double(npass);
    hitMean.push_back(mean);
    hitIntegral.push_back(inte);
  }

  // printf(" FINISHED npass %u nfail %u output file  %s \n", npass, nfail, fout->GetName());
  printf(" finished %i ( %i ) pass %i (%i) fail %i ( frac %0.3f ) output file %s  \n",
         npass + nfail,
         int(hEventPass->GetEntries()),
         npass, int(hEventPass->GetBinContent(0)),
         nfail,
         double(nfail) / double(npass + nfail),
         // int(hEventPass->GetBinContent(8 + 1)),
         fout->GetName());

  // hEventPass->Print("all");
  printf("pass fractions total = %.0f \n", hEventPass->GetEntries());
  for (int ibin = 0; ibin < hEventPass->GetNbinsX(); ++ibin)
  { // include error on poisson probability
    double nbin = hEventPass->GetBinContent(ibin);
    double ntot = hEventPass->GetEntries();
    double prob = nbin / ntot;
    double perror = sqrt(prob * (1. - prob) / ntot);
    printf(" bin %i fail %.f frac %.3f +/- %.3f name %s \n", ibin, hEventPass->GetBinContent(ibin), prob, perror, codeNames[ibin].Data());
  }

  for (int idet = 0; idet < hTotSum.size(); ++idet)
  {
    printf(" \t chan %i means: tot %.2f pre %.2f trig %.2f late %.2f \n", idet,
           hTotSum[idet]->GetMean(),
           hPreSum[idet]->GetMean(),
           hTrigSum[idet]->GetMean(),
           hLateSum[idet]->GetMean());
  }

  for (int idet = 0; idet < hitMean.size(); ++idet)
    printf("chan %i wave integral %.4E average hits per event %.4f \n ", idet, hitIntegral[idet], hitMean[idet]);

  printf("PMT HIT MULTIPLICITY cut %0.f \n", qpeakCosmicCut);
  hCosmicMult->Print("all");
  fout->Write();
  fout->Close();
  printf(" ***** FINISHED ****** %s entries %lld \n", fout->GetName(), nentries);
  return nentries;
}

anaCRun::anaCRun(TString theTag)
{
  // define failure code names

  // pass bit failures
  /*enum FAILURECODES
  {
    PASS = 0,
    GAUSFAIL = 0x1,
    EARLYCUT = 0x2,
    FIRSTTIME = 0x4,
    COSMIC = 0x8,
    GAMMA = 0x16,
    TOTALCODES
  };
  */

  for (unsigned ic = 0; ic < TOTALCODES; ++ic)
    codeNames.push_back(TString("mixed"));
  codeNames[PASS] = TString("pass");
  codeNames[GAUSFAIL] = TString("baseline");
  codeNames[EARLYCUT] = TString("earlycut");
  codeNames[FIRSTTIME] = TString("firsttime");
  codeNames[COSMIC] = TString("cosmic");
  codeNames[GAMMA] = TString("gamma");

  tag = theTag;
  // tbrun = new TBRun(tag);
  cout << " anaCRun::anaCRun instance of anaCRun gamma version  with tag= " << tag << " CHANNELS = " << CHANNELS - 1 << " diffStepSipm= " << diffStepSipm << " diffStepPmt= " << diffStepPmt << endl;

  rawBr.clear();

  for (int ichan = 0; ichan < CHANNELS; ++ichan)
  {
    TBRawEvent *rawEv = new TBRawEvent(ichan);
    rawEv->rdigi.resize(7500);
    rawEv->SetName(Form("rawChan%i", ichan));
    rawBr.push_back(rawEv);
  }
}
