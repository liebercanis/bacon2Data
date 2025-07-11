/***Ths is GAMMA version Sept 25 2024 **/
// revised Jan 15 2025
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
  // add two for summed waveform 0-12 channels, 13 summed
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

  // pass bit failures hex
  enum FAILURECODES
  {
    PASS = 0,
    BASEFAIL = 0x1,
    EARLYCUT = 0x2,
    FIRSTTIME = 0x4,
    COSMIC = 0x8,
    GAMMA = 0x10,
    TRIGFAIL = 0x20,
    TOTALCODES = 2 * TRIGFAIL
  };

  enum
  {
    FAILBITS = 7
  };

  std::vector<TString> bitNames;
  int failCode[FAILBITS];

  std::vector<TString> codeNames;

  int badEvent = 5671;
  int failGamma = 0;
  int failCosmic = 0;

  bool doNotOverWrite = true;
  bool theFirstFile = true;
  bool isSim = false;
  std::vector<TDet *> simDet;
  TH1D *hSimFoundTimeDiff;
  int badEventDirMax = 1000;
  int exampleDirMax = 1000;
  int missedDirMax = 1000;
  bool reportFailures = true;
  double noiseToSignal = 0.04;
  TBRun *tbrun;
  TFile *fout;
  TFile *fin;
  TTree *rawTree;
  TTree *simTree;
  // ntuples to check cuts
  TNtuple *ntBase;
  TNtuple *ntTrig;
  TNtuple *ntNonTrig;
  TH1D *hPreSumCut;
  TH1D *hCosmicCut;
  TH1D *hGammaCut;
  TNtuple *ntHit;
  TNtuple *ntSimMatch;
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
  TNtuple *ntFailures;
  vector<TH1D *> baseHist;
  vector<TH1D *> noiseHist;
  vector<TH1D *> skewHist;
  vector<TH1D *> sumWave;
  vector<TH1D *> sumHitWave;
  vector<TH1D *> sumPeakWave;
  vector<TH1D *> valHist;
  vector<TH1D *> sumWaveA;
  vector<TH1D *> sumWaveB;
  vector<TH1D *> valHistB;
  std::vector<std::vector<TH1D *>> sumWaveFail;

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

  TH1D *hTrigSumNoCut;
  TH1D *hTrigSumCut;
  std::vector<TH1D *> hTrigSumCutRatio;
  TH1D *hPreQpeak;
  TH1D *hLateQpeak;
  TH1D *hCountPre;
  TH1D *hCountLate;
  TH1D *hCountLateTime;
  TH2D *hCountLateTimeQpeak;
  TH2D *hTriangle;
  TH1D *evCount;
  TH1D *histQSum;
  TH1D *hEventPass;
  TH1D *hEventFail;
  TH1D *histHitCount;
  TH1D *hNoPeak;
  TH1D *hSumPMT;
  TH1D *threshHist;
  TH2D *threshValueHist;
  TH1D *crossHist;
  TH1D *hCosmicMult;
  // TH1D *histQPE;
  TH1D *histQPrompt;
  TH1D *hTriggerTimeDiff;
  TH1D *hTriggerTimeAllVal;
  TH1D *hTriggerTimeAllValPmt;
  TH1D *hTriggerHitTimeAll;
  TH1D *hTriggerTime;
  TH1D *hTriggerShift;

  // sim comparison histos by channel
  std::vector<TH1D *> hWaveHitFound;
  std::vector<TH1D *> hWaveHitMissed;
  std::vector<TH1D *> hWaveHitNoise;

  //
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
  // vector<TBWave *> waveList;/
  hitFinder *finder;
  TString tag;
  int currentBuffer;
  Long64_t currentBufferCount;
  Long64_t eventNumber;
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
  bool simTimeMatch(double stime, double ftime);

  /*
  void setTBRun(TBRun *theTBRun)
  {
    tbrun = theTBRun;
  }
  */
  void makeTernary(double a, double b, double c, double &x, double &y);
  std::vector<std::vector<double>> fixedDigi; // all the fixed waveforms
  std::vector<unsigned> trigTimes;
  std::vector<unsigned> sTrigTimes; // after correction
  std::vector<double> adcBin;
  std::vector<double> speCount;
  TDirectory *threshDir;
  TDirectory *earlyPeakDir;
  TDirectory *rawSumDir;
  TDirectory *exampleDir;
  TDirectory *missedDir;
  TDirectory *sumDir;
  TDirectory *anaDir;
  TDirectory *badEventDir;
  TDirectory *pmtDir;
  TDirectory *fftDir;
  TDirectory *templateDir;
  TDirectory *simDir;

  Long64_t nentries;
  double QPEPeak;
  //
  int MaxSPEShape = 4;
  unsigned trigStart = 600;
  int nominalTrigger = 753; // was 729; this is nominal trigger sample

  /**************** define nominal gains ***************/
  double nominalGain = 134.786401;     // 170.;     // was 160.0; set Jue 13 2025
  double nominalTrigGain = 735.688747; //
  double nominalQsumGain = 4940.503519;
  double nominalQsumTrigGain = 32056.789775;
  double nominalPmtGain = 502.;
  double nominalQsumPmtGain = 1713;
  double landauMax = 1.0; // 0.018063;
  double qsumGain[CHANNELS];
  //  227.4; // average
  //   double nominalGain = 160.0; // average
  unsigned firstTime;       // corrected trigger time for event
  unsigned timeOffset = 13; // changed from 17 may 13, 2024
  double passValEarlyCut = 100.0;
  /// double passValEarlyPmtCut = 225.0;
  ULong_t triggerEnd = 800; // 740;
  ULong_t lateTimeStart = 900;
  ULong_t triggerStart = 730; // 740;
  ULong_t timeVeryLateCut = 3500;
  /* need to tune these cuts on data */
  double trigRatioCutLow = 0.2;             // qsum fraction
  double trigRatioCutHigh = 0.8;            // qsum fraction
  double preSumCut = 4. * nominalGain;      ///
  double totCosmicCut = 50. * nominalGain;  //
  double lateGammaCut = 100. * nominalGain; //
  double trigSumCut = 3.0;

  double prePeakCut = 0.5;
  double latePeakCut = 3.5;                 // march 18 2024 2.5;
  double diffStepSipm = 3.;                 // 6 ns steps for SIPM
  double diffStepPmt = 1.;                  // back to one on Oct 15 2024
  double cosmicCut = 3.E3;                  // set Nov 3 2024
  double qpeakCosmicCut = 3. * nominalGain; // 3*SPE
  double hitThresholdPmt = 30.;             // set Nov 13 2024
};
// I do this in two place, so I wanted to be sure to do it the same.
bool anaCRun::simTimeMatch(double stime, double ftime)
{
  bool rc = false;
  double timeDiff = stime - ftime;
  if (timeDiff > 20 && timeDiff < 80)
    rc = true;
  return rc;
}

//// https://mathworld.wolfram.com/TernaryDiagram.html
void anaCRun::makeTernary(double a, double b, double c, double &x, double &y)
{
  double s = a + b + c;
  x = 0.5 * (a + 2. * b) / s;
  y = sqrt(3.) / 2. * a / s;
}

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
    // printf("line299 %i raw %u base %f \n ", j, rawBr[ic]->rdigi[j], idet->base);
    if (val > adc)
    {
      adc = val;
      utime = j + off;
    }
  }
  // printf("line304 ..... event %lld chan %i adc %f utime %u \n", eventNumber, ic, adc, utime);
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
  // start out with array filled with zeros
  for (unsigned ib = 0; ib < NONSUMCHANNELS; ++ib)
  {
    // FIX time shift bug Nov 24 2024!
    std::vector<double> fDigi;
    fDigi.clear();
    fDigi.resize(rawBr[0]->rdigi.size());
    std::fill(fDigi.begin(), fDigi.end(), 0);
    int timeShift = nominalTrigger - firstTime; // nominalTrigger = 730
    // printf("line434 first %u  shift %i off %u chan %u \n", firstTime, timeShift, timeOffset, ib);
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
      // val *= nominalGain / sipmGain[ib];
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
  // get new nominal gains
  double newNominalGain = 0;
  double newNominalTrigGain = 0;

  printf("line466 got %lu gains \n", sipmGain.size());
  for (unsigned long j = 0; j < sipmGain.size(); ++j)
  {
    printf(" %lu  gain %.4f error %.4f   \n", j, sipmGain[j], sipmGainError[j]);
    if (j < 9)
      newNominalGain += sipmGain[j];
    if (j > 8 && j < 12)
      newNominalTrigGain += sipmGain[j];
  }
  newNominalGain /= double(9);
  newNominalTrigGain /= double(3);
  printf(" GGGGGGGG nominal gains %f trig %f  GGGGGGGGGG\n", newNominalGain, newNominalTrigGain);
}

bool anaCRun::readGains(TString fileName)
{
  for (int i = 0; i < CHANNELS; ++i)
    qsumGain[i] = nominalQsumGain;
  qsumGain[9] = nominalQsumTrigGain;
  qsumGain[10] = nominalQsumTrigGain;
  qsumGain[11] = nominalQsumTrigGain;

  /* define nominal */
  sipmGain.clear();
  sipmGainError.clear();
  sipmGain.resize(NONSUMCHANNELS);
  sipmGainError.resize(NONSUMCHANNELS);
  /*    preliminaru gains
        no trig 107
        trig 509
        PMT 395
  */
  for (unsigned long j = 0; j < sipmGain.size(); ++j)
  {
    if (j < 9)
    {
      sipmGain[j] = nominalGain;
      sipmGainError[j] = sqrt(nominalGain);
    }
    else if (j < 12)
    {
      sipmGain[j] = nominalTrigGain;
      sipmGainError[j] = sqrt(nominalTrigGain);
    }
    else
    {
      sipmGain[j] = nominalPmtGain;
      sipmGainError[j] = sqrt(nominalPmtGain);
    }
  }
  /* look for gain file */
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
    printf(" fopen couldnt open template file %s\n", fileName.Data());
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
  fin->GetObject("gains-05_19_2025-05_19_2025", gGain);
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
  sumWaveA.clear();
  sumWaveB.clear();
  sumWaveFail.clear();
  sumWaveFail.resize(FAILBITS);
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
  // updated June 11 2025 file run-05_19_2025-file.root
  // set to 2  times digi sigma
  for (unsigned long j = 0; j < channelSigmaValue.size(); ++j)
  {
    // chanThreshold[j] = 2. * 36.4; from bad data
    chanThreshold[j] = 4. * 36.; //
  }
  // trigger SIPMs do three sigma
  chanThreshold[9] = 4. * 99.;
  chanThreshold[10] = 4. * 99.;
  chanThreshold[11] = 4. * 99.;
  chanThreshold[12] = 3. * 2.4; // based on histogram sigma, had been 5.3;
  chanThreshold[13] = 3. * 33.; // should be same as trigger sipm
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
  rawTree = nullptr;
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

  simTree = nullptr;
  isSim = false;
  fin->GetObject("SimTree", simTree);
  if (simTree)
  {
    isSim = true;
    printf("line 696  anaEvent file %s THIS IS SIMULATION\n", fileName.Data());
  }

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
  eventNumber = entry;
  TTree *tree = NULL;
  fout->GetObject("RunTree", tree);
  if (!tree)
  {
    printf("line 531 ERROR!! anaEvent no tree event %lld \n", entry);
    fout->ls();
  }
  // get sim branches
  simDet.clear();
  if (simTree)
  {
    simTree->GetEntry(entry); // have to load the entry
    // simTree->GetListOfBranches()->ls();
    //    get branch pointers and save in detList
    TIter next(simTree->GetListOfBranches());
    TBranchElement *aBranch = NULL;
    while ((aBranch = (TBranchElement *)next()))
    {
      simDet.push_back((TDet *)aBranch->GetObject());
    }
    // print out sim branches for fist entry
    if (entry == 0)
    {
      printf("getting simDet branches %lu \n", simDet.size());
      for (unsigned idet = 0; idet < simDet.size(); ++idet)
        printf("GOTSIM det %i %s \n", idet, simDet[idet]->GetName());
    }
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

  // also fill chan 13
  TDet *tdet13 = tbrun->getDet(NONSUMCHANNELS); // get channel 13 det
  tdet13->clear();
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
    {
      printf("rdigi bad size  event %lld channel %u %lu \n", entry, ib, rawBr[ib]->rdigi.size());
      continue;
    }
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
    // printf("line761 \t\t\t event %lld base %f \n", entry, base);

    // baseline correction from fitted Gaussian
    hEvGaus[ib]->Reset("ICES");
    hEvRawWave[ib]->Reset("ICES");
    // for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
    for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
    {
      double val = double(rawBr[ib]->rdigi[j]) - base; // base is > digi value!
      if (j < trigStart)
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
        printf("@line804 failed Baseline Cut event %llu chan %i base %f ave %f status %i \n", entry, ib, base, ave, fitStatus);
      if (badEventDir->GetList()->GetEntries() < badEventDirMax)
      {
        badEventDir->cd();
        TH1D *EvRawWave = (TH1D *)hEvRawWave[ib]->Clone(Form("EvRawBaseFailEvent%lld-Ch%i", entry, ib));
        EvRawWave->SetTitle(Form("EvRawBaseFailEvent%lld-Ch%i", entry, ib));
        TH1D *hEvGausClone = (TH1D *)hEvGaus[ib]->Clone(Form("EvGausEv%lldchan%imean%.2fsigma%.2fstatus%i", entry, ib, fitMean, sigma, fitStatus));
      }
      passBit |= BASEFAIL;
    }
    ntBase->Fill(entry, ib, base, base + fitMean, fitMean, sigma, fitStatus);

    // printf("@line652 baseline %lld chan %u base %f ave %f  \n", entry, ib, base, ave);

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
    idet->totSum = 0;
    idet->preSum = 0;
    idet->trigSum = 0;
    idet->lateSum = 0;
    idet->totPeakSum = 0;
    idet->prePeakSum = 0;
    idet->trigPeakSum = 0;
    idet->latePeakSum = 0;

    /*********
     * make sums for cuts
     *********/
    for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
    {
      double val = double(rawBr[ib]->rdigi[j]) - idet->base;
      if (val < 3. * idet->sigma)
        continue;
      idet->totSum += val / qsumGain[ib];   // convert to approximate number of photons
      tdet13->totSum += val / qsumGain[ib]; // convert to approximate number of photons
      if (j < triggerStart)
        idet->preSum += val / qsumGain[ib];
      if (j > lateTimeStart)
      {
        idet->lateSum += val / qsumGain[ib];
        tdet13->lateSum += val / qsumGain[ib];
      }
      // channel sum
      baseHist[ichan]->Fill(val);
      if (hChannelGaus.size() > 0)
        hChannelGaus[ib]->Fill(val);
    }

    /* add maxAdc */
    double maxAdc;
    int maxSample;
    getMaxRawAdc(ib, idet->base, maxAdc, maxSample);
    idet->maxAdc = maxAdc;
    idet->maxSample = maxSample;
    // I have added this to the TDet as  maxSample maxAdc
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

  double preSum = 0;
  /*********  early cut to remove photons before trigger *******/
  for (unsigned ib = 0; ib < NONSUMCHANNELS; ++ib)
  {
    unsigned ichan = ib;
    TDet *idet = tbrun->getDet(ichan);
    preSum += idet->preSum / qsumGain[ib];
  }

  hPreSumCut->Fill(preSum);
  if (preSum > preSumCut)
  {
    printf("line919 fail EARLYCUT cut %f val %f \n", preSumCut, preSum);
    if (badEventDir->GetList()->GetEntries() < badEventDirMax)
    {
      badEventDir->cd();
      TH1D *EvRawWave = (TH1D *)hEvRawWave[13]->Clone(Form("EvRawEarlyCutFailEvent%lld-Ch%i", entry, 13));
      EvRawWave->SetTitle(Form("EvRawBaselineEvent%lld-Ch%i", entry, 13));
    }
    passBit |= EARLYCUT;
  }

  /******* trigger time cut *********/
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

    /* set passBit 1 val cut is different for PMT
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
    }*/

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

  /* ave non trig  sipm times before shift
  double nonTimeAve = 0;
  double nonTimeSigma = 0;
  unsigned chanBad2;
  double dmax2;
  getTriggerTimeStats(&trigTimes[6], nonTimeAve, nonTimeSigma, chanBad2, dmax2);
  */

  hTriggerTime->Fill(double(firstTime));
  // Bug fix-- was and fixed to Or Jan 27 2025
  if (firstTime > triggerEnd || firstTime < triggerStart)
  {
    if (reportFailures)
      printf("@line893 failed triggerEnd event %llu cut %lu time %u \n", entry, triggerEnd, firstTime);

    if (badEventDir->GetList()->GetEntries() < badEventDirMax)
    {
      badEventDir->cd();
      TH1D *EvRawWave = (TH1D *)hEvRawWave[13]->Clone(Form("EvRawFirstTimelEvent%lld-Ch%i", entry, 13));
      EvRawWave->SetTitle(Form("EvRawFirstTimeEvent%lld-Ch%i", entry, 13));
    }
    passBit |= FIRSTTIME;
  }

  /******  cosmic cut based on light in PMT *********/
  TDet *tdetPmt = tbrun->getDet(12);
  hCosmicCut->Fill(tdetPmt->totSum);
  if (tdetPmt->totSum > totCosmicCut)
  {
    passBit |= COSMIC;
    if (badEventDir->GetList()->GetEntries() < badEventDirMax)
    {
      badEventDir->cd();
      TH1D *EvRawWave = (TH1D *)hEvRawWave[12]->Clone(Form("EvRawCosmicEvent%lld-Ch%i", entry, 12));
      EvRawWave->SetTitle(Form("EvRawCosmicEvent%lld-Ch%i", entry, 12));
    }
    ++failCosmic;
    if (reportFailures)
      printf("@line1077 failed cosmic event %llu bit %i cut %E totSum %E  \n", entry, passBit, cosmicCut, tdetPmt->totSum);
  }

  /********** gamma cut *********/
  TDet *idet9 = tbrun->getDet(9);
  TDet *idet10 = tbrun->getDet(10);
  TDet *idet11 = tbrun->getDet(11);
  hGammaCut->Fill(tbrun->getDet(13)->lateSum);
  if (tbrun->getDet(13)->lateSum > lateGammaCut)
  {
    if (badEventDir->GetList()->GetEntries() < badEventDirMax)
    {
      badEventDir->cd();
      TH1D *EvRawWave = (TH1D *)hEvRawWave[13]->Clone(Form("EvRawGammaEvent%lld-Ch%i", entry, 13));
      EvRawWave->SetTitle(Form("EvRawGammaEvent%lld-Ch%i", entry, 13));
    }
    passBit |= GAMMA;
    ++failGamma;
    if (reportFailures)
    {
      printf("@line1090 failed gamma event %llu bit %i cut %E lateSum %E \n", entry, passBit, lateGammaCut, tbrun->getDet(13)->lateSum);
      printf("@line1091  %f %f %f sum %E \n", idet9->lateSum, idet10->lateSum, idet11->lateSum, idet9->lateSum + idet10->lateSum + idet11->lateSum);
    }
  }

  // if (passBit == 0)
  for (unsigned ib = 0; ib < NONSUMCHANNELS; ++ib)
  {
    ntNonTrig->Fill(double(entry), double(ib), tbrun->getDet(ib)->totSum);
  }
  double trigSum = idet9->totSum + idet10->totSum + idet11->totSum;
  hTrigSumNoCut->Fill(trigSum);

  /******   trigger cut ********/
  /* try a cut like TUM */
  double triggerSum = idet9->totSum + idet10->totSum + idet11->totSum;
  double qFraction[3];
  qFraction[0] = idet9->totSum / triggerSum;
  qFraction[1] = idet10->totSum / triggerSum;
  qFraction[2] = idet11->totSum / triggerSum;
  for (unsigned iratio = 0; iratio < hTrigSumCutRatio.size(); ++iratio)
    hTrigSumCutRatio[iratio]->Fill(qFraction[iratio]);

  int failsTrigger = 0;
  if (qFraction[0] < trigRatioCutLow || qFraction[0] > trigRatioCutHigh)
    failsTrigger |= 0x2;
  if (qFraction[1] < trigRatioCutLow || qFraction[1] > trigRatioCutHigh)
    failsTrigger |= 0x4;
  if (qFraction[2] < trigRatioCutLow || qFraction[2] > trigRatioCutHigh)
    failsTrigger |= 0x8;

  double xternQ, yternQ;
  makeTernary(qFraction[0], qFraction[1], qFraction[2], xternQ, yternQ);
  // printf("line1097 %f %f %f %f %f \n", qFraction[0], qFraction[1], qFraction[2], xternQ, yternQ);
  ntTrig->Fill(double(entry), idet9->totSum, idet10->totSum, idet11->totSum, tdet13->totSum, qFraction[0], qFraction[1], qFraction[2], xternQ, yternQ, failsTrigger);
  if (failsTrigger != 0)
  {
    if (badEventDir->GetList()->GetEntries() < badEventDirMax)
    {
      badEventDir->cd();
      TH1D *EvRawWave = (TH1D *)hEvRawWave[13]->Clone(Form("EvRawTrigFailEvent%lld-Ch%i", entry, 13));
      EvRawWave->SetTitle(Form("EvRawTrigFailineEvent%lld-Ch%i", entry, 13));
    }
  }
  // softer trig cut
  if (triggerSum < trigSumCut)
    passBit |= TRIGFAIL;
  // printf("line1054 TRIGFAIL %lld cut %f chan 9 %f,%f chan 10 %f chan 11 %f ratio 9-10 %f ratio 9-11 %f ratio 10-11 %f \n", entry, trigRatioCutLow, trigRatioCutHigh, idet9->totSum, idet10->totSum, idet11->totSum, qFraction[0], qFraction[1], qFraction[2]);
  if (failsTrigger == 0)
    hTrigSumCut->Fill(trigSum);

  // fill triangle plot
  if (passBit == 0)
    hTriangle->Fill(xternQ, yternQ);

  /********************************************************
   * now that we have the firstTime
        align to nominalTrigger
        defined as timeShift>0 shift right
        normalize to nominal gain
   ********************************************************/
  doTimeShiftAndNorm();

  /***********   fill ntuple for threshold setting loop over channels ************/
  for (unsigned long ib = 0; ib < NONSUMCHANNELS; ++ib)
  {
    digi.clear();
    digi = fixedDigi[ib]; // get the time shifted and normed waveform
    TDet *idet = tbrun->getDet(ib);

    double peakMax = 0;
    // line931  recalculate and save the average
    /*
    double theAve = 0;
    for (unsigned j = 0; j < digi.size(); ++j)
    {
      theAve += digi[j];
    }
    theAve /= double(rawBr[ib]->rdigi.size());
    idet->ave = theAve;
    */
    // do digi sums on fixedDigi
    idet->totSum = 0;
    idet->preSum = 0;
    idet->lateSum = 0;
    for (unsigned j = 0; j < digi.size(); ++j)
    {
      idet->totSum += digi[j] / qsumGain[ib];
      if (j < triggerStart)
        idet->preSum += digi[j] / qsumGain[ib];

      if (j > triggerStart && j < triggerEnd)
      {
        idet->trigSum += digi[j] / qsumGain[ib];
        if (digi[j] > peakMax)
          peakMax = digi[j];
      }

      if (j > lateTimeStart)
        idet->lateSum += digi[j] / qsumGain[ib];
    }
    // add some other variables
    idet->pass = passBit;
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

  /********************************************************
        start of second channel loop doing pulse finding
   *******************************************************/

  // also fill chan 13
  // TDet *tdet13 = tbrun->getDet(NONSUMCHANNELS); // get channel 13 det
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

    evCount->Fill(ib);                        // chan 0 from GetBinContent(0)
    double hitThreshold = 0.75 * nominalGain; // 500.0;
    if (trig)
      hitThreshold = 0.75 * nominalTrigGain;
    if (ib == 12)
      hitThreshold = 0.75 * nominalPmtGain; // this is 5*(6 sigma noise)
    double theStep = diffStepSipm;
    if (ib == 12)
    {
      theStep = diffStepPmt;
    }
    /******************************** call hitFinder  *******************************/
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
  // TDet *tdetPmt = tbrun->getDet(12);
  TDet *tdet1 = tbrun->getDet(1);
  TDet *tdet4 = tbrun->getDet(4);
  for (unsigned ihit = 0; ihit < tdetPmt->hits.size(); ++ihit)
  {
    if (tdetPmt->hits[ihit].qpeak > qpeakCosmicCut)
      ++nCosmicHits;
  }

  int nPmtTrigger = 0;
  for (unsigned ihit = 0; ihit < tdetPmt->hits.size(); ++ihit)
  {
    if (tdetPmt->hits[ihit].startTime > triggerStart && tdetPmt->hits[ihit].startTime < triggerEnd)
      ++nPmtTrigger;
  }

  for (unsigned ihit = 0; ihit < tdet1->hits.size(); ++ihit)
  {
    if (tdet1->hits[ihit].startTime > triggerStart && tdet1->hits[ihit].startTime < triggerEnd)
      ++nPmtTrigger;
  }
  for (unsigned ihit = 0; ihit < tdet4->hits.size(); ++ihit)
  {
    if (tdet4->hits[ihit].startTime > triggerStart && tdet4->hits[ihit].startTime < triggerEnd)
      ++nPmtTrigger;
  }

  // if (nPmtTrigger > 0)
  //   printf("line1146 event %llu size hits 1,4,pmt %lu %lu %lu nPmtTrigger %i \n", entry, tdet1->hits.size(), tdet4->hits.size(), tdetPmt->hits.size(), nPmtTrigger);

  // do cosmic cut based on large pulse counting
  // hCosmicMult->Fill(double(nCosmicHits));
  /*
  hCosmicCut->Fill(tdetPmt->totSum);
  if (nCosmicHits > 0 || tdetPmt->totSum > totCosmicCut)
  {
    passBit |= COSMIC;
    ++failCosmic;
    if (reportFailures)
      printf("@line1077 failed cosmic event %llu bit %i cut %E totSum %E nCosmicHits %i \n", entry, passBit, cosmicCut, tdetPmt->totSum, nCosmicHits);
    if (badEventDir->GetList()->GetEntries() < badEventDirMax)
    {
      badEventDir->cd();
      TH1D *EvRawWave = (TH1D *)hEvRawWave[12]->Clone(Form("EvRawCosmicEvent%lldTotSum%.0E-Ch%i", entry, tdetPmt->totSum, 12));
      EvRawWave->SetTitle(Form("EvRawCosmicEvent%lldTotSum%.3E-Ch%i", entry, tbrun->getDet(12)->totSum, 12));
    }
  }
    */

  // look at good PMT events
  if (passBit == 0 && tdetPmt->peakMax > hitThresholdPmt && pmtDir->GetList()->GetEntries() < 1000)
  {
    pmtDir->cd();
    // printf("@line1171 print event %llu peakMax %E \n", entry, tdetPmt->peakMax);
    TH1D *EvRawWave = (TH1D *)hEvRawWave[12]->Clone(Form("EvRawPMTEvent%lldVal%.0E-Ch%i", entry, tdetPmt->peakMax, 12));
    EvRawWave->SetTitle(Form("EvRawPMTEvent%lldVal%.3E-Ch%i", entry, tdetPmt->peakMax, 12));
    if (tdetPmt->hits.size() > 0)
      finder->plotEvent(pmtDir, tdetPmt->channel, entry);
  }
  // for  late hits
  /*
  int nLate9 = 0;
  TDet *tdet9 = tbrun->getDet(9);
  int startLast = 0;
  for (unsigned ihit = 0; ihit < tdet9->hits.size(); ++ihit)
  {
    if (tdet9->hits[ihit].startTime > 7435)
    {
      ++nLate9;
      startLast = tdet9->hits[ihit].startTime;
    }
  }
  */
  /* just collect some events */
  TDet *tdet9 = tbrun->getDet(9);
  TDet *tdet10 = tbrun->getDet(10);
  TDet *tdet11 = tbrun->getDet(11);
  /* refill */
  // baseline correction from fitted Gaussian

  /*
  if (exampleDir->GetList()->GetEntries() < exampleDirMax)
  {
    exampleDir->cd();
    TH1D *EvRawWave = (TH1D *)hEvRawWave[9]->Clone(Form("EvRawEvent%lld-Ch%i-totSum%.0f", entry, 9, tdet9->totSum));
    EvRawWave->SetTitle(Form("EvRawEvent%lld-Ch%i", entry, 9));
    EvRawWave = (TH1D *)hEvRawWave[10]->Clone(Form("EvRawEvent%lld-Ch%i-totSum%.0f", entry, 10, tdet10->totSum));
    EvRawWave->SetTitle(Form("EvRawEvent%lld-Ch%i", entry, 10));

    EvRawWave = (TH1D *)hEvRawWave[11]->Clone(Form("EvRawEvent%lld-Ch%i-totSum%.0f", entry, 11, tdet11->totSum));
    EvRawWave->SetTitle(Form("EvRawEvent%lld-Ch%i", entry, 11));
    // finder->plotEvent(exampleDir, tdet9->channel, entry);
    //  printf("@line1192 print event %llu start %i printed %i \n", entry, startLast, exampleDir->GetList()->GetEntries());
  }
    */

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
    if (hitStartTime < triggerStart && tdet13->hits[ihit].qpeak > prePeakCut)
    {
      ++nPreHits;
      // printf("event preHits %llu cut %lu hitStartTime %lu  qpeak %.2f nPreHits %i \n", entry, triggerStart, hitStartTime, hiti.qpeak, nPreHits);
    }
    if (hitStartTime > triggerEnd)
      hLateQpeak->Fill(tdet13->hits[ihit].qpeak / nominalGain);
    if (hitStartTime > triggerEnd && tdet13->hits[ihit].qpeak / nominalGain > latePeakCut)
    {
      ++nLateHits;
      // printf("event lateHits %llu cut %lu hitStartTime %lu  qpak %.2f nLateHits %i \n", entry, firstTimeCut, hitStartTime, hiti.qpeak, nLateHits);
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
    // update pass bit
    tbrun->getDet(ib)->pass = passBit;
    unsigned totHits = tbrun->getDet(ib)->hits.size();
    hMult[ib]->Fill(double(totHits)); // hit multiplicity
    ntFailures->Fill(entry, ib, totHits, passBit);
    // if (totHits > 0)
    //   printf("line1239 chan %u tot hits %u pass %i \n ", ib, totHits, passBit);

    /* take care here for summed ib=CHANNELS-2 and set appropriate hitThreshold */
    digi.clear();
    digi = fixedDigi[ib];
    if (passBit == 0)
    { // good waves
      for (unsigned j = 0; j < digi.size(); ++j)
      {
        sumWave[ib]->SetBinContent(j + 1, sumWave[ib]->GetBinContent(j + 1) + digi[j]);
        valHist[ib]->Fill(digi[j]);
      }
    } // check cosmic,gamma failure events

    // make this sum All instead of Bad
    for (unsigned j = 0; j < digi.size(); ++j)
    {
      sumWaveA[ib]->SetBinContent(j + 1, sumWaveA[ib]->GetBinContent(j + 1) + digi[j]);
    }

    // make this for Bad
    if (passBit != 0)
      for (unsigned j = 0; j < digi.size(); ++j)
      {
        sumWaveB[ib]->SetBinContent(j + 1, sumWaveB[ib]->GetBinContent(j + 1) + digi[j]);
        valHistB[ib]->Fill(digi[j]);
      }

    // sum wave by failure code
    for (int ic = 0; ic < FAILBITS; ++ic)
    {
      if (passBit & failCode[ic])
        for (unsigned j = 0; j < digi.size(); ++j)
        {
          sumWaveFail[ic][ib]->SetBinContent(j + 1, sumWaveFail[ic][ib]->GetBinContent(j + 1) + digi[j]);
        }
    }

    if (passBit == 0 && totHits > 0)
    { // waves with hits
      for (unsigned j = 0; j < digi.size(); ++j)
      {
        sumHitWave[ib]->SetBinContent(j + 1, sumHitWave[ib]->GetBinContent(j + 1) + digi[j]);
      }
    } // check cosmic,gamma failure events
    // if ((passBit & int(pow(2, 3))) != 0 | (passBit & int(pow(2, 4))) != 0)
    //{

    /* just collect some events */
    if (passBit == 0 && tbrun->getDet(ib)->hits.size() > 0 && ib < 9)
    //&& (tbrun->getDet(ib)->hits[0].qpeak > 200 && tbrun->getDet(ib)->hits[0].qpeak < 250)
    {
      if (exampleDir->GetList()->GetEntries() < exampleDirMax)
      {
        exampleDir->cd();
        TH1D *EvRawWave = (TH1D *)hEvRawWave[ib]->Clone(Form("EvRawEvent%lld-Ch%i-qpeak%0.f", entry, ib, tbrun->getDet(ib)->hits[0].qpeak));
        EvRawWave->SetTitle(Form("EvRawEvent%lld-Ch%i", entry, ib));
        finder->plotEvent(exampleDir, tbrun->getDet(ib)->channel, entry);
        // printf("@line1192 print event %llu start %i printed %i \n", entry, startLast, exampleDir->GetList()->GetEntries());
      }
    }

    //}
  }

  // if (passBit != 0) return passBit;
  // printf("line818  event %lld passbit %i \n",entry,passBit);
  if (passBit != 0)
  {
    /* collect example of failing evnets */
    // printf("@line913 event %lld passBit %i det %i nhits %u \n",
    //        entry, int(passBit), NONSUMCHANNELS, tbrun->detList[NONSUMCHANNELS]->nhits());
    return passBit;
  }

  /***************************************
  **** good events, passBit ==0 ******
  ****************************************/
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

    /*
    ****        loop over all hits for this detector
    */
    for (unsigned ihit = 0; ihit < tdet->hits.size(); ++ihit)
    {
      TDetHit thit = tdet->hits[ihit];

      if (thit.qpeak < 1)
        printf("line822 chan %i ihit %i startTime %i  peak %f\n", tdet->channel, ihit, int(thit.startTime), thit.qpeak);
      // do not scale these June 11 2025
      hQSum[idet]->Fill(thit.qsum);
      hQPeak[idet]->Fill(thit.qpeak);
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
      // sumHitWave[idet]->SetBinContent(thit.firstBin + 1, sumHitWave[idet]->GetBinContent(thit.firstBin + 1) + thit.qsum);
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

    /*
    ********** if simulation double loop for sim comparison
    */
    if (simTree && idet < 13) // no summed det in simulation
    {
      // each found hit is matched only once
      std::vector<int> foundHitUsed;
      foundHitUsed.resize(tdet->hits.size(), 0); // size and set to zero
      TDet *sdet = simDet[idet];                 // get sim det
      std::vector<int> missedHit;
      for (unsigned isim = 0; isim < sdet->hits.size(); ++isim) // loop over sim hits
      {
        missedHit.clear();
        TDetHit simHit = sdet->hits[isim];
        bool hitMatch = false;
        double timeDiff;
        int ihitMatch = -1;
        for (unsigned ihit = 0; ihit < tdet->hits.size(); ++ihit) // loop over hitFinder hits
        {
          TDetHit finderHit = tdet->hits[ihit];
          // printf(" fimder chan %i sim %i \n", tdet->channel, sdet->channel);
          int foundTimeBim = hWaveHitNoise[0]->FindBin(finderHit.startTime);

          // is hitFinder hit a sim hit?
          timeDiff = simHit.startTime - finderHit.startTime;
          hSimFoundTimeDiff->Fill(simHit.startTime - finderHit.startTime);
          if (simTimeMatch(simHit.startTime, finderHit.startTime) && foundHitUsed[ihit] == 0)
          {
            ihitMatch = ihit;
            foundHitUsed[ihit] = 1;
          }
          printf("ev %lld channel %i sim hit %i time %f found %i time %f isused %i \n", entry, sdet->channel, isim, simHit.startTime, ihit, finderHit.startTime, foundHitUsed[ihit]);
        } // loop over finder hits
        int simTimeBin = hWaveHitFound[0]->FindBin(simHit.startTime);
        if (ihitMatch != -1) // match
        {
          hWaveHitFound[idet]->SetBinContent(simTimeBin, hWaveHitFound[idet]->GetBinContent(simTimeBin) + 1);
          ntSimMatch->Fill(double(entry), 1, double(sdet->channel), double(isim), double(ihitMatch),
                           simDet[idet]->hits[isim].startTime - tdet->hits[ihitMatch].startTime, simDet[idet]->hits[isim].startTime, tdet->hits[ihitMatch].startTime,
                           simDet[idet]->hits[isim].qpeak, tdet->hits[ihitMatch].qpeak);
        }
        else // no match
        {
          missedHit.push_back(isim);
          hWaveHitMissed[idet]->SetBinContent(simTimeBin, hWaveHitMissed[idet]->GetBinContent(simTimeBin) + 1);
        }

      } // end loop over sim hits

      // print out events with missed hits
      if (missedHit.size() > 0 && missedDir->GetList()->GetEntries() < missedDirMax)
      {
        missedDir->cd();
        TH1D *EvRawWave = (TH1D *)hEvRawWave[idet]->Clone(Form("EvRawEvent%lld-Ch%i-timeBin%i", entry, idet,
                                                               int(sdet->hits[missedHit[0]].startTime)));
        EvRawWave->SetTitle(Form("EvRawEvent%lld-Ch%i", entry, idet));
        finder->plotEvent(missedDir, tbrun->getDet(idet)->channel, entry);
        // printf("@line1192 print event %llu start %i printed %i \n", entry, startLast, missedDir->GetList()->GetEntries());
      }

      if (missedHit.size() > 0)
      {
        for (unsigned isim = 0; isim < missedHit.size(); ++isim) // loop over missed sim
        {
          int imissed = missedHit[isim];                            // index of missed hit
          for (unsigned ihit = 0; ihit < tdet->hits.size(); ++ihit) // loop over found hits
            ntSimMatch->Fill(double(entry), 0, double(sdet->channel), double(imissed), double(ihit),
                             simDet[idet]->hits[imissed].startTime - tdet->hits[ihit].startTime, simDet[idet]->hits[imissed].startTime, tdet->hits[ihit].startTime,
                             simDet[idet]->hits[imissed].qpeak, tdet->hits[ihit].qpeak);
        }
      }

      // reverse order for noise hits
      for (unsigned ihit = 0; ihit < tdet->hits.size(); ++ihit) // loop over hitFinder hits
      {
        bool hitMatch = false;
        TDetHit finderHit = tdet->hits[ihit];
        for (unsigned isim = 0; isim < sdet->hits.size(); ++isim) // loop over sim hits
        {
          TDetHit simHit = sdet->hits[isim];
          if (simTimeMatch(simHit.startTime, finderHit.startTime))
          {
            hitMatch = true;
          }
        }
        int foundTimeBin = hWaveHitNoise[0]->FindBin(finderHit.startTime);
        if (!hitMatch)
        {
          hWaveHitNoise[idet]->SetBinContent(foundTimeBin, hWaveHitNoise[idet]->GetBinContent(foundTimeBin + 1));
        }
      }
    } // simTree in file

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

  cout << " anaCRunFile  for rootData/ input file shortName= " << theFile << endl;

  if (!openFile(theFile)) // and get branches
  {
    printf("anaCRun no such file %s \n", theFile.Data());
    return -1;
  }

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

  templateDir = fout->mkdir("templateDir");
  rawSumDir = fout->mkdir("rawSumDir");
  badEventDir = fout->mkdir("badEventDir");
  pmtDir = fout->mkdir("pmtDir");
  exampleDir = fout->mkdir("exampleDir");
  missedDir = fout->mkdir("missedDir");
  threshDir = fout->mkdir("threshDir");
  earlyPeakDir = fout->mkdir("earlyPeakDir");
  anaDir = fout->mkdir("anadir");
  sumDir = fout->mkdir("sumDir");
  TDirectory *finderDir = fout->mkdir("finderDir");
  TDirectory *splitDir = fout->mkdir("splitDir");
  // TDirectory *fitSingletDir = fout->mkdir("fitSingletDir");
  TDirectory *sumWaveDir = fout->mkdir("sumWaveDir");
  TDirectory *fftDir = fout->mkdir("fftDir");
  TDirectory *simDir = fout->mkdir("simDir");
  fout->ls();

  currentBuffer = -1;
  currentBufferCount = 0;
  printf(" anaCRun::anaCRunFile starting anaCRun file %s maxEntries %llu firstEntry %llu \n",
         theFile.Data(), maxEntries, firstEntry);

  // new gain file
  TString gainFileName = TString(getenv("BOBJ")) + TString("/gains-2025-06-18-18-48.root");
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

  ntSimMatch = new TNtuple("ntSimMatch", "sim found comparison ntuple", "event:match:chan:shit:fhit:tdiff:stime:ftime:qpeaks:qpeakf");
  // histograms for event cuts
  ntBase = new TNtuple("ntBase", " baseline ntuple ", "event:chan:base0:base1:fitMean:sigma:status"); // Fill(entry, ib, ave, sigma, fitStatus);;
  ntAdc = new TNtuple("ntAdc", " ADC ntuple ", "event:chan:sample:digi");
  ntTrig = new TNtuple("ntTrig", " trigger cut  ntuple ", "event:qsum9:qsum10:qsum11:qsum13:ratio910:ratio911:ratio1011:xternQ:yternQ:fails");
  ntNonTrig = new TNtuple("ntNonTrig", " non trigger ntuple ", "event:chan:qsum");
  hTriggerTime = new TH1D("TriggerTime", " ave of trigger Sipm times ", 1000, 0, 1000);
  hPreSumCut = new TH1D("PreSumCut", " pre trigger sum /nominal gain ", 100, 0, 2 * preSumCut);
  hCosmicCut = new TH1D("CosmicCut", " PMT sum /nominal gain", 1000, 0, 2. * totCosmicCut);
  hGammaCut = new TH1D("GammaCut", "gamma late sum chan 13 /nominal gain ", 1000, 0, 2. * lateGammaCut);
  hTrigSumNoCut = new TH1D("TrigSumNoCut", " before cut qsum9+qsum10+qsum11  in units nominal PE ", 160, 0, 40.);
  hTrigSumCut = new TH1D("TrigSumCut", " qsum9+qsum10+qsum11  in units nominal PE ", 160, 0, 40.);
  hTriangle = new TH2D("Triangle", "ytern vs xtern", 100, 0., 1., 100, 0., 1.);

  TString hName;
  hName.Form("TrigRatio%i-9-10", 0);
  hTrigSumCutRatio.push_back(new TH1D(hName, hName, 50, 0., 10.));
  hName.Form("TrigRatio%i-9-11", 1);
  hTrigSumCutRatio.push_back(new TH1D(hName, hName, 50, 0., 10.));
  hName.Form("TrigRatio%i-10-11", 2);
  hTrigSumCutRatio.push_back(new TH1D(hName, hName, 50, 0., 10.));

  // Fill(entry, ib, ave, sigma, fitStatus);;
  ntFailures = new TNtuple("ntFailures", " failures ntuple ", "event:chan:totHits:pass");
  ntThresholdAll = new TNtuple("ntThresholdAll", "ntThreshold no cuts", "event:chan:sample:ddigi");
  ntThresholdAdc = new TNtuple("ntThresholdAdc", "ntThresholdAdc", "event:chan:sampleLow:sampleHigh:maxBin:adcMax");
  ntThreshold = new TNtuple("ntThreshold", "ntThreshold passing cuts", "event:chan:sampleLow:ddigiLow:sampleHigh:ddigiHigh:maxBin:adcMax");
  ntHit = new TNtuple("ntHit", "hit ntuple", "event:flag:chan:time:peakTime:qpeak");
  ntChan = new TNtuple("ntChan", "channel ntuple", "trig:chan:ave:sigma:skew:base:peakmax:totSum:lateSum:negcrossings:pass");
  ntSpeYield = new TNtuple("ntSpeYield", "spe per sipm",
                           "event:spe0:spe1:spe2:spe3:spe4:spe5:spe6:spe7:spe8:spe9:spe10:spe11");
  ntSetTrigTime = new TNtuple("ntSetTrigTime", " trig time and val", "event:chan:time:val");
  ntTrigTime = new TNtuple("ntTrigTime", "trigger time check ntuple", "entry:chan:firstTime:time:adc:ftime:fadc");
  ntChanSum = new TNtuple("ntchansum", "channel ntuple", "sum0:sum1:sum2:sum3:sum4:sum5:sum6:sum7:sum8:sum9:sum10:sum11:sum12:pass");
  evCount = new TH1D("eventcount", "event count", CHANNELS, 0, CHANNELS);
  hEventPass = new TH1D("EventPass", " event failures", TOTALCODES, 0, TOTALCODES);
  hEventFail = new TH1D("EventFail", " event fail bit", FAILBITS, 0, FAILBITS); // first bin is pass
  hNoPeak = new TH1D("noPeak", "no peak events count by channel", CHANNELS, 0, CHANNELS);
  histHitCount = new TH1D("hitCount", "hit count by channel", CHANNELS, 0, CHANNELS);
  histQSum = new TH1D("histqsum", "qsum by channel", CHANNELS, 0, CHANNELS);
  // nn/histqpe = new th1d("histqpe", "qpe by channel", CHANNELS, 0, CHANNELS);
  histQPrompt = new TH1D("histqprompt", "qprompt by channel", CHANNELS, 0, CHANNELS);
  histQSum->Sumw2();
  histQPrompt->Sumw2();
  // hCosmicMult = new TH1D("CosmicMult", "CosmicMult", 10, 0, 10);

  //
  anaDir->cd();
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
  // hCountLateTime = new TH1D("CountLateTime ", htitle, 30, 0, 7500);
  // hCountLateTime->GetXaxis()->SetTitle("sample time");
  // hCountLateTime->Sumw2();
  hCountLateTimeQpeak = new TH2D("CountLateTimeQpeak", " sum qpeak vs time ", 750, 0, 7500, 80, 0, 20);
  hCountLateTimeQpeak->GetXaxis()->SetTitle("sample time");
  hCountLateTimeQpeak->GetYaxis()->SetTitle("qpeak [SPE]");

  printf("line1724 rawBr.size %lu \n", rawBr.size());
  double qpeakLimit;
  double qsumLimit;
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

    for (int ih = 0; ih < hEvGaus.size(); ++ih)
      hEvGaus[ih]->SetDirectory(nullptr);

    // for summary //
    qpeakLimit = 5. * nominalGain;
    qsumLimit = 5. * nominalQsumGain;

    bool trigger = ichan == 9 || ichan == 10 || ichan == 11;
    if (trigger)
    {
      qpeakLimit = 5. * nominalTrigGain;
      qsumLimit = 5. * nominalQsumTrigGain;
    }
    if (ichan == 12)
    {
      qpeakLimit = 5. * nominalPmtGain;
      qsumLimit = 5. * nominalQsumPmtGain;
    }
    int nbins = 700.;
    hTotSum.push_back(new TH1D(Form("TotPeakSumChan%i", i), Form("tot peak sum chan %i", i), nbins, 0, qpeakLimit));
    hPreSum.push_back(new TH1D(Form("PrePeakSumChan%i", i), Form("pre peak sum chan %i", i), nbins, 0, qpeakLimit));
    hTrigSum.push_back(new TH1D(Form("TrigPeakSumChan%i", i), Form("trig peak sum chan %i", i), nbins, 0, qpeakLimit));
    hLateSum.push_back(new TH1D(Form("LatePeakSumChan%i", i), Form("late peak sum chan %i", i), nbins, 0, qpeakLimit));
  }

  // one more for summed channel 13
  // hEvRawWave.push_back(new TH1D(Form("evRawWave%i", 13), Form("evRawWave%i", 13), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));

  threshValueHist = new TH2D("threshValueHist", " threshold crossings value channels by time  ", 7500, 0, 7500, 1000, 0, 100000);
  threshHist = new TH1D("threshHist", " threshold crossings trig channels ", 20, 0, 20);
  crossHist = new TH1D("crossHist", "  negative crossings non trigger channels", 100, 0, 100);
  sumDir->cd();

  for (unsigned i = 0; i < rawBr.size(); ++i)
  {
    unsigned ichan = i;
    qpeakLimit = 5. * nominalGain;
    qsumLimit = 5. * nominalQsumGain;

    bool trigger = ichan == 9 || ichan == 10 || ichan == 11;
    if (trigger)
    {
      qpeakLimit = 5. * nominalTrigGain;
      qsumLimit = 5. * nominalQsumTrigGain;
    }
    if (ichan == 12)
    {
      qpeakLimit = 5. * nominalPmtGain;
      qsumLimit = 5. * nominalQsumPmtGain;
    }

    hQPeak.push_back(new TH1D(Form("QPeakChan%i", ichan), Form("QPeakChan%i", ichan), 700, 0, qpeakLimit));
    hQSum.push_back(new TH1D(Form("QSumChan%i", ichan), Form("QSumChan%i", ichan), 1000, 0, qsumLimit));
    hQSpe.push_back(new TH1D(Form("QSpeChan%i", ichan), Form("QSpeChan%i", ichan), 9, 0, 9.));
    sumWave.push_back(new TH1D(Form("sumWave%i", ichan), Form("sumWave%i", ichan), rawBr[0]->rdigi.size(), 0, 2 * rawBr[0]->rdigi.size()));
    sumWaveA.push_back(new TH1D(Form("sumWaveAll%i", ichan), Form("sumWaveAll%i", ichan), rawBr[0]->rdigi.size(), 0, 2 * rawBr[0]->rdigi.size()));
    sumWaveB.push_back(new TH1D(Form("sumWaveBad%i", ichan), Form("sumWaveBad%i", ichan), rawBr[0]->rdigi.size(), 0, 2 * rawBr[0]->rdigi.size()));

    for (int ic = 0; ic < FAILBITS; ++ic)
    {
      sumWaveFail[ic].push_back(new TH1D(Form("sumWaveFail%iChan%i", failCode[ic], ichan), Form("sumWaveFail%iChan%i", failCode[ic], ichan), rawBr[0]->rdigi.size(), 0, 2 * rawBr[0]->rdigi.size()));
    }

    sumHitWave.push_back(new TH1D(Form("sumHitWave%i", ichan), Form("sumHitWave%i", ichan), rawBr[0]->rdigi.size(), 0, 2 * rawBr[0]->rdigi.size()));
    sumPeakWave.push_back(new TH1D(Form("sumPeakWave%i", ichan), Form("sumPeakWave%i", ichan), rawBr[0]->rdigi.size(), 0, 2 * rawBr[0]->rdigi.size()));
  }

  for (int ih = 0; ih < sumHitWave.size(); ++ih)
  {
    sumWave[ih]->GetXaxis()->SetTitle("time [ns]");
    sumWaveA[ih]->GetXaxis()->SetTitle("time [ns]");
    sumWaveB[ih]->GetXaxis()->SetTitle("time [ns]");
    sumHitWave[ih]->GetXaxis()->SetTitle("time [ns]");
    sumHitWave[ih]->GetXaxis()->SetTitle("time [ns]");
    sumPeakWave[ih]->GetXaxis()->SetTitle("time [ns]");
  }
  // SPE Shapes
  templateDir->cd();
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

  // maks simulation histogrms
  simDir->cd();
  hSimFoundTimeDiff = new TH1D("SimFoundTimeDiff", "sim - found time", 2000, -1000, 1000);
  for (unsigned i = 0; i < rawBr.size(); ++i)
  {
    unsigned ichan = i;
    hWaveHitFound.push_back(new TH1D(Form("waveHitFoundChan%i", ichan), Form("WaveHitFoundChan%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));
    hWaveHitMissed.push_back(new TH1D(Form("waveHitMissedChan%i", ichan), Form("WaveHitMissedChan%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));

    hWaveHitNoise.push_back(new TH1D(Form("waveHitNoiseChan%i", ichan), Form("WaveHitNoiseChan%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));
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
      // printf("line1990 event %lld passes\n", entry);
      hEventFail->SetBinContent(1, hEventFail->GetBinContent(1) + 1);
    }
    else
    {
      ++nfail;
      for (int ic = 0; ic < FAILBITS; ++ic)
      {
        if (passBit & failCode[ic])
        { // logical and
          hEventFail->SetBinContent(ic + 1, hEventFail->GetBinContent(ic + 1) + 1);
          // printf("line2049 event %lld fails bit %i\n", entry, ic);
        }
      }
      // set pass bit and fill tbrun
      for (int idet = 0; idet < tbrun->detList.size(); ++idet)
      {
        tbrun->detList[idet]->pass = passBit;
      }
    }
    hEventPass->SetBinContent(passBit, hEventPass->GetBinContent(passBit) + 1);
    // failure rate by bit

    // sum wave by failure code

    tbrun->fill();
    if (entry / 100 * 100 == entry)
    {
      printf("... entry %llu pass %u fail %u  failures by bit:\n", entry, npass, nfail);
      // printf(" FINISHED npass %u nfail %u output file  %s \n", npass, nfail, fout->GetName());
      printf(" line 2209 entry %i ( %i ) pass %i (%i) fail %i ( frac %0.3f )  \n",
             npass + nfail,
             int(hEventFail->GetEntries()),
             npass, int(hEventFail->GetBinContent(1)),
             nfail,
             double(nfail) / double(npass + nfail));

      for (int ibin = 1; ibin <= hEventFail->GetNbinsX(); ++ibin)
        printf(" bin %i content %.0f %s \n", ibin, hEventFail->GetBinContent(ibin), bitNames[ibin - 1].Data());
      hEventFail->Print("all");

      // if (npass > 0)
      if (0)
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

  // hEventPass->Print("all");
  printf("pass fractio/ns total = %.0f fail cosmic %i fail gamma %i \n", hEventPass->GetEntries(), failCosmic, failGamma);
  for (int ibin = 0; ibin < hEventPass->GetNbinsX(); ++ibin)
  { // inc/lude error on poisson probability
    double nbin = hEventPass->GetBinContent(ibin);
    double ntot = hEventPass->GetEntries();
    double prob = nbin / ntot;
    double perror = sqrt(prob * (1. - prob) / ntot);
    printf(" bin %i fail %.f frac %.3f +/- %.3f name %s \n", ibin, hEventPass->GetBinContent(ibin), prob, perror, codeNames[ibin].Data());
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
  printf(" fail bit frequency total = %.0f  pass %i fail %i fail cosmic %i fail gamma %i \n", hEventFail->GetEntries(), npass, nfail, failCosmic, failGamma);
  for (int ibin = 1; ibin <= hEventFail->GetNbinsX(); ++ibin)
  { // include error on poisson probability
    double nbin = hEventFail->GetBinContent(ibin);
    double ntot = hEventFail->GetEntries();
    double prob = nbin / ntot;
    double perror = sqrt(prob * (1. - prob) / ntot);
    printf(" bit %i fail %.f frac %.3f +/- %.3f  %s \n", ibin, hEventFail->GetBinContent(ibin), prob, perror, bitNames[ibin - 1].Data());
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
  // print out pulse finding stats
  if (isSim)
  {
    printf("pulse finding stats from simulation:\n");
    for (unsigned ichan = 0; ichan < hWaveHitFound.size() - 1; ++ichan)
    {
      double foundFraction = 1.;
      double tot = hWaveHitFound[ichan]->GetEntries() + hWaveHitMissed[ichan]->GetEntries();
      if (tot > 0)
        foundFraction = hWaveHitFound[ichan]->GetEntries() / tot;
      printf(" \t chan %i found %.0f missed %.0f fake %.0f  found fraction %.3f \n", ichan, hWaveHitFound[ichan]->GetEntries(),
             hWaveHitMissed[ichan]->GetEntries(), hWaveHitNoise[ichan]->GetEntries(), foundFraction);
    }
  }
  // hCosmicMult->Print("all");
  fout->Write();
  fout->Close();
  printf(" ***** FINISHED ****** %s entries %lld \n", fout->GetName(), nentries);
  return nentries;
}

anaCRun::anaCRun(TString theTag)
{

  failCode[0] = PASS;
  failCode[1] = BASEFAIL;
  failCode[2] = EARLYCUT;
  failCode[3] = FIRSTTIME;
  failCode[4] = COSMIC;
  failCode[5] = GAMMA;
  failCode[6] = TRIGFAIL;

  bitNames.resize(FAILBITS);
  bitNames[0] = TString("pass");
  bitNames[1] = TString("baseline");
  bitNames[2] = TString("earlycut");
  bitNames[3] = TString("firsttime");
  bitNames[4] = TString("cosmic");
  bitNames[5] = TString("gamma");
  bitNames[6] = TString("trigger");

  for (unsigned ic = 0; ic < TOTALCODES; ++ic)
    codeNames.push_back(TString("mixed"));
  codeNames[PASS] = TString("pass");
  codeNames[BASEFAIL] = TString("baseline");
  codeNames[TRIGFAIL] = TString("trigger");
  codeNames[EARLYCUT] = TString("earlycut");
  codeNames[FIRSTTIME] = TString("firsttime");
  codeNames[COSMIC] = TString("cosmic");
  codeNames[GAMMA] = TString("gamma");

  printf(" failure codes \n");
  for (unsigned ic = 0; ic < FAILBITS; ++ic)
    printf("bit %i hex %x %s\n", ic, failCode[ic], bitNames[ic].Data());

  printf("TRIGSUMCUT %f,%f nominalGain %f\n", trigRatioCutLow, trigRatioCutHigh, nominalGain);

  tag = theTag;
  // tbrun = new TBRun(tag);
  cout << " anaCRun::anaCRun instance of anaCRun gamma version 2  with tag= " << tag << " CHANNELS = " << CHANNELS - 1 << " diffStepSipm= " << diffStepSipm << " diffStepPmt= " << diffStepPmt << endl;

  rawBr.clear();

  for (int ichan = 0; ichan < CHANNELS; ++ichan)
  {
    TBRawEvent *rawEv = new TBRawEvent(ichan);
    rawEv->rdigi.resize(7500);
    rawEv->SetName(Form("rawChan%i", ichan));
    rawBr.push_back(rawEv);
  }
}