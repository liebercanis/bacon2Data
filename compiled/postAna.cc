/*
       program to analyze RunTree chain from date tag this is second pass after pulse findiing anacRunGamma.cc has been run
           uses TTree RunTree making a chain from date tag xx_xx_yyyy
               M Gold Nov 12 2025

       ........modified impliment event cuts...... Dec 2025 *
*/
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
#include "TReadGains.hxx"

using namespace TMath;

// class to read and store gains
TReadGains *readGains;
TChain *RunTree;
TFile *fout;
TString tag;
TNtuple *ntTrig;
Long64_t totalEntries;
Long64_t maxEntry;
Long64_t totalPass;
std::vector<TString> fileListName;
std::vector<std::vector<double>> vecFail;
TBEventData *eventData;
TDirectory *sumDir;
TDirectory *anaDir; //
TDirectory *cutDir;
// directory for gain plots
TDirectory *gainDir;
// vectors for hist pointers
std::vector<TH1D *> hQPeak;
std::vector<TH1D *> hQSum;

TH1D *hEventPass;
TH1D *eventCount;
TH1D *hEventPassNew;
TH1D *hPassBitNew;
TH2D *hTriangleUn;
TH2D *hTriangle;
TH2D *hTriangleSecondUn;
TH2D *hTriangleSecond;
TH2D *hTrianglePass;
TH1D *hGammaPeak;
TH1D *hGammaPeakPass;
TH1D *hGammaCut;
TH1D *hCosmicCut;

std::vector<TH1D *> hLightCurve;
std::vector<TH1D *> hLightNorm;

std::vector<double> qsumGain; // read from class TReadGain

// cut values
double cosmicCut = 58.; // value normlized to nominalPmtGain
double gammaCut = 10;   // was 140.; // was 150 normalized to nominalGain; //

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
  TRIANGLE = 0x40, // 2^6
  TOTALCODES = 2 * TRIANGLE
};

enum
{
  FAILBITS = 8
};

std::vector<TString> bitNames;
std::vector<TString> codeNames;
std::vector<int> failCode;

//
enum
{
  CHANNELS = 14,
  NONSUMCHANNELS = CHANNELS - 1
};

enum
{
  MAXSAMPLES = 7500
};

std::vector<int> vTotal;
std::vector<int> vPass;

int nFiles;
TString theStartTag;
TString theEndTag;
std::string sdate;
TDatime dateTime;
time_t time0;
time_t time1;
struct tm tmStruct;
double ntotal;
double npass;
vector<int> filePass;
vector<int> fileTotal;

void normalize(int ichan)
{
  // printf("at line174 %s normalize to %d\n", hSave->GetName(), totalPass);
  TH1D *hist = hLightCurve[ichan];
  TH1D *hSave = hLightNorm[ichan];
  for (int ibin = 0; ibin < hist->GetNbinsX(); ++ibin)
  {
    double xbin = hist->GetBinContent(ibin);
    double ebin = sqrt(abs(hist->GetBinContent(ibin)));
    // divide by channel gain and totalPass
    hSave->SetBinContent(ibin, xbin / double(totalPass) / readGains->sipmPeakGain[ichan]);
    hSave->SetBinError(ibin, ebin / double(totalPass) / readGains->sipmPeakGain[ichan]);
  }
}

//// https://mathworld.wolfram.com/TernaryDiagram.html
void makeTernary(double a, double b, double c, double &x, double &y)
{
  double s = a + b + c;
  x = 0.5 * (a + 2. * b) / s;
  y = sqrt(3.) / 2. * a / s;
}

int passEventCuts(Long64_t entry)
{ // there is only 1 passbit but stored in all det. I admit is bad MGold
  RunTree->GetEntry(entry);
  // RunTree->GetListOfBranches()->ls();
  //   get branch pointers and save in detList
  TIter next(RunTree->GetListOfBranches());
  TBranchElement *aBranch = NULL;
  int passBit = 0;
  std::vector<TDet *> detList;
  while ((aBranch = (TBranchElement *)next()))
  {
    int idet = TString(TString(aBranch->GetName())(4, 2)).Atoi();
    if (TString(aBranch->GetName()) == TString("eventData"))
    { // skip this branch
      continue;
    }

    TDet *det = (TDet *)aBranch->GetObject();
    // collect TDet in list
    detList.push_back(det);
    // loop over hits
    if (idet == 9)
    {
      passBit = det->pass;
    }
  }
  // clear old bits
  passBit &= ~(TRIGFAIL); // triangle bit in anaCRunGamma.cc
  passBit &= ~(COSMIC);
  passBit &= ~(GAMMA);
  passBit &= ~(TRIANGLE);

  // if recalculating a cut, would do it here
  // for (int idet = 0; idet < detList.size(); ++idet)
  //{
  ///  printf("det %i totSum %f \n", idet, detList[idet]->totSum);
  //}

  /* be careful to rmove nominal gain used in anaCRunGamma */
  // scale factor to new gain
  double scale[3];
  for (int i = 0; i < 3; ++i)
    scale[i] = readGains->sipmSumGain[9 + i] / readGains->nominalQsumTrigGain;

  if (entry == 0)
    for (int i = 0; i < 3; ++i)
      printf("gain scale factor: channel %i  ratio gain/nominal Qsum %f  \n", i, scale[i]);

  double triggerSum = 0;
  double triggerSumUn = 0;
  for (int i = 0; i < 3; ++i)
  {
    triggerSum += detList[9 + i]->totSum * scale[i];
    triggerSumUn += detList[9 + i]->totSum;
  }

  hGammaPeak->Fill(triggerSum);
  double qFraction[3];
  qFraction[0] = detList[9]->totSum * scale[0] / triggerSum;
  qFraction[1] = detList[10]->totSum * scale[1] / triggerSum;
  qFraction[2] = detList[11]->totSum * scale[2] / triggerSum;

  double qFractionUn[3];
  qFractionUn[0] = detList[9]->totSum / triggerSumUn;
  qFractionUn[1] = detList[10]->totSum / triggerSumUn;
  qFractionUn[2] = detList[11]->totSum / triggerSumUn;

  double xternQ, yternQ;
  makeTernary(qFraction[0], qFraction[1], qFraction[2], xternQ, yternQ);
  hTriangle->Fill(xternQ, yternQ);

  double xternQun, yternQun;
  makeTernary(qFractionUn[0], qFractionUn[1], qFractionUn[2], xternQun, yternQun);
  hTriangleUn->Fill(xternQun, yternQun);

  /* look at second peak triangle*/

  /* look at second peak triangle*/
  if (triggerSumUn > 50.)
    hTriangleSecondUn->Fill(xternQun, yternQun);

  if (triggerSum > 50.)
    hTriangleSecond->Fill(xternQ, yternQ);

  /******   triangle cut ** try a cut like TUM */
  bool passTriangle = true;
  for (int itr = 0; itr < 3; ++itr)
    if (qFraction[itr] < 0.2 || qFraction[itr] > 0.8)
      passTriangle = false;

  // plot passing
  if (passTriangle)
    hTrianglePass->Fill(xternQ, yternQ);

  // set triangle bit
  if (!passTriangle)
    passBit |= TRIANGLE;

  // gamma cut
  double pmtLateSum = detList[12]->lateSum * readGains->sipmSumGain[12] / readGains->nominalQsumPmtGain;
  hGammaCut->Fill(pmtLateSum);
  if (pmtLateSum > gammaCut)
    passBit |= GAMMA;

  // cosmic cut on summed SIPM
  // det 13 is sum of alll SIPMS
  double totSum13 = 0;
  for (int i = 0; i < 9; ++i)
    totSum13 += detList[i]->totSum * readGains->sipmSumGain[i] / readGains->nominalQsumGain;

  for (int i = 9; i < 12; ++i)
    totSum13 += detList[i]->totSum * readGains->sipmSumGain[i] / readGains->nominalQsumTrigGain;

  hCosmicCut->Fill(totSum13);
  if (totSum13 > cosmicCut)
    passBit |= COSMIC;

  // loop over fail bits
  for (int ic = 0; ic < FAILBITS; ++ic)
  {
    if (passBit & failCode[ic])
    {
      // printf("event %lld det %i bin %s pass %i \n", entry, idet, bitNames[ic].Data(), passBit);
      hPassBitNew->SetBinContent(ic, hPassBitNew->GetBinContent(ic) + 1);
    }
  }

  ntTrig->Fill(double(entry), pmtLateSum, totSum13, triggerSum, qFractionUn[0], qFractionUn[1], qFractionUn[2], qFraction[0], qFraction[1], qFraction[2], xternQ, yternQun, double(passBit));

  // fill passing gamma peak
  if (passBit == 0)
    hGammaPeakPass->Fill(triggerSum);

  return passBit;
}

string currentDate()
{
  time_t rawtime;
  struct tm *timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  char output[30];
  strftime(output, 30, "%Y-%m-%d-%H-%M", timeinfo);
  return string(output);
}

// get all pointers we need
bool getPointers(TFile *f)
{
  // printf("line176 getPointers file %s\n", f->GetName());
  bool isGoodFile = true;
  if (!f)
  {
    isGoodFile = false;
  }
  TString name(f->GetName());
  if (f->IsZombie())
  {
    cout << "line922 skipping zombie " << name << endl;
    isGoodFile = false;
  }

  TTree *RunTree = NULL;
  f->GetObject("RunTree", RunTree);
  if (RunTree == NULL)
  {
    cout << "line1215 skipping BAD file no RunTree" << name << endl;
    isGoodFile = false;
  }
  else
  {
    totalEntries += RunTree->GetEntries();
    maxEntry = totalEntries;
    printf("\t\t file %s has %lld RunTree entries total %lld \n", f->GetName(), RunTree->GetEntries(), totalEntries);
  }
  // if (isGoodFile)
  //   printf("good 1 \n");

  sumDir = nullptr;
  f->GetObject("sumDir", sumDir);
  if (sumDir == NULL)
  {
    cout << "line203 skipping BAD file no sumDir" << name << endl;
    isGoodFile = false;
  }

  anaDir = nullptr;
  f->GetObject("anaDir", anaDir);
  if (sumDir == NULL)
  {
    cout << "line214 skipping BAD file no anaDir" << name << endl;
    isGoodFile = false;
  }

  eventCount = nullptr;
  f->GetObject("eventcount", eventCount);
  if (!eventCount)
  {
    cout << "line223 skipping BAD file no eventcount " << name << endl;
    isGoodFile = false;
  }

  TH1D *hEventPass = nullptr;
  f->GetObject("EventPass", hEventPass);
  if (!hEventPass)
  {
    cout << "line1230 skipping BAD file no EventPass " << name << endl;
    isGoodFile = false;
  }

  TH1D *hGammaPeakFile = nullptr;
  f->GetObject("GammaPeak", hGammaPeakFile);
  if (!hGammaPeakFile)
  {
    cout << "line1230 no GammaPeak " << name << endl;
  }

  TH1D *hGammaPeakCut = nullptr;
  f->GetObject("GammaPeakCut", hGammaPeakCut);
  if (!hGammaPeakCut)
  {
    cout << "line1230 no GammaPeakCut " << name << endl;
  }

  eventData = new TBEventData();
  // RunTree->GetListOfBranches()->ls();
  if (RunTree)
    RunTree->SetBranchAddress("eventData", &eventData);
  if (!eventData)
    isGoodFile = false;

  // ***** not fatal if missing *****
  TH1D *hTrigSumNoCut = nullptr;
  f->GetObject("TrigSumNoCut", hTrigSumNoCut);
  TH1D *hTrigSumCut = nullptr;
  f->GetObject("TrigSumCut", hTrigSumCut);
  //
  return isGoodFile;
}

// count subruns and channels
unsigned long countFiles()
{
  totalEntries = 0;
  TString dirName = TString("caenData");
  TString dirNameSlash = TString("caenData/");
  cout << " count files in dir " << dirName << endl;
  TSystemDirectory dir(dirName, dirName); // TSystemDirectory
  TList *files = dir.GetListOfFiles();
  // print list
  TIter next(files); // Create an iterator for the TList
  TSystemFile *file;
  while ((file = (TSystemFile *)next()))
  {
    string name = string(file->GetName());
    TString tname = TString(name.c_str());
    string exten = name.substr(name.find_last_of(".") + 1);
    if (exten != string("root"))
      continue;
    // see if file is correct time
    int month = TString(tname(tname.Last('n') + 2, 2)).Atoi();
    int day = TString(tname(tname.Last('n') + 5, 2)).Atoi();
    int year = TString(tname(tname.Last('n') + 8, 4)).Atoi();
    tmStruct.tm_year = year;
    tmStruct.tm_mon = month;
    tmStruct.tm_mday = day;
    time_t fileTime = mktime(&tmStruct);
    const auto diff0 = std::difftime(fileTime, time0);
    const auto diff1 = std::difftime(fileTime, time1);
    // printf("line914 info : file %s time %s", tname.Data(), asctime(gmtime(&fileTime)));
    // cout << " \t ..... " << diff0 << " " << diff1 << endl;
    bool timetest = diff0 >= 0 && diff1 <= 0;
    if (!timetest)
    {
      // cout << "line1148 skip out of time file " << name << endl;
      continue;
    }
    TString fullName = dirNameSlash + TString(name.c_str());
    TFile *f = new TFile(fullName, "READONLY");
    /* not really needed but checks for file contents */
    if (getPointers(f))
      fileListName.push_back(TString(name.c_str()));
    f->Close();
  }
  return fileListName.size();
}

TDatime getTime(int ifile, Long64_t ievent = 0)
{
  TDatime datime;
  if (RunTree)
  {
    RunTree->GetEntry(ievent);
    printf(" file %i event %lli FileYear = %u , FileMonth = %u , FileDay = %u , FileHour = %u , FileMin = %u , FileSec = %u \n", ifile, ievent, eventData->year + 1900, eventData->mon + 1, eventData->day, eventData->hour, eventData->min, eventData->sec);
    datime.Set(eventData->year, eventData->mon + 1, eventData->day, eventData->hour, eventData->min, eventData->sec);
  }
  return datime;
}

void setTime(TString startTag, TString endTag)
{
  theStartTag = startTag;
  theEndTag = endTag;
  tag = theStartTag + TString("-") + theEndTag;
  int month0 = TString(theStartTag(0, 2)).Atoi();
  int day0 = TString(theStartTag(3, 2)).Atoi();
  int year0 = TString(theStartTag(6, 4)).Atoi();
  int month1 = TString(theEndTag(0, 2)).Atoi();
  int day1 = TString(theEndTag(3, 2)).Atoi();
  int year1 = TString(theEndTag(6, 4)).Atoi();
  // printf(" start %i %i %i ene %i %i %i  \n",month0,day0,year0,month1,day1,year1);

  /* fill in values for 2019-08-22 23:22:26 */
  tmStruct.tm_year = year0;
  tmStruct.tm_mon = month0;
  tmStruct.tm_mday = day0;
  time0 = mktime(&tmStruct);
  printf("set start %s\n", asctime(gmtime(&time0)));
  tmStruct.tm_year = year1;
  tmStruct.tm_mon = month1;
  tmStruct.tm_mday = day1;
  time1 = mktime(&tmStruct);
  printf("set end %s\n", asctime(gmtime(&time1)));
}

void loop()
{
  totalPass = 0;
  printf(" start of entry loop maxEntry=%lld\n", maxEntry);
  // loop over entries
  for (Long64_t entry = 0; entry < maxEntry; ++entry)
  {
    if (entry / 10000 * 10000 == entry)
    {
      printf("line330 .....loop entry %lld \n", entry);
      // force printing to log file
      fflush(stdout);
    }
    int passBit = passEventCuts(entry);
    hEventPassNew->SetBinContent(passBit, hEventPassNew->GetBinContent(passBit) + 1);
    if (passBit != 0)
      continue;

    ++totalPass;

    RunTree->GetEntry(entry);
    // RunTree->GetListOfBranches()->ls();
    //   get branch pointers and save in detList
    TIter next(RunTree->GetListOfBranches());
    TBranchElement *aBranch = NULL;
    // loop over branches

    while ((aBranch = (TBranchElement *)next()))
    {
      // skip eventData branch
      if (TString(aBranch->GetName()) == TString("eventData"))
      { // skip this branch
        continue;
      }
      int idet = TString(TString(aBranch->GetName())(4, 2)).Atoi();
      bool trig = false; // define trigger sipms
      if (idet == 9 || idet == 10 || idet == 11)
        trig = true;

      /* the branch is class TDet so cast it as such */
      TDet *det = (TDet *)aBranch->GetObject();
      // printf("det %i hits %lu \n", idet, det->hits.size());
      //  check if passes eventCuts

      // loop over hits
      for (unsigned ihit = 0; ihit < det->hits.size(); ++ihit)
      {
        TDetHit thit = det->hits[ihit];
        // fill light curve
        hLightCurve[idet]->SetBinContent(thit.firstBin + 1, hLightCurve[idet]->GetBinContent(thit.firstBin + 1) + thit.qpeak);
        /* fill gain histograms */
        hQPeak[idet]->Fill(thit.qpeak);
        hQSum[idet]->Fill(thit.qsum);
      } // end branch loop
    }
  }
}

/* build the TChain and call loop */
void post(TString tag)
{
  vecFail.resize(FAILBITS);
  /*gains-2024-02-01-17-06.root*/

  gStyle->SetOptStat(1001101);
  /* get RunTree */
  RunTree = new TChain("RunTree");
  //** add files  */
  for (unsigned ifile = 0; ifile < fileListName.size(); ++ifile)
  {
    TString fullName = TString("caenData/") + fileListName[ifile];
    printf("RunTree add file %s \n", fullName.Data());
    RunTree->Add(fullName);
  }

  if (!RunTree)
    return;
  printf("files in chain:\n");
  RunTree->GetListOfFiles()->Print();
  Long64_t ntriggers = RunTree->GetEntries();
  printf(" in post: tag %s total triggers in this chain %lld \n", tag.Data(), ntriggers);
  // RunTree->GetListOfBranches()->ls();

  // trigger info ntuple
  // ntTrig->Fill( pmtLateSum , totSum13 , triggerSum ,qFraction[0] ,  qFraction[1] , qFraction[2] , double(passBit) );
  ntTrig = new TNtuple("ntTrig", "trigger info", "event:pmtLateSum:totSum13:triggerSum:qFracUn0:qFracUn1:qFracUn2:qFraction0:qFraction1:qFraction2:xQ:yQ:passBit");

  // make histograms
  hPassBitNew = new TH1D("PassBitNew", "pass bit", FAILBITS, 0, FAILBITS);
  hEventPassNew = new TH1D("EventPassNew", " remade event failures", TOTALCODES, 0, TOTALCODES);
  hGammaCut = new TH1D("GammaCut", "gamma pmt lateSum/nominal gain ", 4000, 0, 5. * gammaCut);
  hCosmicCut = new TH1D("CosmicCut", " cosmic qsum13/nominal gain", 4000, 0, 5. * cosmicCut);
  hTriangleUn = new TH2D("TriangleUn", "ytern vs xtern unscaled", 100, 0., 1., 100, 0., 1.);
  hTriangle = new TH2D("Triangle", "ytern vs xtern", 100, 0., 1., 100, 0., 1.);
  hTriangleSecondUn = new TH2D("TriangleSecondUn", "ytern vs xtern in second gamma peak", 100, 0., 1., 100, 0., 1.);
  hTriangleSecond = new TH2D("TriangleSecond", "ytern vs xtern in second gamma peak gains", 100, 0., 1., 100, 0., 1.);
  hTrianglePass = new TH2D("TrianglePass", "ytern vs xtern", 100, 0., 1., 100, 0., 1.);
  hGammaPeak = new TH1D("GammaPeak", "gamma peak (photons)", 150, 0., 300.);
  hGammaPeakPass = new TH1D("GammaPeakPass", "gamma peak  pass triangle (photons)", 150, 0., 300.);

  for (unsigned i = 0; i < CHANNELS; ++i)
  {
    // normalized to SPE
    hLightCurve.push_back(new TH1D(Form("LightCurveChan%i", i), Form("LightCurveChan%i", i), MAXSAMPLES, 0, 2 * MAXSAMPLES));
    hLightCurve[hLightCurve.size() - 1]->GetXaxis()->SetTitle("time [ns]");
    hLightCurve[hLightCurve.size() - 1]->GetYaxis()->SetTitle("number of photons/2ns");

    hLightNorm.push_back(new TH1D(Form("LightNormChan%i", i), Form("LightNormChan%i", i), MAXSAMPLES, 0, 2 * MAXSAMPLES));
    hLightNorm[hLightNorm.size() - 1]->GetXaxis()->SetTitle("time [ns]");
    hLightNorm[hLightNorm.size() - 1]->GetYaxis()->SetTitle("normalized number of photons per event /2ns");
  }
  /* make gain hisograms*/
  fout->cd("gainDir");
  double qpeakLimit;
  double qsumLimit;
  for (unsigned ichan = 0; ichan < CHANNELS; ++ichan)
  {
    qpeakLimit = 5. * readGains->nominalGain;
    qsumLimit = 5. * readGains->nominalQsumGain;

    bool trigger = ichan == 9 || ichan == 10 || ichan == 11;
    if (trigger)
    {
      qpeakLimit = 5. * readGains->nominalTrigGain;
      qsumLimit = 5. * readGains->nominalQsumTrigGain;
    }
    if (ichan == 12)
    {
      qpeakLimit = 5. * readGains->nominalPmtGain;
      qsumLimit = 5. * readGains->nominalQsumPmtGain;
    }

    hQPeak.push_back(new TH1D(Form("QPeakChan%i", ichan), Form("QPeakChan%i", ichan), 2000, 0, qpeakLimit));
    hQSum.push_back(new TH1D(Form("QSumChan%i", ichan), Form("QSumChan%i", ichan), 2000, 0, qsumLimit));
  }

  /*
   *  loop over events
   */
  loop();

  printf("total %llu pass %llu \n", maxEntry, totalPass);
  // hEventPassNew->Print("all");
  printf("pass fractions total = %.0f  \n", hEventPassNew->GetEntries());
  for (int ibin = 0; ibin < hEventPassNew->GetNbinsX(); ++ibin)
  { // inc/lude error on poisson probability
    double nbin = hEventPassNew->GetBinContent(ibin);
    double ntot = hEventPassNew->GetEntries();
    double prob = nbin / ntot;
    double perror = sqrt(prob * (1. - prob) / ntot);
    if (nbin > 0)
      printf(" bin %i fail %.f frac %.3f +/- %.3f name %s \n", ibin, hEventPassNew->GetBinContent(ibin), prob, perror, codeNames[ibin].Data());
  }

  // do not normilzed summed chan 13
  for (int ich = 0; ich < hLightCurve.size() - 1; ++ich)
    normalize(ich);

  hPassBitNew->Print("all");
  // loop over fail bits
  printf("summary of bit failures %llu pass %llu \n", maxEntry, totalPass);
  for (int ic = 0; ic < FAILBITS; ++ic)
  {
    double prob = hPassBitNew->GetBinContent(ic) / double(maxEntry);
    double perror = sqrt(prob * (1. - prob)) / double(maxEntry);
    printf("bit %i %s val %.0f frac %.3f +/- %.3f \n", ic, bitNames[ic].Data(), hPassBitNew->GetBinContent(ic), prob, perror);
  }

  // fout->ls();
  fout->Write();
  // fout->ls();
}

int main(int argc, char *argv[])
{

  fout = nullptr;
  cout << "executing " << argv[0] << " post hit finding analysis  " << endl;
  printf(" usage:  start date string <stag> end date string <etag> max entries <default all> \n ");
  if (argc < 2)
  {
    printf("require file date start string <stag> args.\n  exit \n");
    exit(0);
  }

  /* read gains from saved file */
  readGains = new TReadGains();

  // store qsumGain[ib];
  for (unsigned ch = 0; ch < readGains->sipmSumGain.size(); ++ch)
    qsumGain.push_back(readGains->sipmSumGain[ch]);

  /* for failure bits */
  failCode.resize(FAILBITS);
  failCode[0] = PASS;
  failCode[1] = BASEFAIL;
  failCode[2] = EARLYCUT;
  failCode[3] = FIRSTTIME;
  failCode[4] = COSMIC;
  failCode[5] = GAMMA;
  failCode[6] = TRIGFAIL;
  failCode[7] = TRIANGLE;

  vecFail.resize(FAILBITS);
  bitNames.resize(FAILBITS);
  bitNames[0] = TString("Pass");
  bitNames[1] = TString("Baseline");
  bitNames[2] = TString("Earlycut");
  bitNames[3] = TString("Firsttime");
  bitNames[4] = TString("Cosmic");
  bitNames[5] = TString("Gamma");
  bitNames[6] = TString("Trigger");
  bitNames[7] = TString("Triangle");

  codeNames.resize(TOTALCODES);
  // build trigger bit pattern names
  for (int ic = 0; ic < TOTALCODES; ++ic)
  {
    for (int ibit = 0; ibit < FAILBITS; ++ibit)
    {
      if (ic & failCode[ibit])
        codeNames[ic] += bitNames[ibit];
    }
  }

  printf("failure codes: \n");
  for (int icode = 0; icode < 8; ++icode)
    printf("bit %i hex value %i name %s \n", icode, failCode[icode], bitNames[icode].Data());

  TString dirName = TString("caenData");
  TString dirNameSlash = TString("caenData/");
  theStartTag = TString(argv[1]);

  printf(" input args %i \n ", argc);
  for (int jarg = 1; jarg < argc; ++jarg)
    printf(" %i= %s ", jarg, argv[jarg]);
  printf("\n");

  theEndTag = TString(argv[1]);
  if (argc == 3)
  {
    theEndTag = TString(argv[2]);
  }

  // convert string to dates
  setTime(theStartTag, theEndTag);

  dateTime = TDatime(2023, 3, 9, 22, 0, 0);

  /* count files between dates */
  unsigned nfiles = countFiles();
  printf("count files from %s to %s total files  %ld \n", theStartTag.Data(), theEndTag.Data(), fileListName.size());
  if (nfiles == 0)
  {
    printf(" >>>> datatype no files found <<<<\n");
    exit(0);
  }

  if (argc > 3)
  {
    maxEntry = atoi(argv[3]);
  }

  printf(" >>>>> analyze %u  files from %s to %s tag %s totalEntries %lld maxEntry %lld <<<<<\n", nfiles, theStartTag.Data(), theEndTag.Data(), tag.Data(), totalEntries, maxEntry);

  sdate = currentDate();
  tag = theStartTag + TString("-") + theEndTag;
  TString sentries;
  sentries.Form("-%llu", maxEntry);

  fout = new TFile(TString("post-") + tag + sentries + TString(".root"), "recreate");
  gainDir = fout->mkdir("gainDir");
  // pick up first pass hEventCount now that fout is open
  for (int i = 0; i < fileListName.size(); ++i)
  {
    if (!fout)
      fout = new TFile(TString("post-") + tag + sentries + TString(".root"), "update");
    cout << i << "  " << fileListName[i] << endl;
    TString dirNameSlash = TString("caenData/");
    TString fullName = dirNameSlash + fileListName[i];
    /* open file in list */
    TFile *f = new TFile(fullName, "READONLY");
    TH1D *hEventPass = nullptr;
    f->GetObject("EventPass", hEventPass);
    if (hEventPass && fout)
    {
      TH1D *hEventPassFile = (TH1D *)hEventPass->Clone(Form("EventPassFile%i", i));
      // cout << "line589 append " << hEventPassFile->GetName() << endl;
      fout->Add(hEventPassFile);
      fout->Write();
      /* have to close output file to save*/
      fout->Close();
      delete fout;
      fout = nullptr;
      f->Close();
    }
  }

  if (!fout)
    fout = new TFile(TString("post-") + tag + sentries + TString(".root"), "update");
  printf("after added evenCount \n");
  fout->ls();

  cout << " starting summary for   " << fileListName.size() << " on " << sdate << " writing to file " << fout->GetName() << endl;
  /* here we make the TCHain and them loop over it */

  post(tag);
}