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
#include "modelAllFit.hh" // for geo ff function and distance levels

using namespace TMath;

// class to read and store gains
TReadGains *readGains;
TChain *RunTree;
TFile *fout;
bool isSimulation = false;
bool isLedRun = false;
double triggerSum;
TString currentFileName = TString("");
int currentFileNumber = -1;

std::vector<double> scaleSum;
std::vector<double> scalePeak;

std::vector<double> scaleLateSum;
std::vector<double> scalePreSum;

std::vector<double> earlyHitCount;
std::vector<double> lateHitCount;
double earlyHitCountFile = 0;
double lateHitCountFile = 0;

TString tag;
Long64_t totalEntries;
TNtuple *ntTrig;
TNtuple *ntGamma;
TNtuple *ntLateSum;
TNtuple *ntPreSum;
TNtuple *ntLateInt;
TNtuple *ntHitCount;
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
TH1D *hGammaPeakHit;
TH1D *hGammaPeakPass;
TH1D *hGammaCut;
TH1D *hCosmicCut;
TH1D *hQsumChannel;
TH1D *hQsumChannelEff;

std::vector<TH1D *> hLateSumChan;
std::vector<TH1D *> hLightCurve;
std::vector<TH1D *> hLightNorm;
std::vector<TH1D *> hLightEff;

std::vector<double> qsumGain; // read from class TReadGain
double aveGain;

// cut values
double cosmicCut = 100.; // value normlized to nominalPmtGain
double gammaCut = 10;    // was 140.; // was 150 normalized to nominalGain; //

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

double gainFunc(int ich)
{
  double gain = readGains->sipmPeakGain[ich];
  /* use same gains in simulation
  if (!isSimulation)
    return gain;
  if (ich < 9)
    gain = readGains->nominalGain;
  else if (ich > 8 && ich < 12)
    gain = readGains->nominalTrigGain;
  else
    gain = readGains->nominalPmtGain;
  */
  return gain;
}

void normalize(int ichan)
{
  double gain = gainFunc(ichan);

  // printf("at line174 %s normalize to %d\n", hSave->GetName(), totalPass);
  TH1D *hist = hLightCurve[ichan];
  TH1D *hSave = hLightNorm[ichan]; /// normalized to total pass and gain
  TH1D *hSave2 = hLightEff[ichan]; /// normalized to total pass and gain
  for (int ibin = 0; ibin < hist->GetNbinsX(); ++ibin)
  {
    double xbin = hist->GetBinContent(ibin);
    double ebin = sqrt(abs(hist->GetBinContent(ibin)));
    // divide by channel gain and totalPass
    {
      hSave->SetBinContent(ibin, xbin / double(totalPass) / gain);
      hSave->SetBinError(ibin, ebin / double(totalPass) / gain);
      hSave2->SetBinContent(ibin, xbin / double(totalPass) / gain);
      hSave2->SetBinError(ibin, ebin / double(totalPass) / gain);
    }
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

  if (entry == 0)
    readGains->printGains();

  triggerSum = 0;
  double triggerSumUn = 0;
  for (int i = 9; i < 12; ++i)
  {
    // triggerSum += detList[9 + i]->totSum * scale[i];
    triggerSum += detList[i]->totSum * scaleSum[i];
    triggerSumUn += detList[i]->totSum;
  }

  // printf("postAna .... %lld totSum %f\n,", entry, detList[9]->totSum);

  /* set trigfail if triggerSum>230 */
  if (triggerSum > 230)
    passBit |= TRIGFAIL;

  hGammaPeak->Fill(triggerSum);
  double qFraction[3];
  qFraction[0] = detList[9]->totSum * scaleSum[9] / triggerSum;
  qFraction[1] = detList[10]->totSum * scaleSum[10] / triggerSum;
  qFraction[2] = detList[11]->totSum * scaleSum[11] / triggerSum;

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
  double pmtLateSum = 0;
  if (!isSimulation)
    pmtLateSum = detList[12]->lateSum * scaleSum[12];
  else
    pmtLateSum = detList[12]->lateSum;
  hGammaCut->Fill(pmtLateSum);
  if (pmtLateSum > gammaCut)
    passBit |= GAMMA;

  // cosmic cut on summed SIPM
  // det 13 is sum of alll SIPMS
  double totSum13 = 0;
  for (int i = 0; i < 12; ++i)
    totSum13 += detList[i]->totSum * scaleSum[i];

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

  ntTrig->Fill(double(entry), pmtLateSum, totSum13, triggerSum, detList[9]->totSum, detList[10]->totSum, detList[11]->totSum,
               detList[9]->totSum * scaleSum[9], detList[10]->totSum * scaleSum[10], detList[11]->totSum * scaleSum[11], xternQ, yternQun, double(passBit));

  // fill passing gamma peak
  if (passBit == 0)
    hGammaPeakPass->Fill(triggerSum);

  for (unsigned i = 0; i < 9; ++i)
  {
    hQsumChannel->Fill(i + 1, detList[i]->totSum);
    hQsumChannelEff->Fill(i + 1, detList[i]->totSum * scaleSum[i] / effGeoFunc(i));
  }

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
  /* nominal gains have been applied in pulse finding step */
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
    if (passBit != 0 && !isLedRun)
      continue;

    ++totalPass;

    RunTree->GetEntry(entry);
    TString theFileName = TString(RunTree->GetCurrentFile()->GetName());

    if (theFileName != currentFileName)
    {
      printf("Processing file %d previous early %.0f late %.0f \n", currentFileNumber, earlyHitCountFile, lateHitCountFile);
      fflush(stdout);
      ++currentFileNumber;
      currentFileName = theFileName;

      if (currentFileNumber >= 0)
      {
        earlyHitCount.push_back(earlyHitCountFile);
        lateHitCount.push_back(lateHitCountFile);
      }
      earlyHitCountFile = 0;
      lateHitCountFile = 0;
    }

    // RunTree->GetListOfBranches()->ls();
    //   get branch pointers and save in detList
    TIter next(RunTree->GetListOfBranches());
    TBranchElement *aBranch = NULL;
    // loop over branches

    double eventTriggerHitQsum = 0;
    double photonSum[CHANNELS];
    double qsumSum[CHANNELS];
    double qsumLate[CHANNELS];

    for (int ich = 0; ich < CHANNELS; ++ich)
    {
      photonSum[ich] = 0;
      qsumSum[ich] = 0;
      qsumLate[ich] = 0;
    }
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
      // want to subtract off noise hits from preSum lateSum 3000-5500 ULong_t triggerStart = 730;
      // double scale = readGains->nominalQsumGain / aveGain;
      qsumLate[idet] = det->lateSum * scaleSum[idet]; // nominal gain applied in pulse finding step
      // printf("det %i hits %lu \n", idet, det->hits.size());
      //  check if passes eventCuts

      hLateSumChan[idet]->Fill(det->lateSum * scaleSum[idet]);
      ntPreSum->Fill(double(entry), double(idet), effGeoFunc(idet), det->preSum);
      ntLateSum->Fill(double(entry), double(idet), effGeoFunc(idet), det->lateSum);
      // loop over hits
      for (unsigned ihit = 0; ihit < det->hits.size(); ++ihit)
      {
        TDetHit thit = det->hits[ihit];
        // count early and late hits
        if (thit.firstBin < 600)
          earlyHitCountFile++;
        if (thit.firstBin >= 7500 - 600)
          lateHitCountFile++;
        // fill light curve
        hLightCurve[idet]->SetBinContent(thit.firstBin + 1, hLightCurve[idet]->GetBinContent(thit.firstBin + 1) + thit.qpeak / readGains->sipmPeakGain[idet]);
        /* fill gain histograms */
        photonSum[idet] += thit.qpeak / readGains->sipmPeakGain[idet];
        qsumSum[idet] += thit.qsum / readGains->sipmSumGain[idet];
        if (trig)
        {
          eventTriggerHitQsum += thit.qsum / readGains->sipmSumGain[idet];
        }
        // for ledData only look after 6000
        if (isLedRun && thit.firstBin < 6000)
          continue;
        hQPeak[idet]->Fill(thit.qpeak);
        hQSum[idet]->Fill(thit.qsum);
        // want to subtract off noise hits from preSum
        // if (idet > 8 && idet < 12)
        //  printf("... idet %i scale %f qsum %f eventTriggerHitQsum %f \n", idet, scale[idet], thit.qsum, eventTriggerHitQsum);
      } // end branch loop
      ntLateInt->Fill(double(entry), double(idet), qsumLate[0], qsumLate[1], qsumLate[2], qsumLate[3], qsumLate[4], qsumLate[5], qsumLate[6], qsumLate[7], qsumLate[8], qsumLate[9], qsumLate[10], qsumLate[11]);
    } // branch

    hGammaPeakHit->Fill(eventTriggerHitQsum);

    // ntGamma = new TNtuple("ntGamma", "gamma peak", "event:ph9:qsum9:ph10:qsum10:ph11:qsum11:peak");
    ntGamma->Fill(double(entry), photonSum[9], qsumSum[9], photonSum[10], qsumSum[10], photonSum[11], qsumSum[11], eventTriggerHitQsum, triggerSum);

  } // entry
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

  // geometric eff
  bool geoVersionOld = false;
  setDistanceLevels(geoVersionOld);
  for (unsigned i = 0; i < 12; ++i)
  {
    printf("chan %i distance %f geo eff %.2E\n", i, distanceLevel[getLevel(i)], effGeoFunc(i));
  }

  aveGain = 0;
  for (int i = 0; i < 9; ++i)
  {
    printf("chan %i gain %f \n", i, readGains->sipmSumGain[i]);
    aveGain += readGains->sipmSumGain[i];
  }
  aveGain /= double(9);
  printf("average gain %f \n", aveGain);

  /* get scale factor */
  // scale factor to new gain
  scaleSum.resize(readGains->sipmSumGain.size());
  scalePeak.resize(readGains->sipmSumGain.size());
  for (unsigned i = 0; i < readGains->sipmSumGain.size(); ++i)
  {
    scalePeak[i] = 1.0;
    scaleSum[i] = 1.0;
    if (!isSimulation)
    {
      scalePeak[i] = readGains->getNominalPeak(i) / readGains->sipmPeakGain[i];
      scaleSum[i] = readGains->getNominalSum(i) / readGains->sipmSumGain[i];
    }
  }

  // trigger info ntuple
  // ntTrig->Fill( pmtLateSum , totSum13 , triggerSum ,qFraction[0] ,  qFraction[1] , qFraction[2] , double(passBit) );
  ntTrig = new TNtuple("ntTrig", "trigger info", "event:pmtLateSum:totSum13:triggerSum:qun0:qun1:qun2:q0:q1:q2:xQ:yQ:passBit");
  ntLateInt = new TNtuple("ntLateInt", "late integral", "event:pmtLateSum:totSum13:triggerSum:qun0:qun1:qun2:q0:q1:q2:xQ:yQ:passBit");
  ntHitCount = new TNtuple("ntHitCount", "hit count", "file:early:late");
  ntGamma = new TNtuple("ntGamma", "gamma peak", "event:ph9:qsum9:ph10:qsum10:ph11:qsum11:eventSum:sum");

  ntLateSum = new TNtuple("ntLateSum", "late sum info", "event:chan:geo:lateSum");
  ntPreSum = new TNtuple("ntPreSum", "pre sum info", "event:chan:geo:preSum");
  ntLateInt = new TNtuple("ntLateInt", "late integral by channel", "event:chan:int0:int1:int2:int3:int4:int5:int6:int7:int8:int9:int10:int11");
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
  hGammaPeak->GetXaxis()->SetTitle("gamma peak (summed ADC)");
  hGammaPeakPass = new TH1D("GammaPeakPass", "gamma peak  pass triangle (photons)", 150, 0., 300.);
  hGammaPeakPass->GetXaxis()->SetTitle("gamma peak (summed ADC)");
  hGammaPeakHit = new TH1D("GammaPeakHit", "gamma peak (photons)", 150, 0., 300.);
  hGammaPeakHit->GetXaxis()->SetTitle("gamma peak (hit area qsum)");
  hQsumChannel = new TH1D("QsumChannel", "qsum channel (photons)", 9, 0., 9.);
  hQsumChannel->GetYaxis()->SetTitle("summed qsum [SPE]");
  hQsumChannel->GetXaxis()->SetTitle("channel");
  hQsumChannelEff = new TH1D("QsumChannelEff", "qsum channel (photons)", 9, 0., 9.);
  hQsumChannelEff->GetYaxis()->SetTitle("summed qsum [geo scaled]");
  hQsumChannelEff->GetXaxis()->SetTitle("channel");

  TDirectory *ledDir = fout->mkdir("ledDir");
  ledDir->cd();
  for (unsigned i = 0; i < CHANNELS; ++i)
  {
    hLateSumChan.push_back(new TH1D(Form("LateSumChan%i", i), Form("LateSumChan%i", i), 600, -10., 50.));
    hLateSumChan[hLateSumChan.size() - 1]->GetXaxis()->SetTitle("summed late photons [SPE]");
    hLateSumChan[hLateSumChan.size() - 1]->GetYaxis()->SetTitle("evemts");
  }

  fout->cd();

  for (unsigned i = 0; i < CHANNELS; ++i)
  {
    // normalized to SPE
    hLightCurve.push_back(new TH1D(Form("LightCurveChan%i", i), Form("LightCurveChan%i", i), MAXSAMPLES, 0, 2 * MAXSAMPLES));
    hLightCurve[hLightCurve.size() - 1]->GetXaxis()->SetTitle("time [ns]");
    hLightCurve[hLightCurve.size() - 1]->GetYaxis()->SetTitle("number of photons/2ns");

    hLightNorm.push_back(new TH1D(Form("LightNormChan%i", i), Form("LightNormChan%i", i), MAXSAMPLES, 0, 2 * MAXSAMPLES));
    hLightNorm[hLightNorm.size() - 1]->GetXaxis()->SetTitle("time [ns]");
    hLightNorm[hLightNorm.size() - 1]->GetYaxis()->SetTitle("normalized number of photons per event /2ns");

    hLightEff.push_back(new TH1D(Form("LightEffChan%i", i), Form("LightEffChan%i eff corrected ", i), MAXSAMPLES, 0, 2 * MAXSAMPLES));
    hLightEff[hLightEff.size() - 1]->GetXaxis()->SetTitle("time [ns]");
    hLightEff[hLightEff.size() - 1]->GetYaxis()->SetTitle("normalized number of photons per event /2ns");
  }
  /* make gain hisograms*/
  fout->cd("gainDir");
  double qpeakLimit;
  double qsumLimit;
  for (unsigned ichan = 0; ichan < CHANNELS; ++ichan)
  {
    qpeakLimit = 50. * (readGains->sipmPeakGain[ichan]);
    qsumLimit = 50. * (readGains->sipmSumGain[ichan]);
    hQPeak.push_back(new TH1D(Form("QPeakChan%i", ichan), Form("QPeakChan%i", ichan), 2000, 0, qpeakLimit));
    hQSum.push_back(new TH1D(Form("QSumChan%i", ichan), Form("QSumChan%i", ichan), 2000, 0, qsumLimit));
  }

  /*
   *  loop over events
   */
  earlyHitCountFile = 0;
  lateHitCountFile = 0;
  loop();

  earlyHitCount.push_back(earlyHitCountFile);
  lateHitCount.push_back(lateHitCountFile);

  for (int i = 0; i < earlyHitCount.size(); ++i)
  {
    ntHitCount->Fill(i, earlyHitCount[i], lateHitCount[i]);
  }

  for (unsigned i = 0; i < earlyHitCount.size(); ++i)
  {
    printf("File %d: Early hits = %.0f, Late hits = %.0f\n", i, earlyHitCount[i], lateHitCount[i]);
  }

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

  isLedRun = false;
  if (theStartTag.Contains("02_26_2026"))
    isLedRun = true;

  printf(" input args %i \n ", argc);
  for (int jarg = 1; jarg < argc; ++jarg)
    printf(" %i= %s ", jarg, argv[jarg]);
  printf("\n");
  if (isLedRun)
  {
    printf("\n\n*************** THIS IS LED RUN ***************\n\n");
  }

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
  if (theStartTag.Contains("btbSim"))
    isSimulation = true;
  if (nfiles == 0)
  {
    printf(" >>>> datatype no files found <<<<\n");
    exit(0);
  }
  if (isSimulation)
    printf("****** this is simulation data **** \n");

  /** set gains read gains from saved file */

  if (isSimulation) // use default gains
    readGains = new TReadGains(false);
  else // use gains from file
    readGains = new TReadGains();

  // store qsumGain[ib];
  for (unsigned ch = 0; ch < readGains->sipmSumGain.size(); ++ch)
    qsumGain.push_back(readGains->sipmSumGain[ch]);

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

  cout << " starting summary for   " << fileListName.size() << " on " << sdate << " writing to file " << fout->GetName() << endl;
  /* here we make the TCHain and them loop over it */

  post(tag);
  printf("end of job \n");
  // fout->ls();
}