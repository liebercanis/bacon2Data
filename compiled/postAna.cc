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
#include <ctime>     // for time() and ctime()
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
#include "failCodes.hh"   // for fail bit names and codes

using namespace TMath;

// class to read and store gains
TReadGains *readGains;
TChain *RunTree;
TFile *fout;
bool isSimulation = false;
bool isLedRun = false;
double triggerSum;
TString currentFileName = TString("");
int currentFileNumber = 0;
int nominalTrigger = 729;
int theMaximumBin = 0;

std::vector<double> scaleSum;
std::vector<double> scalePeak;

std::vector<double> scaleLateSum;
std::vector<double> scalePreSum;

std::vector<double> earlyHitCountFile;
std::vector<double> lateHitCountFile;
int hitCountNev = 0;

TString tag;
Long64_t totalEntries;
TNtuple *ntTrig;
TNtuple *ntTrigChan;
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
TDirectory *crossDir;

// vectors for hist pointers
TH1D *hTrigTimeEvent;
std::vector<int> vTrigTime;
std::vector<TH1D *> hTrigTime;
std::vector<TH1D *> hQPeak;
std::vector<TH1D *> hQSum;
std::vector<TH1D *> hNextHitTime;
std::vector<TH1D *> hNextHitTimeOther;

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
TH1D *hOverlapEvent;
TH1D *hCosmicCut;
TH1D *hQsumChannel;
TH1D *hQsumChannelEff;

std::vector<TH1D *> hMulti;
std::vector<TH1D *> hLateSumChan;
std::vector<TH1D *> hLightCurve;
std::vector<TH1D *> hLightNorm;
std::vector<TH1D *> hLightEff;

std::vector<double> qsumGain; // read from class TReadGain
double aveGain;

// cut values based on 04_16_2026 revised June 9 2026
double overlapEventCut = 70.; // was 100.;           //         // value normlized to nominalPmtGain
double cosmicEventCut = 10.;  // 10; // was 140.; // was 150 normalized to nominalGain; //

// pass bit failures hex
// COSMIC AND GAMMA now refined as cosmicEvent and overlapEvent

enum
{
  TRIGGERSAMPLE = 695
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

bool fileExists(TString theFile)
{
  bool exists = false;
  FILE *aFile;
  aFile = fopen(theFile.Data(), "r");
  if (aFile)
  {
    fclose(aFile);
    exists = true;
  }
  return exists;
}

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
      // already normalized to gain line 653 , so just divide by total pass
      hSave->SetBinContent(ibin, xbin / double(totalPass));
      hSave->SetBinError(ibin, ebin / double(totalPass));
      hSave2->SetBinContent(ibin, xbin / double(totalPass));
      hSave2->SetBinError(ibin, ebin / double(totalPass));
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

  // total events in undeflow
  hPassBitNew->SetBinContent(0, hPassBitNew->GetBinContent(0) + 1);

  for (int i = 9; i < 12; ++i)
  {
    // triggerSum += detList[9 + i]->totSum * scale[i];
    triggerSum += detList[i]->totSum * scaleSum[i];
    triggerSumUn += detList[i]->totSum;
  }

  // printf("postAna .... %lld totSum9 %f triggerSum %f\n,", entry, detList[9]->totSum, triggerSum);

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

  // cosmic cut on PMT late sum fixed to totSum jun 10 2026
  double pmtTotSum = 0;
  if (!isSimulation)
    pmtTotSum = detList[12]->totSum * scaleSum[12];
  else
    pmtTotSum = detList[12]->totSum;
  hCosmicCut->Fill(pmtTotSum);
  if (pmtTotSum > cosmicEventCut)
    passBit |= COSMIC;

  // det 13 is sum of alll SIPMS overlapEvent cut on this
  double totSum13 = 0;
  for (int i = 0; i < 12; ++i)
    totSum13 += detList[i]->totSum * scaleSum[i];

  hOverlapEvent->Fill(totSum13);
  if (totSum13 > overlapEventCut)
    passBit |= GAMMA;

  ntTrig->Fill(double(entry), pmtTotSum, totSum13, triggerSum, detList[9]->totSum, detList[10]->totSum, detList[11]->totSum, xternQ, yternQun, double(passBit));

  // fill passing gamma peak
  if (passBit == 0)
    hGammaPeakPass->Fill(triggerSum);

  for (unsigned i = 0; i < 9; ++i)
  {
    hQsumChannel->Fill(i + 1, detList[i]->totSum);
    hQsumChannelEff->Fill(i + 1, detList[i]->totSum * scaleSum[i] / effGeoFunc(i));
  }

  bool multCut = false;
  // finally a SPE mult cut on not trigger SIPM
  for (unsigned idet = 0; idet < CHANNELS; ++idet)
  {
    // hit loop
    for (unsigned ihit = 0; ihit < detList[idet]->hits.size(); ++ihit)
    {
      TDetHit thit = detList[idet]->hits[ihit];
      hMulti[idet]->Fill(thit.qpeak / readGains->sipmPeakGain[idet]);
      // cut on 1.5 SPE for filling light curve should I remove entire event?
      if (idet < 9 && thit.qpeak / readGains->sipmPeakGain[idet] > 1.5)
        multCut = true;
    }
  }

  if (multCut)
  {
    passBit |= MULTI;
    // printf("MESSAGE event %lld fails multi pass = %i %i \n", entry, passBit, passBit & MULTI);
  }

  // trig time cut
  // reset trig time histos
  for (unsigned idet = 0; idet < hTrigTime.size(); ++idet)
  {
    hTrigTime[idet]->Reset("ICESM");
  }

  // for each event loop over trigger dets to get start time
  for (unsigned idetNumber = 0; idetNumber < 12; ++idetNumber)
  {
    // idetNumber hit loop
    for (unsigned ihit = 0; ihit < detList[idetNumber]->hits.size(); ++ihit)
    {
      TDetHit iDetHit = detList[idetNumber]->hits[ihit];
      hTrigTime[idetNumber]->Fill(iDetHit.firstBin);
    }
  }

  // start time is earliest trigger time of the 3 trigger SIPMS
  vTrigTime.clear();
  for (unsigned idetNumber = 9; idetNumber < 12; ++idetNumber)
    vTrigTime.push_back(hTrigTime[idetNumber]->GetXaxis()->GetBinCenter(hTrigTime[idetNumber]->GetMaximumBin()));

  // printf("line 411 event %lld  9 max bin: %i 10 max bin: %i 11 max bin: %i\n", entry, vTrigTime[0], vTrigTime[1], vTrigTime[2]);
  std::sort(vTrigTime.begin(), vTrigTime.end());
  theMaximumBin = vTrigTime[0];
  if (theMaximumBin < 700)
    theMaximumBin = nominalTrigger;
  hTrigTimeEvent->Fill(vTrigTime[0]);
  if (vTrigTime[0] < 660 || vTrigTime[0] > 740)
  {
    passBit |= TRIGTIME;
    // printf("MESSAGEline 434 event %lld fails trigtime %i pass = %i %i \n", entry, vTrigTime[0], passBit, passBit & TRIGTIME);
  }

  // gainDir->cd();
  if (entry < 100)
    for (unsigned idetNumber = 9; idetNumber < 12; ++idetNumber)
      hTrigTime[idetNumber]->Clone(TString::Format("hTrigTimeDet%i_event%lld", idetNumber, entry));

  // loop over fail bits
  for (int ic = 0; ic < FAILBITS; ++ic)
  {
    if (passBit == 0 && ic == 0)
    {
      hPassBitNew->SetBinContent(ic + 1, hPassBitNew->GetBinContent(ic + 1) + 1);
      // printf("event %lld bin %s passbit %i fail code %i pass  %f \n", entry, bitNames[ic].Data(), passBit, failCode[ic], hPassBitNew->GetBinContent(ic + 1));
    }
    else if (passBit & failCode[ic])
    {
      // printf("MESSAGE line 415 event %lld fails pass = %i bit %i name %s \n", entry, passBit, ic, bitNames[ic + 1].Data());
      hPassBitNew->SetBinContent(ic + 1, hPassBitNew->GetBinContent(ic + 1) + 1);
    }
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
    printf("\t\t MESSAGE file %s has %lld RunTree entries total %lld \n", f->GetName(), RunTree->GetEntries(), totalEntries);
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
  cout << " MESSAGE line 475 count files in dir " << dirName << endl;
  TSystemDirectory dir(dirName, dirName); // TSystemDirectory
  TList *files = dir.GetListOfFiles();
  // files->ls();
  //  print list
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
    tmStruct.tm_year = year - 1900; // struct tm counts from 1900
    tmStruct.tm_mon = month - 1;    // struct tm is 0-based
    tmStruct.tm_mday = day;
    time_t fileTime = mktime(&tmStruct);
    const auto diff0 = std::difftime(fileTime, time0);
    const auto diff1 = std::difftime(fileTime, time1);
    bool timetest = diff0 >= 0 && diff1 <= 0;
    // printf("fileTime: %s  month=%i day=%i year=%i  diff0=%.0f diff1=%.0f pass=%i  file=%s\n",
    //        asctime(localtime(&fileTime)), month, day, year, diff0, diff1, int(timetest), tname.Data());
    if (!timetest)
    {
      cout << "   skip out of time file " << name << endl;
      continue;
    }
    TString fullName = dirNameSlash + TString(name.c_str());
    printf("MESSAGE: openFile name %s \n", fullName.Data());
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
  tmStruct.tm_year = year0 - 1900; // struct tm counts from 1900
  tmStruct.tm_mon = month0 - 1;    // struct tm is 0-based
  tmStruct.tm_mday = day0;
  time0 = mktime(&tmStruct);
  printf("set start %s\n", asctime(localtime(&time0)));
  tmStruct.tm_year = year1 - 1900;
  tmStruct.tm_mon = month1 - 1;
  tmStruct.tm_mday = day1;
  time1 = mktime(&tmStruct);
  printf("set end %s\n", asctime(localtime(&time1)));
}

void loop()
{
  /* nominal gains have been applied in pulse finding step */
  totalPass = 0;
  printf(" start of entry loop maxEntry=%lld\n", maxEntry);
  // loop over entries
  for (Long64_t entry = 0; entry < maxEntry; ++entry)
  {

    ++hitCountNev;

    if (entry / 10000 * 10000 == entry)
    {
      printf("line330 .....loop entry %lld nev this file %i \n", entry, hitCountNev);
      // force printing to log file
      fflush(stdout);
    }
    int passBit = passEventCuts(entry);
    // set to pass for debugging
    // cut on passBit passBit = 0;

    /*
    if (TRIGTIME & passBit)
    {
      printf("MESSAGE line 340 event %lld fails trigtime %i pass = %i %i \n", entry, vTrigTime[0], passBit, passBit & TRIGTIME);
    }
      */
    hEventPassNew->SetBinContent(passBit, hEventPassNew->GetBinContent(passBit) + 1);
    if (passBit != 0 && !isLedRun && !isSimulation)
      continue;

    ++totalPass;

    RunTree->GetEntry(entry);
    TString theFileName = TString(RunTree->GetCurrentFile()->GetName());

    if (theFileName != currentFileName)
    {
      ++currentFileNumber;
      currentFileName = theFileName;
      printf("Processing file %s number %i  \n", currentFileName.Data(), currentFileNumber);
      fflush(stdout);

      // Extract file number from filename string (e.g., "84" from "anaCRun-run-05_14_2026-file_84.root-0.root")
      string stringFileName = string(currentFileName.Data());
      size_t pos = stringFileName.find("file_");
      int fileNum = -1; // Default value if file number is not found
      if (pos != string::npos)
      {
        size_t endPos = stringFileName.find_first_not_of("0123456789", pos + 5);
        string fileNumStr = stringFileName.substr(pos + 5, endPos - (pos + 5));
        fileNum = stoi(fileNumStr);
        cout << "MESSAGE line 595: Extracted file number: " << fileNum << endl;
      }

      if (currentFileNumber > 1)
      {
        printf("Filling file %d nev %i chan7 noise %f \n", currentFileNumber, hitCountNev, earlyHitCountFile[7]);
        fflush(stdout);
        for (int ichan = 0; ichan < earlyHitCountFile.size(); ++ichan)
        {
          ntHitCount->Fill(fileNum, hitCountNev, ichan, earlyHitCountFile[ichan], lateHitCountFile[ichan]);
        }
        // reset hit counts for new file
        std::fill(earlyHitCountFile.begin(), earlyHitCountFile.end(), 0);
        std::fill(lateHitCountFile.begin(), lateHitCountFile.end(), 0);
        hitCountNev = 0;
      }
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

    // get list of branches and save in detList/
    std::vector<TDet *> detList;
    while ((aBranch = (TBranchElement *)next()))
    {
      // skip eventData branch
      if (TString(aBranch->GetName()) == TString("eventData"))
      { // skip this branch
        continue;
      }
      // int idet = TString(TString(aBranch->GetName())(4, 2)).Atoi();

      /* the branch is class TDet so cast it as such */
      TDet *det = (TDet *)aBranch->GetObject();
      detList.push_back(det);
    }
    for (unsigned idet = 0; idet < detList.size(); ++idet)
    {
      TDet *det = detList[idet];
      bool trig = false; // define trigger sipms
      if (idet == 9 || idet == 10 || idet == 11)
        trig = true;

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
        if (idet < 13) // include PMT
        {
          if (thit.firstBin < 600)
            earlyHitCountFile[idet] = earlyHitCountFile[idet] + 1;
          if (thit.firstBin >= 7500 - 600)
            lateHitCountFile[idet] = lateHitCountFile[idet] + 1;
        }

        // trigger time shift to time of first trigger SIPM
        int theHitTime = thit.firstBin + nominalTrigger - theMaximumBin;
        // int theHitTime = thit.firstBin;
        // if (thit.firstBin < nominalTrigger && idet == 9)
        //  printf("line 729 event %lld det %i hit %i firstBin %i theMaximumBin %i theHitTime %i\n", entry, idet, ihit, thit.firstBin, theMaximumBin, theHitTime);
        ntTrigChan->Fill(double(entry), vTrigTime[0], idet, theHitTime, vTrigTime[0], double(passBit));

        // fill light curve
        // printf("line 729 det %i max bin: %i\n", idet, hTrigTime[idet]->GetMaximumBin());
        hLightCurve[idet]
            ->SetBinContent(theHitTime + 1, hLightCurve[idet]->GetBinContent(theHitTime + 1) + thit.qpeak / readGains->sipmPeakGain[idet]);
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
        // get the next hit
        double nextHitStartTime = 7500; // default to end of window startTime is a double
        if (ihit < det->hits.size() - 1)
          nextHitStartTime = det->hits[ihit + 1].startTime;
        hNextHitTime[idet]->Fill(nextHitStartTime);
      } // end branch loop
      ntLateInt->Fill(double(entry), double(idet), qsumLate[0], qsumLate[1], qsumLate[2], qsumLate[3], qsumLate[4], qsumLate[5], qsumLate[6], qsumLate[7], qsumLate[8], qsumLate[9], qsumLate[10], qsumLate[11]);
    } // branch

    hGammaPeakHit->Fill(eventTriggerHitQsum);
    // triangle variables
    double qFraction[3];
    qFraction[0] = qsumSum[9];
    qFraction[1] = qsumSum[10];
    qFraction[2] = qsumSum[11];
    double xternQ, yternQ;
    makeTernary(qFraction[0], qFraction[1], qFraction[2], xternQ, yternQ);

    // ntGamma = new TNtuple("ntGamma", "gamma peak", "event:ph9:qsum9:ph10:qsum10:ph11:qsum11:peak");
    ntGamma->Fill(double(entry), photonSum[9], qsumSum[9], photonSum[10], qsumSum[10], photonSum[11], qsumSum[11], photonSum[12], qsumSum[12], eventTriggerHitQsum, triggerSum, xternQ, yternQ);

    // cross talk loop over dets
    for (unsigned idetNumber = 0; idetNumber < detList.size(); ++idetNumber)
    {
      // idetNumber hit loop
      for (unsigned ihit = 0; ihit < detList[idetNumber]->hits.size(); ++ihit)
      {
        TDetHit iDetHit = detList[idetNumber]->hits[ihit];

        // loop over all other detectors
        for (unsigned jdetNumber = 0; jdetNumber < detList.size(); ++jdetNumber)
        {
          double nextHitStartTime = 7500;
          if (jdetNumber == idetNumber)
            continue;
          // other det hits loop
          for (unsigned jhit = 0; jhit < detList[jdetNumber]->hits.size(); ++jhit)
          {
            TDetHit jDetHit = detList[jdetNumber]->hits[jhit];
            if (jDetHit.startTime <= iDetHit.startTime)
              continue;
            nextHitStartTime = jDetHit.startTime;
            break; // only want the next hit after this one
          } // other hit loop
          hNextHitTimeOther[idetNumber]->Fill(nextHitStartTime);
        } // other det loop
      } // this det hit loop
    } // this det loop
  } // entry
} // end of loop function
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
  printf("MESSAGE line 710 files in chain:\n");
  RunTree->GetListOfFiles()->Print();
  Long64_t ntriggers = RunTree->GetEntries();
  printf("MESSAGE line 715  in post: tag %s total triggers in this chain %lld \n", tag.Data(), ntriggers);
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
  printf("MESSAGE line 745 average gain %f \n", aveGain);

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
  fout->cd();
  // trigger info ntuple
  ntHitCount = new TNtuple("ntHitCount", "hit count", "file:nev:chan:early:late");
  ntTrig = new TNtuple("ntTrig", "trigger info", "event:pmtTotSum:totSum13:triggerSum:qun0:qun1:qun2:q0:q1:q2:xQ:yQ:passBit");
  ntTrigChan = new TNtuple("ntTrigChan", "trigger info by channel", "event:trigTime:chan:hitTime:maxBin:passBit");
  // ntLateInt = new TNtuple("ntLateInt", "late integral", "event:lateTotSum:totSum13:triggerSum:qun0:qun1:qun2:q0:q1:q2:xQ:yQ:passBit");
  ntGamma = new TNtuple("ntGamma", "gamma peak", "event:ph9:qsum9:ph10:qsum10:ph11:qsum11:ph12:qsum12:hitSum:ADCSum:xternQ:yternQ");

  ntLateSum = new TNtuple("ntLateSum", "late sum info", "event:chan:geo:lateSum");
  ntPreSum = new TNtuple("ntPreSum", "pre sum info", "event:chan:geo:preSum");
  ntLateInt = new TNtuple("ntLateInt", "late integral by channel", "event:chan:int0:int1:int2:int3:int4:int5:int6:int7:int8:int9:int10:int11");
  // make histograms
  // upper edge of last bin = 8
  hPassBitNew = new TH1D("PassBitNew", "pass bit", FAILBITS - 1, 0, FAILBITS - 1);
  hEventPassNew = new TH1D("EventPassNew", " remade event failures", TOTALCODES, 0, TOTALCODES);
  hCosmicCut = new TH1D("CosmicCut", "cosmic cut pmtot totSum/nominal gain ", 1500., 0, 1500.);
  hOverlapEvent = new TH1D("OverlapEvent", " overlap event qsum13/nominal gain", 15000, 0., 15000.);
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

  // cross talk plots
  crossDir = fout->mkdir("crossDir");
  crossDir->cd();
  for (unsigned ichan = 0; ichan < CHANNELS; ++ichan)
  {
    hNextHitTime.push_back(new TH1D(Form("NextHitTimeChan%i", ichan), Form("NextHitTimeChan%i samples", ichan), 500, 0, 500));
    hNextHitTimeOther.push_back(new TH1D(Form("NextHitTimeOtherChan%i", ichan), Form("NextHitOtherTimeChan%i samples", ichan), 500, 0, 500));
    hMulti.push_back(new TH1D(Form("MultiChan%i", ichan), Form("MultiChan%i samples", ichan), 100, 0., 10.));
    hMulti.back()->GetXaxis()->SetTitle("number of SPE");
    hMulti.back()->GetYaxis()->SetTitle("number of hits");
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
  hTrigTimeEvent = new TH1D(Form("TrigTimeEvent"), Form("TrigTimeEvent samples"), 1000, 0, 1000);
  for (unsigned ichan = 0; ichan < CHANNELS; ++ichan)
  {
    qpeakLimit = 50. * (readGains->sipmPeakGain[ichan]);
    qsumLimit = 50. * (readGains->sipmSumGain[ichan]);
    hQPeak.push_back(new TH1D(Form("QPeakChan%i", ichan), Form("QPeakChan%i", ichan), 2000, 0, qpeakLimit));
    hQSum.push_back(new TH1D(Form("QSumChan%i", ichan), Form("QSumChan%i", ichan), 2000, 0, qsumLimit));

    hTrigTime.push_back(new TH1D(Form("TrigTimeChan%i", ichan), Form("TrigTimeChan%i samples", ichan), 50, 700, 750));
    hTrigTime.back()->GetXaxis()->SetTitle("time [samples]");
    hTrigTime.back()->GetYaxis()->SetTitle("number of hits");
  }

  printf("MESSAGE line 835 starting loop \n");
  fout->ls();
  /*
   *  loop over events
   */
  earlyHitCountFile.resize(12);
  lateHitCountFile.resize(12);
  loop();

  // store from last file
  printf("MESSAGE line 870 Filling file %s number %d nev %i \n", currentFileName.Data(), currentFileNumber, hitCountNev);
  for (int ichan = 0; ichan < earlyHitCountFile.size(); ++ichan)
  {
    ntHitCount->Fill(currentFileNumber, hitCountNev, ichan, earlyHitCountFile[ichan], lateHitCountFile[ichan]);
  }

  // Print ntHitCount entries
  printf("\n MESSAGE line 885=== ntHitCount entries %lli ===\n", ntHitCount->GetEntries());
  printf("File    Chan    Early   Late\n");
  printf("----    ----    -----   ----\n");
  // ntHitCount->Scan("file:chan:early:late", "", "");

  float fnfile = 0;
  float fchan = 0;
  float fearly = 0;
  float flate = 0;
  float fnev = 0;
  ntHitCount->SetBranchAddress("file", &fnfile);
  ntHitCount->SetBranchAddress("nev", &fnev);
  ntHitCount->SetBranchAddress("chan", &fchan);
  ntHitCount->SetBranchAddress("early", &fearly);
  ntHitCount->SetBranchAddress("late", &flate);
  for (int i = 0; i < ntHitCount->GetEntries(); i++)
  {
    ntHitCount->GetEntry(i);
    printf("MESSAGE line 905 file %.0f  events  %.0f chan   %.0f  early  %.0f  late  %.0f \n", fnfile, fnev, fchan, fearly, flate);
  }

  printf("MESSAGE line 920 total %llu pass %llu \n", maxEntry, totalPass);
  hPassBitNew->Print("all");
  printf("MESSAGE line 925 pass fractions total = %.0f \n", hEventPassNew->GetEntries());
  for (int ibin = 0; ibin < hEventPassNew->GetNbinsX(); ++ibin)
  { // inc/lude error on poisson probability
    double nbin = hEventPassNew->GetBinContent(ibin);
    double ntot = hEventPassNew->GetEntries();
    double prob = nbin / ntot;
    double perror = sqrt(prob * (1. - prob) / ntot);

    /*
    if (nbin > 0)
      printf("MESSAGE line 926 bin %i fail %.f frac %.3f +/- %.3f name %s \n", ibin, hEventPassNew->GetBinContent(ibin), prob, perror, codeNames[ibin].Data());
      */
  }

  // do not normilzed summed chan 13
  for (int ich = 0; ich < hLightCurve.size() - 1; ++ich)
    normalize(ich);

  hPassBitNew->Print("all");
  // loop over fail bits
  printf("MESSAGE line 1050 summary of bit failures %llu pass %llu all %0.f\n", maxEntry, totalPass, hPassBitNew->GetBinContent(0));
  for (int ic = 0; ic < hPassBitNew->GetNbinsX(); ++ic)
  {
    double prob = hPassBitNew->GetBinContent(ic + 1) / double(maxEntry);
    double perror = sqrt(prob * (1. - prob)) / double(maxEntry);
    printf("bit %i %s number %.0f frac %.5f +/- %.5f \n", ic, bitNames[ic + 1].Data(), hPassBitNew->GetBinContent(ic + 1), prob, perror);
  }

  // fout->ls();
  // fout->ls();
}

int main(int argc, char *argv[])
{

  fout = nullptr;
  time_t now = time(0);
  cout << "MESSAGE line 1070 executing " << argv[0] << " post hit finding analysis " << "Date and time: " << ctime(&now) << endl;

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
  failCode[7] = TRIGTIME;
  failCode[8] = TRIANGLE;
  failCode[9] = MULTI;
  failCode[10] = TOTALCODES;

  vecFail.resize(FAILBITS);
  bitNames.resize(FAILBITS);
  // offset by one to include ALL
  bitNames[0] = TString("All");
  bitNames[1] = TString("Pass");
  bitNames[2] = TString("Baseline");
  bitNames[3] = TString("Earlycut");
  bitNames[4] = TString("Firsttime");
  bitNames[5] = TString("Cosmic");
  bitNames[6] = TString("Gamma");
  bitNames[7] = TString("Trigger");
  bitNames[8] = TString("TriggerTime");
  bitNames[9] = TString("Triangle");
  bitNames[10] = TString("Mult");

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

  printf("MESSAGE line 960 failure codes: \n");
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
    printf("\n\nMESSAGE line 975 *************** THIS IS LED RUN ***************\n\n");
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
  printf("MESSAGE line 985 count files from %s to %s total files  %ld \n", theStartTag.Data(), theEndTag.Data(), fileListName.size());

  if (nfiles == 0)
  {
    printf(" >>>> datatype no files found <<<<\n");
    exit(-1);
  }

  if (fileListName[0].Contains("btbSim"))
    isSimulation = true;

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

  printf("MESSAGE line 1000 >>>>>> analyze %u  files from %s to %s tag %s totalEntries %lld maxEntry %lld <<<<<\n", nfiles, theStartTag.Data(), theEndTag.Data(), tag.Data(), totalEntries, maxEntry);

  sdate = currentDate();
  tag = theStartTag + TString("-") + theEndTag;
  TString sentries;
  sentries.Form("-%llu", maxEntry);

  TString outFileName;
  outFileName = TString("post-") + tag + sentries + TString(".root");
  if (isSimulation)
    outFileName = TString("post-btbSim-") + tag + sentries + TString(".root");

  fout = new TFile(outFileName, "recreate");
  gainDir = fout->mkdir("gainDir");
  // pick up first pass hEventCount now that fout is open
  for (int i = 0; i < fileListName.size(); ++i)
  {
    if (!fout)
      fout = new TFile(outFileName, "update");
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
  {
    fout = new TFile(outFileName, "update");
    gainDir = fout->mkdir("gainDir");
  }

  cout << "MESSAGE line 1030 starting summary for   " << fileListName.size() << " on " << sdate << " writing to file " << fout->GetName() << endl;
  /* here we make the TCHain and them loop over it */

  post(tag);
  printf("MESSAGE line 1050 write file at end of job \n");
  fout->Write();
  fout->Close();
  delete RunTree;
  printf("MESSAGE line 1060 end of job \n");
}