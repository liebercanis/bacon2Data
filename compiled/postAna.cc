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
TNtuple *ntSum;
TNtuple *ntHit;
TNtuple *ntTDiff;
Long64_t totalEntries;
Long64_t maxEntry;
std::vector<TString> fileListName;
std::vector<std::vector<double>> vecFail;
TBEventData *eventData;
TDirectory *sumDir;
TDirectory *anaDir;
TDirectory *cutDir;

TH1D *hEventPass;
TH1D *eventCount;
TH1D *hEventPassNew;
TH1D *hPassBitNew;
std::vector<TH1D *> hLightCurve;

std::vector<double> qsumGain; // read from class TReadGain

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
int failCode[FAILBITS];

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
int totalPass;
int totalEvents;

double nominalGain;
double nominalTrigGain;

int passEventCuts(Long64_t entry)
{
  // there is only 1 passbit but stored in all det. I admit is bad MGold
  RunTree->GetEntry(entry);
  // RunTree->GetListOfBranches()->ls();
  //   get branch pointers and save in detList
  TIter next(RunTree->GetListOfBranches());
  TBranchElement *aBranch = NULL;
  int passBit = 0;
  while ((aBranch = (TBranchElement *)next()))
  {
    int idet = TString(TString(aBranch->GetName())(4, 2)).Atoi();
    if (TString(aBranch->GetName()) == TString("eventData"))
    { // skip this branch
      continue;
    }

    TDet *det = (TDet *)aBranch->GetObject();
    // loop over hits
    if (idet == 9)
    {
      passBit = det->pass;
      // loop over fail bits
      for (int ic = 0; ic < FAILBITS; ++ic)
      {
        if (passBit & failCode[ic])
        {
          // printf("event %lld det %i bin %s pass %i \n", entry, idet, bitNames[ic].Data(), passBit);
          hPassBitNew->SetBinContent(ic, hPassBitNew->GetBinContent(ic) + 1);
        }
      }
    }
    // if recalculating a cut, would do it here
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

  TH1D *hGammaPeak = nullptr;
  f->GetObject("GammaPeak", hGammaPeak);
  if (!hGammaPeak)
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
  // loop over entries
  for (Long64_t entry = 0; entry < maxEntry; ++entry)
  {
    if (entry / 1000 * 1000 == entry)
      printf("line330 .....loop entry %lld \n", entry);
    int passBit = passEventCuts(entry);
    hEventPassNew->SetBinContent(passBit, hEventPassNew->GetBinContent(passBit) + 1);
    if (passBit != 0)
      continue;

    RunTree->GetEntry(entry);
    // RunTree->GetListOfBranches()->ls();
    //   get branch pointers and save in detList
    TIter next(RunTree->GetListOfBranches());
    TBranchElement *aBranch = NULL;
    // loop over branches

    while ((aBranch = (TBranchElement *)next()))
    {
      int idet = TString(TString(aBranch->GetName())(4, 2)).Atoi();
      bool trig = false; // define trigger sipms
      if (idet == 9 || idet == 10 || idet == 11)
        trig = true;

      // skip eventData branch
      if (TString(aBranch->GetName()) == TString("eventData"))
      { // skip this branch
        continue;
      }

      TDet *det = (TDet *)aBranch->GetObject();
      // printf("det %i hits %lu \n", idet, det->hits.size());
      //  check if passes eventCuts

      // loop over hits
      for (unsigned ihit = 0; ihit < det->hits.size(); ++ihit)
      {
        TDetHit thit = det->hits[ihit];
        // fill light curve
        hLightCurve[idet]->SetBinContent(thit.firstBin + 1, hLightCurve[idet]->GetBinContent(thit.firstBin + 1) + thit.qpeak);
      } // end branch loop
    }
  }
}

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
  hPassBitNew = new TH1D("PassBitNew", "pass bit", FAILBITS, 0, FAILBITS);

  // make histograms
  hEventPassNew = new TH1D("EventPassNew", " remade event failures", TOTALCODES, 0, TOTALCODES);
  for (unsigned i = 0; i < CHANNELS; ++i)
  {
    // normalized to SPE
    hLightCurve.push_back(new TH1D(Form("LightCurveChan%i", i), Form("LightCurveChan%i", i), MAXSAMPLES, 0, 2 * MAXSAMPLES));
    hLightCurve[hLightCurve.size() - 1]->GetXaxis()->SetTitle("time [ns]");
    hLightCurve[hLightCurve.size() - 1]->GetYaxis()->SetTitle("normilized number of photons/2ns");
  }
  // loop over events
  loop();

  printf("total %llu \n", maxEntry);
  // hEventPassNew->Print("all");
  printf("pass fractions total = %.0f  \n", hEventPassNew->GetEntries());
  for (int ibin = 0; ibin < hEventPassNew->GetNbinsX(); ++ibin)
  { // inc/lude error on poisson probability
    double nbin = hEventPassNew->GetBinContent(ibin);
    double ntot = hEventPassNew->GetEntries();
    double prob = nbin / ntot;
    double perror = sqrt(prob * (1. - prob) / ntot);
    printf(" bin %i fail %.f frac %.3f +/- %.3f name %s \n", ibin, hEventPassNew->GetBinContent(ibin), prob, perror, codeNames[ibin].Data());
  }

  hPassBitNew->Print("all");

  // pick up first pass hEventCount now that fout is open
  for (int i = 0; i < fileListName.size(); ++i)
  {
    cout << "getPointers " << i << "  " << fileListName[i] << endl;
    TString dirNameSlash = TString("caenData/");
    TString fullName = dirNameSlash + fileListName[i];
    TFile *f = new TFile(fullName, "READONLY");

    TH1D *hEventPass = nullptr;
    f->GetObject("EventPass", hEventPass);
    if (hEventPass && fout)
    {
      TH1D *hEventPassFile = (TH1D *)hEventPass->Clone(Form("EventPassFile%i", i));
      cout << "line466 append " << hEventPassFile->GetName() << endl;
      fout->Add(hEventPassFile);
      fout->Write();
    }
    f->Close();
  }

  fout->ls();
  // fout->Write();
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

  readGains = new TReadGains();
  // use nominal gains for now FIXME
  nominalGain = readGains->nominalGain;
  nominalTrigGain = readGains->nominalTrigGain;
  readGains->printGains();

  // store qsumGain[ib];
  for (unsigned ch = 0; ch < readGains->sipmSumGain.size(); ++ch)
    qsumGain.push_back(readGains->sipmSumGain[ch]);

  for (unsigned ic = 0; ic < TOTALCODES; ++ic)
    codeNames.push_back(TString("mixed"));
  codeNames[PASS] = TString("pass");
  codeNames[BASEFAIL] = TString("baseline");
  codeNames[TRIANGLE] = TString("triangle");
  codeNames[EARLYCUT] = TString("earlycut");
  codeNames[FIRSTTIME] = TString("firsttime");
  codeNames[COSMIC] = TString("cosmic");
  codeNames[GAMMA] = TString("gamma");
  codeNames[TRIGFAIL] = TString("trigfail");
  codeNames[TRIANGLE] = TString("traingle");

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
  bitNames[0] = TString("pass");
  bitNames[1] = TString("baseline");
  bitNames[2] = TString("earlycut");
  bitNames[3] = TString("firsttime");
  bitNames[4] = TString("cosmic");
  bitNames[5] = TString("gamma");
  bitNames[6] = TString("trigger");
  bitNames[7] = TString("trianlge");

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

  setTime(theStartTag, theEndTag);

  dateTime = TDatime(2023, 3, 9, 22, 0, 0);

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

  printf(" >>>>> analyze %u  files from %s to %s tag %s totalEntries %lld maxEntry %lldb<<<<<\n", nfiles, theStartTag.Data(), theEndTag.Data(), tag.Data(), totalEntries, maxEntry);

  sdate = currentDate();
  tag = theStartTag + TString("-") + theEndTag;
  TString sentries;
  sentries.Form("-%llu", maxEntry);
  fout = new TFile(TString("post-") + tag + sentries + TString(".root"), "recreate");

  // pick up first pass hEventCount now that fout is open
  for (int i = 0; i < fileListName.size(); ++i)
  {
    cout << i << "  " << fileListName[i] << endl;
    /*
    TString dirNameSlash = TString("caenData/");
    TString fullName = dirNameSlash + fileListName[i];
    TFile *f = new TFile(fullName, "READONLY");
    getPointers(f);
    f->Close();
    */
  }

  cout << " starting summary for   " << fileListName.size() << " on " << sdate << " writing to file " << fout->GetName() << endl;
  post(tag);
}