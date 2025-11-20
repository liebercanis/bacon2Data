/*
    program to analyze RunTree chain from date tag
    this is second pass after pulse findiing anacRunGamma.cc has been run
    uses TTree RunTree making a chain from date tag xx_xx_yyyy
        M Gold Nov 12 2025
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
TChain *RunTree;
TFile *fout;
TString tag;
TNtuple *ntSum;
TNtuple *ntHit;
TNtuple *ntTDiff;
int passBit;

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

//
enum
{
  CHANNELS = 14,
  NONSUMCHANNELS = CHANNELS - 1
};

enum
{
  MAXSAMPLES = 7500
}

TReadGains *readGains;
std::vector<TString> bitNames;

std::vector<int> vTotal;
std::vector<int> vPass;

int nFiles;
TString tag;
TString theStartTag;
TString theEndTag;
std::string sdate;
vector<double> normQsum;
vector<double> normQPE;
TDatime dateTime;
time_t time0;
time_t time1;
struct tm tmStruct;
Long64_t maxFiles;
double ntotal;
double npass;
vector<int> filePass;
vector<int> fileTotal;
int totalPass;
int totalEvents;

// old stuff below
//  hSipmLatePeaksFitograms as vectors
std::vector<TH1D *> hTotSum;
std::vector<TH1D *> hPreSum;
std::vector<TH1D *> hTrigSum;
std::vector<TH1D *> hLateSum;

double y[NCHAN];
double nominalGain;

TH1D *hPassBit;
TH1D *hTrigEventSumArea;
TH1D *hTrigEventSumPeak;
TH2D *hTrigHitPeakTime;
TH2D *hTrigHitPeakTimeCoarse;
TH2D *hAllHitPeakTimeCoarse;
TH2D *hAllHitPeakTime;
TH1D *hHitTimeDiff;

TH1D *hTrigLatePeaks;
TH1D *hSipmLatePeaks;
TH1D *hTrigLatePeaksFit;
TH1D *hSipmLatePeaksFit;
TH1D *hCountPre;
TH1D *hCountLate;
TH1D *hCountLateTime;
TH2D *hCountLateTimeQpeak;

unsigned sipmCut;
unsigned trigCut;
double peakCut;

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
// count subruns and channels
unsigned long countFiles()
{
  dirName = TString("caenData");
  dirNameSlash = TString("caenData/");
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

    // see if file is good
    TString fullName = dirNameSlash + TString(name.c_str());
    // open file and get pointers
    TFile *f = new TFile(fullName, "READONLY");
    if (getPointers(f))
      fileList.push_back(TString(name.c_str()));
    f->Close();
  }
  return fileList.size();
}

TDatime getTime(int ifile, Long64_t ievent = 0)
{
  TDatime datime;
  if (runTree)
  {
    runTree->GetEntry(ievent);
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

void loop(Long64_t maxEntry)
{
  sipmCut = 0;
  trigCut = 0;
  // loop over entries
  for (Long64_t entry = 0; entry < maxEntry; ++entry)
  {
    passBit = 0;
    int nPre = 0;
    int nLate = 0;
    // get entry
    RunTree->GetEntry(entry);
    // get branch pointers and save in detList
    TIter next(RunTree->GetListOfBranches());
    TBranchElement *aBranch = NULL;
    // loop over branches
    double trigPre = 0;
    double trigTrig = 0;
    double trigLate = 0;
    double sipmPre = 0;
    double sipmTrig = 0;
    double sipmLate = 0;
    double trigEventSumArea = 0;
    double trigEventSumPeak = 0;
    double trigTotPeak = 0;

    bool trigIsCut = false;
    bool sipmIsCut = false;

    while ((aBranch = (TBranchElement *)next()))
    {
      int id = TString(TString(aBranch->GetName())(4, 2)).Atoi();
      if (id >= NCHAN) // skip PMT
        continue;
      bool trig = false; // define trigger sipms
      if (id == 9 || id == 10 || id == 11)
        trig = true;
      if (TString(aBranch->GetName()) == TString("eventData")) // skip this branch
        continue;
      TDet *det = (TDet *)aBranch->GetObject();
      // if(id==0) cout << "branch " << aBranch->GetName() << " entry " << entry << " TotSum " << det->totSum << endl;
      // do not normlize to nominal
      double oldNominal = 270.5;
      hTotSum[id]->Fill(det->totPeakSum);
      hPreSum[id]->Fill(det->prePeakSum);
      hTrigSum[id]->Fill(det->trigPeakSum);
      hLateSum[id]->Fill(det->latePeakSum);
      if (trig)
      {
        trigPre += det->prePeakSum / nominalGain;
        trigTrig += det->trigPeakSum / nominalGain;
        trigLate += det->latePeakSum / nominalGain;
        trigEventSumPeak += det->trigPeakSum / nominalGain;
        trigEventSumArea += det->trigSum;
        trigTotPeak += (det->trigPeakSum + det->latePeakSum) / nominalGain;
      }
      else
      {
        sipmPre += det->prePeakSum / nominalGain;
        sipmTrig += det->trigPeakSum / nominalGain;
        sipmLate += det->latePeakSum / nominalGain;
      }
      // loop over hits

      for (unsigned ihit = 0; ihit < det->hits.size(); ++ihit)
      {
        TDetHit hiti = det->hits[ihit];
        // printf(" det %i time %.0f qpeak %f \n",id,det->hits[ihit].startTime, det->hits[ihit].qpeak);
        if (trig)
          hTrigHitPeakTime->Fill(det->hits[ihit].startTime, det->hits[ihit].qpeak / nominalGain);
        if (trig)
          hTrigHitPeakTimeCoarse->Fill(det->hits[ihit].startTime, det->hits[ihit].qpeak / nominalGain);
        hAllHitPeakTime->Fill(det->hits[ihit].startTime, det->hits[ihit].qpeak / nominalGain);
        hAllHitPeakTimeCoarse->Fill(det->hits[ihit].startTime, det->hits[ihit].qpeak / nominalGain);
        ntHit->Fill(double(entry), double(id), det->hits[ihit].startTime, det->hits[ihit].qpeak / nominalGain);

        if (trig && det->hits[ihit].startTime > 800)
          hTrigLatePeaks->Fill(det->hits[ihit].qpeak / nominalGain);
        if (!trig && det->hits[ihit].startTime > 800)
          hSipmLatePeaks->Fill(det->hits[ihit].qpeak / nominalGain);

        if (trig && (det->hits[ihit].startTime > 800 || det->hits[ihit].startTime < 600) && det->hits[ihit].qpeak / nominalGain > 10)
          trigIsCut = true;
        if (!trig && (det->hits[ihit].startTime > 800 || det->hits[ihit].startTime < 600) && det->hits[ihit].qpeak / nominalGain > 10)
          sipmIsCut = true;

        /*
        if ((det->hits[ihit].startTime > 800 || det->hits[ihit].startTime < 600) && det->hits[ihit].qpeak / nominalGain > 10)
          printf("xxxxx event %llu chan %i time %f qpeak %f\n", entry, id, det->hits[ihit].startTime, det->hits[ihit].qpeak / nominalGain);
          */
      }

      // look at time difference between hits.
      if (id < 9)
      {
        for (unsigned ihit = 0; ihit < det->hits.size(); ++ihit)
        {
          TDetHit hiti = det->hits[ihit];
          for (unsigned jhit = ihit + 1; jhit < det->hits.size(); ++jhit)
          {
            TDetHit hitj = det->hits[ihit];
            double tdiff = double(det->hits[jhit].startTime) - double(det->hits[ihit].startTime);
            // printf("%i %i %f\n", ihit, jhit, tdiff);
            hHitTimeDiff->Fill(tdiff);
            ntTDiff->Fill(double(det->hits[ihit].startTime), double(det->hits[jhit].startTime), det->hits[ihit].qpeak, det->hits[jhit].qpeak); //= new TNtuple("ntTDiff","time diff","startTime1:startTime2:qpeak1:qpeak2");
          }
        }
      }

      // pre cut
      int npreHits = 0;
      int nlateHits = 0;
      if (id == 13)
      {
        for (unsigned ihit = 0; ihit < det->hits.size(); ++ihit)
        {
          TDetHit hiti = det->hits[ihit];
          if (det->hits[ihit].startTime < 600)
          {
            ++npreHits;
            ++nPre;
          }
          hCountLateTimeQpeak->Fill(det->hits[ihit].startTime, det->hits[ihit].qpeak / nominalGain);
          if (det->hits[ihit].startTime > 1000 && det->hits[ihit].qpeak / nominalGain > peakCut)
          {
            ++nlateHits;
            ++nLate;
            hCountLateTime->Fill(det->hits[ihit].startTime);
          }
        }
        hCountPre->Fill(npreHits);
        hCountLate->Fill(nlateHits);
      }
    } // end branch loop
    if (trigEventSumPeak < 0)
      trigEventSumPeak = 0;
    ntSum->Fill(trigEventSumArea, trigEventSumPeak, trigTotPeak, trigPre, trigTrig, trigLate, sipmPre, sipmTrig, sipmLate);
    if (entry == entry / 10000 * 10000)
      printf("... event %llu trig sum %.2E %.2E \n", entry, trigEventSumArea, trigEventSumPeak);
    hTrigEventSumArea->Fill(trigEventSumArea);
    hTrigEventSumPeak->Fill(trigEventSumPeak);
    if (trigIsCut)
      ++trigCut;
    if (sipmIsCut)
      ++sipmCut;
    if (nPre > 0)
      passBit |= 0x1;
    if (nLate > 0)
      passBit |= 0x2;
    hPassBit->Fill(passBit);
  }
}

void post(TString tag = TString("10_06_2025"), Long64_t maxEntry = 0)
{
  peakCut = 6.5;
  /*gains-2024-02-01-17-06.root*/
  y[0] = 229.5;
  y[1] = 221;
  y[2] = 236;
  y[3] = 200.6;
  y[4] = 229;
  y[5] = 229; // not fit same as 4
  y[6] = 231.6;
  y[7] = 237.2;
  y[8] = 232.6;
  y[9] = 646.9;
  y[10] = 619.5;
  y[11] = 605;
  nominalGain = 0;
  for (int k = 0; k < 9; ++k)
    nominalGain += y[k];
  nominalGain /= 9.;
  printf(" the tag is %s nominal gain = %f ", tag.Data(), nominalGain);
  gStyle->SetOptStat(1001101);
  /* get RunTree */
  RunTree = new TChain("RunTree");
  TString name;
  name.Form("caenData/anaCRun*%s*.root", tag.Data());
  printf("open chain with %s \n", name.Data());
  RunTree->Add(name);
  if (!RunTree)
    return;
  printf("files in chain:\n");
  RunTree->GetListOfFiles()->Print();
  Long64_t ntriggers = RunTree->GetEntries();
  printf(" total triggers in this chain %lld \n", ntriggers);
  if (maxEntry == 0)
    maxEntry = ntriggers;
  RunTree->GetListOfBranches()->ls();
  TString sentries;
  sentries.Form("-%llu", maxEntry);
  fout = new TFile(TString("post-") + tag + sentries + TString(".root"), "recreate");

  hHitTimeDiff = new TH1D("HitTimeDiff", "samples between hits ", 3000, 0, 3000);
  hTrigEventSumArea = new TH1D("TrigEventSumArea", "trig sipm trigger window sum area ", 600, 0, 6.E5);
  hTrigEventSumPeak = new TH1D("TrigEventSumPeak", "trig sipm trigger window sum peak ", 400, 0, 40);
  hTrigHitPeakTime = new TH2D("TrigHitPeakTime", "trig summed trig peak versus time ", 7500, 0, 7500, 400, 0, 40);
  hAllHitPeakTime = new TH2D("AllHitPeakTime", "trig summed trig peak versus time ", 7500, 0, 7500, 400, 0, 40);
  hTrigHitPeakTimeCoarse = new TH2D("TrigHitPeakTimeCoarse", "trig summed trig peak versus time ", 75, 0, 7500, 400, 0, 40);
  hAllHitPeakTimeCoarse = new TH2D("AllHitPeakTimeCoarse", "trig summed trig peak versus time ", 75, 0, 7500, 400, 0, 40);

  //
  TString htitle;
  hCountPre = new TH1D("CountPre", " hits sample<600 in sum", 20, 0, 20);
  htitle.Form("hits qpeak>%.2f SPE sample>1000 in sum", peakCut);
  hCountLate = new TH1D("CountLate", htitle, 20, 0, 20);
  htitle.Form("umber of late time hits with qpeak>%.2f", peakCut);
  hCountLate->GetXaxis()->SetTitle(htitle);
  htitle.Form("hits qpeak>%.2f SPE sample>1000 in sum", peakCut);
  hCountLateTime = new TH1D("CountLateTime ", htitle, 30, 0, 7500);
  hCountLateTime->GetXaxis()->SetTitle("sample time");
  hCountLateTime->Sumw2();
  hCountLateTimeQpeak = new TH2D("CountLateTimeQpeak", " sample>1000 in sum qpeak vs time ", 30, 0, 7500, 20, 0, 20);
  hCountLateTimeQpeak->GetXaxis()->SetTitle("sample time");
  hCountLateTimeQpeak->GetYaxis()->SetTitle("qpeak [SPE]");

  ntTDiff = new TNtuple("ntTDiff", "time diff", "startTime1:startTime2:qpeak1:qpeak2");
  ntSum = new TNtuple("ntSum", " ADC sums ", "trigSumArea:trigSumPeak:trigTotPeak:trigPre:trigTrig:trigLate:sipmPre:sipmTrig:sipmLate");
  ntHit = new TNtuple("ntHit", " hits ", "event:chan:time:qpeak");

  hTrigLatePeaks = new TH1D("TrigEventSumLatePeak", "trig sipm trigger late peaks ", 1000, 0, 100);
  hTrigLatePeaksFit = new TH1D("TrigEventSumPeakFit", "trig sipm trigger late peaks ", 1000, 0, 100);
  hSipmLatePeaks = new TH1D("SipmEventSumPeak", "non-trig sipm late peaks ", 1000, 0, 100);
  hSipmLatePeaksFit = new TH1D("SipmEventSumPeakFit", "non-trig sipm late peaks ", 1000, 0, 100);
  hPassBit = new TH1D("PassBit", "pass bit", 4, 0, 4);

  // make hSipmLatePeaksFitos
  for (unsigned i = 0; i < NCHAN; ++i)
  {
    double limit = 40;
    int nbins = 400.;
    // normalized to SPE
    hTotSum.push_back(new TH1D(Form("TotPeakSumChan%i", i), Form("tot peak sum chan %i", i), nbins, 0, limit));
    hPreSum.push_back(new TH1D(Form("PrePeakSumChan%i", i), Form("pre peak sum chan %i", i), nbins, 0, limit));
    hTrigSum.push_back(new TH1D(Form("TrigPeakSumChan%i", i), Form("trig peak sum chan %i", i), nbins, 0, limit));
    hLateSum.push_back(new TH1D(Form("LatePeakSumChan%i", i), Form("late peak sum chan %i", i), nbins, 0, limit));
  }

  //
  loop(maxEntry);

  // poisson fit
  TF1 *fpoi1 = new TF1("fpoi1", "[1]*pow([0],x)*Exp(-[0])/Gamma(x+1.)", 0, 10);
  double norm = hSipmLatePeaks->Integral();
  // set non-zero initial values for parameters
  fpoi1->SetParameter(0, 1);
  fpoi1->SetParameter(1, norm);
  hSipmLatePeaks->Fit("fpoi1", "R");
  for (int ib = 1; ib < hSipmLatePeaksFit->GetNbinsX(); ++ib)
  {
    double xbin = hSipmLatePeaksFit->GetBinCenter(ib) - 0.5;
    double fbin = fpoi1->Eval(xbin);
    hSipmLatePeaksFit->SetBinContent(ib, fbin);
    hSipmLatePeaksFit->SetBinError(ib, 0);
    hSipmLatePeaksFit->GetYaxis()->SetTitle("yield");
    hSipmLatePeaksFit->GetXaxis()->SetTitle("SPE");
  }
  hSipmLatePeaksFit->SetLineColor(kRed);

  TCanvas *canSipmLate = new TCanvas("sipmLate", "sipmLate");
  hSipmLatePeaks->Draw();
  hSipmLatePeaksFit->Draw("same");

  // R" = fit between "xmin" and "xmax" of the "f1"

  // poisson fit
  TF1 *fpoi2 = new TF1("fpoi2", "[1]*pow([0],x)*Exp(-[0])/Gamma(x+1.)", 0, 10);
  norm = hTrigLatePeaks->Integral();
  // set non-zero initial values for parameters
  fpoi1->SetParameter(0, 1);
  fpoi1->SetParameter(1, norm);
  hTrigLatePeaks->Fit("fpoi1", "R");
  for (int ib = 1; ib < hTrigLatePeaksFit->GetNbinsX(); ++ib)
  {
    double xbin = hTrigLatePeaksFit->GetBinCenter(ib) - 0.5;
    double fbin = fpoi1->Eval(xbin);
    hTrigLatePeaksFit->SetBinContent(ib, fbin);
    hTrigLatePeaksFit->SetBinError(ib, 0);
    hTrigLatePeaksFit->GetYaxis()->SetTitle("yield");
    hTrigLatePeaksFit->GetXaxis()->SetTitle("SPE");
  }
  hTrigLatePeaksFit->SetLineColor(kRed);

  TCanvas *canTrigLate = new TCanvas("trigLate", "trigLate");
  hTrigLatePeaks->Draw();
  hTrigLatePeaksFit->Draw("same");

  printf("total %llu trig cut %u (%f) sipm cut %u (%f) \n",
         maxEntry, trigCut, double(trigCut) / double(maxEntry), sipmCut, double(sipmCut) / double(maxEntry));

  fout->Write();
  printf(" cut is %.2f \n", peakCut);
  hPassBit->Print("all");
  // fout->ls();
}

int main(int argc, char *argv[])
{
  cout << "executing " << argv[0] << " post hit finding analysis  " << endl;
  printf(" usage:  start date string <stag> end date string <etag> max files <default all> \n ");
  if (argc < 2)
  {
    printf("reguire file date start string <stag> args\n");
    exit(0);
  }
  vecFail.resize(FAILBITS);
  bitNames.resize(FAILBITS);
  bitNames[0] = TString("pass");
  bitNames[1] = TString("baseline");
  bitNames[2] = TString("earlycut");
  bitNames[3] = TString("firsttime");
  bitNames[4] = TString("cosmic");
  bitNames[5] = TString("gamma");
  bitNames[6] = TString("trigger");

  TReadGains *readGains = new TReadGains();
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

  printf("count files from %s to %s \n", theStartTag.Data(), theEndTag.Data());
  unsigned nfiles = countFiles();
  if (nfiles == 0)
  {
    printf(" >>>> datatype no files found <<<<\n");
    exit(0);
  }

  printf(" >>>>> files from %s to %s tag %s <<<<<\n", theStartTag.Data(), theEndTag.Data(), tag.Data());

  for (int i = 0; i < fileList.size(); ++i)
    cout << i << "  " << fileList[i] << endl;

  maxFiles = fileList.size();
  if (argc > 3)
  {
    maxFiles = atoi(argv[3]);
  }
  printf(" for %s found %lu files maxFiles %lli \n", tag.Data(), fileList.size(), maxFiles);

  sdate = currentDate();
  cout << " starting summary for   " << maxFiles << endl;

  fout = new TFile(Form("postAna-%s-nfiles-%lld-created-%s.root", tag.Data(), maxFiles, sdate.c_str()), "recreate");
}