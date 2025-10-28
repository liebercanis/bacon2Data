// removed SumHitWave July 25, 2024
#include <ctime>
#include <iostream>
#include <iterator>
#include <locale>
#include <TString.h>
#include <TROOT.h>
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TDatime.h"
#include "TDirectory.h"
#include "TSystemDirectory.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TList.h"
#include "TF1.h"
#include "TNtuple.h"
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
#include <TLeaf.h>
#include <TFormula.h>
#include <TStyle.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
// bobj classes
#include "TBWave.hxx"
#include "TBRun.hxx"
#include "TBRawEvent.hxx"
#include "hitFinder.hxx"
#include "TBFile.hxx"
#include "TBEventData.hxx"
#include "TReadGains.hxx"

//
enum
{
  CHANNELS = 14,
  NONSUMCHANNELS = CHANNELS - 1
};

enum
{
  FAILBITS = 7
};

TReadGains *readGains;
std::vector<TString> bitNames;

std::vector<int> vTotal;
std::vector<int> vPass;
std::vector<std::vector<double>> vMaxValue;
std::vector<std::vector<double>> vIntegral;

static int waveBins = 7500.;
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
vector<double> sumHits;
static double startTime = 660.; // hWave->GetBinLowEdge(maxBin) + hWave->GetBinWidth(maxBin) / 2.;
static double endTime = 75000.0;

double nominalGain = 134.786401;     // 170.;     // was 160.0; set Jue 13 2025
double nominalTrigGain = 735.688747; //
double nominalQsumGain = 4940.503519;
double nominalQsumTrigGain = 32056.789775;
double nominalPmtGain = 502.;
double nominalQsumPmtGain = 1713;

std::vector<TString> fileList;
std::vector<double> filenum;
std::vector<double> efilenum;
std::vector<double> fileTime;
std::vector<TDatime> fileDatime;
TBFile *bf;
TBEventData *eventData;
TTree *runTree;
TH1D *hGammaPeak;    // on output file
TH1D *hGammaPeakCut; // on output file
TH1D *hEventPass;    // must be in input file
TH1D *hTrigSumCut;
TH1D *hTrigSumNoCut;
TH1D *hRunTrigSumNoCut;
TH1D *hRunTrigSumCut;
std::vector<std::vector<double>> vecFail;
std::vector<double> vecFile;
std::vector<double> QPEMean;
std::vector<double> vecQPEMean;
std::vector<double> vecQPENave;
std::vector<vector<double>> vecQsum;
std::vector<vector<double>> vecEQsum;
std::vector<vector<double>> vecQsumUn;
std::vector<vector<double>> vecEQsumUn;
std::vector<vector<double>> vecPeaksum;
std::vector<vector<double>> vecEPeaksum;
std::vector<vector<double>> vecQPE;
std::vector<vector<double>> vecEQPE;
std::vector<vector<double>> vecPeak;
std::vector<vector<double>> vecEPeak;
std::vector<vector<double>> vecQPESigma;
std::vector<vector<double>> vecEQPESigma;
std::vector<vector<double>> vSlope;
std::vector<vector<double>> vESlope;
// summed waves
std::vector<vector<TH1D *>> vRunQSum;
std::vector<vector<TH1D *>> vRunPeakSum;
std::vector<vector<TH1D *>> vRunPeakWave;
std::vector<TH1D *> hSumWave;
std::vector<TH1D *> hRunPeakWave;
std::vector<TH1D *> hRunSumWave;
std::vector<TH1D *> hUnNormedPeakWave;
std::vector<TH1D *> hUnNormedSumWave;
std::vector<TH1D *> hRunLatePeakSum;
std::vector<TH1D *> vNormByFile;
std::vector<TGraph *> gInte;
std::vector<vector<double>> runSums;
std::vector<vector<TString>> runSumNames;

// these must be in input file
TDirectory *sumDir;
TDirectory *anaDir;
// made on output file
TDirectory *fitSumDir;
TDirectory *waveSumDir;
TDirectory *qpeSumDir;
TDirectory *runSumDir;
TDirectory *gainSumDir;
TFile *fin;
TFile *fout;
bool first = true;
TH1D *hQPEChan;
TH1D *hQPESigmaChan;
TH1D *eventCount;
void makeGraphs();
enum dataType
{
  SIS = 0,
  CAEN = 1
};
int theDataType;
TString dirName;
TString dirNameSlash;

double effQuantum128 = 0.17;
double nPhotons = 50.E3 * 5.486; // norm photon yield 2.743E+05

TString startDate;
TString endDate;

// vectors for gains
std::vector<double> sipmGain;
std::vector<double> sipmGainError;
std::vector<double> sipmSumGain;
std::vector<double> sipmSumGainError;
double xWaveLow = 0;
double xWaveHigh = 7500; // max sample

void saveToOutput(TDirectory *dir, TString inName, TString outName)
{
  TH1D *hIn = NULL;
  /* clone hGammaPeak hGammaPeakCut*/
  fin->GetObject(inName, hIn);
  if (hIn)
  {
    TH1D *hOut = NULL;
    dir->GetObject(outName, hOut);
    if (hOut == NULL)
    {
      hOut = (TH1D *)hIn->Clone(outName);
      hOut->SetTitle(outName);
      // cout << "line613 ... adding  " << hIn->GetName() << " file "  << fin->GetName() << " hit QSum " << hOut->GetEntries() << endl;
      fout->Add(hOut);
    }
    else
    {
      fout->GetObject(outName, hOut);
      hOut->Add(hIn);
    }
  } // if hIn
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

  hEventPass = nullptr;
  f->GetObject("EventPass", hEventPass);
  if (!hEventPass)
  {
    cout << "line1230 skipping BAD file no EventPass " << name << endl;
    isGoodFile = false;
  }

  hGammaPeak = nullptr;
  f->GetObject("GammaPeak", hGammaPeak);
  if (!hGammaPeak)
  {
    cout << "line1230 no GammaPeak " << name << endl;
  }

  hGammaPeakCut = nullptr;
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
  hTrigSumNoCut = nullptr;
  f->GetObject("TrigSumNoCut", hTrigSumNoCut);
  hTrigSumCut = nullptr;
  f->GetObject("TrigSumCut", hTrigSumCut);
  //
  return isGoodFile;
}

/* start of code */

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

// normalize to total pass
void normalizeTotalPass(TString histSet)
{
  printf("line195  \t for set %s  in normalizeTotalPass runSumDir has %d entries total pass %d \n", histSet.Data(), runSumDir->GetList()->GetEntries(), totalPass);

  for (int ichan = 0; ichan < NONSUMCHANNELS; ++ichan)
    printf(" chan %i gain %f \n", ichan, readGains->sipmPeakGain[ichan]);

  TString histName;
  for (int ichan = 0; ichan < NONSUMCHANNELS; ++ichan)
  {
    histName.Form("UnNormed%sChan%i", histSet.Data(), ichan);
    TH1D *hist;
    runSumDir->GetObject(histName, hist);
    if (!hist)
      printf("at line206 %s not found \n", histName.Data());
    TH1D *hSave;
    if (hist)
    {
      if (histSet.Contains("Peak"))
        hSave = hRunPeakWave[ichan];
      else
        hSave = hRunSumWave[ichan];

      sumHits[ichan] = hist->Integral(startTime, endTime);

      // printf("at line174 %s normalize to %d\n", hSave->GetName(), totalPass);
      for (int ibin = 0; ibin < hist->GetNbinsX(); ++ibin)
      {
        double xbin = hist->GetBinContent(ibin);
        // double ebin = abs(hist->GetBinError(ibin));
        double ebin = sqrt(abs(hist->GetBinContent(ibin)));

        // cout << " at line212 " << hist->GetName() << " ibin " << ibin << " xbin " << xbin << " ebin " << ebin << " tot " << totalPass << endl;
        // divide by nominal gain july 22 2024
        hSave->SetBinContent(ibin, xbin / double(totalPass) / readGains->sipmPeakGain[ichan]);
        hSave->SetBinError(ibin, ebin / double(totalPass) / readGains->sipmSumGain[ichan]);
        if (ichan == 7 && ibin > 847 && ibin < 851 && histSet.Contains("Peak"))
          printf("line360XXXXXXX %i %f totalPass %i gain %f corr %f\n", ibin, xbin, totalPass, readGains->sipmPeakGain[ichan], xbin / double(totalPass) / readGains->sipmPeakGain[ichan]);
      }

      cout << "at line171"
           << " file " << filenum.size() << " totalPass " << totalPass << "  " << hist->GetName()
           << " peak value " << hist->GetBinContent(hist->GetMaximumBin()) << " integral " << hist->Integral(startTime, endTime)
           << " normed  " << hSave->GetName() << " integral " << hSave->Integral(startTime, endTime) << endl;
    }
  }
}

void fitQPE()
{

  cout << " ***** fitQPE ***** " << endl;

  // for (unsigned ichan = 0; ichan < vRunQSum.size(); ++ichan)
  for (unsigned ichan = 0; ichan < CHANNELS; ++ichan)
  {

    cout << " fitQPE " << ichan << " num histos " << vRunPeakSum[ichan].size() << endl;
    for (unsigned ihist = 0; ihist < vRunQSum[ichan].size(); ++ihist)
    {
      TString cloneName;
      TH1D *hclone = NULL;
      cloneName.Form("runPeakSumCh%iFile%i", ichan, ihist);
      qpeSumDir->GetObject(cloneName, hclone);
      if (hclone)
      {
        cout << ".... ichan " << ichan << " file " << ihist << " " << hclone->GetName() << endl;
      }
      else
      {
        cout << ".... ichan " << ichan << " file " << ihist << "  not found " << endl;
        continue;
      }

      // fit QPE
      double xlow = 0.;
      double xhigh = 0;
      double par1 = 1;
      double epar1 = 0;
      double par2 = 0;
      double epar2 = 0;
      if (theDataType == SIS)
      { // For SIS data
        if (ichan == 0 || ichan == 1 || ichan == 2 || ichan == 5)
        {
          xlow = 120.;
          xhigh = 500.;
        }
        else if (ichan == 4 || ichan == 6 || ichan == 7 || ichan == 8)
        {
          xlow = 120.;
          xhigh = 500.;
        }
        else if (ichan == 12)
        {
          xlow = 120.;
          xhigh = 500.;
        }
      }
      else
      { // For CAEN data
        xlow = 200.;
        xhigh = 300.;
      }
      // if (ichan == 9 || ichan == 10 || ichan == 11)
      //{
      //   continue;
      // }
      TF1 *gfit = NULL;
      bool trigger = ichan == 9 || ichan == 10 || ichan == 11;
      if (!trigger)
      {
        hclone->Fit("gaus", "Q", " ", xlow, xhigh);
        gfit = (TF1 *)hclone->GetListOfFunctions()->FindObject("gaus");
      }

      if (gfit)
      {
        par1 = gfit->GetParameter(1);
        epar1 = gfit->GetParError(1);
        if (epar1 > par1)
          epar1 = par1;
        par2 = gfit->GetParameter(2);
        epar2 = gfit->GetParError(2);
        if (epar2 > par2)
          epar2 = par2;
      }
      else
        printf("Fit to chan %u file  %u fails\n", ichan, ihist);
      if (isnan(par1) || isinf(par1) || par1 > 1.E9 || par1 <= 0)
      {
        par1 = 1;
        epar1 = 0;
        par2 = 0;
        epar2 = 0;
      }
      hQPEChan->SetBinContent(ichan, par1);
      hQPEChan->SetBinError(ichan, epar1);
      hQPESigmaChan->SetBinContent(ichan, par2);
      hQPESigmaChan->SetBinError(ichan, epar2);
      vecQPE[ichan].push_back(par1);
      vecEQPE[ichan].push_back(epar1);
      vecQPESigma[ichan].push_back(par2);
      vecEQPESigma[ichan].push_back(epar2);
      if (ichan == 12)
        printf(" fit QPE chan %u  file %u mean %f %f  sigma %f %f  \n", ichan, ihist, par1, epar1, par2, epar2);
    }
  }
}

void setTimeGraph(TMultiGraph *mg, TString ylabel)
{
  mg->GetXaxis()->SetTimeDisplay(1);
  mg->GetXaxis()->SetNdivisions(-205);
  mg->GetXaxis()->SetTimeFormat("%m-%d-%H");
  mg->GetXaxis()->SetTimeOffset(0, "gmt");
  mg->GetYaxis()->SetTitle(ylabel);
  mg->GetXaxis()->SetTitle("Month-Day-Hour");
}

Int_t get_month_index(TString name)
{
  int imonth;
  if (name.EqualTo("Jan"))
    imonth = 1;
  if (name.EqualTo("Feb"))
    imonth = 2;
  if (name.EqualTo("Mar"))
    imonth = 3;
  if (name.EqualTo("Apr"))
    imonth = 4;
  if (name.EqualTo("May"))
    imonth = 5;
  if (name.EqualTo("Jun"))
    imonth = 6;
  if (name.EqualTo("Jul"))
    imonth = 7;
  if (name.EqualTo("Aug"))
    imonth = 8;
  if (name.EqualTo("Sep"))
    imonth = 9;
  if (name.EqualTo("Oct"))
    imonth = 10;
  if (name.EqualTo("Nov"))
    imonth = 11;
  if (name.EqualTo("Dec"))
    imonth = 12;
  return imonth;
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

/*
*****  all pointers found using getPointers()
*/
void fileLoop()
{
  // must have output file opened
  if (!fout)
    return;

  totalEvents = 0;
  totalPass = 0;
  nFiles = 0;
  filenum.clear();
  efilenum.clear();
  printf("+++++++ fileLoop over %lld files +++++ \n", maxFiles);
  // DEF would be nice to put in a way to look at the last maxFiles files

  for (unsigned ifile = 0; ifile < maxFiles; ++ifile)
  {
    TString fullName = dirNameSlash + fileList[ifile];
    // printf("line479 %s\n", fullName.Data());
    fin = new TFile(fullName, "readonly");

    /* pick up cut failures */
    TH1D *hEventFail = NULL;
    fin->GetObject("EventFail", hEventFail);
    if (hEventFail)
    {
      vecFile.push_back(double(ifile));
      // printf("line539 file %i %s EventFail %.0f \n", ifile, fin->GetName(), hEventFail->GetEntries());
      for (int ibin = 0; ibin < hEventFail->GetNbinsX(); ++ibin)
      {
        vecFail[ibin].push_back(hEventFail->GetBinContent(ibin + 1));
        // printf("bin %i contents %f \n", ibin + 1, hEventFail->GetBinContent(ibin + 1));
      }
    }
    else
    {
      printf("line543 file %i %s has no EventFail!! \n", ifile, fin->GetName());
    }

    // for summing
    hRunTrigSumNoCut = nullptr;
    hRunTrigSumCut = nullptr;
    // get ponters
    if (!getPointers(fin)) // get pointers for this file
      continue;

    /* get statistics */
    ntotal = hEventPass->GetEntries();
    npass = hEventPass->GetBinContent(0);
    filePass.push_back(int(npass));
    fileTotal.push_back(int(ntotal));

    // 0 = ntriggers, 1 = npass
    // printf("\n\n line485 total events file %i %s %f events passed  %f \n ", ifile, fileList[ifile].Data(), eventCount->GetBinContent(0), eventCount->GetBinContent(1));
    // TH1D *hevcount = (TH1D *)eventCount->Clone(Form("eventCount%i", ifile));
    // fout->Add(hevcount);

    filenum.push_back(double(ifile));
    efilenum.push_back(0);
    TDatime dateTime = getTime(ifile);
    fileDatime.push_back(dateTime);
    fileTime.push_back(dateTime.Convert());
    printf(" \n ***** starting file %i , %lu  %s  pass %.0f *******\n", ifile, filenum.size(), fin->GetName(), npass);

    /* clone EventPass
    TString cloneName;
    TString fileTag(fileList[ifile](13, 9));
    cloneName.Form("EventPass%uDate%s", ifile, fileTag.Data());
    TH1D *hClone = (TH1D *)hEventPass->Clone(cloneName);
    printf("line521 %s \n", hClone->GetName());
    hClone->SetTitle(cloneName);
    fout->Add(hClone);
    */

    /* clone EventPass */
    TString OutName;
    TString InName;
    OutName.Form("EventPassSum");
    InName.Form("EventPsss");
    TH1D *hIn = NULL;
    fin->GetObject(InName, hIn);
    if (hIn)
    {
      TH1D *hOut = NULL;
      fout->GetObject(OutName, hOut);
      if (hOut == NULL)
      {
        hOut = (TH1D *)hIn->Clone(OutName);
        hOut->SetTitle(OutName);
        fout->Add(hOut);
      }
      else
      {
        fout->GetObject(OutName, hOut);
        hOut->Add(hIn);
      }
    } // if hIn

    /* clone hGammaPeak hGammaPeakCut*/
    TString gammaOutName;
    TString gammaInName;
    gammaOutName.Form("GammaPeakSum");
    gammaInName.Form("GammaPeak");
    saveToOutput(fout, gammaInName, gammaOutName);

    /* clone hGammaPeak hGammaPeakCut*/
    gammaOutName.Form("GammaPeakCutSum");
    gammaInName.Form("GammaPeakCut");
    saveToOutput(fout, gammaInName, gammaOutName);

    /* clone hGammaPeak hTrianglePlot*/
    gammaOutName.Form("TriangleSum");
    gammaInName.Form("Triangle");
    saveToOutput(fout, gammaInName, gammaOutName);

    gammaOutName.Form("TriangleCutSum");
    gammaInName.Form("TriangleCut");
    saveToOutput(fout, gammaInName, gammaOutName);
    /*
    fin->GetObject(gammaInName, hIn);
    if (hIn)
    {
      TH1D *hOut = NULL;
      fout->GetObject(gammaOutName, hOut);
      if (hOut == NULL)
      {
        hOut = (TH1D *)hIn->Clone(gammaOutName);
        hOut->SetTitle(gammaOutName);
        // cout << "line613 ... adding  " << hIn->GetName() << " file "  << fin->GetName() << " hit QSum " << hOut->GetEntries() << endl;
        fout->Add(hOut);
      }
      else
      {
        fout->GetObject(gammaOutName, hOut);
        hOut->Add(hIn);
      }
    } // if hIn
     */

    /* clone hGammaPeak hGammaPeakCut*/
    /*
    gammaOutName.Form("GammaPeakCutSum");
    gammaInName.Form("GammaPeakCut");
    fin->GetObject(gammaInName, hIn);
    if (hIn)
    {
      TH1D *hOut = NULL;
      fout->GetObject(gammaOutName, hOut);
      if (hOut == NULL)
      {
        hOut = (TH1D *)hIn->Clone(gammaOutName);
        hOut->SetTitle(gammaOutName);
        // cout << "line613 ... adding  " << hIn->GetName() << " file "  << fin->GetName() << " hit QSum " << hOut->GetEntries() << endl;
        fout->Add(hOut);
      }
      else
      {
        fout->GetObject(gammaOutName, hOut);
        hOut->Add(hIn);
      }
    } // if hIn
    */

    /******
     * add peak and sum gains  to gainSumDir
     * ******/
    /**************** add peak ****************/
    //  loop over channels
    for (int ichan = 0; ichan < NONSUMCHANNELS; ++ichan)
    {
      TString gainInName;
      TString gainOutName;
      gainOutName.Form("GainPeakChan%i", ichan);
      gainInName.Form("QPeakChan%i", ichan);
      TH1D *hIn = NULL;
      sumDir->GetObject(gainInName, hIn);
      if (hIn)
      {
        TH1D *hOut = NULL;
        gainSumDir->GetObject(gainOutName, hOut);
        if (hOut == NULL)
        {
          hOut = (TH1D *)hIn->Clone(gainOutName);
          hOut->SetTitle(gainOutName);
          // cout << "line516 ... adding  " << hIn->GetName() << " file "
          //      << fin->GetName() << " hit QPeak " << hOut->GetEntries() << endl;
          gainSumDir->Add(hOut);
        }
        else
        {
          gainSumDir->GetObject(gainOutName, hOut);
          hOut->Add(hIn);
          // cout << "line526 ... found for " << gainOutName
          //      << " in entries " << hIn->GetEntries() << " hit QPeak " << hOut->GetEntries() << endl;
        }
      } // if hIn
    } // channel loop

    /***************** add qsum ****************/
    //  loop over channels
    for (int ichan = 0; ichan < NONSUMCHANNELS; ++ichan)
    {
      TString gainOutName;
      TString gainInName;
      gainOutName.Form("GainSumChan%i", ichan);
      gainInName.Form("QSumChan%i", ichan);
      TH1D *hIn = NULL;
      sumDir->GetObject(gainInName, hIn);
      if (hIn)
      {
        TH1D *hOut = NULL;
        gainSumDir->GetObject(gainOutName, hOut);
        if (hOut == NULL)
        {
          hOut = (TH1D *)hIn->Clone(gainOutName);
          hOut->SetTitle(gainOutName);
          // cout << "line613 ... adding  " << hIn->GetName() << " file "  << fin->GetName() << " hit QSum " << hOut->GetEntries() << endl;
          gainSumDir->Add(hOut);
        }
        else
        {
          gainSumDir->GetObject(gainOutName, hOut);
          hOut->Add(hIn);
          // cout << "line621 ... found for " << gainOutName
          //      << " in entries " << hIn->GetEntries() << " hit QSum " << hOut->GetEntries() << endl;
        }
      } // if hIn
    } // channel loop

    // TrigSumNoCut
    // make summed histos on output
    if (!hRunTrigSumNoCut && hTrigSumNoCut)
    {
      printf("line550 GOT hTrigSumNoCut IN FILE %s\n", fin->GetName());
      hRunTrigSumNoCut = (TH1D *)hTrigSumNoCut->Clone("RunTrigSumNoCut");
      fout->Add(hRunTrigSumNoCut);
      printf("hRunTrigSumCut IN FILE %s named %s\n", fin->GetName(), hRunTrigSumNoCut->GetName());
    }

    if (hRunTrigSumCut && hTrigSumCut)
    {
      printf("line562 hRunTrigSumCut IN FILE named %s \n", fin->GetName());
      hRunTrigSumCut = (TH1D *)hTrigSumCut->Clone("RunTrigSumCut");
      fout->Add(hRunTrigSumCut);
      printf("hRunTrigSumCut IN FILE %s named %s\n", fin->GetName(), hRunTrigSumCut->GetName());
    }

    // sum with Add
    if (hRunTrigSumNoCut && hTrigSumNoCut)
    {
      // printf("line571 hRunTrigSumNOCut IN FILE %s %s add to %s \n", fin->GetName(), hTrigSumNoCut->GetName(), hRunTrigSumNoCut->GetName());
      hRunTrigSumNoCut->Add(hTrigSumNoCut);
    }

    if (hRunTrigSumCut && hTrigSumCut)
    {
      // printf("line573 hRunTrigSumCut IN FILE %s %s\n", fin->GetName(), hTrigSumCut->GetName());
      hRunTrigSumCut->Add(hTrigSumCut);
    }

    /******  loop over sumDir *****/
    TList *sumList = sumDir->GetListOfKeys();
    TIter next(sumList);
    TKey *key;
    // printf("line714 in fileLoop addsumDirHistos %u \n", sumList->GetEntries());
    //  sumDir->GetListOfKeys()->ls();
    while (TKey *key = (TKey *)next())
    {
      TKey *keyprev = NULL;
      keyprev = (TKey *)sumDir->GetListOfKeys()->Before(key);
      if (keyprev && ((key->GetName(), keyprev->GetName()) == 0))
        continue;

      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TH1D"))
        continue;
      TH1D *h = (TH1D *)key->ReadObj();
      std::string name = string(h->GetName());
      // cout << " name " << name << "  " << name.find("sumWave") << " " << name.find("Bad") << " npos " << std::string::npos << endl;

      /****** get run sum ****/
      TH1D *hClone;
      TH1D *hRunClone;
      if (name.find("sumWave") != std::string::npos && name.find("Bad") == std::string::npos && name.find("All") == std::string::npos && name.find("Fail") == std::string::npos)
      {
        // waveSumDir->cd();
        // cout << "line700  sumWave clone " << name << " file " << ifile << endl;
        string chan = name.substr(name.find_last_of("e") + 1);
        int ichan = stoi(chan);
        TString cloneName;
        cloneName.Form("RunSumWaveFile%uChan%i", ifile, ichan);
        hClone = (TH1D *)h->Clone(cloneName);
        hClone->SetTitle(cloneName);
        hSumWave.push_back(hClone);
        waveSumDir->Add(hClone);
      }
      if (name.find("sumPeakWave") != std::string::npos)
      {
        // waveSumDir->cd();
        string chan = name.substr(name.find_last_of("e") + 1);
        int ichan = stoi(chan);
        // printf("line760 sumPeakWave clone %s file %i chan %i\n", name.c_str(), ifile, ichan);
        TString cloneName;
        cloneName.Form("RunPeakWaveFile%uChan%i", ifile, ichan);
        hClone = (TH1D *)h->Clone(cloneName);
        hClone->SetTitle(cloneName);
        // TString histName;
        //  histName.Form("RunPeakWaveChan%i", ichan);
        //  hRunClone = (TH1D *)h->Clone(histName);
        waveSumDir->Add(hClone);
        vRunPeakWave[ichan].push_back(hClone);
        /* do not sum here */
      }

      /****** get PeakChan by channel ****/
      if (name.find("QPeakChan") != std::string::npos)
      {
        string chan = name.substr(name.find_last_of("n") + 1);
        int ichan = stoi(chan);
        // cout << "line739 " << name << "string chan" << chan << " int " << ichan << endl;

        TString cloneName;
        cloneName.Form("runPeakSumCh%ifile%i", ichan, ifile);
        TH1D *hpAdd = (TH1D *)h->Clone(cloneName);
        // cout << " line744 found RunQSum " << ichan << " " << fileList[ifile] << " " << h->GetName() << " " << hpAdd->GetName() << endl;
        hpAdd->SetMarkerStyle(20);
        hpAdd->SetMarkerSize(0.5);
        vRunPeakSum[ichan].push_back(hpAdd);
        qpeSumDir->Add(hpAdd);
      }
      // get QSumChan by channel
      if (name.find("QSumChan") != std::string::npos)
      {
        string chan = name.substr(name.find_last_of("n") + 1);
        int ichan = stoi(chan);
        TString cloneName;
        cloneName.Form("runQSumCh%iFile%i", ichan, ifile);
        TH1D *hqAdd = (TH1D *)h->Clone(cloneName);
        /*
        cout << " line758 found RunQSum " << ichan << " " << fileList[ifile] << " " << h->GetName() << " " << hqAdd->GetName()
             << " max x  " << h->GetBinLowEdge(h->GetNbinsX() + 1) << endl;
             */
        hqAdd->SetMarkerStyle(20);
        hqAdd->SetMarkerSize(0.5);
        vRunQSum[ichan].push_back(hqAdd);
        qpeSumDir->Add(hqAdd);
      }
      /****** get SingletShape by channel ****/
      if (name.find("SingletShape") != std::string::npos)
      {
        // does sum exist on output file? in directory runSumDir
        string chan = name.substr(name.find_last_of("e") + 1);
        int ichan = stoi(chan);
        TString cloneName;
        cloneName.Form("singletShapeRunSumCh%i", ichan);
        // does sum exist on output file? in directory runSumDir
        TH1D *hpAdd = NULL;
        runSumDir->GetObject(cloneName, hpAdd);
        if (!hpAdd)
        {
          TH1D *hpAdd = (TH1D *)h->Clone(cloneName);
          hpAdd->SetMarkerStyle(20);
          hpAdd->SetMarkerSize(0.5);
          qpeSumDir->Add(hpAdd);
        }
        else
        { // already exits so sum this histo
          hpAdd->Add(h);
        }
      }
    } // sum over keys

    fout->Write();
    // if(ifile==0) fout->ls();
    // waveSumDir->Write();
    totalEvents += int(ntotal);
    totalPass += int(npass);
    vTotal.push_back(ntotal);
    vPass.push_back(npass);

    ++nFiles;
    // close file
    fin->Close();
    printf("line812 end loop over sumDir keys file %i of %i named %s events %i file pass %i totalPass %i WaveSumDir keys %i \n", ifile, nFiles, fin->GetName(), int(ntotal), int(npass), totalPass, waveSumDir->GetNkeys());
  } // end loop over files

  printf("line865 end of fileLoop finished over %lld good files %d  totalEvents %i totalPass %i waveSumDir has %d \n", maxFiles, nFiles, totalEvents, totalPass, waveSumDir->GetNkeys()); // DEF would be nice to put in a way to look at the last maxFiles files
}

void sumHistosChannel(int ichan, TString histSet)
{
  // printf("line921 chan %i set %s nFiles %i \n", ichan, histSet.Data(), nFiles);
  //  sum over files
  for (int ih = 0; ih < nFiles; ++ih)
  {
    // printf("line925 in sumHistosChannel ichan %i set %s ifile %i \n", ichan, histSet.Data(), ih);
    //  waveToSum already in output file by run and channel
    TString histName;
    histName.Form("Run%sFile%uChan%i", histSet.Data(), ih, ichan);
    TH1D *waveToSum = NULL;
    waveSumDir->GetObject(histName, waveToSum);
    if (waveToSum == NULL)
    {
      printf("line856 skipping %s  %s chan %i file %i \n", histSet.Data(), histName.Data(), ichan, ih);
      continue;
    }
    printf("line936 at %s  %s chan %i file %i \n", histSet.Data(), histName.Data(), ichan, ih);
    cout << "line937  waveToSum channel " << ih << " file "
         << waveToSum->GetName() << " passing for file  " << filePass[ih] << endl;

    // new histogram
    int nbinsx = waveToSum->GetNbinsX();
    double xlow = waveToSum->GetXaxis()->GetBinLowEdge(0);
    double xup = waveToSum->GetXaxis()->GetBinUpEdge(nbinsx);

    // normalize to qpe.
    fitSumDir->cd();
    TH1D *hWaveToFit;
    TH1D *hWaveToFitNotNormed;
    TString fitwaveName;
    fitwaveName.Form("fitwave%sChan%iFile%i", histSet.Data(), ichan, ih);
    TString notNormedName;
    notNormedName.Form("notNormedWave%sChan%iFile%i", histSet.Data(), ichan, ih);
    hWaveToFit = (TH1D *)waveToSum->Clone(fitwaveName);
    hWaveToFitNotNormed = (TH1D *)waveToSum->Clone(notNormedName);
    // fitSumDir->ls();

    // cout << hWaveToFit->GetName() << endl;
    hWaveToFit->SetMarkerStyle(21);
    hWaveToFit->SetMarkerSize(0.2);
    hWaveToFit->GetListOfFunctions()->Clear();
    hWaveToFitNotNormed->SetMarkerStyle(21);
    hWaveToFitNotNormed->SetMarkerSize(0.2);
    hWaveToFitNotNormed->GetListOfFunctions()->Clear();

    for (int ibin = 0; ibin < waveToSum->GetNbinsX(); ++ibin)
    {
      double xbin = waveToSum->GetBinContent(ibin);
      double ebin = waveToSum->GetBinError(ibin);
      hWaveToFitNotNormed->SetBinContent(ibin, xbin);
      hWaveToFitNotNormed->SetBinError(ibin, ebin);
      // norm to total pass in file and to gain
      hWaveToFit->SetBinContent(ibin, xbin / double(filePass[ih]) / readGains->sipmPeakGain[ichan]);
      hWaveToFit->SetBinError(ibin, ebin / double(filePass[ih]) / readGains->sipmPeakGain[ichan]);
    }
    if (histSet.Contains("Peak"))
    {
      vNormByFile.push_back(hWaveToFit);
      // printf("line927 chan %i file %i push back vNormByFile size %lu \n \n", ichan,ih,vNormByFile.size());
    }

    // add and save in output file runSumDir;
    if (histSet == TString("PeakWave"))
    {
      runSumDir->cd();
      if (hUnNormedPeakWave[ichan] == NULL)
      {
        histName.Form("UnNormed%sChan%i", histSet.Data(), ichan);
        hUnNormedPeakWave[ichan] = (TH1D *)hWaveToFitNotNormed->Clone(histName);
        hUnNormedPeakWave[ichan]->SetTitle(histName);
        runSumDir->Add(hUnNormedPeakWave[ichan]);
        vIntegral[ichan].push_back(hWaveToFit->Integral(startTime, endTime));
        vMaxValue[ichan].push_back(hWaveToFit->GetBinContent(hWaveToFit->GetMaximumBin()));
        // printf("line993 file %i %s Int %f peak %f \n", ih, hWaveToFit->GetName(), vIntegral[ichan][vIntegral[ichan].size() - 1], vMaxValue[ichan][vMaxValue[ichan].size() - 1]);
      }
      else
      {
        histName.Form("UnNormed%sChan%i", histSet.Data(), ichan);
        // hRunHitWave[ichan]->Add(hWaveToFit);
        runSumDir->GetObject(histName, hUnNormedPeakWave[ichan]);
        if (hUnNormedPeakWave[ichan] == NULL)
        {
          // printf("line 951 NULL chan %i file %i %s \n", ichan, ih, histName.Data());
          runSumDir->ls();
        }
        hUnNormedPeakWave[ichan]->Add(hWaveToFit);
        // printf("line1003 file %i Int %f peak %f \n", ih, hWaveToFit->Integral(startTime, endTime), hWaveToFit->GetBinContent(hWaveToFit->GetMaximumBin()));
        vIntegral[ichan].push_back(hWaveToFit->Integral(startTime, endTime));
        vMaxValue[ichan].push_back(hWaveToFit->GetBinContent(hWaveToFit->GetMaximumBin()));
        // printf("line1006 file %i %s Int %f peak %f \n", ih, hWaveToFitNotNormed->GetName(), vIntegral[ichan][vIntegral[ichan].size() - 1], vMaxValue[ichan][vMaxValue[ichan].size() - 1]);
      }
      printf("line1011 sumHistos file %i chan %i hist %s integral %f \n", ih, ichan, hUnNormedPeakWave[ichan]->GetName(), hUnNormedPeakWave[ichan]->Integral(startTime, endTime));
    }
    else if (histSet == TString("SumWave"))
    {
      runSumDir->cd();
      histName.Form("UnNormed%sChan%i", histSet.Data(), ichan);
      if (hUnNormedSumWave[ichan] == NULL)
      {
        hUnNormedSumWave[ichan] = (TH1D *)waveToSum->Clone(histName);
        hUnNormedSumWave[ichan]->SetTitle(histName);
      }
      else
      {
        histName.Form("UnNormed%sChan%i", histSet.Data(), ichan);
        runSumDir->GetObject(histName, hUnNormedSumWave[ichan]);
        if (hUnNormedSumWave[ichan] == NULL)
          printf("line950 no RunSumWave %s\n", histName.Data());
        hUnNormedSumWave[ichan]->Add(hWaveToFit);
      }
    }
  } // sum over files
}

void sumHistos()
{
  // histograms time in ns
  // make sum histograms on output file in directory runSumDir
  runSumDir->cd();
  TString histName;
  for (unsigned ichan = 0; ichan < CHANNELS; ++ichan)
  {
    histName.Form("Run%sChan%i", "PeakWave", ichan);
    hRunPeakWave[ichan] = new TH1D(histName, histName, waveBins, 0, 2. * waveBins);
    hRunPeakWave[ichan]->GetXaxis()->SetTitle("time [ns]");
    hRunPeakWave[ichan]->GetYaxis()->SetTitle("yield [SPE] ");
    histName.Form("Run%sChan%i", "SumWave", ichan);
    hRunSumWave[ichan] = new TH1D(histName, histName, waveBins, 0, 2. * waveBins);
    hRunSumWave[ichan]->GetXaxis()->SetTitle("time [ns]");
    hRunSumWave[ichan]->GetYaxis()->SetTitle("yield [SPE] ");
  }

  printf("line1013 sumHistos: Nfiles %d Number of waveSumDir %d\n", nFiles, waveSumDir->GetList()->GetEntries());
  // loop over channels
  for (int ichan = 0; ichan < CHANNELS; ++ichan)
    sumHistosChannel(ichan, TString("PeakWave"));
  for (int ichan = 0; ichan < CHANNELS; ++ichan)
    sumHistosChannel(ichan, TString("SumWave"));
  // normalize each channel
  // waveSumDir->ls();
  normalizeTotalPass(TString("PeakWave"));
  normalizeTotalPass(TString("SumWave"));
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

int main(int argc, char *argv[])
{
  cout << "executing " << argv[0] << " make summary plots  " << endl;
  printf(" usage: summary start date string <stag> end date string <etag> max files <default all> \n ");
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

  readGains = new TReadGains();
  dirName = TString("caenData");
  dirNameSlash = TString("caenData/");
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

  fout = new TFile(Form("summary-%s-nfiles-%lld-created-%s.root", tag.Data(), maxFiles, sdate.c_str()), "recreate");
  fitSumDir = fout->mkdir("fitSumDir");
  waveSumDir = fout->mkdir("waveSumDir");
  qpeSumDir = fout->mkdir("qpeSumDir");
  runSumDir = fout->mkdir("runSumDir");
  gainSumDir = fout->mkdir("gainSumDir");
  fout->cd();

  hQPEChan = new TH1D("QPEChan", "QPE  by channel", 12, 0, 12);
  hQPESigmaChan = new TH1D("QPESigmaChan", "QPE  by channel", 12, 0, 12);
  hQPEChan->Sumw2();
  hQPESigmaChan->Sumw2();
  vMaxValue.resize(CHANNELS);
  vIntegral.resize(CHANNELS);
  vecQsum.resize(CHANNELS);
  vecEQsum.resize(CHANNELS);
  vecQsumUn.resize(CHANNELS);
  vecEQsumUn.resize(CHANNELS);
  vecPeaksum.resize(CHANNELS);
  vecEPeaksum.resize(CHANNELS);
  vecQPE.resize(CHANNELS);
  vecQPEMean.resize(CHANNELS);
  vecQPENave.resize(CHANNELS);
  vecEQPE.resize(CHANNELS);
  vecPeak.resize(CHANNELS);
  vecEPeak.resize(CHANNELS);
  vecQPESigma.resize(CHANNELS);
  vecEQPESigma.resize(CHANNELS);
  vSlope.resize(CHANNELS);
  vESlope.resize(CHANNELS);
  vRunQSum.resize(CHANNELS);
  vRunPeakSum.resize(CHANNELS);
  vRunPeakWave.resize(CHANNELS);

  hRunPeakWave.resize(CHANNELS);
  hRunSumWave.resize(CHANNELS);
  hUnNormedPeakWave.resize(CHANNELS);
  hUnNormedSumWave.resize(CHANNELS);
  sumHits.resize(CHANNELS);
  hRunLatePeakSum.resize(CHANNELS);

  fout->cd();

  fileLoop();
  printf("line1284 \t\t >>> after fileLoop processed << %li  total pass %i channels %lu <<<<< \n", filenum.size(), totalPass, vRunPeakWave.size());

  for (unsigned jfile = 0; jfile < filenum.size(); ++jfile)
  {
    printf(" file %i %s total %i pass %i pass frac %.3f \n", int(filenum[jfile]), fileList[jfile].Data(),
           fileTotal[jfile], filePass[jfile], double(filePass[jfile]) / double(fileTotal[jfile]));
  }
  fout->Write();

  // call function to fit slopes and fill vSlope, vESlope
  if (filenum.size() > 0)
  {
    cout << "line1357 call sumHistos files " << filenum.size() << " total pass " << totalPass << endl;
    sumHistos();
  }
  // endNow:

  // print totalHits
  printf("line1367 \t >>> end of job: files processed << %li  total pass %i <<<<< \n", filenum.size(), totalPass);

  /* runSum integrals of waveSuDir */
  printf("line1416 from Directory %s with %i keys do integrals: \n", waveSumDir->GetName(), waveSumDir->GetNkeys());
  // waveSumDir->ls();
  runSums.resize(NONSUMCHANNELS);
  runSumNames.resize(NONSUMCHANNELS);
  TList *fitList = waveSumDir->GetListOfKeys();
  TIter next2(fitList);
  while (TKey *key = (TKey *)next2())
  {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1D"))
      continue;
    TH1D *h = (TH1D *)key->ReadObj();

    TString sname = h->GetName();
    int chanNumber = TString(sname(sname.Last('n') + 1, sname.Length())).Atoi();
    if (chanNumber > 12)
      continue;
    int fileNumber = TString(sname(sname.Last('e') + 1, sname.Length())).Atoi();
    double inte = h->Integral();
    // printf("line1435 file %i hist %s chan %i %s pass %i integral %.3E \n", fileNumber, h->GetName(), chanNumber, h->GetName(), filePass[fileNumber], inte);
    runSums[chanNumber].push_back(inte);
    runSumNames[chanNumber].push_back(h->GetName());
  }

  // calculate mean hits from waveforms
  /*
  printf("line1461 calculate mean hits from hRunPeakWave %lu from runSumDir \n", hRunPeakWave.size());
  for (int idet = 0; idet < hRunPeakWave.size(); ++idet)
  {
    TString histName;
    histName.Form("UnNormedPeakWaveChan%i", idet);
    TH1D *hist;
    runSumDir->GetObject(histName, hist);
    if (hist)
      printf(" idet %i unnormed %.3E inte %.3E \n", idet, hist->Integral(), hRunPeakWave[idet]->Integral());
    // else printf(" did not find %s \n",histName.Data());
  }
    */

  /* report cleanup cut failures */
  printf("line1256 vecFail size %lu \n", vecFail[0].size());
  std::vector<std::vector<double>> normFailures; // normalized to number of events
  normFailures.resize(FAILBITS);
  // std::vector<std::vector<double>> normFailuresError;
  for (unsigned ifile = 0; ifile < vecFail[0].size(); ++ifile)
  {
    printf("file %i PASS %.0f BASEFAIL %.0f EARLYCUT %.0f FIRSTTIME %.0f COSMIC %.0f GAMMA %.0f TRIGFAIL %.0f \n ", ifile, vecFail[0][ifile], vecFail[1][ifile], vecFail[2][ifile], vecFail[3][ifile], vecFail[4][ifile], vecFail[5][ifile], vecFail[6][ifile]);
    double sum = 0;
    for (int icode = 0; icode < vecFail.size(); ++icode)
      sum += vecFail[icode][ifile];
    printf("line1277 totals for file %i %.0f total pass %i \n", ifile, sum, totalPass);

    for (int icode = 0; icode < FAILBITS; ++icode)
    {
      normFailures[icode].push_back(vecFail[icode][ifile] / sum);
    }
  }

  printf("line1284 number of files = normFailureds[0] size %lu \n", normFailures[0].size());

  /*
  for (unsigned ifile = 0; ifile < vecFail[0].size(); ++ifile)
  {
    // printf("file %i NORMED PASS %.0f BASEFAIL %.0f EARLYCUT %.0f FIRSTTIME %.0f COSMIC %.0f GAMMA %.0f TRIGFAIL %.0f \n ", ifile, normFailures[0][ifile], normFailures[1][ifile], normFailures[2][ifile], normFailures[3][ifile], normFailures[4][ifile], normFailures[5][ifile], normFailures[6][ifile]);
    for (int icode = 0; icode < FAILBITS; ++icode)
    {
      printf("line1265 ifile %i code %i  %s normed %f \n", ifile, icode, bitNames[icode].Data(), normFailures[icode][ifile]);
    }
  }
    */

  /* make vector normed with errors */
  TGraph *gFailures[FAILBITS];
  for (unsigned icode = 0; icode < FAILBITS; ++icode)
  {
    gFailures[icode] = new TGraph(normFailures[icode].size(), &vecFile[0], &normFailures[icode][0]);
    gFailures[icode]->SetName(Form("%s-%s", bitNames[icode].Data(), tag.Data()));
    gFailures[icode]->SetTitle(Form("%s-%s", bitNames[icode].Data(), tag.Data()));
    gFailures[icode]->GetHistogram()->GetXaxis()->SetTitle("file number");
    gFailures[icode]->GetHistogram()->GetYaxis()->SetTitle("failure fraction");
    gFailures[icode]->SetMarkerColor(kGreen);
    gFailures[icode]->SetMarkerStyle(20);
    fout->Add(gFailures[icode]);
  }

  printf(" nfiles %lu \n", vMaxValue[0].size());
  for (int ifile = 0; ifile < vTotal.size(); ++ifile)
  {
    printf(" file %i total %i pass %i \n", ifile, vTotal[ifile], vPass[ifile]);

    for (int ich = 0; ich < NONSUMCHANNELS; ++ich)
      printf("\t file %i  channel %i MaxValue %f integral %f \n", ifile, ich, vMaxValue[ich][ifile], vIntegral[ich][ifile]);
  }

  /* make MaxValue graphs */
  printf("line1345 make MaxValue file %lu %lu\n", vecFile.size(), vMaxValue[0].size());
  TGraph *gMaxValue[NONSUMCHANNELS];
  for (unsigned ich = 0; ich < NONSUMCHANNELS; ++ich)
  {
    gMaxValue[ich] = new TGraph(vMaxValue[ich].size(), &vecFile[0], &vMaxValue[ich][0]);
    gMaxValue[ich]->SetName(Form("MaxValueChan%i-%s", ich, tag.Data()));
    gMaxValue[ich]->SetTitle(Form("MaxValueChan%i-%s", ich, tag.Data()));
    gMaxValue[ich]->GetHistogram()->GetXaxis()->SetTitle("file number");
    gMaxValue[ich]->GetHistogram()->GetYaxis()->SetTitle("Peak MaxValue [SPE]");
    gMaxValue[ich]->SetMarkerColor(kRed);
    gMaxValue[ich]->SetMarkerStyle(21);
    fout->Add(gMaxValue[ich]);
  }

  /* make Integral graphs */
  TGraph *gIntegral[NONSUMCHANNELS];
  for (unsigned ich = 0; ich < NONSUMCHANNELS; ++ich)
  {
    gIntegral[ich] = new TGraph(vIntegral[ich].size(), &vecFile[0], &vIntegral[ich][0]);
    gIntegral[ich]->SetName(Form("IntegralChan%i-%s", ich, tag.Data()));
    gIntegral[ich]->SetTitle(Form("IntegralChan%i-%s", ich, tag.Data()));
    gIntegral[ich]->GetHistogram()->GetXaxis()->SetTitle("file number");
    gIntegral[ich]->GetHistogram()->GetYaxis()->SetTitle("Peak Integral [SPE]");
    gIntegral[ich]->SetMarkerColor(kBlue);
    gIntegral[ich]->SetMarkerStyle(22);
    fout->Add(gIntegral[ich]);
  }

  cout << "line1474 summary finished "
       << " total pass " << totalPass << " maxFiles  " << maxFiles << " good file " << nfiles << " files written to " << fout->GetName() << endl;

  fout->Purge(1);
  fout->Write();
  fout->Close();

  exit(0);
}
/*  put this all at bottom */
void makeGraphs()
{
  printf(" \t\t make graphs %lu \n", filenum.size());
  cout << " \n\t ******* makeGraphs *****  " << endl;
  cout << " \n\t vecQsum " << vecQsum.size() << " vecQPE  " << vecQPE.size() << endl;
  if (vecQsum.size() == 0)
    return;
  int myColor[13] = {41, 42, 43, 44, 45, 46, 2, 3, 4, 31, 32, 33, 34};
  int myStyle[13] = {21, 22, 23, 24, 25, 26, 21, 22, 23, 31, 32, 33, 34};

  // normalize to first file

  normQsum.resize(vecQsum.size());
  for (unsigned ic = 0; ic < vecQsum.size(); ++ic)
  {
    // printf(" vecQsum %i %lu \n", ic, vecQsum[ic].size());
    normQsum[ic] = 1;
    if (vecQsum[ic].size() > 0)
    { // ave over before doping
      double beforeSum = 0;
      int normCount = 0;
      if (vecQsum[ic].size() > 20)
      {
        for (unsigned jt = 0; jt < 20; ++jt)
        {
          if (!isinf(vecQsum[ic][jt]) && vecQsum[ic][jt] > 0 && fileDatime[jt].Convert() < dateTime.Convert())
          {
            beforeSum += vecQsum[ic][jt];
            ++normCount;
          }
        }
      }
      if (normCount > 0)
        normQsum[ic] = beforeSum / double(normCount);
    }
    // printf("\t  normQsum =  %f  \n", normQsum[ic]);
  }

  normQPE.resize(vecQPE.size());
  for (unsigned ic = 0; ic < vecQPE.size(); ++ic)
  {
    // printf(" vecQPE %i %lu \n", ic, vecQPE[ic].size());
    normQPE[ic] = 1;
    if (vecQPE[ic].size() > 0)
    {
      printf("\t vecQPE =  %f  \n", vecQPE[ic][0]);
      if (!isinf(vecQPE[ic][0]) && vecQPE[ic][0] > 0)
        normQPE[ic] = vecQPE[ic][0];
    }
    // printf("\t normQPE =  %f  \n", normQPE[ic]);
  }

  // slope graphs
  // one graph per channel
  TString ylabel;
  vector<TGraphErrors *> graphSlope;
  TMultiGraph *mgslope = new TMultiGraph();
  for (unsigned ic = 0; ic < 9; ++ic)
  {
    cout << " add " << ic << " size " << vSlope[ic].size() << " size " << vESlope[ic].size() << endl;
    cout << "     " << ic << " size " << fileTime.size() << " size " << efilenum.size() << endl;
    graphSlope.push_back(new TGraphErrors(filenum.size(), &fileTime[0], &(vSlope[ic][0]), &efilenum[0], &(vESlope[ic][0])));
    unsigned ilast = graphSlope.size() - 1;
    graphSlope[ilast]->SetName(Form("slopeChan%i", ic));
    graphSlope[ilast]->SetTitle(Form("slope-chan-%i", ic));
    graphSlope[ilast]->SetMarkerSize(1);
    graphSlope[ilast]->SetMarkerColor(myColor[ic]);
    graphSlope[ilast]->SetMarkerStyle(myStyle[ic]);
    fout->Add(graphSlope[ilast]);
    if (ic != 5 && ic != 6 && ic != 3 && ic < 9)
      mgslope->Add(graphSlope[ilast]);
  }
  ylabel.Form(" fitted Lifetime [microsec] ");
  setTimeGraph(mgslope, ylabel);
  TCanvas *canSlope = new TCanvas(Form("SlopeFit-%s", sdate.c_str()), Form("SlopeFit-%s", sdate.c_str()));
  mgslope->Draw("ap");
  mgslope->GetYaxis()->SetRangeUser(0.5, 2.5);
  gPad->Update();
  canSlope->BuildLegend();
  canSlope->SetGrid();
  canSlope->Print(".png");
  fout->Append(canSlope);

  // one graph per channel
  for (unsigned ic = 0; ic < vecQsum.size(); ++ic)
  {
    printf("QSUM ch %i size %lu \n", ic, vecQsum[ic].size());
    for (unsigned ih = 0; ih < vecQsum[ic].size(); ++ih)
    {
      double qpe = vecQPE[ic][ih];
      if (qpe <= 1. || isnan(qpe) || isinf(qpe))
        qpe = 1.;
      printf(" \t\t QSUM chan %i file%i qpe %.3E  qsum %f  new  %f \n",
             ic, ih, vecQPE[ic][ih], vecQsum[ic][ih], vecQsum[ic][ih] / qpe);
      vecQsum[ic][ih] = vecQsum[ic][ih] / qpe;
      vecEQsum[ic][ih] = vecEQsum[ic][ih] / qpe;
    }
  }

  cout << " graph normalized vecQsum " << endl;

  TNtuple *ntQsum = new TNtuple("ntQsum", " normalied qsum by run ", "run:q0:q1:q2:q4:q5:q7:q8:q12");
  fout->Append(ntQsum);

  for (unsigned ih = 0; ih < vecQsum[0].size(); ++ih)
    ntQsum->Fill(float(ih), vecQsum[0][ih], vecQsum[1][ih], vecQsum[2][ih], vecQsum[4][ih], vecQsum[5][ih], vecQsum[7][ih], vecQsum[8][ih], vecQsum[12][ih]);

  vector<TGraphErrors *> gqsum;
  TMultiGraph *mgsum = new TMultiGraph();
  for (unsigned ic = 0; ic < vecQsum.size(); ++ic)
  {
    // cout << " add " << ic << endl;
    gqsum.push_back(new TGraphErrors(filenum.size(), &fileTime[0], &(vecQsum[ic][0]), &efilenum[0], &(vecEQsum[ic][0])));
    gqsum[ic]->SetName(Form("qsumChanNorm%i", ic));
    gqsum[ic]->SetTitle(Form("qsum-chan-%i", ic));
    gqsum[ic]->SetMarkerSize(1);
    gqsum[ic]->SetMarkerColor(myColor[ic]);
    gqsum[ic]->SetMarkerStyle(myStyle[ic]);
    fout->Add(gqsum[ic]);
    // if (ic == 6 || ic == 7 || ic == 8)
    if (ic != 5 && ic != 6 && ic != 3 && ic != 9 && ic != 10 && ic != 11)
      mgsum->Add(gqsum[ic]);
  }
  // overlay all channel graphs on canvas
  int ndiv = 10 + 100 * 5 + 10000 * 3;
  ylabel.Form("integrated charge (qpe) / effOther  ");
  setTimeGraph(mgsum, ylabel);
  TCanvas *can = new TCanvas(Form("Qsummary-%s", sdate.c_str()), Form("Qsummary-%s", sdate.c_str()));
  mgsum->Draw("ap");
  gPad->Update();
  can->BuildLegend();
  can->SetGrid();
  can->Print(".png");
  fout->Append(can);
  /**/

  vector<TGraphErrors *> gqsum2;
  TMultiGraph *mgsum2 = new TMultiGraph();
  for (unsigned ic = 0; ic < vecQsum.size(); ++ic)
  {
    // cout << " add " << ic << endl;
    gqsum2.push_back(new TGraphErrors(filenum.size(), &filenum[0], &(vecQsum[ic][0]), &efilenum[0], &(vecEQsum[ic][0])));
    gqsum2[ic]->SetName(Form("qsumChanNorm%i", ic));
    gqsum2[ic]->SetTitle(Form("qsum-chan-%i", ic));
    gqsum2[ic]->SetMarkerSize(1);
    gqsum2[ic]->SetMarkerColor(myColor[ic]);
    gqsum2[ic]->SetMarkerStyle(myStyle[ic]);
    fout->Add(gqsum2[ic]);
    // if (ic == 6 || ic == 7 || ic == 8)
    if (ic != 6 && ic != 3 && ic != 9 && ic != 10 && ic != 11)
      mgsum2->Add(gqsum2[ic]);
  }
  // overlay all channel graphs on canvas
  TCanvas *can2 = new TCanvas(Form("QsummaryByRun-%s", sdate.c_str()), Form("QsummaryByRun-%s", sdate.c_str()));
  mgsum2->Draw("ap");
  gPad->Update();
  can2->BuildLegend();
  can2->SetGrid();
  can2->Print(".png");
  fout->Append(can2);

  // QPE graphs one graph per channel
  vector<TGraphErrors *> gqpe;
  TMultiGraph *mgQPE = new TMultiGraph();
  for (unsigned ic = 0; ic < vecQPE.size(); ++ic)
  {
    // for (int ifile = 0; ifile < vecQPE[ic].size(); ++ifile)
    //   printf(" chan %u file %i %f\n" , ic,ifile,vecQPE[ic][ifile]);
    //  cout << " add " << ic << endl;
    gqpe.push_back(new TGraphErrors(vecQPE[ic].size(), &fileTime[0], &(vecQPE[ic][0]), &efilenum[0], &(vecEQPE[ic][0])));
    gqpe[ic]->SetName(Form("GraphQPEChan%i", ic));
    gqpe[ic]->SetTitle(Form("Graph-QPE-chan-%i", ic));
    gqpe[ic]->SetMarkerSize(1);
    gqpe[ic]->SetMarkerColor(myColor[ic]);
    gqpe[ic]->SetMarkerStyle(myStyle[ic]);
    fout->Add(gqpe[ic]);
    // if (ic == 6 || ic == 7 || ic == 8)
    if (ic != 5 && ic != 6 && ic != 3 && ic < 9)
      mgQPE->Add(gqpe[ic]);
  }
  // overlay all channel graphs on canvas
  // mg->GetXaxis()->SetNdivisions(1010);
  //    n = n1 + 100 * n2 + 10000 * n3 Where n1 is the number of primary divisions, n2 is the number of second order divisions and n3 is the number of third order divisions. n < 0, the axis will be forced to use exactly n divisions.
  ylabel.Form("single photon charge (normed)");
  setTimeGraph(mgQPE, ylabel);
  TCanvas *canqpe = new TCanvas(Form("Graph-QPE-%s", sdate.c_str()), Form("Graph-QPE-%s", sdate.c_str()));
  mgQPE->Draw("ap");
  canqpe->BuildLegend();
  canqpe->SetGrid();
  canqpe->Print(".png");
  fout->Append(canqpe);

  // QPE sigma graphs one graph per channel
  vector<TGraphErrors *> gqpeSigma;
  TMultiGraph *mgQPESigma = new TMultiGraph();
  for (unsigned ic = 0; ic < vecQPESigma.size(); ++ic)
  {
    // for (int ifile = 0; ifile < vecQPE[ic].size(); ++ifile)
    //   printf(" chan %u file %i %f\n" , ic,ifile,vecQPE[ic][ifile]);
    //  cout << " add " << ic << endl;
    gqpeSigma.push_back(new TGraphErrors(vecQPESigma[ic].size(), &fileTime[0], &(vecQPESigma[ic][0]), &efilenum[0], &(vecEQPESigma[ic][0])));
    gqpeSigma[ic]->SetName(Form("GraphQPESigmaChan%i", ic));
    gqpeSigma[ic]->SetTitle(Form("Graph-QPE-Sigma-chan-%i", ic));
    gqpeSigma[ic]->SetMarkerSize(1);
    gqpeSigma[ic]->SetMarkerColor(myColor[ic]);
    gqpeSigma[ic]->SetMarkerStyle(myStyle[ic]);
    fout->Add(gqpeSigma[ic]);
    if (ic == 6 || ic == 7 || ic == 8)
      mgQPESigma->Add(gqpeSigma[ic]);
  }
  // overlay all channel graphs on canvas
  // mg->GetXaxis()->SetNdivisions(1010);

  //    n = n1 + 100 * n2 + 10000 * n3 Where n1 is the number of primary divisions, n2 is the number of second order divisions and n3 is the number of third order divisions. n < 0, the axis will be forced to use exactly n divisions.
  ylabel.Form("single photon sigma ");
  setTimeGraph(mgQPESigma, ylabel);
  TCanvas *canqpesigma = new TCanvas(Form("Graph-QPESigma-%s", sdate.c_str()), Form("Graph-QPESigma-%s", sdate.c_str()));
  mgQPESigma->Draw("ap");
  canqpesigma->BuildLegend();
  canqpesigma->SetGrid();
  canqpesigma->Print(".png");
  fout->Append(canqpesigma);

  // graphs without norm  one graph per channel

  vector<TGraphErrors *> gqsumUn;
  for (unsigned ic = 0; ic < CHANNELS; ++ic)
  {
    // cout << " add " << ic << endl;
    gqsumUn.push_back(new TGraphErrors(filenum.size(), &fileTime[0], &(vecQsumUn[ic][0]), &efilenum[0], &(vecEQsumUn[ic][0])));
    gqsumUn[ic]->SetName(Form("qsumChanUn%i", ic));
    gqsumUn[ic]->SetTitle(Form("qsum-unnormalized-chan-%i", ic));
    gqsumUn[ic]->SetMarkerSize(1);
    gqsumUn[ic]->SetMarkerColor(myColor[ic]);
    gqsumUn[ic]->SetMarkerStyle(myStyle[ic]);
    fout->Add(gqsumUn[ic]);
  }
  ylabel.Form("integrated charge");
  TMultiGraph *mgL1 = new TMultiGraph();
  // mgL1->Add(gqsumUn[6]);
  mgL1->Add(gqsumUn[7]);
  mgL1->Add(gqsumUn[8]);
  setTimeGraph(mgL1, ylabel);
  TCanvas *canL1 = new TCanvas(Form("QsummaryL1-%s", sdate.c_str()), Form("QsummaryL1-%s", sdate.c_str()));
  mgL1->Draw("ap");
  gPad->Update();
  canL1->BuildLegend();
  canL1->SetGrid();
  fout->Append(canL1);

  TMultiGraph *mgL2 = new TMultiGraph();
  mgL2->Add(gqsumUn[4]);
  // mgL2->Add(gqsumUn[4]);
  mgL2->Add(gqsumUn[5]);
  setTimeGraph(mgL2, ylabel);
  TCanvas *canL2 = new TCanvas(Form("QsummaryL2-%s", sdate.c_str()), Form("QsummaryL2-%s", sdate.c_str()));
  mgL2->Draw("ap");
  gPad->Update();
  canL2->BuildLegend();
  canL2->SetGrid();
  fout->Append(canL2);

  TMultiGraph *mgL3 = new TMultiGraph();
  mgL3->Add(gqsumUn[0]);
  mgL3->Add(gqsumUn[1]);
  mgL3->Add(gqsumUn[2]);
  setTimeGraph(mgL3, ylabel);
  TCanvas *canL3 = new TCanvas(Form("QsummaryL3-%s", sdate.c_str()), Form("QsummaryL3-%s", sdate.c_str()));
  mgL3->Draw("ap");
  gPad->Update();
  canL3->BuildLegend();
  canL3->SetGrid();
  fout->Append(canL3);
}