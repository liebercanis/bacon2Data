// revised April 11 2024 for new gains and new anaCRun.cc
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

//
enum
{
  CHANNELS = 14,
  NONSUMCHANNELS = CHANNELS - 1
};

static int waveBins =7500.;
int nFiles;
TString tag;
TString theStartTag;
TString theEndTag;
std::string sdate;
vector<double> normQsum;
vector<double> normQPE;
TDatime dopeTime;
time_t time0;
time_t time1;
struct tm tmStruct;
Long64_t maxFiles;
double ntotal;
double npass;
vector<double> filePass;
int totalPass;
vector<double> effOther;
vector<double> sumHits;
static double startTime = 660.; // hWave->GetBinLowEdge(maxBin) + hWave->GetBinWidth(maxBin) / 2.;
static double endTime=75000.0;
static double nominalGain = 227.4;

std::vector<TString> fileList;
std::vector<double> filenum;
std::vector<double> efilenum;
std::vector<double> fileTime;
std::vector<TDatime> fileDatime;
TBFile *bf;
TBEventData *eventData;
TTree *runTree;
TH1D *hEventPass;
TH1D *hThreshHist;
TH1D *hCrossHist;
TH1D *hCosmicCut1;
TH1D *hCosmicCut2;
TH1D *hRunEventPass;
TH1D *hRunThreshHist;
TH1D *hRunCrossHist;
TH1D *hRunCosmicCut1;
TH1D *hRunCosmicCut2;
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
std::vector<vector<TH1D *>> vRunHitWave;
std::vector<vector<TH1D *>> vRunPeakWave;
std::vector<TH1D *> hSumWave;
std::vector<TH1D *> hRunHitWave;
std::vector<TH1D *> hRunSumWave;
std::vector<TH1D *> hUnNormedHitWave;
std::vector<TH1D *> hUnNormedSumWave;
std::vector<TH1D *> hRunLatePeakSum;
TDirectory *sumDir;
TDirectory *fitSumDir;
TDirectory *waveSumDir;
TDirectory *qpeSumDir;
TDirectory *peakSumDir;
TDirectory *runSumDir;
TDirectory *gainSumDir;
TFile *fin;
TFile *fout;
bool first = true;
TH1D *hQPEChan;
TH1D *hQPESigmaChan;
TH1D *eventCount;
void makeGraphs();
TNtuple *ntQhit8;
enum dataType
{
  SIS = 0,
  CAEN = 1
};
int theDataType;
TString dirName;
TString dirNameSlash;

double effQuantum128 = 0.17;
double nPhotons = 50.E3 * 5.486; //norm photon yield 2.743E+05 

TString startDate;
TString endDate;

// vectors for gains
std::vector<double> sipmGain;
std::vector<double> sipmGainError;
double xWaveLow = 0;
double xWaveHigh = 7500; // max sample

/* start of code */

void setTime(TString startTag, TString endTag)
{
    theStartTag = startTag;
    theEndTag=endTag;
    tag = theStartTag+TString("-")+theEndTag;
    int month0=TString(theStartTag(0,2)).Atoi(); 
    int day0  =TString(theStartTag(3,2)).Atoi(); 
    int year0 =TString(theStartTag(6,4)).Atoi(); 
    int month1=TString(theEndTag(0,2)).Atoi(); 
    int day1  =TString(theEndTag(3,2)).Atoi(); 
    int year1 =TString(theEndTag(6,4)).Atoi(); 
    //printf(" start %i %i %i ene %i %i %i  \n",month0,day0,year0,month1,day1,year1);

    /* fill in values for 2019-08-22 23:22:26 */
    tmStruct.tm_year = year0;
    tmStruct.tm_mon = month0;
    tmStruct.tm_mday = day0;
    time0  = mktime(&tmStruct);
    printf("set start %s\n", asctime(gmtime(&time0)));
    tmStruct.tm_year = year1;
    tmStruct.tm_mon = month1;
    tmStruct.tm_mday = day1;
    time1  = mktime(&tmStruct);
    printf("set end %s\n", asctime(gmtime(&time1))); 
}

// normalize to total pass
void normalizeTotalPass(TString histSet)
{
  TString histName;
  for (int ichan = 0; ichan < NONSUMCHANNELS; ++ichan)
  {
    if (ichan == 12) // skip broken PMT
      continue;

    histName.Form("UnNormed%sChan%i", histSet.Data(), ichan);
    TH1D *hist;
    runSumDir->GetObject(histName, hist);
    if (!hist)
      printf("at line165 %s not found \n", histName.Data());
    TH1D* hSave;
    if (hist)
    {
      if(histSet.Contains("Hit")) hSave =  hRunHitWave[ichan];
      else hSave =  hRunSumWave[ichan];

      sumHits[ichan] = hist->Integral(startTime,endTime)/nominalGain;
      //printf("at line174 %s normalize to %d\n", hSave->GetName(), totalPass);
      for (int ibin = 0; ibin < hist->GetNbinsX(); ++ibin)
      {
        double xbin = hist->GetBinContent(ibin);
        // double ebin = abs(hist->GetBinError(ibin));
        double ebin = sqrt(abs(hist->GetBinContent(ibin)));
        // if (xbin < 0)
        //cout << " at line212 " << hist->GetName() << " ibin " << ibin << " xbin " << xbin << " ebin " << ebin << " tot " << totalPass << endl;
        // divide by nominal gain july 22 2024
        hSave->SetBinContent(ibin, xbin / double(totalPass));
        hSave->SetBinError(ibin, ebin / double(totalPass));
      }
      cout << "at line171"
           << " totalPass " << totalPass << "  " << hist->GetName() << " integral " << hist->Integral(startTime,endTime)
           << " normed  " << hSave->GetName() << " integral " << hSave->Integral(startTime,endTime) << endl;
    }
  }
}

bool readGains(TString fileName)
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

  printf("\t\t\t stored gains %lu \n", sipmGain.size());
  for (unsigned long j = 0; j < sipmGain.size(); ++j)
  {
    printf(" %lu  gain %.4f error %.4f   \n", j, sipmGain[j], sipmGainError[j]);
  }
  return true;
}

double effGeo(int ichan)
{
  int ilevel = -1;
  if (ichan == 6 || ichan == 7 || ichan == 8)
    ilevel = 0;
  else if (ichan == 3 || ichan == 4 || ichan == 5)
    ilevel = 1;
  else if (ichan == 0 || ichan == 1 || ichan == 2)
    ilevel = 2;

  double e = 1.0;
  if (ilevel < 0)
    return e;
  /*
  Area of SiPMs is 6.0mm x 6.0mm

      Channels 6, 7, and 8 are at 11.6 cm
      from the source Channels 3, 4, and 5 are at 23.2 cm
      from the source Channels 0, 1, and 2 are at 34.8 cm from the source Channel 12 is at 36 cm from the source.
      */
  double a = pow(0.6, 2.);
  double distance2[3];
  double b = 4.0 * TMath::Pi();
  distance2[0] = pow(11.6, 2.);
  distance2[1] = pow(23.2, 2.);
  distance2[2] = pow(34.8, 2.);

  e = a / b / distance2[ilevel];
  return e;
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
  if (theDataType == SIS)
  {
    bf = NULL;
    fin->GetObject("tbfile", bf);
    if (!bf)
    {
      printf(" no TBfile file %s \n", fileList[ifile].Data());
      return datime;
    }
    TString strfileTime = bf->modified;
    TString tempstring = TString(strfileTime(20, 4));
    int fileYear = tempstring.Atoi();
    tempstring = TString(strfileTime(4, 3));
    int fileMonth = get_month_index(tempstring);
    tempstring = TString(strfileTime(8, 2));
    int fileDay = tempstring.Atoi();
    tempstring = TString(strfileTime(11, 2));
    int fileHour = tempstring.Atoi();
    tempstring = TString(strfileTime(14, 2));
    int fileMin = tempstring.Atoi();
    tempstring = TString(strfileTime(17, 2));
    int fileSec = tempstring.Atoi();
    datime.Set(fileYear, fileMonth, fileDay, fileHour, fileMin, fileSec);
    printf("FileYear = %u , FileMonth = %u , FileDay = %u , FileHour = %u , FileMin = %u , FileSec = %u \n", fileYear, fileMonth, fileDay, fileHour, fileMin, fileSec);
  }
  else
  {
    if (!eventData)
    {
      printf(" no eventData in file %s \n", fileList[ifile].Data());
      return datime;
    }

    runTree->GetEntry(ievent);

    printf(" file %i event %lli FileYear = %u , FileMonth = %u , FileDay = %u , FileHour = %u , FileMin = %u , FileSec = %u \n", ifile, ievent, eventData->year + 1900, eventData->mon + 1, eventData->day, eventData->hour, eventData->min, eventData->sec);
    datime.Set(eventData->year, eventData->mon + 1, eventData->day, eventData->hour, eventData->min, eventData->sec);
  }

  return datime;
}

void fileLoop()
{
  totalPass = 0;
  nFiles = 0;
  printf("\n \n \n +++++++ fileLoop over %lld files +++++ \n", maxFiles);
  // DEF would be nice to put in a way to look at the last maxFiles files
  for (unsigned ifile = 0; ifile < maxFiles; ++ifile)
  {
    TString fullName = dirNameSlash + fileList[ifile];
    fin = new TFile(fullName, "readonly");

    printf(" \n ***** starting file %i  %s  *******\n", ifile, fin->GetName());

    eventCount = NULL;

    cout << " get event Count " << ifile << endl;
    fin->GetObject("eventcount", eventCount);
    if (eventCount)
    {
      ntotal = eventCount->GetBinContent(0);
      npass = eventCount->GetBinContent(1);
      filePass.push_back(npass);
      // 0 = ntriggers, 1 = npass
      //printf("\n\n line485 total events file %i %s %f events passed  %f \n ", ifile, fileList[ifile].Data(), eventCount->GetBinContent(0), eventCount->GetBinContent(1));
      // TH1D *hevcount = (TH1D *)eventCount->Clone(Form("eventCount%i", ifile));
      // fout->Add(hevcount);
    }
    else
    {
      printf("NO EVENT COUNT IN FILE %i  %s\n", ifile, fin->GetName());
      continue;
    }
    /* add gains */
    TDirectory *anaDir = NULL;
    fin->GetObject("anadir", anaDir); // typeO in directory name
    if (!anaDir)
    {
      //printf("line489 anaDir not found in file %s \n", fin->GetName());
      continue;
    }
    else
    {
      //printf("line502  anaDir is found in file %s \n", fin->GetName());
      // loop over channels
      TString gainInName;
      TString gainOutName;
      for (int ichan = 0; ichan < NONSUMCHANNELS; ++ichan)
      {
        gainOutName.Form("GainSumChan%i", ichan);
        gainInName.Form("LatePeakSumChan%i", ichan);
        TH1D *hIn = NULL;
        anaDir->GetObject(gainInName, hIn);
        if (hIn)
        {
          TH1D *hOut = NULL;
          gainSumDir->GetObject(gainOutName, hOut);
          if (hOut == NULL)
          {
            hOut = (TH1D *)hIn->Clone(gainOutName);
            hOut->SetTitle(gainOutName);
            cout << "line516 ... addding  " << hIn->GetName() << " file "
                 << fin->GetName() << " summed entries " << hOut->GetEntries() << endl;
            gainSumDir->Add(hOut);
          }
          else
          {
            gainSumDir->GetObject(gainOutName, hOut);
            hOut->Add(hIn);
            cout << "line526 ... found for " << gainOutName
                 << " in entries " << hIn->GetEntries() << " summed entries " << hOut->GetEntries() << endl;
          }
        } // if hIn
      } // channel loop
    } // if anaDir

    // get hist of cleanup cut pass bit
    fin->GetObject("EventPass", hEventPass);
    if (!hEventPass)
    {
      printf("NO EVENT PASS IN FILE %s\n", fin->GetName());
      continue;
    }
    if (ifile == 0)
    {
      hRunEventPass = (TH1D *)hEventPass->Clone("RunEventPass");
      fout->Add(hRunEventPass);
    }
    else
    {
      fout->GetObject("RunEventPass", hRunEventPass);
      hRunEventPass->Add(hEventPass);
    }
    // thresh
    anaDir->GetObject("threshHist", hThreshHist);
    if (!hThreshHist)
      cout << "  thresh" << endl;
    if (ifile == 0)
    {
      hRunThreshHist = (TH1D *)hThreshHist->Clone("RunThreshHist");
      fout->Add(hRunThreshHist);
    }
    else
    {
      fout->GetObject("RunThreshHist", hRunThreshHist);
      hRunThreshHist->Add(hThreshHist);
    }

    // crossings
    anaDir->GetObject("crossHist", hCrossHist);
    if (!hCrossHist)
      cout << " cross" << endl;
    if (ifile == 0)
    {
      hRunCrossHist = (TH1D *)hCrossHist->Clone("RunCrossHist");
      fout->Add(hRunCrossHist);
    }
    else
    {
      fout->GetObject("RunCrossHist", hRunCrossHist);
      hRunCrossHist->Add(hCrossHist);
    }
    // cosmic 1
    anaDir->GetObject("cosmicCut1", hCosmicCut1);
    if (!hCosmicCut1)
      cout << " cosmic 1" << endl;
    if (ifile == 0)
    {
      hRunCosmicCut1 = (TH1D *)hCosmicCut1->Clone("RunCosmicCut1");
      fout->Add(hRunCosmicCut1);
    }
    else
    {
      fout->GetObject("RunCosmicCut1", hRunCosmicCut1);
      hRunCosmicCut1->Add(hCosmicCut1);
    }
    // cosmic 2
    anaDir->GetObject("cosmicCut2", hCosmicCut2);
    if (!hCosmicCut2)
      cout << " cosmic 2" << endl;
    if (ifile == 0)
    {
      hRunCosmicCut2 = (TH1D *)hCosmicCut2->Clone("RunCosmicCut2");
      fout->Add(hRunCosmicCut2);
    }
    else
    {
      fout->GetObject("RunCosmicCut2", hRunCosmicCut2);
      hRunCosmicCut2->Add(hCosmicCut2);
    }

    //

    sumDir = NULL;
    fin->GetObject("sumDir", sumDir);
    if (!sumDir)
    {
      printf(" sumDir not found in file %s \n", fin->GetName());
      continue;
    }
    printf(" sumDir for file %s \n", fin->GetName());

    TH1D *hqsum = NULL;
    TH1D *hqprompt = NULL;
    if (theDataType == SIS)
    {
      fin->GetObject("histQsum", hqsum);
      fin->GetObject("histQprompt", hqprompt);
    }
    else
    {
      fin->GetObject("histqsum", hqsum);
      fin->GetObject("histqprompt", hqprompt);
    }
    runTree = NULL;
    eventData = NULL;
    if (theDataType == CAEN)
    {
      eventData = new TBEventData();
      fin->GetObject("RunTree", runTree);
      if (!runTree)
        printf("file has not RunTree %s \n", fin->GetName());
      else
        runTree->SetBranchAddress("eventData", &eventData);
    }

    if (hqsum && hqprompt)
    {
      TH1D *hq = (TH1D *)hqsum->Clone(Form("hqsum%i", ifile));
      peakSumDir->Add(hq);
      TH1D *hp = (TH1D *)hqsum->Clone(Form("hqprompt%i", ifile));
      peakSumDir->Add(hp);
    }

    filenum.push_back(double(ifile));
    efilenum.push_back(0);
    TDatime datetime = getTime(ifile);
    fileDatime.push_back(datetime);
    fileTime.push_back(datetime.Convert());

    cout << " TIME FILE  " << ifile << " " << fileList[ifile] << " modified "
         << " TDatime " << datetime.AsString() << " as int " << datetime.Convert() << endl;

    if (hqsum && hqprompt)
    {
      for (int i = 0; i < hqsum->GetNbinsX() - 1; ++i)
      {
        double eval = 0;
        double val = 0;
        if (!isnan(hqsum->GetBinContent(i + 1)) && !isinf(hqsum->GetBinContent(i + 1)))
        {
          val = hqsum->GetBinContent(i + 1);
          eval = hqsum->GetBinError(i + 1);
        }
        vecQsum[i].push_back(val);
        vecEQsum[i].push_back(eval);
        vecQsumUn[i].push_back(val);
        vecEQsumUn[i].push_back(eval);
        cout << "chan " << i << " qsum " << hqsum->GetBinContent(i + 1) << " size " << vecQsum[i].size() << endl;
      }
    }

    /******  loop over sumDir *****/
    TList *sumList = sumDir->GetListOfKeys();
    TIter next(sumList);
    TKey *key;
    printf("addsumDirHistos %u \n", sumList->GetEntries());
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

      /****** get run sum ****/
      TH1D *hClone;
      TH1D *hRunClone;
      if (name.find("sumWave") != std::string::npos && name.find("Bad") == std::string::npos)
      {
        // waveSumDir->cd();
        cout << "line700  sumWave clone " << name << " file " << ifile << endl;
        string chan = name.substr(name.find_last_of("e") + 1);
        int ichan = stoi(chan);
        TString cloneName;
        cloneName.Form("RunSumWaveFile%uChan%i", ifile, ichan);
        hClone = (TH1D *)h->Clone(cloneName);
        hClone->SetTitle(cloneName);
        hSumWave.push_back(hClone);
        waveSumDir->Add(hClone);
      }
      if (name.find("sumHitWave") != std::string::npos)
      {
        // waveSumDir->cd();
        cout << "line709  sumHitWave clone " << name << " file " << ifile << endl;
        string chan = name.substr(name.find_last_of("e") + 1);
        int ichan = stoi(chan);
        TString cloneName;
        cloneName.Form("RunHitWaveFile%uChan%i", ifile, ichan);
        hClone = (TH1D *)h->Clone(cloneName);
        hClone->SetTitle(cloneName);
        //TString histName;
        // histName.Form("RunHitWaveChan%i", ichan);
        // hRunClone = (TH1D *)h->Clone(histName);
        waveSumDir->Add(hClone);

        cout << "line718  sumHitWave clone " << hClone->GetName()
             << " file " << ifile << " vRunHitWave " << vRunHitWave[ichan].size() << endl;
        /*  do not sum here */
      }
      if (name.find("sumPeakWave") != std::string::npos)
      {
        // waveSumDir->cd();
        cout << "line722 sumPeakWave clone " << name << " file " << ifile << endl;
        string chan = name.substr(name.find_last_of("e") + 1);
        int ichan = stoi(chan);
        TString cloneName;
        cloneName.Form("RunPeakWaveFile%uChan%i", ifile, ichan);
        hClone = (TH1D *)h->Clone(cloneName);
        hClone->SetTitle(cloneName);
        //TString histName;
        // histName.Form("RunPeakWaveChan%i", ichan);
        // hRunClone = (TH1D *)h->Clone(histName);
        waveSumDir->Add(hClone);
        vRunPeakWave[ichan].push_back(hClone);
        /* do not sum here */
      }

      /****** get PeakChan by channel ****/
      if (name.find("QPeakChan") != std::string::npos)
      {
        string chan = name.substr(name.find_last_of("n") + 1);
        int ichan = stoi(chan);
        cout << "line739 " << name << "string chan" << chan << " int " << ichan << endl;

        TString cloneName;
        cloneName.Form("runPeakSumCh%ifile%i", ichan, ifile);
        TH1D *hpAdd = (TH1D *)h->Clone(cloneName);
        cout << " line744 found RunQSum " << ichan << " " << fileList[ifile] << " " << h->GetName() << " " << hpAdd->GetName() << endl;
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
        cout << " line758 found RunQSum " << ichan << " " << fileList[ifile] << " " << h->GetName() << " " << hqAdd->GetName()
             << " max x  " << h->GetBinLowEdge(h->GetNbinsX() + 1) << endl;
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
        TH1D* hpAdd=NULL;
        runSumDir->GetObject(cloneName,hpAdd);
        if(!hpAdd) {
          TH1D *hpAdd = (TH1D *)h->Clone(cloneName);
          hpAdd->SetMarkerStyle(20);
          hpAdd->SetMarkerSize(0.5);
          qpeSumDir->Add(hpAdd);
        } else { // already exits so sum this histo
          hpAdd->Add(h);
        }
      }
    } // sum over keys
    cout << " end loop over sumDir keys " << endl; // loop over sumDir keys
    // waveSumDir->ls();
    fout->Write();
    // if(ifile==0) fout->ls();
    // waveSumDir->Write();
    fin->Close();
    totalPass += int(npass);
    ++nFiles;
    cout << "fileLoop finished file " << ifile << "  " << nFiles << endl;
  } // end loop over files
  printf("fileLoop finished over %lld good files %d  totalPass %i \n", maxFiles, nFiles, totalPass); // DEF would be nice to put in a way to look at the last maxFiles files
}

void sumHistosChannel(int ichan, TString histSet)
{
  if (ichan == 12) // skip PMT
    return;
  vecQPEMean[ichan] = QPEMean[ichan]; // starting value
  /* if (ichan > 8 && ichan < 12)
    return; // skip trigger sipms
  */
  //printf("line792 sumHistosChannel:: sum histos total pass %i ichan %i set %s QPEMean %f  \n", totalPass, ichan, histSet.Data(), QPEMean[ichan]);

  // sum over files
  for (int ih = 0; ih < nFiles; ++ih)
  {
    // waveToSum already in output file by run and channel
    TString histName;
    histName.Form("Run%sFile%uChan%i", histSet.Data(), ih, ichan);
    TH1D *waveToSum = NULL;
    waveSumDir->GetObject(histName, waveToSum);
    if (waveToSum == NULL)
    {
      //printf("line799 skipping %s  %s chan %i file %i \n", histSet.Data(), histName.Data(), ichan, ih);
      continue;
    }

    // cout << "at line 765 waveToSum " << waveToSum->GetName() << endl;

    // new histogram
    int nbinsx = waveToSum->GetNbinsX();
    double xlow = waveToSum->GetXaxis()->GetBinLowEdge(0);
    double xup = waveToSum->GetXaxis()->GetBinUpEdge(nbinsx);

    // normalize to qpe.
    double qpeChan = 227.4; // nominal gain vecQPE[ichan][ih]; come back to thhis
                            /*
                            if (qpeChan <= 0 || isnan(qpeChan)) // bad value
                              qpeChan = vecQPEMean[ichan];
                            else if (qpeChan < 0.75 * vecQPEMean[ichan] || qpeChan > 1.50 * vecQPEMean[ichan]) // check if too far from mean;
                              qpeChan = vecQPEMean[ichan];
                            else // update mean
                            {
                              vecQPEMean[ichan] = (vecQPEMean[ichan] * double(vecQPENave[ichan]) + qpeChan) / double(vecQPENave[ichan] + 1);
                              vecQPENave[ichan] += 1;
                            }
                            */

    fitSumDir->cd();
    TH1D *hWaveToFit;
    TH1D *hWaveToFitNotNormed;
    TString fitwaveName;
    fitwaveName.Form("fitwave%sChan%iFile%i", histSet.Data(), ichan, ih);
    TString notNormedName;
    notNormedName.Form("notNormedWave%sChan%iFile%i", histSet.Data(), ichan, ih);
    hWaveToFit = (TH1D *)waveToSum->Clone(fitwaveName);
    hWaveToFitNotNormed = (TH1D *)waveToSum->Clone(notNormedName);
    // cout << hWaveToFit->GetName() << endl;
    hWaveToFit->SetMarkerStyle(21);
    hWaveToFit->SetMarkerSize(0.2);
    hWaveToFit->GetListOfFunctions()->Clear();
    hWaveToFitNotNormed->SetMarkerStyle(21);
    hWaveToFitNotNormed->SetMarkerSize(0.2);
    hWaveToFitNotNormed->GetListOfFunctions()->Clear();

    //printf(" line 850 NNNNNN normalize to qpe chan %i nbins %i file %i qpe %f !!!!!\n", ichan, waveToSum->GetNbinsX(), ih, qpeChan);
    /* loop over histogram bins */
    //printf("line862 waveToSum %s  chan %i file %i entries %f total event %f  \n",
           //waveToSum->GetName(), ichan, ih, waveToSum->GetEntries(), double(filePass[ih]));
    for (int ibin = 0; ibin < waveToSum->GetNbinsX(); ++ibin)
    {
      // apply effOther correction
      double xbin = waveToSum->GetBinContent(ibin) / effOther[ichan];
      if (qpeChan <= 0 || qpeChan > 1.0E9)
      {
        printf("BBBBBBBB bad  QPE %f chan %i file %i \n", qpeChan, ichan, ih);
        qpeChan = 1.0;
      }
      // apply QPE norm to hit wave
      if (histSet == TString("HitWave"))
        xbin /= qpeChan;
      if (isnan(xbin))
      {
        printf("BBBBBBBB bad xbin value xbin  %f chan %i file %i qpe %f \n", xbin, ichan, ih, qpeChan);
        xbin = 0.0;
      }
      double ebin;
      if (histSet == TString("HitWave"))
        ebin = sqrt(xbin);
      else
        ebin = waveToSum->GetBinError(ibin) / effOther[ichan];

      // norm to number of triggers
      hWaveToFitNotNormed->SetBinContent(ibin, xbin);
      hWaveToFitNotNormed->SetBinError(ibin, ebin);
      // norm to total pass in file
      hWaveToFit->SetBinContent(ibin, xbin / double(filePass[ih]));
      hWaveToFit->SetBinError(ibin, ebin / double(filePass[ih]));
      if (ichan == 8)
        ntQhit8->Fill(float(ih), float(ibin), qpeChan, hWaveToFit->GetBinContent(ibin));
    }
    /*
    cout << "at line894 "
         << " filePass " << filePass[ih] << "  " << hWaveToFit->GetName()
         << " sum entries " << waveToSum->GetEntries()
         << " sum integral " << waveToSum->Integral()
         << " fit integral  " << hWaveToFit->Integral() << endl;
         */

    // check histo
    bool addIt = true;
    /*
    for (int ibin = 1; ibin < hWaveToFit->GetNbinsX(); ++ibin)
      if (isnan(hWaveToFit->GetBinContent(ibin)) || isinf(hWaveToFit->GetBinContent(ibin)))
      {
        printf("XXXXXXXXXXXXX WARNING bin %i chan %i hist %s file %i %s XXXXXXX ", ibin, ichan, histName.Data(), ih, fileList[ih].Data());
        addIt = false;
      }
    if (!addIt)
      cout << endl;
      */

    /*
    if (addIt)
      printf("addIt line917 totalPass %i set %s hist %s hist %s file %i %s \n",
             totalPass, histSet.Data(), hWaveToFit->GetName(), waveToSum->GetName(), ih, fileList[ih].Data());
    else
      printf("NOT addIt %s %s \n", hWaveToFit->GetName(), fileList[ih].Data());
      */

    // add and save in output file runSumDir;
    if (histSet == TString("HitWave"))
    {
      runSumDir->cd();
      if (hUnNormedHitWave[ichan] == NULL && addIt)
      {
        // histName.Form("Run%sChan%i", histSet.Data(), ichan);
        //  do norm to tatal pass after sum over files
        // hRunHitWave[ichan] = (TH1D *)hWaveToFitNotNormed->Clone(histName);
        histName.Form("UnNormed%sChan%i", histSet.Data(), ichan);
        //cout << " line915 NNNNNNNNNN new clone " << histName << " of " << hWaveToFitNotNormed->GetName() << endl;
        hUnNormedHitWave[ichan] = (TH1D *)hWaveToFitNotNormed->Clone(histName);
        hUnNormedHitWave[ichan]->SetTitle(histName);
        // runSumDir->Add(hRunHitWave[ichan]);
        runSumDir->Add(hUnNormedHitWave[ichan]);
      }
      else if (addIt)
      {
        histName.Form("UnNormed%sChan%i", histSet.Data(), ichan);
        // hRunHitWave[ichan]->Add(hWaveToFitNotNormed);
        runSumDir->GetObject(histName, hUnNormedHitWave[ichan]);
        if (hUnNormedHitWave[ichan] == NULL)
        {
          //printf("line 951 NULL chan %i file %i %s \n", ichan, ih, histName.Data());
          runSumDir->ls();
        }
        hUnNormedHitWave[ichan]->Add(hWaveToFitNotNormed);
      }
      //cout << " at line943 " << fileList[ih].Data() << " totalPass " << totalPass << "  " << hUnNormedHitWave[ichan]->GetName() << " integral " << hUnNormedHitWave[ichan]->Integral() << endl;

      // end add RunHitWave
    }
    else if (histSet == TString("SumWave"))
    {
      runSumDir->cd();
      histName.Form("UnNormed%sChan%i", histSet.Data(), ichan);
      // cout << " line896 sum to " << histName << " summing " << waveToSum->GetName() << endl;
      if (hUnNormedSumWave[ichan] == NULL && addIt)
      {
        //cout << " line939 NNNNNNNNNN new clone " << histName << " " << waveToSum->GetName() << endl;
        hUnNormedSumWave[ichan] = (TH1D *)waveToSum->Clone(histName);
        hUnNormedSumWave[ichan]->SetTitle(histName);
      }
      else if (addIt)
      {
        histName.Form("UnNormed%sChan%i", histSet.Data(), ichan);
        runSumDir->GetObject(histName, hUnNormedSumWave[ichan]);
        if (hUnNormedSumWave[ichan] == NULL)
          printf("line950 no RunSumWave %s\n", histName.Data());
        hUnNormedSumWave[ichan]->Add(hWaveToFit);
      }
    }
  } // sum over files
    // end add RunHitWave
}

void sumHistos()
{
  // histograms time in ns
  runSumDir->cd();
  TString histName;
  for (unsigned ichan = 0; ichan < CHANNELS; ++ichan)
  {
    histName.Form("Run%sChan%i","HitWave", ichan);
    hRunHitWave[ichan] = new TH1D(histName,histName,waveBins,0,2.*waveBins);
    hRunHitWave[ichan]->GetXaxis()->SetTitle("time [ns]");
    hRunHitWave[ichan]->GetYaxis()->SetTitle("yield [SPE] ");
    histName.Form("Run%sChan%i","SumWave", ichan);
    hRunSumWave[ichan] = new TH1D(histName,histName,waveBins,0,2.*waveBins);
    hRunSumWave[ichan]->GetXaxis()->SetTitle("time [ns]");
    hRunSumWave[ichan]->GetYaxis()->SetTitle("yield [SPE] ");
  }

  runSumDir->ls();

  printf("line963 sumHistos: Number of files  %d vRunHitWave size %lu\n", nFiles, vRunHitWave.size());
  // loop over channels
  for (int ichan = 0; ichan < CHANNELS; ++ichan)
    sumHistosChannel(ichan, TString("HitWave"));
  for (int ichan = 0; ichan < CHANNELS; ++ichan)
    sumHistosChannel(ichan, TString("SumWave"));
  // normalize each channel
  // waveSumDir->ls();
  normalizeTotalPass(TString("HitWave"));
  normalizeTotalPass(TString("SumWave"));
}

void fitSlopes()
{
  printf("Make slope graph.  Number of files = %lu Number of channels = %lu \n", filenum.size(), vRunPeakWave.size());
  double flow = 0.;
  double fhigh = 0;
  if (theDataType == SIS)
  { // For SIS data
    flow = 8 * 120.;
    fhigh = 8 * 500.;
  }
  else
  { // For CAEN data
    flow = 2. * 1000.;
    fhigh = 2. * 3000.;
  }

  for (unsigned ichan = 0; ichan < vRunHitWave.size(); ++ichan)
  {
    if (ichan > 8 && ichan < 12)
      continue; // skip trigger sipms
    TString histName;
    TH1D *hWaveToFit = NULL;
    for (unsigned ih = 0; ih < vRunPeakWave[ichan].size(); ++ih)
    {
      histName.Form("fitwaveChan%iFlile%i", ichan, ih);
      fitSumDir->GetObject(histName, hWaveToFit);
      if (hWaveToFit == NULL)
        continue;
      printf("slope graph ichan %i %s \n", ichan, hWaveToFit->GetName());
      TF1 *gslopefit = NULL;
      /* chisq fit def, L for likelihood
        L Uses a log likelihood method(default is chi - square method).
        To be used when the histogram represents counts.
      */

      hWaveToFit->Fit("expo", "QLQ", " ", flow, fhigh);
      gslopefit = (TF1 *)hWaveToFit->GetListOfFunctions()->FindObject("expo");
      double mfit = 0;
      double emfit = 0;
      double time = 0;
      double etime = 0;
      double MicroSecPerNs = 1. / 1000.; // DEF added in switch for different digitizer MicroSecPerDac
      if (gslopefit)
      {
        mfit = -1. * gslopefit->GetParameter(1);
        emfit = gslopefit->GetParError(1);
        time = 0;
        etime = 0;
        if (mfit != 0)
        {
          time = MicroSecPerNs / mfit;
          etime = gslopefit->GetParError(1) / abs(gslopefit->GetParameter(1)) * abs(time);
        }
      }
      /* correct but confusing MG
        if (time > 40) time = 0;
        etime = gslopefit->GetParError(1) *(-1)*time / mfit;
        if (etime > time) etime = time;
      */
      printf(" Slope for file %u, channel %i = %f tau = %f tau_error = %f \n", ih, ichan, mfit, time, etime);
      vSlope[ichan].push_back(time);
      vESlope[ichan].push_back(etime);
      // printf(" %u %s %s chan %i file %i time %f %f\n", ih, hSumWave[ih]->GetName(), hSumHitWave[ih]->GetName(), ichan, ifile,time,etime);
    }
  }
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
  cout << " count files in dir " << dirName << endl;
  TSystemDirectory dir(dirName, dirName); // TSystemDirectory
  TList *files = dir.GetListOfFiles();    //
  TIter next(files);
  TSystemFile *file;
  while ((file = (TSystemFile *)next()))
  {
    string name = string(file->GetName());
    TString tname = TString(name.c_str());
    string exten = name.substr(name.find_last_of(".") + 1);
    if (exten != string("root"))
      continue;
      // see if file is correct time
      int month = TString(tname(tname.Last('n')+2,2)).Atoi();
      int day = TString(tname(tname.Last('n')+5,2)).Atoi();
      int year = TString(tname(tname.Last('n')+8,4)).Atoi();
      tmStruct.tm_year = year;
      tmStruct.tm_mon = month;
      tmStruct.tm_mday = day;
      time_t fileTime  = mktime(&tmStruct);
       const auto diff0 = std::difftime( fileTime, time0 );
       const auto diff1 = std::difftime( fileTime, time1 );
      printf("line1139 info : file %s time %s",tname.Data(),asctime(gmtime(&fileTime))); 
      cout << " \t ..... " << diff0 << " "<< diff1 <<endl;
      bool timetest = diff0>=0 && diff1<=0;
      if(!timetest) {
        cout <<"line1148 skip out of time file " << name << endl;
        continue;
      }

      // see if file is good
      TString fullName = dirNameSlash + TString(name.c_str());
      TFile *f = new TFile(fullName, "READONLY");
      TTree *RunTree = NULL;
      f->GetObject("RunTree", RunTree);
      if (RunTree == NULL)
      {
        cout << "line1093 skipping BAD file " << name << endl;
        continue;
      }
      f->Close();
      // good file add to list
      fileList.push_back(TString(name.c_str()));
      cout << "line1099 add file " << name << " nFiles= " << fileList.size() << endl;
  }
  return fileList.size();
}

int main(int argc, char *argv[])
{
  cout << "executing " << argv[0] << " make summary plots  " << endl;
  printf(" usage: summary start date string <stag> end date string <etag> max files <default all> \n ");
  if (argc < 1)
  {
    printf("reguire file date start string <stag> args\n");
    exit(0);
  }
  dirName = TString("caenData");
  dirNameSlash = TString("caenData/");
  theStartTag=TString(argv[1]); 

  printf(" input args %i \n ", argc);
  for (int jarg = 1; jarg < argc; ++jarg)
    printf(" %i= %s ", jarg, argv[jarg]);
  printf("\n");

  theEndTag=TString(argv[1]);
  if (argc ==3)
  {
    theEndTag=TString(argv[2]); 
  }

  setTime(theStartTag,theEndTag);

  dopeTime = TDatime(2023, 3, 9, 22, 0, 0);

  printf("count files from %s to %s \n",theStartTag.Data(),theEndTag.Data());
  unsigned nfiles = countFiles();
  if (nfiles == 0)
  {
    printf(" >>>> datatype no files found <<<<\n");
    exit(0);
  }

  printf(" >>>>> files from %s to %s tag %s <<<<<\n",theStartTag.Data(),theEndTag.Data(),tag.Data());

  for (int i = 0; i < fileList.size(); ++i)
    cout << i << "  " << fileList[i] << endl;

  maxFiles = fileList.size();
  if (argc > 3)
  {
    maxFiles = atoi(argv[3]);
  }
  printf(" for %s found %lu files maxFiles %lli \n", tag.Data(), fileList.size(),maxFiles);
  // cleanup cuts;
  hEventPass = NULL;
  hThreshHist = NULL;
  hCrossHist = NULL;
  hCosmicCut1 = NULL;
  hCosmicCut2 = NULL;
  hRunEventPass = NULL;
  hRunThreshHist = NULL;
  hRunCrossHist = NULL;
  hRunCosmicCut1 = NULL;
  hRunCosmicCut2 = NULL;

  sdate = currentDate();
  cout << " starting summary for   " << maxFiles << endl;

  fout = new TFile(Form("summary-%s-nfiles-%lld-created-%s.root", tag.Data(), maxFiles, sdate.c_str()), "recreate");
  fitSumDir = fout->mkdir("fitSumDir");
  waveSumDir = fout->mkdir("waveSumDir");
  qpeSumDir = fout->mkdir("qpeSumDir");
  peakSumDir = fout->mkdir("peakSumDir");
  runSumDir = fout->mkdir("runSumDir");
  gainSumDir = fout->mkdir("gainSumDir");
  fout->cd();

  hQPEChan = new TH1D("QPEChan", "QPE  by channel", 12, 0, 12);
  hQPESigmaChan = new TH1D("QPESigmaChan", "QPE  by channel", 12, 0, 12);
  hQPEChan->Sumw2();
  hQPESigmaChan->Sumw2();
  effOther.resize(CHANNELS);
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
  vRunHitWave.resize(CHANNELS);
  vRunPeakWave.resize(CHANNELS);

  hRunHitWave.resize(CHANNELS);
  hRunSumWave.resize(CHANNELS);
  hUnNormedHitWave.resize(CHANNELS);
  hUnNormedSumWave.resize(CHANNELS);
  sumHits.resize(CHANNELS);
  hRunLatePeakSum.resize(CHANNELS);

  for (unsigned ichan = 0; ichan < CHANNELS; ++ichan)
  {
    hUnNormedHitWave[ichan] = NULL;
    hUnNormedSumWave[ichan] =  NULL;
    effOther[ichan] = 1.;
    sumHits[ichan] = 0.;
  }


  
  fout->cd();


  // nominal values
  QPEMean.resize(CHANNELS);
  QPEMean[0] = 5.576E+03;
  QPEMean[1] = 5.206E+03;
  QPEMean[2] = 5.488E+03;
  QPEMean[3] = 1.000E+00;
  QPEMean[4] = 1.139E+04;
  QPEMean[5] = 5.103E+03;
  QPEMean[6] = 8.700E+03;
  QPEMean[7] = 1.143E+04;
  QPEMean[8] = 1.203E+04;
  QPEMean[9] = 1.000E+00;
  QPEMean[10] = 1.000E+00;
  QPEMean[11] = 1.000E+00;
  QPEMean[12] = 8.187E+02;

  // from integral
  /*
  effOther[0] = 1.720E-01;
  effOther[1] = 1.558E-01;
  effOther[2] = 1.986E-01;
  effOther[4] = 4.757E-01;
  effOther[5] = 5.218E+00;
  effOther[7] = 2.077E-01;
  effOther[8] = 8.284E-02;
  effOther[12] = 2.801E-01;

  // from fits
  effOther[0] = 1.729E-01;
  effOther[1] = 2.604E-01;
  effOther[2] = 2.176E-01;
  effOther[4] = 4.852E-01;
  effOther[5] = 6.848E+00;
  effOther[7] = 2.127E-01;
  effOther[8] = 8.533E-02;
  effOther[12] = 2.779E-01;
  */

  // read old gain file
  TString gainFileName = TString(getenv("BOBJ"))+TString("/gains-2024-02-15-17-26-save.root");
  cout << "read gains from file " << gainFileName << endl;
  if(!readGains(gainFileName)) {
    printf("no gain file %s so exit \n",gainFileName.Data());
    exit(0);
  }


  // use a relative normalization
  /*
  for (unsigned ichan = 0; ichan < CHANNELS; ++ichan)
  {
    effOther[ichan] = effOther[ichan] / effOther[7];
  }
  */
  cout << "@ fileLoop" << endl;
  fileLoop();

  waveSumDir->ls();
  for (int ichan = 0; ichan < CHANNELS; ++ichan)
    printf("%i vRunHitWave size %lu\n", ichan, vRunHitWave[ichan].size());

  printf("line1270 \t\t >>> after fileLoop processed << %li  total pass %i channels %lu <<<<< \n", filenum.size(), totalPass, vRunPeakWave.size());

  for (unsigned jfile = 0; jfile < filenum.size(); ++jfile)
  {
    printf(" %i %s \n", int(filenum[jfile]), fileList[jfile].Data());
  }

  cout << "@ fitQPE" << endl;
  fitQPE();

  // call function to fit slopes and fill vSlope, vESlope
  ntQhit8 = new TNtuple("ntQhit8", " chan 8  ", "run:bin:qpe:qhit");
  if (filenum.size() > 0)
  {
    cout << "@ line1228 sumHistos files " << filenum.size() << " total pass " << totalPass << endl;
    sumHistos();
    // cout << "@  line1230 fitSlopes" << endl;
    //  fitSlopes();
    // cout << "@ line1232 makeGraphs" << endl;
    //  makeGraphs();
  }

  // report
  /*for (unsigned it = 0; it < fileTime.size(); ++it)
    if (fileDatime[it].Convert() < dopeTime.Convert())
      cout << " before " << fileDatime[it].AsString() << "  " << fileList[it] << endl;
    else
      cout << " after  " << fileDatime[it].AsString() << "  " << fileList[it] << endl;
      */
  // fout->ls();

  cout << "Average QPE" << endl;

  for (int ichan = 0; ichan < CHANNELS; ++ichan)
  {
    double sum = 0;
    for (int ifile = 0; ifile < vecQPE[ichan].size(); ++ifile)
      sum += vecQPE[ichan][ifile];
    printf(" chan %i average over %li files, QPE %f \n", ichan, vecQPE[ichan].size(), sum / double(vecQPE[ichan].size()));
  }

  printf(" total triggers %i \n", totalPass);

  // print totalHits
  printf("\t\t >>> end of job: files processed << %li  total pass %i <<<<< \n", filenum.size(), totalPass);

  double fitStart = 12000;
  double fitEnd = 14000;
  for (unsigned ichan = 0; ichan < CHANNELS; ++ichan)
  {
    if (ichan == 9 || ichan == 10 || ichan == 11 || ichan == 3 || ichan == 6)
      continue;
    TString histName;
    histName.Form("RunHitWaveChan%i", ichan);
    runSumDir->GetObject(histName, hRunHitWave[ichan]);
    if (hRunHitWave[ichan] == NULL)
    {
      cout << " cannot find " << histName << endl;
      continue;
    }
    int nbins = hRunHitWave[ichan]->GetNbinsX();
    double sum = hRunHitWave[ichan]->Integral(1, nbins);
    int startBins = hRunHitWave[ichan]->FindBin(fitStart);
    int endBins = hRunHitWave[ichan]->FindBin(fitEnd);
    double ave = hRunHitWave[ichan]->Integral(startBins, endBins) / double(endBins - startBins);

    TF1 *gfit = NULL;
    double ave2 = 0;
    hRunHitWave[ichan]->Fit("pol0", "LF", " ", fitStart, fitEnd);
    gfit = (TF1 *)hRunHitWave[ichan]->GetListOfFunctions()->FindObject("pol0");
    if (gfit)
    {
      // this shouldnot get divided
      ave2 = gfit->GetParameter(0); // / double(backBins);
      printf("......pol0 fit to chan %i value %f \n", ichan, ave2);
    }
    else
      printf("WARNING pol0 Fit to chan %u fails \n", ichan);

    if (isnan(ave2) || ave2 < 0)
    {
      printf("WARNING pol0 Fit to chan %u fails \n", ichan);
      ave2 = ave;
    }
    double back = ave * double(nbins);
    double back2 = ave2 * double(nbins);
    // from fit to hist
    // if (ichan == 12)
    //  back = 6.45;
    printf("RunHitWave chan %i nbins %i entries %f sum %E  bins (%f,%f)  ave (%.4f,%.4f)  back (%E,%E)\n", ichan, nbins, hRunHitWave[ichan]->GetEntries(), sum, fitStart, fitEnd, ave, ave2, back, back2);
    // sumHits[ichan] = sum - back;
    if (sum - back < 0)
      back = back2;
      sum /= sipmGain[ichan];
    //printf("totalHits[%i]=%E;\n", ichan, sum - back);
    printf("totalHits[%i]=%E SPE;\n", ichan, sum);
  }
  printf("\n totalHits in %i passed triggers \n", totalPass);
  for (unsigned ichan = 0; ichan < CHANNELS; ++ichan)
    if (sumHits[ichan] != 0)
      printf("totalHits[%i]=%E ;\n", ichan, sumHits[ichan]);

  printf("\n qpeMean in %i passed triggers \n", totalPass);
  for (unsigned ichan = 0; ichan < CHANNELS; ++ichan)
    printf("QPEMean[%i]=%.3E ;\n", ichan, vecQPEMean[ichan]);

  printf("\n effOther\n");
  for (unsigned ichan = 0; ichan < CHANNELS; ++ichan)
    printf("effOther[%i]=%f ;\n", ichan, effOther[ichan]);

  //    fout->ls();
  // get rid of first bin
  for (int ihist = 0; ihist < hRunLatePeakSum.size(); ++ihist)
  {
    TString runHistName;
    runHistName.Form("RunLatePeakSumChan%i", ihist);
    fout->GetObject(runHistName, hRunLatePeakSum[ihist]);
    if (hRunLatePeakSum[ihist])
    {
      hRunLatePeakSum[ihist]->SetBinContent(1, 0);
      cout << " set to zero " << runHistName << " " << hRunLatePeakSum[ihist]->GetBinContent(1) << endl;
    }
    // dont need this fout->Add(hRunLatePeakSum[ihist]);
  }

  // get rid of first bin in gaina
  gainSumDir->Purge(1);
  TList *gainList = gainSumDir->GetListOfKeys();
  TIter next(gainList);
  TKey *key;
  while (TKey *key = (TKey *)next())
  {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1D"))
      continue;
    TH1D *h = (TH1D *)key->ReadObj();
    for (int ibin = 1; ibin < 5; ++ibin)
      h->SetBinContent(ibin, 0);
  }
  gainSumDir->Purge();
  fout->Purge(1);
  runSumDir->ls();
  
  cout << "summary finished "
       << " total pass " << totalPass << " maxFiles  " << maxFiles << " files written to " << fout->GetName() << endl;

 // calculate mean hits from waveforms 
  for( int idet=0; idet< hRunHitWave.size(); ++idet )  {
    TString histName;
    histName.Form("UnNormedHitWaveChan%i",idet);
    TH1D *hist;
    runSumDir->GetObject(histName, hist);
    if(hist) printf(" idet %i unnormed %.3E inte %.3E \n",idet,  hist->Integral(), hRunHitWave[idet]->Integral());
    else printf(" did not find %s \n",histName.Data());
  }

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
          if (!isinf(vecQsum[ic][jt]) && vecQsum[ic][jt] > 0 && fileDatime[jt].Convert() < dopeTime.Convert())
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
  // normalize to qpe and effOther
  for (unsigned ic = 0; ic < vecQsum.size(); ++ic)
  {
    printf("QSUM ch %i size %lu \n", ic, vecQsum[ic].size());
    for (unsigned ih = 0; ih < vecQsum[ic].size(); ++ih)
    {
      double qpe = vecQPE[ic][ih];
      if (qpe <= 1. || isnan(qpe) || isinf(qpe))
        qpe = 1.;
      printf(" \t\t QSUM chan %i file%i qpe %.3E effOther %.3E  qsum %f  new  %f \n",
             ic, ih, vecQPE[ic][ih], effOther[ic], vecQsum[ic][ih], vecQsum[ic][ih] / qpe / effOther[ic]);
      vecQsum[ic][ih] = vecQsum[ic][ih] / qpe / effOther[ic];
      vecEQsum[ic][ih] = vecEQsum[ic][ih] / qpe / effOther[ic];
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
/*
// apply norms
for (unsigned ic = 0; ic < vecQsum.size(); ++ic)
{
  for (unsigned jf = 0; jf < vecQsum[ic].size(); ++jf)
  {
    vecQsum[ic][jf] /= normQsum[ic];
    vecEQsum[ic][jf] /= normQsum[ic];
  }
}
// apply norms
for (unsigned ic = 0; ic < vecQPE.size(); ++ic)
{
  for (unsigned jf = 0; jf < vecQPE[ic].size(); ++jf)
  {
    vecQPE[ic][jf] /= normQPE[ic];
    vecEQPE[ic][jf] /= normQPE[ic];
  }
}
*/

// get slopes from graph
/*
unsigned nslope = gslope->GetN();
printf(" nslope %u \n", nslope);
// the graph starts with ichn = 3
for (unsigned ic = 0; ic < CHANNELS;  ++ic)
{
  vSlope[ic].push_back(0.);
  vESlope[ic].push_back(0);
}
for (unsigned ic = 0; ic < nslope ; ++ic)
{
  unsigned ilast = vSlope[ic].size()-1;
  unsigned ichan = unsigned(gslope->GetPointX(ic));
  double yval = gslope->GetPointY(ic);
  double yerr = gslope->GetErrorY(ic);
  printf("filling chan %u file %u %f %f \n", ichan, ilast , yval, yerr);

  vSlope[ichan][ilast] = gslope->GetPointY(ic);
  vESlope[ichan][ilast] = gslope->GetErrorY(ic);
}

printf("\t slopes %lu \n", vSlope.size());
for (unsigned j = 0; j < vSlope.size(); ++j)
{
  for (unsigned k = 0; k < vSlope[j].size(); ++k)
    printf("chan %u file %u %f %f \n", j, k, vSlope[j][k], vESlope[j][k]);
}
*/

// get slopes from graph
/*
unsigned nslope = gslope->GetN();
printf(" nslope %u \n", nslope);
// the graph starts with ichn = 3
for (unsigned ic = 0; ic < CHANNELS;  ++ic)
{
  vSlope[ic].push_back(0.);
  vESlope[ic].push_back(0);
}
for (unsigned ic = 0; ic < nslope ; ++ic)
{
  unsigned ilast = vSlope[ic].size()-1;
  unsigned ichan = unsigned(gslope->GetPointX(ic));
  double yval = gslope->GetPointY(ic);
  double yerr = gslope->GetErrorY(ic);
  printf("filling chan %u file %u %f %f \n", ichan, ilast , yval, yerr);

  vSlope[ichan][ilast] = gslope->GetPointY(ic);
  vESlope[ichan][ilast] = gslope->GetErrorY(ic);
}

printf("\t slopes %lu \n", vSlope.size());
for (unsigned j = 0; j < vSlope.size(); ++j)
{
  for (unsigned k = 0; k < vSlope[j].size(); ++k)
    printf("chan %u file %u %f %f \n", j, k, vSlope[j][k], vESlope[j][k]);
}
*/
