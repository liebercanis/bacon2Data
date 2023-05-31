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
std::string sdate;
unsigned nchan = 13;
vector<double> normQsum;
vector<double> normQPE;
TDatime dopeTime;
Long64_t maxFiles;
double runNorm;

std::vector<TString> fileList;
std::vector<double> filenum;
std::vector<double> efilenum;
std::vector<double> fileTime;
std::vector<TDatime> fileDatime;
TBFile *bf;
TBEventData *eventData;
TTree *runTree;
std::vector<vector<double>> vecQsum;
std::vector<vector<double>> vecEQsum;
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
std::vector<TH1D *> hSumWave;
std::vector<vector<TH1D *>> vRunQSum;
std::vector<vector<TH1D *>> vRunPeakSum;
std::vector<vector<TH1D *>> vRunHitWave;
std::vector<vector<TH1D *>> vRunPeakWave;
std::vector<TH1D *> hRunSumPeakWave;
std::map<int, TH1D *> histMap;
std::map<int, TH1D *> histqMap;
std::map<int, TH1D *> histpMap;
TDirectory *sumDir;
TDirectory *waveSumDir;
TDirectory *qpeSumDir;
TDirectory *peakSumDir;
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

void addRunSumHistos()
{
  cout << "in addRunSumHistos" << endl;
  int nbinsx = 0;
  bool first = true;
  for (unsigned ifile = 0; ifile < maxFiles; ++ifile)
  {

    TString fullName = dirNameSlash + fileList[ifile];
    fin = new TFile(fullName);
    fin->GetObject("sumDir", sumDir);
    if (!sumDir)
    {
      printf("addRunSumHistos file %s no sumDir!!! \n", fileList[ifile].Data());
      continue;
    }

    TList *sumList = sumDir->GetListOfKeys();
    // cout << " >>>>> addRunSumHistos file  " << fileList[ifile] << " sumDir size " << sumList->GetSize() << endl;
    //  sumDir->ls();

    TIter next(sumList);
    TKey *key;
    // printf("addsumDirHistos %u \n", sumList->GetEntries());
    while (TKey *key = (TKey *)next())
    {
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TH1D"))
        continue;
      TH1D *h = (TH1D *)key->ReadObj();
      std::string name = string(h->GetName());

      // save summed waves
      // summed waves
      TH1D *hClone;
      if (name.find("sumWave") != std::string::npos && name.find("Bad") == std::string::npos)
      {
        // cout << " sumWave clone " << name << " file " << ifile << endl;
        hClone = (TH1D *)h->Clone(Form("%s-file%i", h->GetName(), ifile));
        hSumWave.push_back(hClone);
        waveSumDir->Add(hClone);
      }
      if (name.find("sumHitWave") != std::string::npos)
      {
//        cout << " \t\t AAAAAAA "  << endl;
        string chan = name.substr(name.find_last_of("e") + 1);
        int ichan = stoi(chan);
        hClone = (TH1D *)h->Clone(Form("RunHitWave-file%u-chan%i", ifile, ichan));
        waveSumDir->Add(hClone);
        vRunHitWave[ichan].push_back(hClone);
      }
      if (name.find("sumPeakWave") != std::string::npos) 
      {
        string chan = name.substr(name.find_last_of("e") + 1);
        int ichan = stoi(chan);
        hClone = (TH1D *)h->Clone(Form("RunPeakWave-file%u-chan%i", ifile, ichan));
        waveSumDir->Add(hClone);
        vRunPeakWave[ichan].push_back(hClone);
//        for (int ibin = 0; ibin < hClone->GetNbinsX(); ++ibin)
//        {
        if (ifile == 0)
        {
          printf("1st ichan = %i, size of hRunSumPeakWave = %lu \n", ichan, hRunSumPeakWave.size());
          hRunSumPeakWave.push_back(hClone);
        }
        else
        {
          //printf("else ichan = %i, size of hRunSumPeakWave = %lu \n", ichan, hRunSumPeakWave.size());
          hRunSumPeakWave[ichan]->Add(hClone);
        }
//        }
      }
      // get PeakChan by channel
      if (name.find("QPeakChan") != std::string::npos) {        string chan = name.substr(name.find_last_of("n") + 1);
        int ichan = stoi(chan);
        // cout << name << "string chan" << chan << " int " << ichan << endl;

        TString cloneName;
        cloneName.Form("runPeakSumCh%i-file%i", ichan, ifile);
        TH1D *hpAdd = (TH1D *)h->Clone(cloneName);
        // cout << " +vRunQSum " << ichan << " " << fileList[ifile] << " " << h->GetName()  << " " << hAdd->GetName() << endl;
        hpAdd->SetMarkerStyle(20);
        hpAdd->SetMarkerSize(0.5);
        vRunPeakSum[ichan].push_back(hpAdd);
        peakSumDir->Add(hpAdd);
      }
      // get QSumChan by channel
      if (name.find("QSumChan") != std::string::npos) {
        string chan = name.substr(name.find_last_of("n") + 1);
        int ichan = stoi(chan);
        // cout << name << "string chan" << chan << " int " << ichan << endl;

        TString cloneName;
        cloneName.Form("runQSumCh%i-file%i", ichan, ifile);
        TH1D *hqAdd = (TH1D *)h->Clone(cloneName);
        // cout << " +vRunQSum " << ichan << " " << fileList[ifile] << " " << h->GetName()  << " " << hAdd->GetName() << endl;
        hqAdd->SetMarkerStyle(20);
        hqAdd->SetMarkerSize(0.5);
        vRunQSum[ichan].push_back(hqAdd);
        qpeSumDir->Add(hqAdd);
      }
        // fin->Close();
      }
    }
  }

  void QPEFits()
  {

    cout << " ***** QPEFits ***** " << endl;

    // for (unsigned ichan = 0; ichan < vRunQSum.size(); ++ichan)
    for (unsigned ichan = 0; ichan < 9; ++ichan)
    {
      cout << " QPEFits " << ichan << " num histos " << vRunQSum[ichan].size() << endl;
      for (unsigned ihist = 0; ihist < vRunQSum[ichan].size(); ++ihist)
      {
        TH1D *hclone = vRunQSum[ichan][ihist];

        // fit QPE   
        double xlow = 0.;
        double xhigh = 0;
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
          if (ichan == 0 || ichan == 1 || ichan == 2 || ichan == 5)
          {
            xlow = 4000.;
            xhigh = 6000.;
          }
          else if (ichan == 4 || ichan == 6 || ichan == 7 || ichan == 8)
          {
            xlow = 10000.;
            xhigh = 14000.;
          }
          else if (ichan == 12)
          {
            xlow = 400.;
            xhigh = 1200.;
          }
        }
        if (ichan == 9 || ichan == 10 || ichan == 11)
        {
          return;
        }
        hclone->Fit("gaus", "Q", " ", xlow, xhigh);
        TF1 *gfit = (TF1 *)hclone->GetListOfFunctions()->FindObject("gaus");
        double par1 = 0;
        double epar1 = 0;
        double par2 = 0;
        double epar2 = 0;
        if (gfit)
        {
        par1 = gfit->GetParameter(1);
        epar1 = gfit->GetParError(1);
        if (epar1 > par1) epar1 = par1;
        par2 = gfit->GetParameter(2);
        epar2 = gfit->GetParError(2);
        if (epar2 > par2) epar2 = par2;
        }

        printf(" fit QPE %u  mean %f %f  sigma %f %f  \n", ichan, par1, epar1, par2, epar2);
        hQPEChan->SetBinContent(ichan, par1);
        hQPEChan->SetBinError(ichan, epar1);
        hQPESigmaChan->SetBinContent(ichan, par2);
        hQPESigmaChan->SetBinError(ichan, epar2);
        vecQPE[ichan].push_back(par1);
        vecEQPE[ichan].push_back(epar1);
        vecQPESigma[ichan].push_back(par2);
        vecEQPESigma[ichan].push_back(epar2);
      }
    }
  }

  void setTimeGraph(TMultiGraph * mg, TString ylabel)
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

  TDatime getTime(int ifile, Long64_t ievent=0)
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

      printf(" file %i event %lli FileYear = %u , FileMonth = %u , FileDay = %u , FileHour = %u , FileMin = %u , FileSec = %u \n", ifile,ievent,eventData->year+1900, eventData->mon+1, eventData->day, eventData->hour, eventData->min, eventData->sec);
      datime.Set(eventData->year, eventData->mon+1, eventData->day, eventData->hour, eventData->min, eventData->sec);
   }

    return datime;
  }

  void fileLoop()
  {
    printf("fileLoop over %lld files \n", maxFiles);  // DEF would be nice to put in a way to look at the last maxFiles files
    for (unsigned ifile = 0; ifile < maxFiles; ++ifile)
    {
      TString fullName =  dirNameSlash + fileList[ifile];
      fin = new TFile(fullName);

      fin->GetObject("sumDir", sumDir);
      if (!sumDir)
      {
        printf(" sumDir not found in file %s \n", fin->GetName());
        continue;
      }
      printf(" sumDir for file %s \n", fin->GetName());
      // sumDir->ls();

      fin->GetObject("RunTree", runTree); 

      // use bf->modified a std string
      /* migrate fit from post.C:        TF1 *g = (TF1*)hHitSum[i]->GetListOfFunctions()->FindObject("expo");*/

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
        if(!runTree)
          printf("file has not RunTree %s \n", fin->GetName());
        else
          runTree->SetBranchAddress("eventData",&eventData);
      }

      if (hqsum && hqprompt) {
        TH1D *hq = (TH1D *)hqsum->Clone(Form("hqsum%i", ifile));
        fout->Add(hq);
        TH1D *hp = (TH1D *)hqsum->Clone(Form("hqprompt%i", ifile));
        fout->Add(hp);
      }
      eventCount = NULL;
      cout << " get event Count " << endl;
      fin->GetObject("eventcount", eventCount);
      runNorm = 1;
      if (eventCount)
      {
        // 0 = ntriggers, 1 = npass
        runNorm = eventCount->GetBinContent(0) / eventCount->GetBinContent(1);
        printf("normalize bin 0 %f bin 1 %f runNorm %f \n ",eventCount->GetBinContent(0), eventCount->GetBinContent(1),runNorm);
        printf(" eventCount for file %s  entries %E \n", fin->GetName(), eventCount->GetEntries());
        TH1D *hevcount = (TH1D *)eventCount->Clone(Form("eventCount%i", ifile));
        fout->Add(hevcount);
      }
      else
        printf("NO EVENT COUNT IN FILE %s\n", fin->GetName());

      filenum.push_back(double(ifile));
      efilenum.push_back(0);
      TDatime datetime = getTime(ifile);
      fileDatime.push_back(datetime);
      fileTime.push_back(datetime.Convert());

      cout << " TIME FILE  " << ifile << " " << fileList[ifile] << " modified "
           << " TDatime " << datetime.AsString() << " as int " << datetime.Convert() << endl;

      /*
      for (unsigned ic = 0; ic < vecQPE.size() ; ++ic)
      {
        for (int ifile = 0; ifile < vecQPE[ic].size(); ++ifile)
          printf(" chan %u file %i %f\n", ic, ifile, vecQPE[ic][ifile]);
      }
      */

      if (eventCount)
      {
        for (int i = 0; i <= eventCount->GetNbinsX(); ++i)
        {
          double norm = eventCount->GetBinContent(i + 1) / eventCount->GetBinContent(0);
          printf(" %i count %0f norm %f \n", i, eventCount->GetBinContent(i + 1), norm);
        }
      }
      if (hqsum && hqprompt) {
        for (int i = 0; i < hqsum->GetNbinsX() - 1; ++i)
        {
          double norm = 1.0;
          // now overflow is every event and qsum bin is
          if (eventCount)
          norm = eventCount->GetBinContent(i) / eventCount->GetBinContent(0);
        double val = 0;
        double eval = 0;
        if (!isnan(hqsum->GetBinContent(i + 1)) && !isinf(hqsum->GetBinContent(i + 1)))
          {
            val = hqsum->GetBinContent(i + 1) / norm;
            eval = hqsum->GetBinError(i + 1) / norm;
          }
          vecQsum[i].push_back(val);
          vecEQsum[i].push_back(eval);
          cout << "chan " << i << " qsum " << hqsum->GetBinContent(i + 1) << " norm " << norm << " size " << vecQsum[i].size() << endl;
        }
        printf("  summary file %u  %s chan 6 %f  \n", ifile, fileList[ifile].Data(), vecQsum[6][ifile]);
      
      }
      fin->Close();
    } // end loop over files
  }

  void QPESumFits()
  {

    cout << " ***** QPEFits ***** " << histqMap.size() << endl;

    for (std::map<int, TH1D *>::iterator it = histqMap.begin(); it != histqMap.end(); ++it)
    {
      std::cout << it->first << " => " << it->second->GetName() << '\n';
      int ichan = it->first;
      TH1D *hclone = (TH1D *)it->second;

      // fit QPE   
      double xlow = 0.;
      double xhigh = 0;
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
        if (ichan == 0 || ichan == 1 || ichan == 2 || ichan == 5)
        {
          xlow = 4000.;
          xhigh = 6000.;
        }
        else if (ichan == 4 || ichan == 6 || ichan == 7 || ichan == 8)
        {
          xlow = 10000.;
          xhigh = 14000.;
        }
        else if (ichan == 12)
        {
          xlow = 400.;
          xhigh = 1200.;
        }
      }
      if (ichan == 9 || ichan == 10 || ichan == 11)
      {
        return;
      }
      hclone->Fit("gaus", "Q", " ", xlow, xhigh);
      TF1 *gfit = (TF1 *)hclone->GetListOfFunctions()->FindObject("gaus");
      double par1 = 0;
      double epar1 = 0;
      double par2 = 0;
      double epar2 = 0;
      if (gfit)
      {
        par1 = gfit->GetParameter(1);
        epar1 = gfit->GetParError(1);
        if (epar1 > par1) epar1 = par1;
        par2 = gfit->GetParameter(2);
        epar2 = gfit->GetParError(2);
        if (epar2 > par2) epar2 = par2;
      }

      printf(" QPE %u  mean %f %f  sigma %f %f  \n", ichan, par1, epar1, par2, epar2);
      hQPEChan->SetBinContent(ichan, par1);
      hQPEChan->SetBinError(ichan, epar1);
      hQPESigmaChan->SetBinContent(ichan, par2);
      hQPESigmaChan->SetBinError(ichan, epar2);
      vecQPE[ichan].push_back(par1);
      vecEQPE[ichan].push_back(epar1);
      vecQPESigma[ichan].push_back(par2);
      vecEQPESigma[ichan].push_back(epar2);
    }
  }

  void PeakSumFits()
  {

    cout << " ***** PeakFits ***** " << histpMap.size() << endl;

    for (std::map<int, TH1D *>::iterator it = histpMap.begin(); it != histpMap.end(); ++it)
    {
      std::cout << it->first << " => " << it->second->GetName() << '\n';
      int ichan = it->first;
      TH1D *hclone = (TH1D *)it->second;

      // fit QPE   
      double xlow = 0.;
      double xhigh = 0;
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
        if (ichan == 0 || ichan == 1 || ichan == 2 || ichan == 5)
        {
          xlow = 100.;
          xhigh = 175.;
        }
        else if (ichan == 4 || ichan == 6 || ichan == 7 || ichan == 8)
        {
          xlow = 200.;
          xhigh = 275.;
        }
        else if (ichan == 12)
        {
          xlow = 100.;
          xhigh = 250.;
        }
      }
      if (ichan == 9 || ichan == 10 || ichan == 11)
      {
        return;
      }
      hclone->Fit("gaus", "Q", " ", xlow, xhigh);
      TF1 *gfit = (TF1 *)hclone->GetListOfFunctions()->FindObject("gaus");
      double par1 = 0;
      double epar1 = 0;
      double par2 = 0;
      double epar2 = 0;
      if (gfit)
      {
        par1 = gfit->GetParameter(1);
        epar1 = gfit->GetParError(1);
        if (epar1 > par1) epar1 = par1;
        par2 = gfit->GetParameter(2);
        epar2 = gfit->GetParError(2);
        if (epar2 > par2) epar2 = par2;
      }

      printf(" Peak %u  mean %f %f  sigma %f %f  \n", ichan, par1, epar1, par2, epar2);
      hQPEChan->SetBinContent(ichan, par1);
      hQPEChan->SetBinError(ichan, epar1);
      hQPESigmaChan->SetBinContent(ichan, par2);
      hQPESigmaChan->SetBinError(ichan, epar2);
      vecQPE[ichan].push_back(par1);
      vecEQPE[ichan].push_back(epar1);
      vecQPESigma[ichan].push_back(par2);
      vecEQPESigma[ichan].push_back(epar2);
    }
  }

  void addsumDirHistos()
  {
    cout << "in addsumDirHistos" << endl;
    int nbinsx = 0;
    bool first = true;
    for (unsigned ifile = 0; ifile < maxFiles; ++ifile)
    {

      TString fullName = dirNameSlash + fileList[ifile];
      fin = new TFile(fullName);
      fin->GetObject("sumDir", sumDir);
      if (!sumDir)
      {
        printf("addsumDirHistos file %s no sumDir!!! \n", fileList[ifile].Data());
        continue;
      }

      TList *sumList = sumDir->GetListOfKeys();
      //cout << " >>>>> addsumDirHistos file  " << fileList[ifile] << " sumDir size " << sumList->GetSize() << endl;
      // sumDir->ls();

      TIter next(sumList);
      TKey *key;
      // printf("addsumDirHistos %u \n", sumList->GetEntries());
      while (TKey *key = (TKey *)next())
      {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TH1D"))
          continue;
        TH1D *h = (TH1D *)key->ReadObj();
        std::string name = string(h->GetName());

        // Add QSumChan by channel
        if (name.find("QSumChan") != std::string::npos)
        {
        // cout << fileList[ifile] << " " << h->GetName() << endl;
          string chan = name.substr(name.find_last_of("n") + 1);
          int ichan = stoi(chan);
          TString cloneqName;
          cloneqName.Form("AddQSumCh%i", ichan);
          TH1D *hqAdd = (TH1D *)h->Clone(cloneqName);
          if (first) // create summed histos in runQSum and runPeakSum arrays with maps
          {
            if (nbinsx == 0)
              nbinsx = hqAdd->GetNbinsX();
            qpeSumDir->Add(hqAdd);
            TH1D *fileqHist = NULL;
            qpeSumDir->GetObject(cloneqName, fileqHist);
            histqMap.insert(std::pair<int, TH1D *>(ichan, fileqHist));
            // cout << " QSumChan clone " << name << " chan " << ichan << " NbinsX " << hAdd->GetNbinsX() << " integral " << histMap.at(ichan)->Integral() << endl;
          }
          else if (hqAdd->GetNbinsX() == nbinsx) // add
          {
            histqMap.at(ichan)->Add(hqAdd);
            // cout << " add QSumChan " << name << " chan " << ichan << " NbinsX " << hAdd->GetNbinsX() << " integral " << histMap.at(ichan)->Integral() << endl;
          }
          else
          {
            cout << " Warning! skip QSumChan  " << name << " chan " << ichan << " NbinsX " << nbinsx << " , " << hqAdd->GetNbinsX() << endl;
          }
        }
        // Add QSumChan by channel
        else if (name.find("QPeakChan") != std::string::npos)
        {
          // cout << fileList[ifile] << " " << h->GetName() << endl;
          string chan = name.substr(name.find_last_of("n") + 1);
          int ichan = stoi(chan);
          TString clonepName;
          clonepName.Form("AddPeakSumCh%i", ichan);
          TH1D *hpAdd = (TH1D *)h->Clone(clonepName);
          if (first) // create summed histos in runQSum and runPeakSum arrays with maps
          {
            if (nbinsx == 0)
              nbinsx = hpAdd->GetNbinsX();
            peakSumDir->Add(hpAdd);
            TH1D *filepHist = NULL;
            peakSumDir->GetObject(clonepName, filepHist);
            histpMap.insert(std::pair<int, TH1D *>(ichan, filepHist));
            // cout << " QSumChan clone " << name << " chan " << ichan << " NbinsX " << hAdd->GetNbinsX() << " integral " << histMap.at(ichan)->Integral() << endl;
          }
          else if (hpAdd->GetNbinsX() == nbinsx) // add
          {
            histpMap.at(ichan)->Add(hpAdd);
            // cout << " add QSumChan " << name << " chan " << ichan << " NbinsX " << hAdd->GetNbinsX() << " integral " << histMap.at(ichan)->Integral() << endl;
          }
          else
          {
            cout << " Warning! skip QSumChan  " << name << " chan " << ichan << " NbinsX " << nbinsx << " , " << hpAdd->GetNbinsX() << endl;
          }
        }

      } // loop over sumDir
      if (first)
      {
        for (std::map<int, TH1D *>::iterator it = histqMap.begin(); it != histqMap.end(); ++it)
          std::cout << it->first << " => " << it->second->GetName() << '\n';
        for (std::map<int, TH1D *>::iterator it = histpMap.begin(); it != histpMap.end(); ++it)
          std::cout << it->first << " => " << it->second->GetName() << '\n';
      }
      first = false;
      // fin->Close();
    }
  }

  void fitSlopes()
  {
    for (unsigned ichan = 0; ichan < 13; ++ichan) 
    {
      printf("Make slope graph.  Number of files = %lu Number of channels = %lu \n", filenum.size(), vRunPeakWave.size());

      for (unsigned ih = 0; ih < vRunPeakWave[ichan].size(); ++ih)
      {
        if(vRunPeakWave[ichan][ih]==NULL){
          printf("skipping vRunPeakWave chan %i file %i \n", ichan, ih);
          continue;
        }
        double xlow = 0.;
        double xhigh = 0;
        double MicroSecPerDac = 0;  // DEF added in switch for different digitizer MicroSecPerDac
        if (theDataType == SIS)
        { // For SIS data
          xlow = 120.;
          xhigh = 500.;
          MicroSecPerDac = 8. / 1000.; // 8 ns bins
        }
        else
        { // For CAEN data
          xlow = 1000.;
          xhigh = 3000.;
          MicroSecPerDac = 2. / 1000.; // 2 ns bins
        }

        TH1D *hfitwave = (TH1D *)vRunPeakWave[ichan][ih]->Clone("fitwave");
        hfitwave->GetListOfFunctions()->Clear();
        hfitwave->SetDirectory(nullptr);
        TF1 *gslopefit = NULL;
        hfitwave->Fit("expo", "LQ", " ", xlow, xhigh);
        gslopefit = (TF1 *)hfitwave->GetListOfFunctions()->FindObject("expo");
        double mfit = 0;
        double time = 0;
        double etime = 0;
        if (gslopefit)
        {
          mfit = gslopefit->GetParameter(1); 
          time = -1.0 * MicroSecPerDac / mfit;
          if (time > 40) time = 0;
          etime = gslopefit->GetParError(1) *(-1)*time / mfit;
          if (etime > time) etime = time;
          printf(" Slope for file %u, channel %i = %f tau = %f tau_error = %f \n", ih, ichan, mfit, time, etime);
        }
        vSlope[ichan].push_back(time);
        vESlope[ichan].push_back(etime);
        // printf(" %u %s %s chan %i file %i time %f %f\n", ih, hSumWave[ih]->GetName(), hSumHitWave[ih]->GetName(), ichan, ifile,time,etime);
        delete hfitwave;
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
    TSystemDirectory dir(dirName, dirName);  // TSystemDirectory
    TList *files = dir.GetListOfFiles();     //
    TIter next(files);
    TSystemFile *file;
    while ((file = (TSystemFile *)next()))
    {
      string name = string(file->GetName());
      // cout << name << endl;
      string exten = name.substr(name.find_last_of(".") + 1);
      if (exten != string("root"))
        continue;
      if (name.find("root") != std::string::npos)
        fileList.push_back(TString(name.c_str()));
    }
    return fileList.size();
  }

  int main(int argc, char *argv[])
  {
    cout << "executing " << argv[0] << " make summary plots  " << endl;
    printf(" usage: summary  <type=caen/sis> <dir> <max files 0=all, negative=last>  \n ");
    if (argc < 2)
    {
      printf("reguire <type> <dir> args\n");
      exit(0);
    }
    TString tag("run");

    printf(" input args %i \n ",argc);
    for (int jarg = 1; jarg < argc; ++jarg)
      printf(" %i= %s ", jarg, argv[jarg]);
    printf("\n");

    if (TString(argv[1]).Contains("caen") || TString(argv[1]).Contains("CAEN"))
      theDataType = CAEN;
    else
      theDataType = SIS;

    dirName = TString(argv[2]);
    dirNameSlash = TString(argv[2]) + TString("/");

    printf(" >>>> datatype %i (0=sis,1=caen) dir %s <<<<\n", theDataType, dirName.Data());

    dopeTime = TDatime(2023, 3, 9, 22, 0, 0);

    unsigned nfiles = countFiles();
    if (countFiles() == 0)
      exit(0);

    printf(" for %s found %lu files \n", tag.Data(), fileList.size());
    maxFiles = fileList.size();
    if (argc > 3)
    {
      maxFiles = atoi(argv[3]);
    }

    sdate = currentDate();
    cout << " starting summary for  " << maxFiles << " files on " << sdate << endl;

    fout = new TFile(Form("summary-type-%i-dir-%s-%s.root",theDataType,dirName.Data(),sdate.c_str()), "recreate");
    waveSumDir = fout->mkdir("WaveSumDir");
    qpeSumDir = fout->mkdir("qpeSumDir");
    peakSumDir = fout->mkdir("peakSumDir");
    fout->cd();

    hQPEChan = new TH1D("QPEChan", "QPE  by channel", 12, 0, 12);
    hQPESigmaChan = new TH1D("QPESigmaChan", "QPE  by channel", 12, 0, 12);
    hQPEChan->Sumw2();
    hQPESigmaChan->Sumw2();

    vecQsum.resize(nchan);
    vecEQsum.resize(nchan);
    vecPeaksum.resize(nchan);
    vecEPeaksum.resize(nchan);
    vecQPE.resize(nchan);
    vecEQPE.resize(nchan);
    vecPeak.resize(nchan);
    vecEPeak.resize(nchan);
    vecQPESigma.resize(nchan);
    vecEQPESigma.resize(nchan);
    vSlope.resize(nchan);
    vESlope.resize(nchan);
    vRunQSum.resize(nchan);
    vRunPeakSum.resize(nchan);
    vRunHitWave.resize(nchan);
    vRunPeakWave.resize(nchan);
//    hRunSumPeakWave.resize(nchan);

//    for (unsigned ichan = 0; ichan < nchan; ++ichan)
//      hRunSumPeakWave[ichan] = NULL;

    fileLoop();
    // for (unsigned ic = 0; ic < nchan; ++ic)
    printf(" >>> files processed << %li \n", filenum.size());
    for (unsigned jfile = 0; jfile < filenum.size(); ++jfile)
    {
      int ifile = int(filenum[jfile]);
      printf(" %i %s \n", ifile, fileList[ifile].Data());
    }

    addsumDirHistos();
    addRunSumHistos();
    // for (std::map<int, TH1D *>::iterator it = histMap.begin(); it != histMap.end(); ++it)
    //   std::cout << it->first << " => " << it->second->GetName() << '\n';
    QPEFits();
    // PeakSumFits();

    //normalize
    if (eventCount)
    {
    
      // assumes I had normalized to triggers in previous step
      printf("\t\t\t  run norm %f %lu \n", runNorm,hRunSumPeakWave.size());
      for (int ihist = 0; ihist < int(hRunSumPeakWave.size()); ++ihist)
      {
        for (int ibin = 0; ibin < hRunSumPeakWave[ihist]->GetNbinsX(); ++ibin)
          hRunSumPeakWave[ihist]->SetBinContent(ibin, hRunSumPeakWave[ihist]->GetBinContent(ibin) * runNorm);
        fout->Add(hRunSumPeakWave[ihist]);
        }
    }

        // call function to fit slopes and fill vSlope, vESlope
      if (filenum.size() > 0) {
        fitSlopes();
        makeGraphs();
    }

    // report
    /*for (unsigned it = 0; it < fileTime.size(); ++it)
      if (fileDatime[it].Convert() < dopeTime.Convert())
        cout << " before " << fileDatime[it].AsString() << "  " << fileList[it] << endl;
      else
        cout << " after  " << fileDatime[it].AsString() << "  " << fileList[it] << endl;
        */
    // fout->ls();


    cout << "summary closing .... " << endl;
//    fout->ls();
    fout->Write();
    fout->Close();
    cout << "summary finished " << maxFiles << " " << fout->GetName() << endl;
    exit(0);
  }

  /*  put this all at bottom */
  void makeGraphs()
  {
    printf(" \t\t make graphs %lu \n", filenum.size());
    cout << " \n\t ******* makeGraphs *****  "<< endl;
    cout << " \n\t vecQsum "<< vecQsum.size() << " vecQPE  "<< vecQPE.size()<< endl;
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
    for (unsigned ic = 7; ic < 9; ++ic)
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
      mgslope->Add(graphSlope[ilast]);
    }
    ylabel.Form(" fitted Lifetime [microsec] ");
    setTimeGraph(mgslope, ylabel);
    TCanvas *canSlope = new TCanvas(Form("SlopeFit-%s", sdate.c_str()), Form("SlopeFit-%s", sdate.c_str()));
    mgslope->Draw("ap");
    mgslope->GetYaxis()->SetRangeUser(0.5,2.5);
    gPad->Update();
    canSlope->BuildLegend();
    canSlope->SetGrid();
    canSlope->Print(".png");
    fout->Append(canSlope);

    // one graph per channel
    vector<TGraphErrors *> gqsum;
    TMultiGraph *mgsum = new TMultiGraph();
    for (unsigned ic = 0; ic < vecQsum.size(); ++ic)
    {
      // cout << " add " << ic << endl;
      gqsum.push_back(new TGraphErrors(filenum.size(), &fileTime[0], &(vecQsum[ic][0]), &efilenum[0], &(vecEQsum[ic][0])));
      gqsum[ic]->SetName(Form("qsumChan%i", ic));
      gqsum[ic]->SetTitle(Form("qsum-chan-%i", ic));
      gqsum[ic]->SetMarkerSize(1);
      gqsum[ic]->SetMarkerColor(myColor[ic]);
      gqsum[ic]->SetMarkerStyle(myStyle[ic]);
      fout->Add(gqsum[ic]);
      if (ic == 6 || ic == 7 || ic == 8)
        mgsum->Add(gqsum[ic]);
    }
    // overlay all channel graphs on canvas
    int ndiv = 10 + 100 * 5 + 10000 * 3;
    ylabel.Form("integrated charge (normed) ");
    setTimeGraph(mgsum, ylabel);
    TCanvas *can = new TCanvas(Form("Qsummary-%s", sdate.c_str()), Form("Qsummary-%s", sdate.c_str()));
    mgsum->Draw("ap");
    gPad->Update();
    can->BuildLegend();
    can->SetGrid();
    can->Print(".png");
    fout->Append(can);

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
    for (unsigned ic = 0; ic < nchan; ++ic)
    {
      // cout << " add " << ic << endl;
      gqsumUn.push_back(new TGraphErrors(filenum.size(), &fileTime[0], &(vecQsum[ic][0]), &efilenum[0], &(vecEQsum[ic][0])));
      gqsumUn[ic]->SetName(Form("qsumChanUn%i", ic));
      gqsumUn[ic]->SetTitle(Form("qsum-unnormalized-chan-%i", ic));
      gqsumUn[ic]->SetMarkerSize(1);
      gqsumUn[ic]->SetMarkerColor(myColor[ic]);
      gqsumUn[ic]->SetMarkerStyle(myStyle[ic]);
      fout->Add(gqsumUn[ic]);
    }
    ylabel.Form("integrated charge");
    TMultiGraph *mgL1 = new TMultiGraph();
    mgL1->Add(gqsumUn[6]);
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
    mgL2->Add(gqsumUn[3]);
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
  }

  // get slopes from graph
  /*
  unsigned nslope = gslope->GetN();
  printf(" nslope %u \n", nslope);
  // the graph starts with ichn = 3
  for (unsigned ic = 0; ic < nchan;  ++ic)
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
  for (unsigned ic = 0; ic < nchan;  ++ic)
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
