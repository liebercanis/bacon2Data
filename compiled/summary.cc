// revised June 7 to close files after reading
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
double ntotal;
double npass;
vector<double> filePass;
int totalPass;
vector<double> effOther;

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
std::vector<vector<TH1D *>> vRunQSum;
std::vector<vector<TH1D *>> vRunPeakSum;
std::vector<vector<TH1D *>> vRunHitWave;
std::vector<vector<TH1D *>> vRunPeakWave;
std::vector<TH1D *> hSumWave;
std::vector<TH1D *> hRunSumPeakWave;
std::vector<TH1D *> hRunSumHitWave;
TDirectory *sumDir;
TDirectory *fitSumDir;
TDirectory *waveSumDir;
TDirectory *qpeSumDir;
TDirectory *peakSumDir;
TDirectory *runSumDir;
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
double nPhotons = 50.E3 * 5.486;

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

void QPEFits()
{

  cout << " ***** QPEFits ***** " << endl;

  // for (unsigned ichan = 0; ichan < vRunQSum.size(); ++ichan)
  for (unsigned ichan = 0; ichan < nchan; ++ichan)
  {

    cout << " QPEFits " << ichan << " num histos " << vRunQSum[ichan].size() << endl;
    for (unsigned ihist = 0; ihist < vRunQSum[ichan].size(); ++ihist)
    {
      TString cloneName;
      TH1D *hclone = NULL;
      cloneName.Form("runQSumCh%iFile%i", ichan, ihist);
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
        continue;
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
        if (epar1 > par1)
          epar1 = par1;
        par2 = gfit->GetParameter(2);
        epar2 = gfit->GetParError(2);
        if (epar2 > par2)
          epar2 = par2;
      }
      else
        printf("Fit to chan %u file  %u fails\n", ichan, ihist);

      hQPEChan->SetBinContent(ichan, par1);
      hQPEChan->SetBinError(ichan, epar1);
      hQPESigmaChan->SetBinContent(ichan, par2);
      hQPESigmaChan->SetBinError(ichan, epar2);
      vecQPE[ichan].push_back(par1);
      vecEQPE[ichan].push_back(epar1);
      vecQPESigma[ichan].push_back(par2);
      vecEQPESigma[ichan].push_back(epar2);
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
  printf("fileLoop over %lld files \n", maxFiles); // DEF would be nice to put in a way to look at the last maxFiles files
  for (unsigned ifile = 0; ifile < maxFiles; ++ifile)
  {
    TString fullName = dirNameSlash + fileList[ifile];
    fin = new TFile(fullName);

    eventCount = NULL;
    cout << " get event Count " << ifile << endl;
    fin->GetObject("eventcount", eventCount);
    if (eventCount)
    {
      ntotal = eventCount->GetBinContent(0);
      npass = eventCount->GetBinContent(1);
      filePass.push_back(npass);
      // 0 = ntriggers, 1 = npass
      printf("total events %f events passed  %f \n ", eventCount->GetBinContent(0), eventCount->GetBinContent(1));
      // TH1D *hevcount = (TH1D *)eventCount->Clone(Form("eventCount%i", ifile));
      // fout->Add(hevcount);
    }
    else
    {
      printf("NO EVENT COUNT IN FILE %s\n", fin->GetName());
      continue;
    }
    TDirectory *anaDir = NULL;
    fin->GetObject("anadir", anaDir); // typeO in directory name
    if (!anaDir)
    {
      printf(" anaDir not found in file %s \n", fin->GetName());
      continue;
    }
    printf(" anaDir for file %s \n", fin->GetName());
    // get hist of cleanup cut pass bit
    fin->GetObject("EventPass", hEventPass);
    if (hEventPass)
      cout << " event pass " << endl;
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
        cout << "chan " << i << " qsum " << hqsum->GetBinContent(i + 1) << " size " << vecQsum[i].size() << endl;
      }
    }

    /******  loop over sumDir *****/
    TList *sumList = sumDir->GetListOfKeys();
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

      /****** get run sum ****/
      TH1D *hClone;
      TH1D *hRunClone;
      if (name.find("sumWave") != std::string::npos && name.find("Bad") == std::string::npos)
      {
        cout << " sumWave clone " << name << " file " << ifile << endl;
        hClone = (TH1D *)h->Clone(Form("%s-file%i", h->GetName(), ifile));
        hSumWave.push_back(hClone);
        waveSumDir->Add(hClone);
      }
      if (name.find("sumHitWave") != std::string::npos)
      {
        cout << " sumHitWave clone " << name << " file " << ifile << endl;
        string chan = name.substr(name.find_last_of("e") + 1);
        int ichan = stoi(chan);
        hClone = (TH1D *)h->Clone(Form("RunHitWaveFile%uChan%i", ifile, ichan));
        TString histName;
        // histName.Form("RunSumHitWaveChan%i", ichan);
        // hRunClone = (TH1D *)h->Clone(histName);
        waveSumDir->Add(hClone);
        vRunHitWave[ichan].push_back(hClone);
        /*  do not sum here
        if (ifile == 0)
        {
            printf("1st ichan = %i, size of hRunSumHitWave = %lu \n", ichan, hRunSumHitWave.size());
            hRunClone->Sumw2();
            hRunSumHitWave.push_back(hRunClone);
            runSumDir->Add(hRunSumHitWave[hRunSumHitWave.size()-1]);
        }
        else
        {
            printf("file %i  ichan = %i, size of hRunSumHitWave = %lu \n", ifile, ichan, hRunSumHitWave.size());
            runSumDir->GetObject(histName, hRunSumHitWave[ichan]);
            hRunSumHitWave[ichan]->Add(hClone);
        }
        */
      }
      if (name.find("sumPeakWave") != std::string::npos)
      {
        cout << " sumPeakWave clone " << name << " file " << ifile << endl;
        string chan = name.substr(name.find_last_of("e") + 1);
        int ichan = stoi(chan);
        hClone = (TH1D *)h->Clone(Form("RunSumPeakWaveFile%uChan%i", ifile, ichan));
        TString histName;
        // histName.Form("RunSumPeakWaveChan%i", ichan);
        // hRunClone = (TH1D *)h->Clone(histName);
        waveSumDir->Add(hClone);
        vRunPeakWave[ichan].push_back(hClone);
        /* do not sum here
        if (ifile == 0)
        {
            printf("1st ichan = %i, size of hRunSumPeakWave = %lu %s \n", ichan, hRunSumPeakWave.size(), hClone->GetName() );
            hRunSumPeakWave.push_back(hRunClone);
            runSumDir->Add(hRunSumPeakWave[hRunSumPeakWave.size() - 1]);
          }
          else
          {
            printf("file %i  ichan = %i, size of hRunSumPeakWave = %lu %s \n", ifile, ichan, hRunSumPeakWave.size(), hClone->GetName());
            //cout << "adding " << histName << endl;
            runSumDir->GetObject(histName, hRunSumPeakWave[ichan]);
            hRunSumPeakWave[ichan]->Add(hClone);
          }
          */
      }

      /****** get PeakChan by channel ****/
      if (name.find("QPeakChan") != std::string::npos)
      {
        string chan = name.substr(name.find_last_of("n") + 1);
        int ichan = stoi(chan);
        // cout << name << "string chan" << chan << " int " << ichan << endl;

        TString cloneName;
        cloneName.Form("runPeakSumCh%ifile%i", ichan, ifile);
        TH1D *hpAdd = (TH1D *)h->Clone(cloneName);
        // cout << " +vRunQSum " << ichan << " " << fileList[ifile] << " " << h->GetName()  << " " << hAdd->GetName() << endl;
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
        // cout << " +vRunQSum " << ichan << " " << fileList[ifile] << " " << h->GetName()  << " " << hAdd->GetName() << endl;
        hqAdd->SetMarkerStyle(20);
        hqAdd->SetMarkerSize(0.5);
        vRunQSum[ichan].push_back(hqAdd);
        qpeSumDir->Add(hqAdd);
      }

    } // loop over sumDir keys
    if (ifile == 0)
      fout->Write();
    // if(ifile==0) fout->ls();
    qpeSumDir->Write();
    waveSumDir->Write();
    fin->Close();
    totalPass += int(npass);
  } // end loop over files
}

void fitSlopes()
{
  printf("Make slope graph.  Number of files = %lu Number of channels = %lu \n", filenum.size(), vRunPeakWave.size());
  for (unsigned ichan = 0; ichan < vRunHitWave.size(); ++ichan)
  {
    if (ichan > 8 && ichan < 12)
      continue; // skip trigger sipms
    printf("Make slope graph ichan %i \n", ichan);
    // cout << " number of files " << vRunPeakWave[ichan].size() << endl;
    fitSumDir->cd();

    double qpeLast = 1.0;
    /* loop over ih = file number */
    for (unsigned ih = 0; ih < vRunPeakWave[ichan].size(); ++ih)
    {
      TString histName;
      histName.Form("RunHitWaveFile%uChan%i", ih, ichan);
      // histName.Form("RunSumPeakWaveFile%uChan%i", ih, ichan);
      // cout << histName << endl;
      waveSumDir->GetObject(histName, vRunPeakWave[ichan][ih]);
      if (vRunPeakWave[ichan][ih] == NULL)
      {
        printf("skipping vRunPeakWave chan %i file %i \n", ichan, ih);
        continue;
      }
      // cout << vRunPeakWave[ichan][ih]->GetName() << endl;
      vRunPeakWave[ichan][ih]->SetXTitle(" digi count 2 ns per bin ");
      vRunPeakWave[ichan][ih]->SetYTitle("summed yield in QPE");
      // new histogram
      int nbinsx = vRunPeakWave[ichan][ih]->GetNbinsX();
      double xlow = vRunPeakWave[ichan][ih]->GetXaxis()->GetBinLowEdge(0);
      double xup = vRunPeakWave[ichan][ih]->GetXaxis()->GetBinUpEdge(nbinsx);
      TH1D *hfitwave = NULL;
      if (theDataType == SIS)
        hfitwave = new TH1D(Form("fitwaveChan%iFlile%i", ichan, ih), Form("fitwaveChan%iFlile%i", ichan, ih), nbinsx, xlow, 8. * xup);
      else
        hfitwave = new TH1D(Form("fitwaveChan%iFlile%i", ichan, ih), Form("fitwaveChan%iFlile%i", ichan, ih), nbinsx, xlow, 2. * xup);
      // cout << hfitwave->GetName() << endl;
      hfitwave->SetMarkerStyle(21);
      hfitwave->SetMarkerSize(0.2);
      hfitwave->GetListOfFunctions()->Clear();

      hfitwave->SetXTitle(" time ns ");
      hfitwave->SetYTitle("summed yield in QPE");

      // normalize to qpe.
      double qpeChan = vecQPE[ichan][ih];
      if (qpeChan <= 0)
        qpeChan = qpeLast;
      else
        qpeLast = vecQPE[ichan][ih];

      printf(" NNNNNN normalize to qpe chan %i file %i qpe %f !!!!!\n", ichan, ih, qpeChan);
      /* loop over histogram bins */
      for (int ibin = 0; ibin < vRunPeakWave[ichan][ih]->GetNbinsX(); ++ibin)
      {
        double xbin = vRunPeakWave[ichan][ih]->GetBinContent(ibin) / qpeChan;
        if (xbin < 0 || xbin > 1.0E9)
          printf("BBBBBBBB bad xbin value QPE %f chan %i file %i \n", qpeChan, ichan, ih);
        vRunPeakWave[ichan][ih]->SetBinContent(ibin, xbin);
        vRunPeakWave[ichan][ih]->SetBinError(ibin, sqrt(xbin));
        hfitwave->SetBinContent(ibin, xbin);
        hfitwave->SetBinError(ibin, sqrt(xbin));
      }

      // add and save in output file runSumDir;
      if (hRunSumHitWave[ichan] == NULL)
      {
        TString histName;
        histName.Form("RunSumHitWaveChan%i", ichan);
        cout << " new clone " << histName << endl;
        hRunSumHitWave[ichan] = (TH1D *)hfitwave->Clone(histName);
        runSumDir->Add(hRunSumHitWave[ichan]);
      }
      else
      {
        // cout << "adding  " << endl;
        hRunSumHitWave[ichan]->Add(hfitwave);
      }

      /// hfitwave->SetDirectory(nullptr);

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

      TF1 *gslopefit = NULL;
      /* chisq fit def, L for likelihood
        L Uses a log likelihood method(default is chi - square method).
        To be used when the histogram represents counts.
      */
      hfitwave->Fit("expo", "LQ", " ", flow, fhigh);
      gslopefit = (TF1 *)hfitwave->GetListOfFunctions()->FindObject("expo");
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
      }
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

  printf(" input args %i \n ", argc);
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
  if (nfiles == 0)
    exit(0);

  printf(" for %s found %lu files \n", tag.Data(), fileList.size());
  maxFiles = fileList.size();
  if (argc > 3)
  {
    maxFiles = atoi(argv[3]);
  }

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
  cout << " starting summary for  " << maxFiles << " files on " << sdate << endl;

  fout = new TFile(Form("summary-type-%i-dir-%s-%s.root", theDataType, dirName.Data(), sdate.c_str()), "recreate");
  fitSumDir = fout->mkdir("fitSumDir");
  waveSumDir = fout->mkdir("WaveSumDir");
  qpeSumDir = fout->mkdir("qpeSumDir");
  peakSumDir = fout->mkdir("peakSumDir");
  runSumDir = fout->mkdir("runSumDir");
  fout->cd();

  hQPEChan = new TH1D("QPEChan", "QPE  by channel", 12, 0, 12);
  hQPESigmaChan = new TH1D("QPESigmaChan", "QPE  by channel", 12, 0, 12);
  hQPEChan->Sumw2();
  hQPESigmaChan->Sumw2();
  effOther.resize(nchan);
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

  hRunSumHitWave.resize(nchan);
  hRunSumPeakWave.resize(nchan);

  /* from 04_14_2023 */
  effOther[0] = 0.215430;
  effOther[1] = 0.292408;
  effOther[2] = 0.246295;
  effOther[3] = 1.;
  effOther[4] = 0.413158;
  effOther[5] = 0.166633;
  effOther[6] = 0.000492;
  effOther[7] = 0.166086;
  effOther[8] = 0.070555;
  effOther[9] = 1.000000;
  effOther[10] = 1.000000;
  effOther[11] = 1.000000;
  effOther[12] = 1.000000;

  for (int ichan = 0; ichan < nchan; ++ichan)
  {
    hRunSumHitWave[ichan] = NULL;
    hRunSumPeakWave[ichan] = NULL;
  }

  fileLoop();

  printf("\t\t >>> files processed << %li  total pass %i <<<<< \n", filenum.size(), totalPass);
  for (unsigned jfile = 0; jfile < filenum.size(); ++jfile)
  {
    int ifile = int(filenum[jfile]);
    printf(" %i %s \n", ifile, fileList[ifile].Data());
  }

  QPEFits();

  // call function to fit slopes and fill vSlope, vESlope
  if (filenum.size() > 0)
  {
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

  cout << "Average QPE" << endl;

  for (int ichan = 0; ichan < nchan; ++ichan)
  {
    double sum = 0;
    for (int ifile = 0; ifile < vecQPE[ichan].size(); ++ifile)
      sum += vecQPE[ichan][ifile];
    printf(" chan %i average over %li QPE %f \n", ichan, vecQPE[ichan].size(), sum / double(vecQPE[ichan].size()));
  }

  printf(" total triggers %i \n", totalPass);

  for (unsigned ichan = 0; ichan < nchan; ++ichan)
    effOther[ichan] = 1.;

  /*
  for (unsigned ichan = 0; ichan < nchan; ++ichan)
  {
    TString histName;
    histName.Form("RunSumHitWaveChan%i", ichan);
    // histName.Form("RunSumPeakWaveFile%uChan%i", ih, ichan);
    runSumDir->GetObject(histName, hRunSumHitWave[ichan]);
    if (hRunSumHitWave[ichan] == NULL)
      continue;
    double eg = effGeo(ichan);
    double nexpect = nPhotons * (totalPass)*eg * effQuantum128;
    double totalHits = hRunSumHitWave[ichan]->Integral();
    effOther[ichan] = totalHits / nexpect;
    printf(" chan %i effGeo %.3E  nexpect %.3E detected  %.3E  effOther %.3f \n", ichan, eg, nexpect, totalHits, totalHits / nexpect);
  }
  */

  // print totalHits
  printf("\t\t >>> files processed << %li  total pass %i <<<<< \n", filenum.size(), totalPass);
  printf("TTTTTTTTTTTTTTTTTT totalHits\n");
  for (unsigned ichan = 0; ichan < nchan; ++ichan)
    printf("totalHits[%i]=%E;\n", ichan, hRunSumHitWave[ichan]->Integral());

  for (unsigned ichan = 0; ichan < nchan; ++ichan)
    printf("effOther[%i]=%f ;\n", ichan, effOther[ichan]);

  //    fout->ls();
  fout->Purge(1);
  fout->Write();
  fout->Close();
  cout << "summary finished " << maxFiles << " files written to " << fout->GetName() << endl;
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
    if (ic != 6 && ic != 3 && ic < 9)
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
  /*
  for (unsigned ic = 0; ic < 9; ++ic)
  {
    for (unsigned ih = 0; ih < vecQsum[ic].size(); ++ih)
    {
      double qpe = vecQPE[ic][ih];
      if (qpe <=0 )
        continue;
      printf("QSUM chan %i file%i qsum %f qpe%f  ratio %f \n", ic, ih, vecQsum[ic][ih], qpe, vecQsum[ic][ih] / qpe);
      vecQsum[ic][ih] = vecQsum[ic][ih] / qpe / effOther[ic];
      vecEQsum[ic][ih] = vecEQsum[ic][ih] / qpe / effOther[ic];
    }
  }
  */

  cout << " graph normalized vecQsum " << endl;

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
    // if (ic == 6 || ic == 7 || ic == 8)
    if (ic != 6 && ic != 3 && ic < 9)
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
    if (ic != 6 && ic != 3 && ic < 9)
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
