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

//
#include "anaRun.cc"

std::vector<TString> fileList;
std::vector<double> filenum;
std::vector<double> efilenum;
std::vector<double> fileTime;
std::vector<TDatime> fileDatime;
TBFile *bf;
std::vector<vector<double>> vecQsum;
std::vector<vector<double>> vecEQsum;
std::vector<vector<double>> vecQPE;
std::vector<vector<double>> vecEQPE;
std::vector<TH1D *> runQSum;
std::map<int, TH1D *> histMap;
TDirectory *sumDir;
TFile *fin;
TFile *fout;
bool first = true;
TH1D *hQPEChan;
TH1D *eventCount;

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

TDatime getTime()
{
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
  TDatime datime;
  datime.Set(fileYear, fileMonth, fileDay, fileHour, fileMin, fileSec);
  printf("FileYear = %u , FileMonth = %u , FileDay = %u , FileHour = %u , FileMin = %u , FileSec = %u \n", fileYear, fileMonth, fileDay, fileHour, fileMin, fileSec);
  return datime;
}
void QPEFits()
{
  if (!sumDir)
  {
    printf("addSumHistos no sumDir!!! \n");
    return;
  }

  // sumDir->ls();

  TList *sumList = sumDir->GetListOfKeys();
  TIter next(sumList);
  TKey *key;
  printf("QPEFIT %u \n", sumList->GetEntries());
  while (TKey *key = (TKey *)next())
  {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1D"))
      continue;
    TH1D *h = (TH1D *)key->ReadObj();
    std::string name = string(h->GetName());
    if (name.find("QSumChan") == std::string::npos)
      continue;
    string chan = name.substr(name.find_last_of("n") + 1);
    int ichan = stoi(chan);
    cout << " addSumHistos clone " << name << " chan " << ichan << endl;
    TH1D *hclone = (TH1D *)h->Clone(Form("fitQSumCh%i", ichan));
    double xlow = 2000.;
    double xhigh = 5000;
    if (ichan == 9 || ichan == 10 || ichan == 11)
    {
      xlow = 4.E4;
      xhigh = 8.E4;
    }
    hclone->Fit("gaus", " ", " ", xlow, xhigh);
    TF1 *gfit = (TF1 *)hclone->GetListOfFunctions()->FindObject("gaus");
    if (gfit)
    {
      hQPEChan->SetBinContent(ichan, gfit->GetParameter(1));
      hQPEChan->SetBinError(ichan, gfit->GetParError(1));
      printf(" ****** chan %i low %E hight %E QPE %E\n", ichan, xlow, xhigh, gfit->GetParameter(1));
      vecQPE[ichan].push_back(gfit->GetParameter(1));
      vecEQPE[ichan].push_back(gfit->GetParError(1));
    }
    else
    {
      vecQPE[ichan].push_back(0);
      vecEQPE[ichan].push_back(0);
    }
    // fout->Add(hclone);
  }
}

void addSumHistos()
{
  if (!sumDir)
  {
    printf("addSumHistos no sumDir!!! \n");
    return;
  }

  // sumDir->ls();

  TList *sumList = sumDir->GetListOfKeys();
  TIter next(sumList);
  TKey *key;
  printf("addSumHistos %u \n", sumList->GetEntries());
  while (TKey *key = (TKey *)next())
  {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1D"))
      continue;
    TH1D *h = (TH1D *)key->ReadObj();
    std::string name = string(h->GetName());
    if (name.find("QSumChan") == std::string::npos)
      continue;
    string chan = name.substr(name.find_last_of("n") + 1);
    int ichan = stoi(chan);
    if (first)
    {
      cout << " addSumHistos clone " << name << " chan " << ichan << endl;
      TH1D *hclone = (TH1D *)h->Clone(Form("runQSumCh%i", ichan));
      runQSum.push_back(hclone);
      histMap.insert(std::pair<int, TH1D *>(ichan, hclone));
      fout->Add(hclone);
    }
    else
    {
      cout << " addSumHistos add  " << name << " chan " << ichan << endl;
      histMap.at(ichan)->Add(h);
    }
    // name.Form("SumWave-%s-%s", h->GetName(), tag.Data());
    // TH1D *hsave = (TH1D *)h->Clone(name);
    //  rawSumDir->Add(hsave);
  }
  first = false;
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
  TString dirName = TString("myData/");
  TSystemDirectory dir("mydata", dirName); // TSystemDirectory
  TList *files = dir.GetListOfFiles();     //
  TIter next(files);
  TSystemFile *file;
  while ((file = (TSystemFile *)next()))
  {
    string name = string(file->GetName());
    cout << name << endl;
    string exten = name.substr(name.find_last_of(".") + 1);
    if (exten != string("root"))
      continue;
    if (name.find("anaRun-run") != std::string::npos)
      fileList.push_back(TString(name.c_str()));
  }
  return fileList.size();
}

int main(int argc, char *argv[])
{
  cout << "executing " << argv[0] << " make summary plots  " << endl;
  printf(" usage: summary  <max files 0=all>  \n ");
  if (argc < 1)
    exit(0);
  TString tag("run");

  TDatime dopeTime(2023, 3, 9, 22, 0, 0);

  countFiles();
  unsigned nchan = 13;
  vecQsum.resize(nchan);
  vecEQsum.resize(nchan);
  vecQPE.resize(nchan);
  vecEQPE.resize(nchan);
  // for (unsigned ic = 0; ic < nchan; ++ic)
  //   printf(" chan %u vecQsum %lu \n",ic,vecQsum[ic].size());

  printf(" for %s found %lu files \n", tag.Data(), fileList.size());
  Long64_t maxFiles = fileList.size();
  if (argc > 1)
  {
    maxFiles = atoi(argv[2]);
  }

  std::string sdate = currentDate();
  cout << " starting summary for  " << maxFiles << " files on " << sdate << endl;

  fout = new TFile(Form("summary-%s.root", sdate.c_str()), "recreate");

  hQPEChan = new TH1D("QPEChan", "QPE  by channel", 12, 0, 12);
  hQPEChan->Sumw2();

  for (unsigned ifile = 0; ifile < maxFiles; ++ifile)
  {
    cout << " starting anaRunFile " << fileList[ifile] << endl;
    TString fullName = TString("myData/") + fileList[ifile];
    fin = new TFile(fullName);
    TGraph *gslope;
    fin->GetObject("tbfile", bf);
    fin->GetObject("slope-graph", gslope);
    if (!bf)
    {
      // printf(" no timestamp in file %s \n",fileList[i].Data());
      continue;
    }
    if (!gslope)
    {
      // printf(" no timestamp in file %s \n",fileList[i].Data());
      continue;
    }
    cout << bf->GetTitle() << " modified " << bf->modified << " slope graph " << gslope->GetTitle() << endl;

    fin->GetObject("sumDir", sumDir);
    if (!sumDir)
    {
      printf(" sumDir not found in file %s \n", fin->GetName());
      continue;
    }
    printf(" sumDir for file %s \n", fin->GetName());

    // use bf->modified a std string
    /* migrate fit from post.C:        TF1 *g = (TF1*)hHitSum[i]->GetListOfFunctions()->FindObject("expo");*/

    TH1D *hqsum;
    TH1D *hqprompt;
    fin->GetObject("histQsum", hqsum);
    fin->GetObject("histQprompt", hqprompt);

    TH1D *hq = (TH1D *)hqsum->Clone(Form("hqsum%i", ifile));
    TH1D *hp = (TH1D *)hqsum->Clone(Form("hqprompt%i", ifile));
    fout->Add(hq);
    fout->Add(hp);
    eventCount = NULL;
    fin->GetObject("eventCount", eventCount);
    if (eventCount)
    {
      printf(" eventCount for file %s  entries %E \n", fin->GetName(), eventCount->GetEntries());
      TH1D *hevcount = (TH1D *)eventCount->Clone(Form("eventCount%i", ifile));
      fout->Add(hevcount);
    }

    // addSumHistos();
    QPEFits();

    if (!hqsum || !hqprompt)
      continue;

    cout << "for file  " << ifile << " " << fileList[ifile] << " " << hqsum->GetName() << " " << hqsum->GetNbinsX() << endl;

    filenum.push_back(double(ifile));
    efilenum.push_back(0);
    TDatime datetime = getTime();
    fileDatime.push_back(datetime);
    fileTime.push_back(datetime.Convert());
    efilenum.push_back(0);

    cout << " TIME FILE  " << ifile << " " << fileList[ifile] << " modified " << bf->modified << " TDatime " << datetime.AsString() << " as int " << datetime.Convert() << endl;

    if (eventCount)
    {
      for (int i = 0; i <= eventCount->GetNbinsX(); ++i)
      {
        double norm = eventCount->GetBinContent(i + 1) / eventCount->GetBinContent(1);
        printf(" %i count %0f norm %f \n", i - 1, eventCount->GetBinContent(i + 1), norm);
      }
    }

    for (int i = 0; i < hqsum->GetNbinsX() - 1; ++i)
    {
      double norm = 1.0;
      // now overflow is every event and qsum bin is
      if (eventCount)
        norm = eventCount->GetBinContent(i) / eventCount->GetBinContent(0);
      vecQsum[i].push_back(hqsum->GetBinContent(i + 1) / norm);
      vecEQsum[i].push_back(hqsum->GetBinError(i + 1) / norm);
      cout << "chan " << i + 1 << " qsum " << hqsum->GetBinContent(i + 1) << " norm " << norm << " size " << vecQsum[i].size() << endl;
    }
  } // end loop over files
  printf(" files %lu \n", filenum.size());
  // for (unsigned ic = 0; ic < nchan; ++ic)
  //  printf(" chan %u vecQsum %lu  \n", ic, vecQsum[ic].size());
  for (unsigned jfile = 0; jfile < filenum.size(); ++jfile)
  {
    int ifile = int(filenum[jfile]);
    printf("  summary file %u  %s chan6 %f  \n", ifile, fileList[ifile].Data(), vecQsum[7][ifile]);
  }
  int myColor[13] = {41, 42, 43, 44, 45, 46, 2, 3, 4, 31, 32, 33, 34};
  int myStyle[13] = {21, 22, 23, 24, 25, 26, 21, 22, 23, 31, 32, 33, 34};

  // normalize to first file
  vector<double> normQsum;
  normQsum.resize(vecQsum.size());
  for (unsigned ic = 0; ic < vecQsum.size(); ++ic)
  {
    printf(" vecQsum %i %lu \n", ic, vecQsum[ic].size());
    normQsum[ic] = 1;
    if (vecQsum[ic].size() > 0)
    { // ave over before doping
      double beforeSum = 0;
      int normCount = 0;
      for (unsigned jt = 0; jt < 20; ++jt)
      {
        if (!isinf(vecQsum[ic][jt]) && vecQsum[ic][jt] > 0 && fileDatime[jt].Convert() < dopeTime.Convert())
        {
          beforeSum += vecQsum[ic][jt];
          ++normCount;
        }
      }
      normQsum[ic] = beforeSum / double(normCount);
    }
    printf("\t  normQsum =  %f  \n", normQsum[ic]);
  }

  vector<double> normQPE;
  normQPE.resize(vecQPE.size());
  for (unsigned ic = 0; ic < vecQPE.size(); ++ic)
  {
    printf(" vecQPE %i %lu \n", ic, vecQPE[ic].size());
    normQPE[ic] = 1;
    if (vecQPE[ic].size() > 0)
    {
      printf("\t vecQPE =  %f  \n", vecQPE[ic][0]);
      if (!isinf(vecQPE[ic][0]) && vecQPE[ic][0] > 0)
        normQPE[ic] = vecQPE[ic][0];
    }
    printf("\t normQPE =  %f  \n", normQPE[ic]);
  }

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

  // one graph per channel
  vector<TGraphErrors *> gqsum;
  TMultiGraph *mg = new TMultiGraph();
  for (unsigned ic = 0; ic < nchan; ++ic)
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
      mg->Add(gqsum[ic]);
  }
  // overlay all channel graphs on canvas

  // QPE graphs one graph per channel
  vector<TGraphErrors *> gqpe;
  TMultiGraph *mgQPE = new TMultiGraph();

  for (unsigned ic = 0; ic < nchan; ++ic)
  {
    // cout << " add " << ic << endl;
    gqpe.push_back(new TGraphErrors(filenum.size(), &fileTime[0], &(vecQPE[ic][0]), &efilenum[0], &(vecEQPE[ic][0])));
    gqpe[ic]->SetName(Form("QPEChan%i", ic));
    gqpe[ic]->SetTitle(Form("QPE-chan-%i", ic));
    gqpe[ic]->SetMarkerSize(1);
    gqpe[ic]->SetMarkerColor(myColor[ic]);
    gqpe[ic]->SetMarkerStyle(myStyle[ic]);
    fout->Add(gqpe[ic]);
    if (ic == 6 || ic == 7 || ic == 8)
      mgQPE->Add(gqpe[ic]);
  }
  // overlay all channel graphs on canvas
  mg->GetXaxis()->SetTimeDisplay(1);
  // mg->GetXaxis()->SetNdivisions(1010);
  /*
    n = n1 + 100 * n2 + 10000 * n3 Where n1 is the number of primary divisions, n2 is the number of second order divisions and n3 is the number of third order divisions. n < 0, the axis will be forced to use exactly n divisions.
  */
  int ndiv = 10 + 100 * 5 + 10000 * 3;
  mg->GetXaxis()->SetNdivisions(220);
  mg->GetXaxis()->SetTimeFormat("%d:%H");
  mg->GetXaxis()->SetTimeOffset(0, "gmt");
  mg->GetYaxis()->SetTitle("integrated charge (normed)");
  TCanvas *can = new TCanvas("Qsummary", "Qsummary");
  mg->Draw("ap");
  gPad->Update();
  can->BuildLegend();
  can->SetGrid();
  can->Print(".png");
  fout->Append(can);

  mgQPE->GetXaxis()->SetTimeDisplay(1);
  mgQPE->GetXaxis()->SetNdivisions(ndiv);
  mgQPE->GetXaxis()->SetNdivisions(220);
  mgQPE->GetXaxis()->SetTimeFormat("%d:%H");
  mgQPE->GetXaxis()->SetTimeOffset(0, "gmt");
  mgQPE->GetYaxis()->SetTitle(" single photon charge (normed)");
  TCanvas *canqpe = new TCanvas("QPE", "QPE");
  mgQPE->Draw("ap");
  canqpe->BuildLegend();
  canqpe->SetGrid();
  fout->Append(canqpe);
  fout->ls();
  fout->Write();

  for (unsigned it = 0; it < fileTime.size(); ++it)
    if (fileDatime[it].Convert() < dopeTime.Convert())
      cout << " before " << fileDatime[it].AsString() << "  " << fileList[it] << endl;
    else
      cout << " after  " << fileDatime[it].AsString() << "  " << fileList[it] << endl;
  cout << "summary finished " << maxFiles << " " << fout->GetName() << endl;
  exit(0);
}