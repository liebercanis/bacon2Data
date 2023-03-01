#include <ctime>
#include <iostream>
#include <iterator>
#include <locale>
#include <TString.h>
#include <TROOT.h>
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TSystemDirectory.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TList.h"
//
#include "anaRun.cc"

std::vector<TString> fileList;
std::vector<double> filenum;
std::vector<double> efilenum;
std::vector<vector<double>> vecQsum;
std::vector<vector<double>> vecEqsum;
std::vector<TH1D*> runQSum;
std::map<int, TH1D*> histMap;
TDirectory *sumDir;
TFile *fin;
TFile *fout;
bool first = true;

void addSumHistos() {
  if(!sumDir){
    printf("addSumHistos no sumDir!!! \n");
    return;
  }

  //sumDir->ls();

  TList *sumList = sumDir->GetListOfKeys();
  TIter next(sumList);
  TKey *key;
  printf("addSumHistos %u \n",sumList->GetEntries());
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
   if(first){
    cout << " addSumHistos clone "   << name << " chan " << ichan << endl;
    TH1D *hclone = (TH1D *)h->Clone(Form("runQSumCh%i",ichan));
    runQSum.push_back(hclone);
    histMap.insert(std::pair<int,TH1D*>(ichan,hclone));
    fout->Add(hclone);
   } else {
     cout << " addSumHistos add  "   << name << " chan " << ichan << endl;
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
unsigned long  countFiles()
{
  TString dirName = TString("myData/");
  TSystemDirectory dir("mydata", dirName); // TSystemDirectory
  TList* files = dir.GetListOfFiles();            //
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
  if(argc<1)
    exit(0);
  TString tag("run");

  countFiles();
  unsigned nchan = 14;
  vecQsum.resize(nchan);
  vecEqsum.resize(nchan);
  //for (unsigned ic = 0; ic < nchan; ++ic)
  //  printf(" chan %u vecQsum %lu \n",ic,vecQsum[ic].size());

  printf(" for %s found %lu files \n", tag.Data(), fileList.size());
  Long64_t maxFiles = fileList.size();
  if (argc > 1)
  {
    maxFiles = atoi(argv[2]);
  }

  std::string sdate = currentDate();
  cout << " starting summary for  " << maxFiles << " files on " << sdate << endl;

  fout = new TFile(Form("summary-%s.root", sdate.c_str()), "recreate");

  for (unsigned ifile = 0; ifile < maxFiles; ++ifile){
    cout << " starting anaRunFile " << fileList[ifile] << endl;
    TString fullName = TString("myData/") + fileList[ifile];
    fin = new TFile(fullName);
    TBFile *bf;
    TGraph *gslope;
    fin->GetObject("tbfile", bf);
    fin->GetObject("slope-graph", gslope);
    if(!bf) {
      //printf(" no timestamp in file %s \n",fileList[i].Data());
      continue;
    }
    if (!gslope)
    {
      // printf(" no timestamp in file %s \n",fileList[i].Data());
      continue;
    }
    cout << bf->GetTitle() << " modified " << bf->modified << " slope graph " << gslope->GetTitle() << endl;

    fin->GetObject("sumDir", sumDir);
    if(!sumDir) {
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

    TH1D *hq = (TH1D *)hqsum->Clone(Form("hqsum%i",ifile));
    TH1D *hp = (TH1D *)hqsum->Clone(Form("hqprompt%i",ifile));
    fout->Add(hq);
    fout->Add(hp);
    
    addSumHistos();

    if (!hqsum || !hqprompt) continue;

    cout << "for file  " << ifile <<" " << fileList[ifile] << " " << hqsum->GetName() << " " << hqsum->GetNbinsX() << endl;

    filenum.push_back(double(ifile));
    efilenum.push_back(0);


    for (int i = 0; i < hqsum->GetNbinsX()-1; ++i)
    {
      vecQsum[i].push_back(hqsum->GetBinContent(i + 1));
      vecEqsum[i].push_back(hqsum->GetBinError(i + 1));
      cout << "chan " << i + 1 << " qsum " << hqsum->GetBinContent(i + 1) << " size "  << vecQsum[i].size() <<  endl;
    }
  } // end loop over files
  printf(" files %lu \n", filenum.size());
  //for (unsigned ic = 0; ic < nchan; ++ic)
   // printf(" chan %u vecQsum %lu  \n", ic, vecQsum[ic].size());
  for (unsigned jfile = 0; jfile < filenum.size() ; ++jfile ) {
    int ifile = int(filenum[jfile]);
    printf("  summary file %u  %s chan6 %f  \n", ifile, fileList[ifile].Data(), vecQsum[7][ifile]);
  }

  // one graph per channel
  vector<TGraphErrors *> gqsum;
  TMultiGraph *mg = new TMultiGraph();
  for (unsigned ic = 0; ic < nchan; ++ic )
  {
    //cout << " add " << ic << endl; 
    gqsum.push_back(new TGraphErrors(filenum.size(), &filenum[0], &(vecQsum[ic][0]), &efilenum[0], &(vecEqsum[ic][0])));
    gqsum[ic]->SetName(Form("qsumChan%i", ic));
    gqsum[ic]->SetTitle(Form("qsum-chan-%i", ic));
    gqsum[ic]->SetMarkerSize(1);
    gqsum[ic]->SetMarkerStyle(3);
    fout->Add(gqsum[ic]);
    if(ic==6||ic==7||ic==8)
      mg->Add(gqsum[ic]);
  }
  // overlay all channel graphs on canvas

  
  TCanvas *can = new TCanvas("Qsummary","Qsummary");
  mg->Draw("ap");
  can->BuildLegend();
  fout->Append(can);
  fout->ls();
  fout->Write();
  cout << "summary finished " << maxFiles <<  " " << fout->GetName() << endl;
  exit(0);
}
