#include <ctime>
#include <iostream>
#include <iterator>
#include <locale>
#include <TString.h>
#include <TROOT.h>
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
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
  for (unsigned ic = 0; ic < nchan; ++ic)
    printf(" chan %u vecQsum %lu \n",ic,vecQsum[ic].size());

  printf(" for %s found %lu files \n", tag.Data(), fileList.size());
  Long64_t maxFiles = fileList.size();
  if (argc > 1)
  {
    maxFiles = atoi(argv[2]);
  }

  std::string sdate = currentDate();
  cout << " starting summary for  " << maxFiles << " files on " << sdate << endl;

  TFile *fout = new TFile(Form("summary-%s.root", sdate.c_str()), "recreate");

  for (unsigned ifile = 0; ifile < maxFiles; ++ifile){
    cout << " starting anaRunFile " << fileList[ifile] << endl;
    TString fullName = TString("myData/") + fileList[ifile];
    TFile *fin = new TFile(fullName);
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
    
    // use bf->modified a std string
    /* migrate fit from post.C:        TF1 *g = (TF1*)hHitSum[i]->GetListOfFunctions()->FindObject("expo");*/

    TH1D *hqsum;
    TH1D *hqprompt;
    fin->GetObject("histQsum", hqsum);
    fin->GetObject("histQprompt", hqprompt);

    if (!hqsum || !hqprompt)
      continue;

    cout << "for file  " << ifile << " " << hqsum->GetName() << " " << hqsum->GetNbinsX()   << endl;

    filenum.push_back(double(ifile));
    efilenum.push_back(0);


    for (int i = 0; i < hqsum->GetNbinsX(); ++i)
    {
      cout << "chan " << i + 1 << " qsum " << hqsum->GetBinContent(i + 1) << endl;
      vecQsum[i].push_back(hqsum->GetBinContent(i + 1));
      vecEqsum[i].push_back(hqsum->GetBinError(i + 1));
    }
  } // end loop over files
  printf(" files %lu \n", filenum.size());
  for (unsigned ic = 0; ic < nchan; ++ic)
    printf(" chan %u vecQsum %lu chan 8 \n", ic, vecQsum[ic].size());
  for (unsigned ifile = 0; ifile < vecQsum[8+1].size(); ++ifile )
    printf(" %u %f  ",ifile, vecQsum[8+1][ifile]);
  cout << endl;

  // one graph per channel
  vector<TGraphErrors *> gqsum;
  TMultiGraph *mg = new TMultiGraph();
  for (unsigned ic = 0; ic < nchan; ++ic )
  {
      gqsum.push_back(new TGraphErrors(filenum.size(), &filenum[0], &(vecQsum[ic][0]), &efilenum[0], &(vecEqsum[ic][0])));
      gqsum[ic]->SetName(Form("qsumChan%i",ic));
      gqsum[ic]->SetTitle(Form("qsum-chan-%i",ic));
      mg->Add(gqsum[ic]);
  }
  // overlay all channel graphs on canvas

  
  TCanvas *can = new TCanvas("Qsummary","Qsummary");
  mg->Draw("ap");
  fout->Append(can);
  fout->ls();
  fout->Write();
  cout << "summary finished " << maxFiles <<  " " << fout->GetName() << endl;
  exit(0);
}
